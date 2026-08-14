from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Iterable
import numpy as np
from scipy.optimize import least_squares

from .kepler import BinaryParams
from .likelihoods import gaussian_1d_loglike, gaussian_2d_loglike
from .model import (
    GaiaResponseConfig,
    predict_gaia_orbital_al,
    predict_relative_astrometry,
    predict_sb2_rv,
)
from .synthetic import JointData


ALL_PARAMETER_NAMES = (
    "period_yr", "t_peri_yr", "eccentricity", "inclination_deg",
    "omega_deg", "node_deg", "m1_msun", "m2_msun", "parallax_mas",
    "gamma_kms", "beta_g",
)

# Instrument-response parameters. These do not belong to BinaryParams because
# they describe the measurement, not the binary. They may be included in
# ``free_names`` to fit or marginalise over an imperfectly known Gaia response
# instead of asserting that its width is known exactly.
RESPONSE_PARAMETER_NAMES = ("sigma_response_mas",)

# Absolute astrometric parameters of the barycentre. They are deliberately NOT
# in the default free set: including them would change every existing
# experiment. Request them explicitly when the Gaia channel carries a sky
# position and epoch, so that the orbit is fitted jointly with parallax and
# proper motion rather than in isolation.
ASTROMETRIC_PARAMETER_NAMES = (
    "pmra_mas_yr", "pmdec_mas_yr", "delta_alpha_star_mas", "delta_delta_mas",
)

ALL_FREE_NAMES = (
    ALL_PARAMETER_NAMES + ASTROMETRIC_PARAMETER_NAMES + RESPONSE_PARAMETER_NAMES
)


@dataclass(frozen=True)
class JointFitResult:
    params: BinaryParams
    success: bool
    message: str
    chi2: float
    n_residuals: int
    n_free: int
    nfev: int
    cost: float
    # Populated only when "sigma_response_mas" was fitted rather than asserted.
    fitted_sigma_response_mas: float | None = None

    @property
    def dof(self) -> int:
        return self.n_residuals - self.n_free

    @property
    def reduced_chi2(self) -> float:
        return self.chi2 / self.dof if self.dof > 0 else np.nan


def _bounds_for(
    name: str,
    initial: BinaryParams,
    sigma_response_mas: float | None = None,
) -> tuple[float, float]:
    if name == "sigma_response_mas":
        if not sigma_response_mas or sigma_response_mas <= 0:
            raise ValueError("fitting sigma_response_mas needs a positive starting width")
        return 0.2 * sigma_response_mas, 5.0 * sigma_response_mas
    if name == "period_yr":
        return max(1e-5, 0.2 * initial.period_yr), 5.0 * initial.period_yr
    if name == "t_peri_yr":
        return initial.t_peri_yr - 2.0 * initial.period_yr, initial.t_peri_yr + 2.0 * initial.period_yr
    if name == "eccentricity":
        return 0.0, 0.98
    if name == "inclination_deg":
        return 0.01, 179.99
    if name in ("omega_deg", "node_deg"):
        return -720.0, 720.0
    if name in ("m1_msun", "m2_msun"):
        return 0.02, 20.0
    if name == "parallax_mas":
        return 0.01, 1000.0
    if name == "gamma_kms":
        return -1000.0, 1000.0
    if name == "beta_g":
        return 0.0, 0.5
    if name in ("pmra_mas_yr", "pmdec_mas_yr"):
        return -1000.0, 1000.0
    if name in ("delta_alpha_star_mas", "delta_delta_mas"):
        return -1000.0, 1000.0
    raise KeyError(name)


def _pack(params: BinaryParams, names: tuple[str, ...]) -> np.ndarray:
    return np.array([float(getattr(params, name)) for name in names], dtype=float)


def _unpack(base: BinaryParams, names: tuple[str, ...], values: np.ndarray) -> BinaryParams:
    return replace(base, **{name: float(value) for name, value in zip(names, values)})


def joint_residuals(params: BinaryParams, data: JointData, gaia_response: GaiaResponseConfig) -> np.ndarray:
    """Return whitened residuals for all three data channels."""
    params.validate()
    pieces: list[np.ndarray] = []

    ast = data.relative_astrometry
    if len(ast.times_yr):
        pred = predict_relative_astrometry(ast.times_yr, params)
        raw = ast.values_mas - pred
        L = np.linalg.cholesky(ast.covariance_mas2)
        white = np.linalg.solve(L, raw[..., None])[..., 0]
        pieces.append(white.ravel())

    rv = data.sb2_rv
    if len(rv.times_yr):
        pred = predict_sb2_rv(rv.times_yr, params)
        pieces.append(((rv.values_kms - pred) / rv.sigma_kms).ravel())

    gaia = data.gaia_al
    if len(gaia.times_yr):
        pred = predict_gaia_orbital_al(
            gaia.times_yr, gaia.scan_angle_deg, params, gaia_response, gaia.astrometry
        )
        pieces.append((gaia.values_mas - pred) / gaia.sigma_mas)

    if not pieces:
        raise ValueError("joint dataset contains no observations")
    return np.concatenate(pieces)


def joint_loglike(
    params: BinaryParams,
    data: JointData,
    gaia_response: GaiaResponseConfig,
    *,
    gaia_jitter_mas: float = 0.0,
    rv_jitter_kms: float = 0.0,
) -> float:
    """Gaussian log-likelihood of all three channels.

    ``joint_residuals`` returns whitened residuals, whose sum of squares is a
    chi-square. That is sufficient for least squares, where the normalisation
    is constant, but not for sampling with free jitter terms, where the
    ``log(variance)`` penalty is what stops the jitter running away. This
    routine therefore builds the full likelihood, including normalisation.
    """
    params.validate()
    if gaia_jitter_mas < 0 or rv_jitter_kms < 0:
        raise ValueError("jitter terms must be >= 0")
    total = 0.0

    ast = data.relative_astrometry
    if len(ast.times_yr):
        total += gaussian_2d_loglike(
            ast.values_mas,
            predict_relative_astrometry(ast.times_yr, params),
            ast.covariance_mas2,
        )

    rv = data.sb2_rv
    if len(rv.times_yr):
        total += gaussian_1d_loglike(
            rv.values_kms,
            predict_sb2_rv(rv.times_yr, params),
            rv.sigma_kms,
            jitter=rv_jitter_kms,
        )

    gaia = data.gaia_al
    if len(gaia.times_yr):
        total += gaussian_1d_loglike(
            gaia.values_mas,
            predict_gaia_orbital_al(
                gaia.times_yr, gaia.scan_angle_deg, params, gaia_response, gaia.astrometry
            ),
            gaia.sigma_mas,
            jitter=gaia_jitter_mas,
        )
    return float(total)


def fit_joint(
    data: JointData,
    initial: BinaryParams,
    gaia_response: GaiaResponseConfig,
    *,
    free_names: Iterable[str] = ALL_PARAMETER_NAMES,
    max_nfev: int = 3000,
) -> JointFitResult:
    """Bounded deterministic fit used for controlled injection/recovery tests.

    This optimizer is deliberately not the final posterior engine. Its role is
    to isolate forward-model bias while keeping inference mechanics simple.
    """
    initial.validate()
    names = tuple(free_names)
    unknown = set(names) - set(ALL_FREE_NAMES)
    if unknown:
        raise ValueError(f"unknown free parameter(s): {sorted(unknown)}")
    if len(set(names)) != len(names):
        raise ValueError("free_names contains duplicates")

    binary_field_names = ALL_PARAMETER_NAMES + ASTROMETRIC_PARAMETER_NAMES
    physical_names = tuple(n for n in names if n in binary_field_names)
    response_names = tuple(n for n in names if n in RESPONSE_PARAMETER_NAMES)
    if response_names and gaia_response.mode == "photocentre":
        raise ValueError("the photocentre response has no width to fit")

    x0 = np.concatenate([
        _pack(initial, physical_names),
        np.array([float(gaia_response.sigma_mas) for _ in response_names], dtype=float),
    ])
    if not names:
        r = joint_residuals(initial, data, gaia_response)
        return JointFitResult(initial, True, "no free parameters", float(r @ r), len(r), 0, 0, 0.5 * float(r @ r))

    bounds = np.array(
        [_bounds_for(name, initial, gaia_response.sigma_mas) for name in physical_names + response_names],
        dtype=float,
    )
    lb, ub = bounds[:, 0], bounds[:, 1]
    n_physical = len(physical_names)

    def residual_vector(x):
        p = _unpack(initial, physical_names, x[:n_physical])
        response = gaia_response
        if response_names:
            response = replace(response, sigma_mas=float(x[n_physical]))
        return joint_residuals(p, data, response)

    result = least_squares(
        residual_vector,
        x0,
        bounds=(lb, ub),
        method="trf",
        x_scale="jac",
        ftol=1e-11,
        xtol=1e-11,
        gtol=1e-11,
        max_nfev=max_nfev,
    )
    fitted = _unpack(initial, physical_names, result.x[:n_physical])
    r = result.fun
    return JointFitResult(
        params=fitted,
        success=bool(result.success),
        message=str(result.message),
        chi2=float(r @ r),
        n_residuals=len(r),
        n_free=len(names),
        nfev=int(result.nfev),
        cost=float(result.cost),
        fitted_sigma_response_mas=(float(result.x[n_physical]) if response_names else None),
    )
