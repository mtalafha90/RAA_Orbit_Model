from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Iterable
import numpy as np
from scipy.optimize import least_squares

from .kepler import BinaryParams
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

    @property
    def dof(self) -> int:
        return self.n_residuals - self.n_free

    @property
    def reduced_chi2(self) -> float:
        return self.chi2 / self.dof if self.dof > 0 else np.nan


def _bounds_for(name: str, initial: BinaryParams) -> tuple[float, float]:
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
        pred = predict_gaia_orbital_al(gaia.times_yr, gaia.scan_angle_deg, params, gaia_response)
        pieces.append((gaia.values_mas - pred) / gaia.sigma_mas)

    if not pieces:
        raise ValueError("joint dataset contains no observations")
    return np.concatenate(pieces)


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
    unknown = set(names) - set(ALL_PARAMETER_NAMES)
    if unknown:
        raise ValueError(f"unknown free parameter(s): {sorted(unknown)}")
    if len(set(names)) != len(names):
        raise ValueError("free_names contains duplicates")

    x0 = _pack(initial, names)
    if not names:
        r = joint_residuals(initial, data, gaia_response)
        return JointFitResult(initial, True, "no free parameters", float(r @ r), len(r), 0, 0, 0.5 * float(r @ r))

    bounds = np.array([_bounds_for(name, initial) for name in names], dtype=float)
    lb, ub = bounds[:, 0], bounds[:, 1]

    def residual_vector(x):
        p = _unpack(initial, names, x)
        return joint_residuals(p, data, gaia_response)

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
    fitted = _unpack(initial, names, result.x)
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
    )
