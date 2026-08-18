from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

import numpy as np
from scipy.optimize import least_squares

from .kepler import BinaryParams, radial_velocities_kms, relative_astrometry_mas

LEGACY_EPOCH_ZERO = 15020.31352
LEGACY_DAYS_PER_YEAR = 365.242198781


@dataclass(frozen=True)
class RealBinaryData:
    object_name: str
    metadata: Mapping[str, str]
    visual_time_yr: np.ndarray
    visual_east_mas: np.ndarray
    visual_north_mas: np.ndarray
    visual_sigma_mas: np.ndarray
    rv1_time_yr: np.ndarray
    rv1_kms: np.ndarray
    rv1_sigma_kms: np.ndarray
    rv2_time_yr: np.ndarray
    rv2_kms: np.ndarray
    rv2_sigma_kms: np.ndarray

    @property
    def n_scalar_constraints(self) -> int:
        return 2 * len(self.visual_time_yr) + len(self.rv1_time_yr) + len(self.rv2_time_yr)


@dataclass(frozen=True)
class RealBinaryFit:
    params: BinaryParams
    free_names: tuple[str, ...]
    chi2: float
    dof: int
    covariance: np.ndarray
    uncertainties: Mapping[str, float]
    chi2_astrometry: float
    chi2_rv1: float
    chi2_rv2: float
    success: bool
    message: str

    @property
    def reduced_chi2(self) -> float:
        return self.chi2 / self.dof


def legacy_rv_epoch_to_year(epoch: float | np.ndarray) -> np.ndarray:
    """Convert the numerical RV epochs used by the legacy PySVOrbit file to decimal year."""
    return 1900.0 + (np.asarray(epoch, dtype=float) - LEGACY_EPOCH_ZERO) / LEGACY_DAYS_PER_YEAR


def _normalise_metadata_key(key: str) -> str:
    k = key.strip().lower()
    if k in {"par", "plx", "parallax"}:
        return "parallax"
    return k


def parse_legacy_binary_csv(path: str | Path) -> RealBinaryData:
    """Parse the legacy visual/SB2 CSV format used by GL765_Test1.csv.

    Metadata accepts ``par``, ``plx``, and ``parallax`` as aliases. Visual rows
    use ``epoch, theta_deg, rho_arcsec, sigma_arcsec, I1, ...`` and are converted
    to tangent-plane East/North coordinates with an isotropic sigma. RV rows use
    ``epoch, velocity_kms, sigma_kms, Va|Vb, ...``.
    """
    metadata: dict[str, str] = {}
    visual: list[tuple[float, float, float, float]] = []
    rv1: list[tuple[float, float, float]] = []
    rv2: list[tuple[float, float, float]] = []

    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or line.startswith("C"):
            continue
        parts = [p.strip() for p in line.split(",")]
        if len(parts) == 2:
            metadata[_normalise_metadata_key(parts[0])] = parts[1]
            continue
        if len(parts) >= 5 and parts[3] in {"Va", "Vb"}:
            row = (float(parts[0]), float(parts[1]), float(parts[2]))
            (rv1 if parts[3] == "Va" else rv2).append(row)
            continue
        if len(parts) >= 6 and parts[4] == "I1":
            visual.append((float(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])))

    if not visual or not rv1 or not rv2:
        raise ValueError("legacy file must contain visual I1 rows and both Va/Vb RV rows")

    vis = np.asarray(visual, dtype=float)
    theta = np.deg2rad(vis[:, 1])
    rho_mas = 1000.0 * vis[:, 2]
    sigma_mas = 1000.0 * vis[:, 3]
    east_mas = rho_mas * np.sin(theta)
    north_mas = rho_mas * np.cos(theta)

    r1 = np.asarray(rv1, dtype=float)
    r2 = np.asarray(rv2, dtype=float)
    return RealBinaryData(
        object_name=metadata.get("object", Path(path).stem),
        metadata=metadata,
        visual_time_yr=vis[:, 0],
        visual_east_mas=east_mas,
        visual_north_mas=north_mas,
        visual_sigma_mas=sigma_mas,
        rv1_time_yr=legacy_rv_epoch_to_year(r1[:, 0]),
        rv1_kms=r1[:, 1],
        rv1_sigma_kms=r1[:, 2],
        rv2_time_yr=legacy_rv_epoch_to_year(r2[:, 0]),
        rv2_kms=r2[:, 1],
        rv2_sigma_kms=r2[:, 2],
    )


FREE_NAMES = (
    "period_yr",
    "t_peri_yr",
    "eccentricity",
    "inclination_deg",
    "omega_deg",
    "node_deg",
    "m1_msun",
    "m2_msun",
    "parallax_mas",
    "gamma_kms",
)


def _params_from_vector(x: np.ndarray, *, fixed_parallax_mas: float | None = None) -> BinaryParams:
    vals = list(map(float, x))
    if fixed_parallax_mas is None:
        P, T, e, inc, omega, node, m1, m2, plx, gamma = vals
    else:
        P, T, e, inc, omega, node, m1, m2, gamma = vals
        plx = float(fixed_parallax_mas)
    return BinaryParams(P, T, e, inc, omega, node, m1, m2, plx, gamma, beta_g=0.0)


def _residual_blocks(data: RealBinaryData, p: BinaryParams) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    ast = relative_astrometry_mas(data.visual_time_yr, p)
    obs = np.column_stack((data.visual_east_mas, data.visual_north_mas))
    r_ast = ((ast - obs) / data.visual_sigma_mas[:, None]).ravel()
    rv1_model = radial_velocities_kms(data.rv1_time_yr, p)[:, 0]
    rv2_model = radial_velocities_kms(data.rv2_time_yr, p)[:, 1]
    r_rv1 = (rv1_model - data.rv1_kms) / data.rv1_sigma_kms
    r_rv2 = (rv2_model - data.rv2_kms) / data.rv2_sigma_kms
    return r_ast, r_rv1, r_rv2


def fit_visual_sb2(
    data: RealBinaryData,
    *,
    fixed_parallax_mas: float | None = None,
    initial: BinaryParams | None = None,
) -> RealBinaryFit:
    """Bounded least-squares fit of visual relative astrometry plus both SB2 curves.

    No Gaia/light-ratio parameter is included because this data set has no Gaia
    measurements; this avoids the unconstrained ``beta_g`` direction that would
    otherwise be present in a generic joint fit.
    """
    if initial is None:
        initial = BinaryParams(
            period_yr=11.8,
            t_peri_yr=1993.2,
            eccentricity=0.24,
            inclination_deg=81.0,
            omega_deg=72.0,
            node_deg=289.0,
            m1_msun=0.78,
            m2_msun=0.81,
            parallax_mas=35.0 if fixed_parallax_mas is None else fixed_parallax_mas,
            gamma_kms=-4.1,
            beta_g=0.0,
        )

    full = np.array([
        initial.period_yr, initial.t_peri_yr, initial.eccentricity,
        initial.inclination_deg, initial.omega_deg, initial.node_deg,
        initial.m1_msun, initial.m2_msun, initial.parallax_mas, initial.gamma_kms,
    ], dtype=float)
    lower = np.array([5.0, 1980.0, 0.0, 0.0, -360.0, 0.0, 0.1, 0.1, 1.0, -30.0])
    upper = np.array([20.0, 2005.0, 0.95, 180.0, 360.0, 360.0, 3.0, 3.0, 100.0, 30.0])

    if fixed_parallax_mas is None:
        x0, lo, hi = full, lower, upper
        free_names = FREE_NAMES
    else:
        keep = [0, 1, 2, 3, 4, 5, 6, 7, 9]
        x0, lo, hi = full[keep], lower[keep], upper[keep]
        free_names = tuple(n for n in FREE_NAMES if n != "parallax_mas")

    def residual(x: np.ndarray) -> np.ndarray:
        p = _params_from_vector(x, fixed_parallax_mas=fixed_parallax_mas)
        return np.concatenate(_residual_blocks(data, p))

    result = least_squares(
        residual, x0, bounds=(lo, hi), method="trf",
        xtol=1e-13, ftol=1e-13, gtol=1e-13, max_nfev=100000,
    )
    p = _params_from_vector(result.x, fixed_parallax_mas=fixed_parallax_mas)
    r_ast, r_rv1, r_rv2 = _residual_blocks(data, p)
    chi2_ast = float(r_ast @ r_ast)
    chi2_rv1 = float(r_rv1 @ r_rv1)
    chi2_rv2 = float(r_rv2 @ r_rv2)
    chi2 = chi2_ast + chi2_rv1 + chi2_rv2
    dof = data.n_scalar_constraints - len(result.x)

    jtj = result.jac.T @ result.jac
    covariance = np.linalg.inv(jtj)
    sigma = np.sqrt(np.diag(covariance))
    uncertainties = {name: float(s) for name, s in zip(free_names, sigma)}
    if fixed_parallax_mas is not None:
        uncertainties["parallax_mas"] = 0.0

    return RealBinaryFit(
        params=p,
        free_names=free_names,
        chi2=chi2,
        dof=dof,
        covariance=covariance,
        uncertainties=uncertainties,
        chi2_astrometry=chi2_ast,
        chi2_rv1=chi2_rv1,
        chi2_rv2=chi2_rv2,
        success=bool(result.success),
        message=str(result.message),
    )
