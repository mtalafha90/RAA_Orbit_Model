from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .kepler import BinaryParams
from .model import (
    GaiaResponseConfig,
    predict_gaia_orbital_al,
    predict_relative_astrometry,
    predict_sb2_rv,
)
from .scanning import GaiaScanSchedule


@dataclass(frozen=True)
class RelativeAstrometryData:
    times_yr: np.ndarray
    values_mas: np.ndarray
    covariance_mas2: np.ndarray


@dataclass(frozen=True)
class SB2RVData:
    times_yr: np.ndarray
    values_kms: np.ndarray
    sigma_kms: np.ndarray


@dataclass(frozen=True)
class GaiaALData:
    times_yr: np.ndarray
    scan_angle_deg: np.ndarray
    values_mas: np.ndarray
    sigma_mas: np.ndarray
    schedule_source: str = ""
    ra_deg: float | None = None
    dec_deg: float | None = None
    release: str = "custom"


@dataclass(frozen=True)
class JointData:
    relative_astrometry: RelativeAstrometryData
    sb2_rv: SB2RVData
    gaia_al: GaiaALData


def _sorted_uniform_times(rng: np.random.Generator, n: int, baseline_yr: float) -> np.ndarray:
    if n <= 0:
        return np.empty(0, dtype=float)
    return np.sort(rng.uniform(0.0, baseline_yr, size=n))


def simulate_joint_data(
    params: BinaryParams,
    gaia_response: GaiaResponseConfig,
    gaia_schedule: GaiaScanSchedule,
    *,
    seed: int = 0,
    n_ast: int = 24,
    n_rv: int = 48,
    baseline_yr: float | None = None,
    ast_sigma_mas: float = 0.20,
    rv_sigma_kms: float = 0.10,
    gaia_sigma_mas: float = 0.10,
) -> JointData:
    """Generate a controlled synthetic joint dataset on an explicit Gaia schedule.

    Gaia epochs and scan angles are never generated internally.  They must be
    supplied through ``GaiaScanSchedule`` so scientific experiments cannot
    silently fall back to uniform random scan angles.  The resolved-astrometry
    and RV epochs remain controlled random samples; by default their time span
    is the mission-relative span through the final supplied Gaia transit.
    """
    params.validate()
    gaia_schedule.validate()
    for name, value in (
        ("ast_sigma_mas", ast_sigma_mas),
        ("rv_sigma_kms", rv_sigma_kms),
        ("gaia_sigma_mas", gaia_sigma_mas),
    ):
        if value <= 0:
            raise ValueError(f"{name} must be > 0")

    if baseline_yr is None:
        baseline = gaia_schedule.mission_span_yr
    else:
        baseline = float(baseline_yr)
    if baseline <= 0 and (n_ast > 0 or n_rv > 0):
        raise ValueError("baseline_yr must be > 0 when astrometry or RV data are simulated")

    rng = np.random.default_rng(seed)

    t_ast = _sorted_uniform_times(rng, n_ast, baseline)
    ast_true = predict_relative_astrometry(t_ast, params)
    ast_cov = np.repeat(np.eye(2)[None, :, :] * ast_sigma_mas**2, n_ast, axis=0)
    ast_obs = ast_true + rng.normal(0.0, ast_sigma_mas, size=ast_true.shape)

    t_rv = _sorted_uniform_times(rng, n_rv, baseline)
    rv_true = predict_sb2_rv(t_rv, params)
    rv_sig = np.full_like(rv_true, rv_sigma_kms, dtype=float)
    rv_obs = rv_true + rng.normal(0.0, rv_sig)

    t_gaia = np.asarray(gaia_schedule.times_yr, dtype=float).copy()
    scan = np.asarray(gaia_schedule.scan_angle_deg, dtype=float).copy()
    gaia_true = predict_gaia_orbital_al(t_gaia, scan, params, gaia_response)
    gaia_sig = np.full(len(t_gaia), gaia_sigma_mas, dtype=float)
    gaia_obs = gaia_true + rng.normal(0.0, gaia_sig)

    return JointData(
        relative_astrometry=RelativeAstrometryData(t_ast, ast_obs, ast_cov),
        sb2_rv=SB2RVData(t_rv, rv_obs, rv_sig),
        gaia_al=GaiaALData(
            t_gaia,
            scan,
            gaia_obs,
            gaia_sig,
            schedule_source=gaia_schedule.source,
            ra_deg=gaia_schedule.ra_deg,
            dec_deg=gaia_schedule.dec_deg,
            release=gaia_schedule.release,
        ),
    )
