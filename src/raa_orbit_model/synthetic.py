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
    *,
    seed: int = 0,
    n_ast: int = 24,
    n_rv: int = 48,
    n_gaia: int = 72,
    baseline_periods: float = 2.5,
    ast_sigma_mas: float = 0.20,
    rv_sigma_kms: float = 0.10,
    gaia_sigma_mas: float = 0.10,
) -> JointData:
    """Generate a controlled synthetic joint dataset.

    Scan angles are sampled uniformly on [0, 180) degrees. This is deliberately
    *not* a Gaia scanning-law simulator; it isolates measurement-model bias from
    scanning-law selection effects in the first injection/recovery experiments.
    """
    params.validate()
    for name, value in (
        ("ast_sigma_mas", ast_sigma_mas),
        ("rv_sigma_kms", rv_sigma_kms),
        ("gaia_sigma_mas", gaia_sigma_mas),
    ):
        if value <= 0:
            raise ValueError(f"{name} must be > 0")
    if baseline_periods <= 0:
        raise ValueError("baseline_periods must be > 0")

    rng = np.random.default_rng(seed)
    baseline = baseline_periods * params.period_yr

    t_ast = _sorted_uniform_times(rng, n_ast, baseline)
    ast_true = predict_relative_astrometry(t_ast, params)
    ast_cov = np.repeat(np.eye(2)[None, :, :] * ast_sigma_mas**2, n_ast, axis=0)
    ast_obs = ast_true + rng.normal(0.0, ast_sigma_mas, size=ast_true.shape)

    t_rv = _sorted_uniform_times(rng, n_rv, baseline)
    rv_true = predict_sb2_rv(t_rv, params)
    rv_sig = np.full_like(rv_true, rv_sigma_kms, dtype=float)
    rv_obs = rv_true + rng.normal(0.0, rv_sig)

    t_gaia = _sorted_uniform_times(rng, n_gaia, baseline)
    scan = rng.uniform(0.0, 180.0, size=n_gaia)
    gaia_true = predict_gaia_orbital_al(t_gaia, scan, params, gaia_response)
    gaia_sig = np.full(n_gaia, gaia_sigma_mas, dtype=float)
    gaia_obs = gaia_true + rng.normal(0.0, gaia_sig)

    return JointData(
        relative_astrometry=RelativeAstrometryData(t_ast, ast_obs, ast_cov),
        sb2_rv=SB2RVData(t_rv, rv_obs, rv_sig),
        gaia_al=GaiaALData(t_gaia, scan, gaia_obs, gaia_sig),
    )
