from __future__ import annotations

from dataclasses import dataclass

from .kepler import BinaryParams, relative_astrometry_mas, radial_velocities_kms
from .gaia import project_along_scan, photocentre_along_scan, blended_gaussian_peak


@dataclass(frozen=True)
class GaiaResponseConfig:
    mode: str = "photocentre"  # "photocentre" or "blended_gaussian_peak"
    sigma_mas: float | None = None


def predict_relative_astrometry(times_yr, params: BinaryParams):
    return relative_astrometry_mas(times_yr, params)


def predict_sb2_rv(times_yr, params: BinaryParams):
    return radial_velocities_kms(times_yr, params)


def predict_gaia_orbital_al(times_yr, scan_angle_deg, params: BinaryParams,
                            response: GaiaResponseConfig = GaiaResponseConfig()):
    rel = relative_astrometry_mas(times_yr, params)
    d_al = project_along_scan(rel[:, 0], rel[:, 1], scan_angle_deg)
    B = params.mass_fraction_secondary
    if response.mode == "photocentre":
        return photocentre_along_scan(d_al, B, params.beta_g)
    if response.mode == "blended_gaussian_peak":
        if response.sigma_mas is None:
            raise ValueError("sigma_mas is required for blended_gaussian_peak mode")
        return blended_gaussian_peak(d_al, B, params.beta_g, response.sigma_mas)
    raise ValueError(f"unknown Gaia response mode: {response.mode}")
