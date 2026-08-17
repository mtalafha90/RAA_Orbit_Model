from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .astrometry import absolute_offsets_mas
from .kepler import BinaryParams, relative_astrometry_mas, radial_velocities_kms
from .gaia import (
    _global_peak_coordinate,
    _unequal_width_peak_and_modes,
    blended_gaussian_peak,
    blended_gaussian_response,
    photocentre_along_scan,
    project_along_scan,
)
from .penoyre import (
    penoyre_oriented_gaussian_peak,
    penoyre_oriented_gaussian_response,
)


@dataclass(frozen=True)
class GaiaResponseConfig:
    mode: str = "photocentre"
    # For ``blended_gaussian_peak`` this is the 1-D research width.  For
    # ``penoyre_gaussian_peak`` it is Penoyre's idealised along-scan/narrow-axis
    # width alpha.  It is never a calibrated Gaia PLSF width in this project.
    sigma_mas: float | None = None
    allow_multi_peak_continuation: bool = False
    # Optional independent width for the secondary component. ``None`` keeps the
    # equal-width 1-D surrogate used for all previously frozen results. Setting
    # it moves the response outside that family for controlled profile-shape
    # misspecification experiments.
    sigma_secondary_mas: float | None = None
    # Across-scan/long-axis width beta of the idealised elongated Gaussian in
    # ``penoyre_gaussian_peak`` mode.  beta=+inf is the 1-D Lindegren/gaiamock
    # limit; finite values introduce orientation dependence.  This is a research
    # parameter, not a Gaia calibration value.
    sigma_ac_mas: float | None = None


@dataclass(frozen=True)
class GaiaALPrediction:
    """AL prediction plus validity metadata for the single-coordinate channel."""

    values_mas: np.ndarray
    single_peak_mask: np.ndarray
    n_peaks: np.ndarray
    projected_separation_mas: np.ndarray
    critical_separation_sigma: float

    @property
    def n_single_peak(self) -> int:
        return int(np.count_nonzero(self.single_peak_mask))

    @property
    def n_multi_peak(self) -> int:
        return int(len(self.single_peak_mask) - self.n_single_peak)


def predict_relative_astrometry(times_yr, params: BinaryParams):
    return relative_astrometry_mas(times_yr, params)


def predict_sb2_rv(times_yr, params: BinaryParams):
    return radial_velocities_kms(times_yr, params)


@dataclass(frozen=True)
class AbsoluteAstrometryConfig:
    """Sky position and epoch needed to evaluate the parallactic ellipse.

    Supplying this to the Gaia channel switches on the five-parameter absolute
    motion: position offset, proper motion and the parallax factors. Without it
    the channel models the orbital wobble alone, which is the configuration all
    previously frozen results were produced with.
    """

    ra_deg: float
    dec_deg: float
    mission_start_decimalyear: float
    reference_time_yr: float = 0.0


def _relative_al(times_yr, scan_angle_deg, params: BinaryParams):
    rel = relative_astrometry_mas(times_yr, params)
    return project_along_scan(rel[:, 0], rel[:, 1], scan_angle_deg)


def absolute_al(times_yr, scan_angle_deg, params: BinaryParams, astrometry):
    """Along-scan projection of the barycentre's absolute motion."""
    if astrometry is None:
        return 0.0
    alpha_star, delta = absolute_offsets_mas(
        times_yr,
        params,
        ra_deg=astrometry.ra_deg,
        dec_deg=astrometry.dec_deg,
        mission_start_decimalyear=astrometry.mission_start_decimalyear,
        reference_time_yr=astrometry.reference_time_yr,
    )
    return project_along_scan(alpha_star, delta, scan_angle_deg)


def predict_gaia_orbital_response(
    times_yr,
    scan_angle_deg,
    params: BinaryParams,
    response: GaiaResponseConfig = GaiaResponseConfig(),
) -> GaiaALPrediction:
    """Predict Gaia-like AL response while preserving single/multi-peak validity.

    For the photocentre model every epoch belongs to the single-coordinate
    channel by construction. For either blended-Gaussian research surrogate,
    multi-peak epochs are returned with NaN values and a false validity mask;
    callers must not reinterpret those epochs as unique centroids.
    """
    d_al = np.asarray(_relative_al(times_yr, scan_angle_deg, params), dtype=float)
    B = params.mass_fraction_secondary

    if response.mode == "photocentre":
        return GaiaALPrediction(
            values_mas=np.asarray(photocentre_along_scan(d_al, B, params.beta_g), dtype=float),
            single_peak_mask=np.ones(len(d_al), dtype=bool),
            n_peaks=np.ones(len(d_al), dtype=int),
            projected_separation_mas=np.abs(d_al),
            critical_separation_sigma=math.inf,
        )

    if response.mode == "blended_gaussian_peak":
        if response.sigma_mas is None:
            raise ValueError("sigma_mas is required for blended_gaussian_peak mode")
        result = blended_gaussian_response(
            d_al, B, params.beta_g, response.sigma_mas,
            sigma_secondary_mas=response.sigma_secondary_mas,
        )
        return GaiaALPrediction(
            values_mas=np.asarray(result.al_mas, dtype=float),
            single_peak_mask=np.asarray(result.single_peak_mask, dtype=bool),
            n_peaks=np.asarray(result.n_peaks, dtype=int),
            projected_separation_mas=np.abs(d_al),
            critical_separation_sigma=float(result.critical_separation_sigma),
        )

    if response.mode == "penoyre_gaussian_peak":
        if response.sigma_mas is None:
            raise ValueError("sigma_mas (alpha) is required for penoyre_gaussian_peak mode")
        if response.sigma_ac_mas is None:
            raise ValueError("sigma_ac_mas (beta) is required for penoyre_gaussian_peak mode")
        rel = relative_astrometry_mas(times_yr, params)
        result = penoyre_oriented_gaussian_response(
            rel[:, 0],
            rel[:, 1],
            scan_angle_deg,
            B,
            params.beta_g,
            response.sigma_mas,
            response.sigma_ac_mas,
        )
        return GaiaALPrediction(
            values_mas=np.asarray(result.al_mas, dtype=float),
            single_peak_mask=np.asarray(result.single_peak_mask, dtype=bool),
            n_peaks=np.asarray(result.n_peaks, dtype=int),
            projected_separation_mas=np.abs(d_al),
            critical_separation_sigma=float(result.critical_separation_sigma),
        )

    raise ValueError(f"unknown Gaia response mode: {response.mode}")


def predict_gaia_orbital_al(
    times_yr,
    scan_angle_deg,
    params: BinaryParams,
    response: GaiaResponseConfig = GaiaResponseConfig(),
    astrometry: AbsoluteAstrometryConfig | None = None,
):
    """Predict a unique AL coordinate in a surrogate single-peak domain.

    ``allow_multi_peak_continuation`` exists only so the deterministic
    least-squares prototype can traverse parameter space continuously. Such an
    intermediate coordinate is not a scientific prediction: final solutions
    must be checked with :func:`predict_gaia_orbital_response` and accepted only
    when every retained epoch is single-peaked.
    """
    d_al = _relative_al(times_yr, scan_angle_deg, params)
    B = params.mass_fraction_secondary
    absolute = absolute_al(times_yr, scan_angle_deg, params, astrometry)
    if response.mode == "photocentre":
        return photocentre_along_scan(d_al, B, params.beta_g) + absolute
    if response.mode == "blended_gaussian_peak":
        if response.sigma_mas is None:
            raise ValueError("sigma_mas is required for blended_gaussian_peak mode")
        if response.allow_multi_peak_continuation:
            scalar = np.ndim(d_al) == 0
            if response.sigma_secondary_mas is None:
                peak = _global_peak_coordinate(d_al, B, params.beta_g, response.sigma_mas)
            else:
                peak, _ = _unequal_width_peak_and_modes(
                    d_al, B, params.beta_g,
                    float(response.sigma_mas), float(response.sigma_secondary_mas),
                )
            peak = peak + absolute
            return float(peak[0]) if scalar else peak
        return blended_gaussian_peak(
            d_al, B, params.beta_g, response.sigma_mas,
            sigma_secondary_mas=response.sigma_secondary_mas,
        ) + absolute
    if response.mode == "penoyre_gaussian_peak":
        if response.sigma_mas is None:
            raise ValueError("sigma_mas (alpha) is required for penoyre_gaussian_peak mode")
        if response.sigma_ac_mas is None:
            raise ValueError("sigma_ac_mas (beta) is required for penoyre_gaussian_peak mode")
        rel = relative_astrometry_mas(times_yr, params)
        return penoyre_oriented_gaussian_peak(
            rel[:, 0],
            rel[:, 1],
            scan_angle_deg,
            B,
            params.beta_g,
            response.sigma_mas,
            response.sigma_ac_mas,
            allow_multi_peak_continuation=response.allow_multi_peak_continuation,
        ) + absolute
    raise ValueError(f"unknown Gaia response mode: {response.mode}")
