"""Idealised orientation-dependent blended-source response after Penoyre (2026).

This module is deliberately a research surrogate, not a calibrated Gaia PLSF.
Penoyre (2026, RASTI 5, rzaf062; corrected by rzag016) shows that two sources
observed with an elongated Gaussian PSF reduce to the same one-dimensional
blend problem as the circular case when lengths are normalised by the
orientation-dependent effective width

    gamma(phi) = alpha / sqrt(1 + [(alpha/beta)^2 - 1] sin(phi)^2).

Here ``alpha_mas`` is the width along the scan axis and ``beta_mas`` is the
width across the scan axis of that idealised Gaussian.  ``phi`` is the angle
between the binary separation vector and the scan axis.  The observed maximum
lies on the line connecting the sources, so its along-scan coordinate is the
maximum along that line multiplied by cos(phi).

No low-order series expansion from Penoyre is used here.  The exact effective
width is fed into the repository's numerical equal-width peak solver.  This
also avoids dependence on the sign errors corrected in rzag016.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .gaia import (
    MultiPeakProfileError,
    _global_peak_coordinate,
    blended_gaussian_peak_count,
    blended_gaussian_response,
    critical_blended_separation_sigma,
    project_along_scan,
)


@dataclass(frozen=True)
class PenoyreGaussianResponse:
    """Validity-aware AL response of the idealised elongated-Gaussian model."""

    al_mas: float | np.ndarray
    n_peaks: int | np.ndarray
    separation_sigma: float | np.ndarray
    effective_width_mas: float | np.ndarray
    pair_scan_angle_rad: float | np.ndarray
    projected_al_separation_mas: float | np.ndarray
    critical_separation_sigma: float

    @property
    def single_peak_mask(self):
        return np.asarray(self.n_peaks) == 1

    @property
    def multi_peak_mask(self):
        return np.asarray(self.n_peaks) > 1


def penoyre_effective_width(alpha_mas: float, beta_mas: float, phi_rad):
    """Return Penoyre's orientation-dependent Gaussian width ``gamma(phi)``.

    ``alpha_mas`` is the along-axis/narrow width and ``beta_mas`` the
    across-axis/long width in the idealised model.  ``beta_mas=np.inf`` is
    allowed because Penoyre notes that this limit reproduces the effectively
    unconstrained across-scan assumption used by the 1-D Lindegren-style
    treatment.

    The widths are research parameters.  They are not Gaia calibration values.
    """
    alpha = float(alpha_mas)
    beta = float(beta_mas)
    if not math.isfinite(alpha) or alpha <= 0.0:
        raise ValueError("alpha_mas must be finite and > 0")
    if math.isnan(beta) or beta <= 0.0:
        raise ValueError("beta_mas must be > 0 or +inf")
    if beta < alpha:
        raise ValueError("the surrogate convention requires beta_mas >= alpha_mas")

    phi = np.asarray(phi_rad, dtype=float)
    inv_beta2 = 0.0 if math.isinf(beta) else 1.0 / beta**2
    inv_gamma2 = np.cos(phi) ** 2 / alpha**2 + np.sin(phi) ** 2 * inv_beta2
    with np.errstate(divide="ignore", invalid="ignore"):
        gamma = 1.0 / np.sqrt(inv_gamma2)
    return float(gamma) if np.ndim(phi_rad) == 0 else gamma


def _oriented_raw_response(
    relative_east_mas,
    relative_north_mas,
    scan_angle_deg,
    mass_fraction_secondary: float,
    beta_g: float,
    alpha_mas: float,
    beta_mas: float,
):
    """Return the raw global maximum and validity metadata before NaN masking."""
    east, north, scan = np.broadcast_arrays(
        np.asarray(relative_east_mas, dtype=float),
        np.asarray(relative_north_mas, dtype=float),
        np.asarray(scan_angle_deg, dtype=float),
    )
    scalar = east.ndim == 0
    shape = east.shape
    east_f = np.atleast_1d(east).ravel()
    north_f = np.atleast_1d(north).ravel()
    scan_f = np.atleast_1d(scan).ravel()

    psi = np.deg2rad(scan_f)
    r = np.hypot(east_f, north_f)
    theta = np.arctan2(east_f, north_f)  # position angle: North through East
    phi = theta - psi
    phi = np.where(r == 0.0, 0.0, phi)
    d_al = np.asarray(project_along_scan(east_f, north_f, scan_f), dtype=float)
    critical = critical_blended_separation_sigma(beta_g)

    # Penoyre notes that beta -> infinity is the limit implicit in the 1-D
    # Lindegren-style AL calculation.  Use the existing implementation exactly
    # in that limit so M2 has a regression-tested reduction to M1.
    if math.isinf(float(beta_mas)):
        baseline = blended_gaussian_response(
            d_al,
            mass_fraction_secondary,
            beta_g,
            alpha_mas,
        )
        raw_peak = np.asarray(
            _global_peak_coordinate(
                d_al,
                mass_fraction_secondary,
                beta_g,
                alpha_mas,
            ),
            dtype=float,
        )
        gamma = np.asarray(penoyre_effective_width(alpha_mas, beta_mas, phi), dtype=float)
        rho = np.abs(d_al) / float(alpha_mas)
        n_peaks = np.asarray(baseline.n_peaks, dtype=int)
    else:
        gamma = np.asarray(penoyre_effective_width(alpha_mas, beta_mas, phi), dtype=float)
        rho = r / gamma
        n_peaks = np.where(rho > critical, 2, 1).astype(int)

        # The extrema lie on the source-separation line.  Solve the same
        # dimensionless equal-width problem at sigma=1, scale back by gamma,
        # then project that peak position onto the Gaia-like AL axis.
        peak_pair_sigma = np.asarray(
            _global_peak_coordinate(
                rho,
                mass_fraction_secondary,
                beta_g,
                1.0,
            ),
            dtype=float,
        )
        raw_peak = gamma * peak_pair_sigma * np.cos(phi)

    def shaped(values):
        arr = np.asarray(values)
        if scalar:
            return arr.ravel()[0]
        return arr.reshape(shape)

    return (
        shaped(raw_peak),
        shaped(n_peaks),
        shaped(rho),
        shaped(gamma),
        shaped(phi),
        shaped(d_al),
        float(critical),
    )


def penoyre_oriented_gaussian_response(
    relative_east_mas,
    relative_north_mas,
    scan_angle_deg,
    mass_fraction_secondary: float,
    beta_g: float,
    alpha_mas: float,
    beta_mas: float,
) -> PenoyreGaussianResponse:
    """Return the single/multi-peak response for an elongated Gaussian PSF.

    For finite ``beta_mas``, mode splitting is governed by ``r/gamma(phi)`` and
    the same light-ratio-dependent critical separation as the equal-width 1-D
    blend.  Multi-peak epochs receive NaN AL coordinates because a unique
    single-coordinate measurement is no longer defined by this surrogate.
    """
    raw, n_peaks, rho, gamma, phi, d_al, critical = _oriented_raw_response(
        relative_east_mas,
        relative_north_mas,
        scan_angle_deg,
        mass_fraction_secondary,
        beta_g,
        alpha_mas,
        beta_mas,
    )
    multi = np.asarray(n_peaks) > 1
    al = np.where(multi, np.nan, raw)
    if np.ndim(raw) == 0:
        al = float(al)
    return PenoyreGaussianResponse(
        al_mas=al,
        n_peaks=n_peaks,
        separation_sigma=rho,
        effective_width_mas=gamma,
        pair_scan_angle_rad=phi,
        projected_al_separation_mas=d_al,
        critical_separation_sigma=critical,
    )


def penoyre_oriented_gaussian_peak(
    relative_east_mas,
    relative_north_mas,
    scan_angle_deg,
    mass_fraction_secondary: float,
    beta_g: float,
    alpha_mas: float,
    beta_mas: float,
    *,
    allow_multi_peak_continuation: bool = False,
):
    """Return the oriented peak, optionally allowing optimizer continuation.

    Continuation through a multi-peak state is a numerical device only.  Final
    scientific solutions must be checked with
    :func:`penoyre_oriented_gaussian_response` and accepted only where the
    retained epochs are single-peaked.
    """
    raw, n_peaks, *_ = _oriented_raw_response(
        relative_east_mas,
        relative_north_mas,
        scan_angle_deg,
        mass_fraction_secondary,
        beta_g,
        alpha_mas,
        beta_mas,
    )
    multi = np.asarray(n_peaks) > 1
    if np.any(multi) and not allow_multi_peak_continuation:
        raise MultiPeakProfileError(int(np.count_nonzero(multi)))
    return raw
