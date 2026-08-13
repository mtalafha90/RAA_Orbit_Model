from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.optimize import brentq


class MultiPeakProfileError(ValueError):
    """Raised when the single-coordinate surrogate is asked to model a multi-peak profile."""

    def __init__(self, n_multi_peak: int):
        self.n_multi_peak = int(n_multi_peak)
        super().__init__(
            f"blended Gaussian profile is multi-peaked for {self.n_multi_peak} transit(s); "
            "a unique single-coordinate response is undefined"
        )


@dataclass(frozen=True)
class BlendedGaussianResponse:
    """Validity-aware response of the equal-width two-Gaussian research surrogate.

    ``al_mas`` is finite only when the blended 1-D profile has one mode.  A
    multi-peak transit is deliberately represented by NaN rather than by an
    arbitrary choice of one of the modes.
    """

    al_mas: float | np.ndarray
    n_peaks: int | np.ndarray
    projected_separation_sigma: float | np.ndarray
    critical_separation_sigma: float

    @property
    def single_peak_mask(self):
        return np.asarray(self.n_peaks) == 1

    @property
    def multi_peak_mask(self):
        return np.asarray(self.n_peaks) > 1


def project_along_scan(delta_alpha_star_mas, delta_delta_mas, scan_angle_deg):
    """Project tangent-plane East/North offsets onto Gaia-like along-scan axes.

    Convention: scan angle psi=0 points North and psi=90 deg points East.
    """
    psi = np.deg2rad(np.asarray(scan_angle_deg, dtype=float))
    east = np.asarray(delta_alpha_star_mas, dtype=float)
    north = np.asarray(delta_delta_mas, dtype=float)
    return east * np.sin(psi) + north * np.cos(psi)


def photocentre_along_scan(relative_al_mas, mass_fraction_secondary, beta_g):
    """Unresolved photocentre relative to the barycentre."""
    return (beta_g - mass_fraction_secondary) * np.asarray(relative_al_mas, dtype=float)


def component_al_positions(relative_al_mas, mass_fraction_secondary):
    """Primary and secondary along-scan positions relative to barycentre."""
    d = np.asarray(relative_al_mas, dtype=float)
    B = float(mass_fraction_secondary)
    return -B * d, (1.0 - B) * d


def critical_blended_separation_sigma(beta_g: float) -> float:
    """Critical separation (in sigma) at which the two-Gaussian mixture becomes bimodal.

    The mixture has equal component widths and flux fractions ``1-beta_g`` and
    ``beta_g`` with ``0 <= beta_g <= 0.5``.  At the mode-splitting boundary,
    the stationary point is simultaneously a root of the first and second
    derivatives.  If the normalized distances from that point to the two
    component centres are ``a`` and ``b``, those conditions imply ``a*b=1``.
    Writing ``a=exp(u)``, ``b=exp(-u)`` gives

        sinh(2u) - 2u = log((1-beta_g)/beta_g),
        d_crit/sigma = 2 cosh(u).

    This is an exact criterion for the present equal-width Gaussian surrogate,
    not a Gaia instrumental resolution threshold.
    """
    beta = float(beta_g)
    if not (0.0 <= beta <= 0.5):
        raise ValueError("prototype assumes component 1 is brighter: 0 <= beta_g <= 0.5")
    if beta == 0.0:
        return math.inf
    if beta == 0.5:
        return 2.0

    log_ratio = math.log((1.0 - beta) / beta)

    def equation(u: float) -> float:
        return math.sinh(2.0 * u) - 2.0 * u - log_ratio

    upper = 1.0
    while equation(upper) < 0.0:
        upper *= 2.0
    u = brentq(equation, 0.0, upper, xtol=1e-13, rtol=1e-13)
    return 2.0 * math.cosh(u)


def blended_gaussian_peak_count(relative_al_mas, beta_g: float, sigma_mas: float):
    """Return one or two modes for each projected separation in the surrogate profile."""
    if sigma_mas <= 0:
        raise ValueError("sigma_mas must be > 0")
    critical = critical_blended_separation_sigma(beta_g)
    scalar = np.ndim(relative_al_mas) == 0
    separation_sigma = np.abs(np.atleast_1d(np.asarray(relative_al_mas, dtype=float))) / float(sigma_mas)
    # At the exact critical separation the two modes have not yet split; the
    # profile has a single degenerate stationary maximum.  Only strictly larger
    # separations are classified as multi-peak.
    n_peaks = np.where(separation_sigma > critical, 2, 1).astype(int)
    return int(n_peaks[0]) if scalar else n_peaks


def _global_peak_coordinate(
    relative_al_mas,
    mass_fraction_secondary,
    beta_g,
    sigma_mas,
    *,
    grid_size: int = 129,
    newton_steps: int = 12,
):
    """Numerically locate the global maximum; caller must ensure it is unique."""
    d = np.atleast_1d(np.asarray(relative_al_mas, dtype=float))
    x1, x2 = component_al_positions(d, mass_fraction_secondary)
    f1, f2 = 1.0 - beta_g, beta_g
    sigma = float(sigma_mas)

    # Exact symmetry limits avoid unnecessary numerical noise.
    if beta_g == 0.0:
        return np.asarray(x1, dtype=float)
    if beta_g == 0.5:
        return 0.5 * (np.asarray(x1, dtype=float) + np.asarray(x2, dtype=float))

    lo = np.minimum(x1, x2) - 5.0 * sigma
    hi = np.maximum(x1, x2) + 5.0 * sigma
    frac = np.linspace(0.0, 1.0, grid_size)
    grid = lo[:, None] + (hi - lo)[:, None] * frac[None, :]
    z1 = (grid - x1[:, None]) / sigma
    z2 = (grid - x2[:, None]) / sigma
    profile = f1 * np.exp(-0.5 * z1**2) + f2 * np.exp(-0.5 * z2**2)
    idx = np.argmax(profile, axis=1)
    x = grid[np.arange(len(d)), idx].copy()

    for _ in range(newton_steps):
        dx1 = x - x1
        dx2 = x - x2
        e1 = np.exp(-0.5 * (dx1 / sigma)**2)
        e2 = np.exp(-0.5 * (dx2 / sigma)**2)
        dI = -(f1 * dx1 * e1 + f2 * dx2 * e2) / sigma**2
        d2I = (
            f1 * ((dx1**2 / sigma**4) - 1.0 / sigma**2) * e1
            + f2 * ((dx2**2 / sigma**4) - 1.0 / sigma**2) * e2
        )
        safe = np.abs(d2I) > 1e-18
        step = np.zeros_like(x)
        step[safe] = dI[safe] / d2I[safe]
        x_new = np.clip(x - step, lo, hi)
        if np.max(np.abs(x_new - x)) < 1e-12:
            x = x_new
            break
        x = x_new
    return x


def blended_gaussian_response(
    relative_al_mas,
    mass_fraction_secondary,
    beta_g,
    sigma_mas,
    *,
    grid_size: int = 129,
    newton_steps: int = 12,
) -> BlendedGaussianResponse:
    """Return a validity-aware single/multi-peak response for each AL separation.

    Multi-peak profiles are flagged and assigned ``NaN`` in ``al_mas``.  The
    function therefore never chooses an arbitrary peak when the single-centroid
    surrogate is no longer well defined.
    """
    if sigma_mas <= 0:
        raise ValueError("sigma_mas must be > 0")
    if grid_size < 17:
        raise ValueError("grid_size must be >= 17")
    critical = critical_blended_separation_sigma(beta_g)

    scalar = np.ndim(relative_al_mas) == 0
    d = np.atleast_1d(np.asarray(relative_al_mas, dtype=float))
    separation_sigma = np.abs(d) / float(sigma_mas)
    n_peaks = np.where(separation_sigma > critical, 2, 1).astype(int)
    valid = n_peaks == 1
    al = np.full(len(d), np.nan, dtype=float)
    if np.any(valid):
        al[valid] = _global_peak_coordinate(
            d[valid],
            mass_fraction_secondary,
            beta_g,
            sigma_mas,
            grid_size=grid_size,
            newton_steps=newton_steps,
        )

    if scalar:
        return BlendedGaussianResponse(
            float(al[0]),
            int(n_peaks[0]),
            float(separation_sigma[0]),
            float(critical),
        )
    return BlendedGaussianResponse(al, n_peaks, separation_sigma, float(critical))


def blended_gaussian_peak(
    relative_al_mas,
    mass_fraction_secondary,
    beta_g,
    sigma_mas,
    grid_size: int = 129,
    newton_steps: int = 12,
):
    """Return the unique blended-profile maximum, rejecting multi-peak inputs.

    This function now represents only the single-peak validity domain of the
    surrogate.  Use :func:`blended_gaussian_response` when per-transit regime
    flags are required.
    """
    response = blended_gaussian_response(
        relative_al_mas,
        mass_fraction_secondary,
        beta_g,
        sigma_mas,
        grid_size=grid_size,
        newton_steps=newton_steps,
    )
    multi = np.asarray(response.n_peaks) > 1
    if np.any(multi):
        raise MultiPeakProfileError(int(np.count_nonzero(multi)))
    return response.al_mas
