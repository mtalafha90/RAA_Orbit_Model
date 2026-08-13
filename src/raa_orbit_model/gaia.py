from __future__ import annotations

import numpy as np
from scipy.optimize import minimize_scalar


def project_along_scan(delta_alpha_star_mas, delta_delta_mas, scan_angle_deg):
    """Gaia-convention tangent-plane projection: psi=0 points North, 90 deg East."""
    psi = np.deg2rad(np.asarray(scan_angle_deg, dtype=float))
    east = np.asarray(delta_alpha_star_mas, dtype=float)
    north = np.asarray(delta_delta_mas, dtype=float)
    return east * np.sin(psi) + north * np.cos(psi)


def photocentre_along_scan(relative_al_mas, mass_fraction_secondary, beta_g):
    """Unresolved photocentre relative to the barycentre."""
    return (beta_g - mass_fraction_secondary) * np.asarray(relative_al_mas, dtype=float)


def component_al_positions(relative_al_mas, mass_fraction_secondary):
    """Primary and secondary AL locations relative to the barycentre."""
    d = np.asarray(relative_al_mas, dtype=float)
    B = float(mass_fraction_secondary)
    return -B * d, (1.0 - B) * d


def blended_gaussian_peak(relative_al_mas, mass_fraction_secondary, beta_g, sigma_mas):
    """Prototype marginal-resolution response from two equal-width 1-D profiles.

    Returns the location of the global maximum of
        F1 exp[-(x-x1)^2/(2 sigma^2)] + F2 exp[-(x-x2)^2/(2 sigma^2)].

    This is a research surrogate, not Gaia's calibrated LSF. `sigma_mas` is
    deliberately explicit and must be calibrated/replaced for real Gaia data.
    Component 1 is assumed to be the brighter/equal source (beta_g <= 0.5).
    """
    if sigma_mas <= 0:
        raise ValueError("sigma_mas must be > 0")
    if not (0.0 <= beta_g <= 0.5):
        raise ValueError("prototype assumes component 1 is brighter: 0 <= beta_g <= 0.5")
    d = np.atleast_1d(np.asarray(relative_al_mas, dtype=float))
    x1, x2 = component_al_positions(d, mass_fraction_secondary)
    f1, f2 = 1.0 - beta_g, beta_g
    out = np.empty_like(d)
    for k, (a, b) in enumerate(zip(np.atleast_1d(x1), np.atleast_1d(x2))):
        lo = min(a, b) - 5.0 * sigma_mas
        hi = max(a, b) + 5.0 * sigma_mas

        def neg_profile(x):
            return -(f1 * np.exp(-0.5 * ((x-a)/sigma_mas)**2)
                     + f2 * np.exp(-0.5 * ((x-b)/sigma_mas)**2))

        grid = np.linspace(lo, hi, 257)
        vals = np.array([neg_profile(x) for x in grid])
        j = int(np.argmin(vals))
        left = grid[max(0, j-1)]
        right = grid[min(len(grid)-1, j+1)]
        res = minimize_scalar(neg_profile, bounds=(left, right), method="bounded")
        out[k] = res.x
    return out if np.ndim(relative_al_mas) else float(out[0])
