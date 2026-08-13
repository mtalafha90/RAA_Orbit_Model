from __future__ import annotations

import numpy as np


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


def blended_gaussian_peak(relative_al_mas, mass_fraction_secondary, beta_g, sigma_mas,
                          grid_size: int = 129, newton_steps: int = 12):
    """Prototype resolution-aware coordinate from two equal-width 1-D profiles.

    The measured coordinate is defined as the global maximum of

        I(x) = F1 exp[-(x-x1)^2/(2 sigma^2)]
             + F2 exp[-(x-x2)^2/(2 sigma^2)].

    A coarse vectorized grid identifies the correct maximum basin, followed by
    Newton refinement of dI/dx=0. This is a research surrogate, not Gaia's
    calibrated LSF. Component 1 is assumed to be brighter/equal (beta_g <= 0.5).
    """
    if sigma_mas <= 0:
        raise ValueError("sigma_mas must be > 0")
    if not (0.0 <= beta_g <= 0.5):
        raise ValueError("prototype assumes component 1 is brighter: 0 <= beta_g <= 0.5")
    if grid_size < 17:
        raise ValueError("grid_size must be >= 17")

    scalar = np.ndim(relative_al_mas) == 0
    d = np.atleast_1d(np.asarray(relative_al_mas, dtype=float))
    x1, x2 = component_al_positions(d, mass_fraction_secondary)
    f1, f2 = 1.0 - beta_g, beta_g
    sigma = float(sigma_mas)

    lo = np.minimum(x1, x2) - 5.0 * sigma
    hi = np.maximum(x1, x2) + 5.0 * sigma
    frac = np.linspace(0.0, 1.0, grid_size)
    grid = lo[:, None] + (hi - lo)[:, None] * frac[None, :]
    z1 = (grid - x1[:, None]) / sigma
    z2 = (grid - x2[:, None]) / sigma
    profile = f1 * np.exp(-0.5 * z1**2) + f2 * np.exp(-0.5 * z2**2)
    idx = np.argmax(profile, axis=1)
    x = grid[np.arange(len(d)), idx].copy()

    # Newton refinement. The basin selection above protects against converging
    # to the weaker local maximum when the pair becomes bimodal.
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

    return float(x[0]) if scalar else x
