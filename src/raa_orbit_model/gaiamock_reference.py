"""Reference implementation of the published gaiamock along-scan binary response.

El-Badry et al. (2024), *Open Journal of Astrophysics* **7**, DOI
10.33232/001c.125461, released ``gaiamock``, whose ``al_bias_binary`` predicts
the Gaia along-scan location of a binary "assuming that the 1D centroid is at
the peak of the combined AL flux profile, following the model from
Lindegren+2022".

This module reproduces that function as published so that the surrogate in
:mod:`raa_orbit_model.gaia` can be scored against it. It is a *reference for
comparison only* and is deliberately not used anywhere in the inference path.

The published routine solves

    x = f xi / (f + exp(xi^2 / 2 - xi x))

by fixed-point iteration, where ``xi`` is the projected along-scan separation
in units of the effective angular resolution ``u`` and ``f = F2 / F1``. That is
the stationary-point condition for the maximum of

    I(x) = exp(-x^2 / 2) + f exp(-(x - xi)^2 / 2),

i.e. an **equal-width two-Gaussian blend** whose measured coordinate is the
profile peak. It is therefore the same physical model as this project's
surrogate, with ``u`` playing the role of the surrogate width. The published
version differs in three respects:

1. it linearises to the photocentre below ``0.1 u`` for numerical stability;
2. it assigns the primary's position above ``(3 - f) u`` rather than solving;
3. its fixed-point iteration starts at the primary and so tracks that mode,
   whereas this project locates the global maximum and refuses to return a
   single coordinate once the profile is genuinely bimodal.
"""

from __future__ import annotations

import numpy as np


def _solve_for_x(ff: float, xi: float, tol: float = 1e-6, niter_max: int = 100) -> float:
    """Fixed-point solve for the blended-profile peak, as published."""
    x = 0.0
    for _ in range(niter_max):
        thisx = ff * xi / (ff + np.exp(0.5 * xi**2 - xi * x))
        if abs(thisx - x) < tol:
            break
        x = thisx
    return x


def al_bias_binary(delta_eta, q: float, f: float, u: float = 90.0):
    """Along-scan offset from the barycentre, following gaiamock.

    ``delta_eta`` is the projected along-scan separation in mas, ``q = M2 / M1``
    is the mass ratio, ``f = F2 / F1`` is the light ratio, and ``u`` is the
    effective angular resolution in mas.

    Note that the published docstring labels ``q`` a flux ratio, but it enters
    as ``q / (1 + q)``, which is the secondary mass fraction; it is the mass
    ratio. Scalar in, scalar out, matching the published signature.
    """
    ratio = np.abs(delta_eta / u)
    if ratio <= 0.1:
        return (f / (1 + f) - q / (1 + q)) * delta_eta
    if ratio <= 3 - f:
        B = _solve_for_x(ff=f, xi=delta_eta / u)
        return u * B - q / (1 + q) * delta_eta
    return -q / (1 + q) * delta_eta


def al_bias_binary_array(delta_eta, q: float, f: float, u: float = 90.0) -> np.ndarray:
    """Vectorised convenience wrapper over :func:`al_bias_binary`."""
    return np.array(
        [al_bias_binary(float(d), q, f, u) for d in np.atleast_1d(delta_eta)],
        dtype=float,
    )


def light_ratio_from_beta(beta_g: float) -> float:
    """Convert this project's flux fraction beta = F2/(F1+F2) to gaiamock's f = F2/F1."""
    beta = float(beta_g)
    if not (0.0 <= beta < 1.0):
        raise ValueError("beta_g must satisfy 0 <= beta_g < 1")
    return beta / (1.0 - beta)


def mass_ratio_from_fraction(mass_fraction_secondary: float) -> float:
    """Convert B = M2/(M1+M2) to gaiamock's q = M2/M1."""
    B = float(mass_fraction_secondary)
    if not (0.0 <= B < 1.0):
        raise ValueError("mass fraction must satisfy 0 <= B < 1")
    return B / (1.0 - B)
