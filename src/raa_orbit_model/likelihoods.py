from __future__ import annotations

import numpy as np

LOG2PI = np.log(2.0 * np.pi)


def gaussian_1d_loglike(observed, model, sigma, jitter=0.0):
    y = np.asarray(observed, dtype=float)
    m = np.asarray(model, dtype=float)
    s = np.asarray(sigma, dtype=float)
    if np.any(s <= 0) or jitter < 0:
        raise ValueError("sigma must be > 0 and jitter >= 0")
    var = s * s + jitter * jitter
    r = y - m
    return float(-0.5 * np.sum(r * r / var + np.log(var) + LOG2PI))


def gaussian_2d_loglike(observed_xy, model_xy, covariance):
    """Vectorised over epochs so this is usable inside a sampler."""
    y = np.asarray(observed_xy, dtype=float)
    m = np.asarray(model_xy, dtype=float)
    C = np.asarray(covariance, dtype=float)
    if y.shape != m.shape or y.ndim != 2 or y.shape[1] != 2:
        raise ValueError("observed_xy and model_xy must have shape (N,2)")
    if C.shape != (len(y), 2, 2):
        raise ValueError("covariance must have shape (N,2,2)")
    if len(y) == 0:
        return 0.0
    sign, logdet = np.linalg.slogdet(C)
    if np.any(sign <= 0):
        raise ValueError("each covariance matrix must be positive definite")
    r = y - m
    quadratic = np.einsum("ni,ni->n", r, np.linalg.solve(C, r[..., None])[..., 0])
    return float(-0.5 * np.sum(quadratic + logdet + 2.0 * LOG2PI))
