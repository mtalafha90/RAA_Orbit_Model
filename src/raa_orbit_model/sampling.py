"""Posterior sampling for the joint orbit model.

`docs/methodology.md` section 10 states that the target is a posterior and that
the deterministic least-squares solution is only a validation backend. Until
now the project produced point estimates alone, so it could report a bias but
not an uncertainty, and could not show which parameters the measurement-model
choice actually degrades through their covariance.

This module implements the affine-invariant ensemble sampler of Goodman & Weare
(2010) directly, in about sixty lines of NumPy, so the package acquires no new
dependency and the algorithm stays inspectable. It is a validation-grade
sampler: adequate for mapping degeneracies in a controlled synthetic study,
and not a substitute for a production sampler on real data.

Priors are uniform over the same bounds the least-squares fitter uses, so the
two backends describe the same parameter space. Anything ``fit_joint`` can
free, this can sample: the orbit, the absolute astrometric parameters and the
instrument response width. It additionally samples per-channel jitter terms,
which the least-squares backend cannot represent because they change the
variance rather than the model, and so need the normalisation that only the
full likelihood carries.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from .fit import RESPONSE_PARAMETER_NAMES, _bounds_for, joint_loglike
from .kepler import BinaryParams
from .model import GaiaResponseConfig
from .synthetic import JointData


@dataclass(frozen=True)
class PosteriorSamples:
    """Chain from the ensemble sampler, shaped (steps, walkers, parameters)."""

    names: tuple[str, ...]
    chain: np.ndarray
    log_prob: np.ndarray
    acceptance_fraction: float
    n_burn: int

    def flat(self, thin: int = 1) -> np.ndarray:
        """Post-burn-in samples, flattened across walkers."""
        if thin < 1:
            raise ValueError("thin must be >= 1")
        kept = self.chain[self.n_burn::thin]
        return kept.reshape(-1, kept.shape[-1])

    def summary(self) -> dict[str, dict[str, float]]:
        """Median and 16th/84th percentiles for each free parameter."""
        flat = self.flat()
        out: dict[str, dict[str, float]] = {}
        for index, name in enumerate(self.names):
            column = flat[:, index]
            q16, median, q84 = np.percentile(column, [16.0, 50.0, 84.0])
            out[name] = {
                "median": float(median),
                "q16": float(q16),
                "q84": float(q84),
                "minus": float(median - q16),
                "plus": float(q84 - median),
                "std": float(np.std(column)),
            }
        return out


#: Nuisance terms that inflate a channel's variance rather than shift its model.
#: ``joint_loglike`` carries the ``log(variance)`` penalty that bounds them, which
#: is the reason the sampler uses the likelihood rather than the chi-square.
JITTER_PARAMETER_NAMES = ("gaia_jitter_mas", "rv_jitter_kms")


def _log_posterior(theta, names, base, data, response, lower, upper) -> float:
    if np.any(theta < lower) or np.any(theta > upper):
        return -np.inf

    binary_updates = {}
    sigma_response = None
    jitter = {"gaia_jitter_mas": 0.0, "rv_jitter_kms": 0.0}
    for name, value in zip(names, theta):
        if name in JITTER_PARAMETER_NAMES:
            jitter[name] = float(value)
        elif name in RESPONSE_PARAMETER_NAMES:
            sigma_response = float(value)
        else:
            binary_updates[name] = float(value)

    params = replace(base, **binary_updates)
    if sigma_response is not None:
        response = replace(response, sigma_mas=sigma_response)
    try:
        params.validate()
        return joint_loglike(
            params, data, response,
            gaia_jitter_mas=jitter["gaia_jitter_mas"],
            rv_jitter_kms=jitter["rv_jitter_kms"],
        )
    except (ValueError, RuntimeError):
        return -np.inf


def sample_posterior(
    data: JointData,
    initial: BinaryParams,
    gaia_response: GaiaResponseConfig,
    *,
    free_names,
    n_walkers: int = 48,
    n_steps: int = 2000,
    n_burn: int = 500,
    stretch_a: float = 2.0,
    ball_fraction: float = 1e-4,
    seed: int = 0,
) -> PosteriorSamples:
    """Sample the joint posterior with an affine-invariant ensemble sampler."""
    names = tuple(free_names)
    if not names:
        raise ValueError("free_names is empty")
    n_dim = len(names)
    if n_walkers < 2 * n_dim:
        raise ValueError(f"n_walkers must be at least 2*n_dim = {2 * n_dim}")
    if n_walkers % 2:
        raise ValueError("n_walkers must be even")
    if not (0 <= n_burn < n_steps):
        raise ValueError("require 0 <= n_burn < n_steps")
    if stretch_a <= 1.0:
        raise ValueError("stretch_a must be > 1")

    if any(n in RESPONSE_PARAMETER_NAMES for n in names) and gaia_response.mode == "photocentre":
        raise ValueError("the photocentre response has no width to sample")

    def bounds_of(name):
        if name == "gaia_jitter_mas":
            return 0.0, 10.0 * float(np.median(data.gaia_al.sigma_mas))
        if name == "rv_jitter_kms":
            return 0.0, 10.0 * float(np.median(data.sb2_rv.sigma_kms))
        return _bounds_for(name, initial, gaia_response.sigma_mas)

    def centre_of(name):
        if name in JITTER_PARAMETER_NAMES:
            # Start just off zero so the walker ball is not pinned to the bound.
            return 0.01 * bounds_of(name)[1]
        if name in RESPONSE_PARAMETER_NAMES:
            return float(gaia_response.sigma_mas)
        return float(getattr(initial, name))

    bounds = np.array([bounds_of(n) for n in names], dtype=float)
    lower, upper = bounds[:, 0], bounds[:, 1]
    centre = np.array([centre_of(n) for n in names], dtype=float)

    rng = np.random.default_rng(seed)
    scale = np.where(np.abs(centre) > 0, np.abs(centre), 1.0) * ball_fraction
    position = centre + scale * rng.standard_normal((n_walkers, n_dim))
    position = np.clip(position, lower, upper)

    def log_prob_of(x):
        return _log_posterior(x, names, initial, data, gaia_response, lower, upper)

    log_prob = np.array([log_prob_of(x) for x in position])
    if not np.any(np.isfinite(log_prob)):
        raise ValueError("no starting walker has finite posterior probability")

    chain = np.empty((n_steps, n_walkers, n_dim), dtype=float)
    log_prob_chain = np.empty((n_steps, n_walkers), dtype=float)
    accepted = 0
    half = n_walkers // 2

    for step in range(n_steps):
        for first in (True, False):
            active = np.arange(0, half) if first else np.arange(half, n_walkers)
            complement = np.arange(half, n_walkers) if first else np.arange(0, half)
            partners = position[rng.choice(complement, size=len(active))]

            # Goodman & Weare stretch move: z ~ g(z) on [1/a, a].
            u = rng.random(len(active))
            z = ((stretch_a - 1.0) * u + 1.0) ** 2 / stretch_a
            proposal = partners + z[:, None] * (position[active] - partners)

            for k, walker in enumerate(active):
                candidate_log_prob = log_prob_of(proposal[k])
                log_accept = (n_dim - 1) * np.log(z[k]) + candidate_log_prob - log_prob[walker]
                if np.log(rng.random()) < log_accept:
                    position[walker] = proposal[k]
                    log_prob[walker] = candidate_log_prob
                    accepted += 1

        chain[step] = position
        log_prob_chain[step] = log_prob

    return PosteriorSamples(
        names=names,
        chain=chain,
        log_prob=log_prob_chain,
        acceptance_fraction=accepted / float(n_steps * n_walkers),
        n_burn=n_burn,
    )
