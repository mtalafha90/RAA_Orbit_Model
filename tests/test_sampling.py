"""Posterior sampling, and the likelihood it is built on.

methodology.md section 10 asks for a posterior; the project previously produced
point estimates only. These tests check the likelihood is consistent with the
existing chi-square path, and that the sampler recovers known parameters with
credible intervals that cover the truth.
"""

import numpy as np
import pytest
from dataclasses import replace

from raa_orbit_model.experiments import single_peak_schedule_for_response
from raa_orbit_model.fit import fit_joint, joint_loglike, joint_residuals
from raa_orbit_model.kepler import BinaryParams
from raa_orbit_model.likelihoods import gaussian_1d_loglike, gaussian_2d_loglike
from raa_orbit_model.model import GaiaResponseConfig
from raa_orbit_model.sampling import sample_posterior
from raa_orbit_model.scanning import schedule_from_arrays
from raa_orbit_model.synthetic import simulate_joint_data

SIGMA = 50.0


def _case(n_gaia=40):
    rng = np.random.default_rng(3)
    schedule = schedule_from_arrays(
        np.sort(rng.uniform(0.0, 5.0, n_gaia)), rng.uniform(0.0, 360.0, n_gaia),
        ra_deg=120.0, dec_deg=30.0,
    )
    base = BinaryParams(
        period_yr=2.0, t_peri_yr=0.15, eccentricity=0.25, inclination_deg=72.0,
        omega_deg=55.0, node_deg=120.0, m1_msun=1.25, m2_msun=0.85,
        parallax_mas=20.0, gamma_kms=7.0, beta_g=0.25,
    )
    truth = replace(base, parallax_mas=SIGMA / base.a_rel_au)
    injected = GaiaResponseConfig("blended_gaussian_peak", SIGMA)
    fit_response = GaiaResponseConfig(
        "blended_gaussian_peak", SIGMA, allow_multi_peak_continuation=True
    )
    selection = single_peak_schedule_for_response(truth, schedule, injected)
    data = simulate_joint_data(truth, injected, selection.schedule, seed=0,
                               baseline_yr=schedule.mission_span_yr,
                               n_ast=12, n_rv=18)
    return truth, data, fit_response


def test_two_dimensional_likelihood_is_vectorised_and_matches_a_manual_sum():
    rng = np.random.default_rng(0)
    y = rng.normal(size=(7, 2))
    m = rng.normal(size=(7, 2))
    cov = np.repeat(np.eye(2)[None] * 0.25, 7, axis=0)
    manual = 0.0
    for r, c in zip(y - m, cov):
        sign, logdet = np.linalg.slogdet(c)
        manual += -0.5 * (r @ np.linalg.solve(c, r) + logdet + 2.0 * np.log(2.0 * np.pi))
    assert gaussian_2d_loglike(y, m, cov) == pytest.approx(manual)


def test_two_dimensional_likelihood_rejects_a_singular_covariance():
    y = np.zeros((2, 2))
    cov = np.repeat(np.zeros((2, 2))[None], 2, axis=0)
    with pytest.raises(ValueError, match="positive definite"):
        gaussian_2d_loglike(y, y, cov)


def test_jitter_lowers_the_penalty_on_a_large_residual():
    observed, model, sigma = np.array([5.0]), np.array([0.0]), np.array([1.0])
    assert gaussian_1d_loglike(observed, model, sigma, jitter=4.0) > gaussian_1d_loglike(
        observed, model, sigma
    )


def test_loglike_and_chi2_differ_by_a_constant():
    """The likelihood must describe the same surface the fitter minimises."""
    truth, data, response = _case()
    offsets = []
    for scale in (1.0, 1.01, 0.98):
        params = replace(truth, period_yr=truth.period_yr * scale)
        residuals = joint_residuals(params, data, response)
        offsets.append(-2.0 * joint_loglike(params, data, response) - float(residuals @ residuals))
    assert max(offsets) - min(offsets) < 1e-6


def test_sampler_validates_its_configuration():
    truth, data, response = _case()
    with pytest.raises(ValueError, match="free_names is empty"):
        sample_posterior(data, truth, response, free_names=())
    with pytest.raises(ValueError, match="at least"):
        sample_posterior(data, truth, response, free_names=("m1_msun", "m2_msun"), n_walkers=2)
    with pytest.raises(ValueError, match="even"):
        sample_posterior(data, truth, response, free_names=("m1_msun",), n_walkers=7)
    with pytest.raises(ValueError, match="n_burn"):
        sample_posterior(data, truth, response, free_names=("m1_msun",),
                         n_walkers=6, n_steps=10, n_burn=10)
    with pytest.raises(ValueError, match="stretch_a"):
        sample_posterior(data, truth, response, free_names=("m1_msun",),
                         n_walkers=6, n_steps=10, n_burn=2, stretch_a=1.0)


def test_posterior_recovers_truth_with_covering_credible_intervals():
    truth, data, response = _case()
    free = ("m1_msun", "m2_msun")
    start = replace(truth, m1_msun=1.28, m2_msun=0.83)
    best = fit_joint(data, start, response, free_names=free)
    assert best.success

    posterior = sample_posterior(data, best.params, response, free_names=free,
                                 n_walkers=12, n_steps=400, n_burn=150, seed=1)
    assert 0.1 < posterior.acceptance_fraction < 0.9
    assert posterior.chain.shape == (400, 12, 2)

    summary = posterior.summary()
    for name in free:
        stats = summary[name]
        true_value = getattr(truth, name)
        assert stats["median"] - 3.0 * stats["minus"] <= true_value
        assert true_value <= stats["median"] + 3.0 * stats["plus"]


def test_flat_respects_burn_in_and_thinning():
    truth, data, response = _case()
    posterior = sample_posterior(data, truth, response, free_names=("m1_msun",),
                                 n_walkers=6, n_steps=60, n_burn=20, seed=0)
    assert posterior.flat().shape == (40 * 6, 1)
    assert posterior.flat(thin=4).shape == (10 * 6, 1)
    with pytest.raises(ValueError, match="thin"):
        posterior.flat(thin=0)


def test_response_width_can_be_sampled_not_only_fitted():
    """fit_joint can free the width, so the sampler must be able to as well."""
    truth, data, response = _case()
    posterior = sample_posterior(data, truth, response,
                                 free_names=("m1_msun", "sigma_response_mas"),
                                 n_walkers=8, n_steps=200, n_burn=80, seed=0)
    width = posterior.summary()["sigma_response_mas"]
    assert width["median"] == pytest.approx(SIGMA, rel=0.1)


def test_sampling_a_width_is_rejected_for_the_photocentre_model():
    truth, data, _ = _case()
    with pytest.raises(ValueError, match="no width to sample"):
        sample_posterior(data, truth, GaiaResponseConfig("photocentre"),
                         free_names=("m1_msun", "sigma_response_mas"),
                         n_walkers=8, n_steps=20, n_burn=5)


def test_jitter_terms_are_sampleable_and_stay_small_without_excess_noise():
    """The log-variance penalty in the likelihood is what bounds these."""
    truth, data, response = _case()
    posterior = sample_posterior(data, truth, response,
                                 free_names=("m1_msun", "gaia_jitter_mas", "rv_jitter_kms"),
                                 n_walkers=8, n_steps=200, n_burn=80, seed=0)
    summary = posterior.summary()
    # The synthetic data carry no excess noise, so both jitters should be
    # driven well below the quoted per-epoch uncertainties rather than running away.
    assert 0.0 <= summary["gaia_jitter_mas"]["median"] < float(np.median(data.gaia_al.sigma_mas))
    assert 0.0 <= summary["rv_jitter_kms"]["median"] < float(np.median(data.sb2_rv.sigma_kms))


def test_masses_and_parallax_are_correlated_in_the_posterior():
    """The degeneracy structure a point estimate cannot report."""
    truth, data, response = _case()
    free = ("m1_msun", "m2_msun", "parallax_mas")
    start = replace(truth, m1_msun=1.28, m2_msun=0.83)
    best = fit_joint(data, start, response, free_names=free)
    posterior = sample_posterior(data, best.params, response, free_names=free,
                                 n_walkers=12, n_steps=400, n_burn=150, seed=2)
    correlation = np.corrcoef(posterior.flat().T)
    assert correlation[0, 1] > 0.3        # the two masses track each other
    assert correlation[0, 2] < 0.0        # and both trade against distance
