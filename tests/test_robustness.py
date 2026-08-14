import numpy as np
import pytest

from raa_orbit_model.fit import ALL_PARAMETER_NAMES, fit_joint
from raa_orbit_model.gaia import blended_gaussian_peak, blended_gaussian_response
from raa_orbit_model.kepler import BinaryParams
from raa_orbit_model.model import GaiaResponseConfig
from raa_orbit_model.robustness import (
    compare_under_misspecification,
    shape_misspecification_scan,
    summarise_by,
)
from raa_orbit_model.scanning import schedule_from_arrays
from raa_orbit_model.synthetic import simulate_joint_data
from raa_orbit_model.experiments import perturbed_start, single_peak_schedule_for_response

SIGMA = 50.0


def _truth(a_over_sigma=1.0, beta_g=0.25):
    base = BinaryParams(
        period_yr=2.0, t_peri_yr=0.15, eccentricity=0.25, inclination_deg=72.0,
        omega_deg=55.0, node_deg=120.0, m1_msun=1.25, m2_msun=0.85,
        parallax_mas=20.0, gamma_kms=7.0, beta_g=beta_g,
    )
    return BinaryParams(**{**base.__dict__, "parallax_mas": a_over_sigma * SIGMA / base.a_rel_au})


def _schedule(n=60, span=5.0, seed=3):
    rng = np.random.default_rng(seed)
    return schedule_from_arrays(
        np.sort(rng.uniform(0.0, span, n)), rng.uniform(0.0, 360.0, n),
        ra_deg=120.0, dec_deg=30.0, source="robustness test",
    )


def test_equal_width_path_is_unchanged_by_the_new_option():
    """The frozen results must remain reproducible bit-for-bit."""
    implicit = blended_gaussian_peak(30.0, 0.4, 0.25, SIGMA)
    explicit = blended_gaussian_peak(30.0, 0.4, 0.25, SIGMA, sigma_secondary_mas=SIGMA)
    assert implicit == pytest.approx(explicit, abs=1e-12)


def test_unequal_width_moves_the_peak_towards_the_narrower_component():
    """A wider secondary has a lower peak, pulling the maximum to the primary."""
    d, B, beta = 40.0, 0.4, 0.25
    equal = blended_gaussian_response(d, B, beta, SIGMA)
    wide = blended_gaussian_response(d, B, beta, SIGMA, sigma_secondary_mas=1.6 * SIGMA)
    primary = -B * d
    assert abs(wide.al_mas - primary) < abs(equal.al_mas - primary)


def test_unequal_width_response_still_counts_modes_and_flags_splitting():
    B, beta = 0.4, 0.45
    close = blended_gaussian_response(10.0, B, beta, 10.0, sigma_secondary_mas=11.0)
    split = blended_gaussian_response(400.0, B, beta, 10.0, sigma_secondary_mas=11.0)
    assert close.n_peaks == 1 and np.isfinite(close.al_mas)
    assert split.n_peaks == 2 and np.isnan(split.al_mas)


def test_response_width_is_recovered_when_fitted_from_a_wrong_start():
    truth = _truth()
    schedule = _schedule()
    injected = GaiaResponseConfig("blended_gaussian_peak", SIGMA)
    selection = single_peak_schedule_for_response(truth, schedule, injected)
    data = simulate_joint_data(truth, injected, selection.schedule, seed=0,
                               baseline_yr=schedule.mission_span_yr)
    start = GaiaResponseConfig("blended_gaussian_peak", 1.2 * SIGMA,
                               allow_multi_peak_continuation=True)
    result = fit_joint(data, perturbed_start(truth), start,
                       free_names=tuple(ALL_PARAMETER_NAMES) + ("sigma_response_mas",))
    assert result.success
    assert result.fitted_sigma_response_mas == pytest.approx(SIGMA, rel=0.05)


def test_fitting_a_width_is_rejected_for_the_photocentre_model():
    truth = _truth()
    schedule = _schedule()
    injected = GaiaResponseConfig("blended_gaussian_peak", SIGMA)
    selection = single_peak_schedule_for_response(truth, schedule, injected)
    data = simulate_joint_data(truth, injected, selection.schedule, seed=0,
                               baseline_yr=schedule.mission_span_yr)
    with pytest.raises(ValueError, match="no width to fit"):
        fit_joint(data, perturbed_start(truth), GaiaResponseConfig("photocentre"),
                  free_names=tuple(ALL_PARAMETER_NAMES) + ("sigma_response_mas",))


def test_in_family_injection_reproduces_the_published_raa_advantage():
    row = compare_under_misspecification(
        _truth(), _schedule(), sigma_response_mas=SIGMA, secondary_width_ratio=1.0, seed=0
    )
    assert row["raa_favoured"]
    assert row["delta_chi2"] > 100.0
    assert row["gaia_n_multi_peak_flagged"] == 0


def test_shape_misspecification_erodes_the_raa_advantage():
    """The decisive control: a wrong profile shape must reduce the advantage."""
    kw = dict(sigma_response_mas=SIGMA, seed=0)
    in_family = compare_under_misspecification(_truth(), _schedule(),
                                               secondary_width_ratio=1.0, **kw)
    out_of_family = compare_under_misspecification(_truth(), _schedule(),
                                                   secondary_width_ratio=1.4, **kw)
    assert out_of_family["delta_chi2"] < in_family["delta_chi2"]


def test_retaining_multi_peak_epochs_reaches_beyond_the_splitting_boundary():
    truth = _truth(a_over_sigma=4.0, beta_g=0.45)
    row = compare_under_misspecification(
        truth, _schedule(), sigma_response_mas=SIGMA, retain_multi_peak=True, seed=0
    )
    assert row["gaia_n_multi_peak_flagged"] > 0
    assert np.isfinite(row["delta_chi2"])


def test_scan_and_summary_group_by_width_ratio():
    rows = shape_misspecification_scan(
        _truth(), _schedule(), secondary_width_ratios=(1.0, 1.4), seeds=(0, 1)
    )
    assert len(rows) == 4
    summary = summarise_by(rows, "secondary_width_ratio")
    assert [s["secondary_width_ratio"] for s in summary] == [1.0, 1.4]
    assert all(s["n"] == 2 for s in summary)
    assert all(0.0 <= s["fraction_raa_favoured"] <= 1.0 for s in summary)
