import numpy as np

from raa_orbit_model.kepler import BinaryParams
from raa_orbit_model.response_controls import (
    ExternalInformationLevel,
    external_information_strength_scan,
    matched_m2_scan,
    summarise_external_information,
    summarise_matched_m2,
    truth_at_response_point,
)
from raa_orbit_model.scanning import schedule_from_arrays


def _truth():
    return BinaryParams(
        period_yr=2.0,
        t_peri_yr=0.15,
        eccentricity=0.25,
        inclination_deg=72.0,
        omega_deg=55.0,
        node_deg=120.0,
        m1_msun=1.25,
        m2_msun=0.85,
        parallax_mas=20.0,
        gamma_kms=7.0,
        beta_g=0.20,
    )


def _schedule(n=36):
    times = np.linspace(0.0, 5.0, n)
    angles = np.mod(np.arange(n) * 137.5, 360.0)
    return schedule_from_arrays(
        times,
        angles,
        ra_deg=120.0,
        dec_deg=30.0,
        source="response control test",
    )


def test_truth_at_response_point_sets_requested_angular_scale_and_light_fraction():
    truth = truth_at_response_point(_truth(), alpha_mas=50.0, a_over_alpha=1.0, beta_g=0.25)
    assert truth.beta_g == 0.25
    assert truth.a_rel_au * truth.parallax_mas / 50.0 == np.testing.assert_allclose(
        truth.a_rel_au * truth.parallax_mas / 50.0, 1.0, rtol=0, atol=1e-12
    )


def test_external_information_control_is_paired_and_tagged():
    levels = (
        ExternalInformationLevel("strong_test", 0.2, 0.1, n_ast=8, n_rv=12),
        ExternalInformationLevel("weak_test", 1.0, 0.5, n_ast=8, n_rv=12),
    )
    rows = external_information_strength_scan(
        _truth(),
        _schedule(),
        levels=levels,
        seeds=(0,),
        alpha_mas=50.0,
        a_over_alpha=0.6,
        beta_g=0.25,
        beta_over_alpha=1.5,
        gaia_sigma_mas=0.1,
    )
    assert len(rows) == 6
    assert {row["model"] for row in rows} == {"photocentre", "equal_width", "penoyre_oriented"}
    assert {row["external_level"] for row in rows} == {"strong_test", "weak_test"}
    assert all(row["seed"] == 0 for row in rows)
    summary = summarise_external_information(rows)
    assert len(summary) == 6
    assert all(item["n"] == 1 for item in summary)


def test_matched_m2_control_fits_only_the_injected_model_and_summarises():
    rows = matched_m2_scan(
        _truth(),
        _schedule(),
        seeds=(0, 1),
        alpha_mas=50.0,
        a_over_alpha=0.6,
        beta_g=0.25,
        beta_over_alpha=1.5,
        n_ast=8,
        n_rv=12,
        ast_sigma_mas=0.2,
        rv_sigma_kms=0.1,
        gaia_sigma_mas=0.1,
    )
    assert len(rows) == 2
    assert all(row["model"] == "penoyre_oriented" for row in rows)
    assert all(row["injection_model"] == "penoyre_oriented" for row in rows)
    assert all(row["success"] for row in rows)
    summary = summarise_matched_m2(rows)
    assert len(summary) == 1
    assert summary[0]["n"] == 2
    assert summary[0]["success_fraction"] == 1.0
