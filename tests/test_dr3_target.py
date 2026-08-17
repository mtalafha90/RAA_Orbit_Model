import math

from raa_orbit_model.dr3_target import (
    GJ765_DEC_DEG,
    GJ765_HIP,
    GJ765_RA_DEG,
    flux_fraction_from_delta_mag,
    gj765_cone_query,
    gj765_photocentre_benchmark,
    gj765_target_query,
    photocentre_axis_mas,
    secondary_mass_fraction,
    summarize_target_rows,
)
from raa_orbit_model.dr3_validation import campbell_to_thiele_innes


def test_gj765_query_uses_hipparcos_crossmatch_and_left_join():
    query = gj765_target_query()
    assert "gaiadr3.hipparcos2_best_neighbour" in query
    assert f"original_ext_source_id = {GJ765_HIP}" in query
    assert "LEFT OUTER JOIN gaiadr3.nss_two_body_orbit" in query
    assert "ipd_frac_multi_peak" in query
    assert "ipd_gof_harmonic_amplitude" in query
    assert "scan_direction_mean_k2" in query


def test_gj765_cone_query_retains_simbad_position_only_as_fallback():
    query = gj765_cone_query(10.0)
    assert f"{GJ765_RA_DEG:.12f}" in query
    assert f"{GJ765_DEC_DEG:.12f}" in query
    assert "CIRCLE('ICRS'" in query


def test_gj765_m0_photocentre_benchmark():
    beta = flux_fraction_from_delta_mag(0.65)
    B = secondary_mass_fraction(0.831, 0.763)
    axis = photocentre_axis_mas(189.0, 0.831, 0.763, beta)
    assert math.isclose(beta, 0.3546475495606176, rel_tol=0, abs_tol=1e-12)
    assert math.isclose(B, 0.4786700125470515, rel_tol=0, abs_tol=1e-12)
    assert math.isclose(axis, 23.440245504436003, rel_tol=0, abs_tol=1e-10)
    bench = gj765_photocentre_benchmark()
    assert math.isclose(bench["predicted_M0_photocentre_axis_mas"], axis)


def test_target_summary_selects_nearest_source_for_cone_fallback_and_converts_nss():
    A, B, F, G = campbell_to_thiele_innes(23.0, 80.0, 250.0, 293.0)
    rows = [
        {
            "source_id": "far",
            "separation_arcsec": "2.0",
            "parallax": "5.0",
            "parallax_error": "1.0",
            "ruwe": "1.0",
            "nss_solution_type": "",
        },
        {
            "source_id": "target",
            "separation_arcsec": "0.02",
            "parallax": "32.0",
            "parallax_error": "0.5",
            "ruwe": "4.2",
            "ipd_frac_multi_peak": "12",
            "ipd_gof_harmonic_amplitude": "0.3",
            "ipd_gof_harmonic_phase": "108",
            "non_single_star": "1",
            "nss_solution_type": "Orbital",
            "nss_period_days": "4350",
            "nss_eccentricity": "0.24",
            "a_thiele_innes": str(A),
            "b_thiele_innes": str(B),
            "f_thiele_innes": str(F),
            "g_thiele_innes": str(G),
        },
    ]
    summary = summarize_target_rows(rows)
    assert summary.source_id == "target"
    assert summary.nss_solution_types == ("Orbital",)
    assert math.isclose(summary.nearest_nss_photocentre_axis_mas, 23.0, abs_tol=1e-10)
    assert math.isclose(summary.nearest_nss_inclination_deg, 80.0, abs_tol=1e-10)


def test_target_summary_preserves_hip2_metadata_and_handles_no_nss():
    rows = [{
        "source_id": "target",
        "separation_arcsec": "0.01",
        "hipparcos2_number_of_neighbours": "1",
        "hipparcos2_gaia_astrometric_params": "5",
        "parallax": "33.0",
        "parallax_error": "0.4",
        "ruwe": "1.7",
        "ipd_frac_multi_peak": "0",
        "ipd_gof_harmonic_amplitude": "0.2",
        "ipd_gof_harmonic_phase": "100",
        "non_single_star": "0",
        "nss_solution_type": "",
    }]
    summary = summarize_target_rows(rows)
    assert summary.hipparcos2_number_of_neighbours == 1
    assert summary.hipparcos2_gaia_astrometric_params == 5
    assert summary.nss_orbit_count == 0
    assert summary.nss_solution_types == ()
    assert math.isnan(summary.nearest_nss_photocentre_axis_mas)
