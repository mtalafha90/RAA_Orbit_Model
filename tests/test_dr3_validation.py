"""Validation step V6 machinery: Thiele-Innes to Campbell, and NSS loading.

DR3 publishes astrometric orbits as Thiele-Innes constants, so comparing them
against this project requires the conversion below. These tests exercise it
offline by round-tripping through the project's own convention, which is the
same one asserted in tests/test_orbit_conventions.py.
"""

import math

import numpy as np
import pytest

from raa_orbit_model.dr3_validation import (
    NSS_ASTROMETRIC_ORBIT_QUERY,
    campbell_table,
    campbell_to_thiele_innes,
    load_nss_csv,
    marginally_resolved_candidates,
    thiele_innes_to_campbell,
)
from raa_orbit_model.kepler import BinaryParams, anomalies, relative_astrometry_mas


def _angle_gap(a, b, modulus=360.0):
    diff = abs(a - b) % modulus
    return min(diff, modulus - diff)


@pytest.mark.parametrize("elements", [
    (10.0, 72.0, 55.0, 120.0),
    (5.0, 30.0, 200.0, 300.0),
    (2.5, 145.0, 10.0, 45.0),
    (8.0, 90.0, 270.0, 90.0),
    (1.0, 5.0, 350.0, 5.0),
])
def test_campbell_thiele_innes_round_trip(elements):
    axis, inclination, omega, node = elements
    constants = campbell_to_thiele_innes(axis, inclination, omega, node)
    recovered = thiele_innes_to_campbell(*constants)
    assert recovered.semi_major_axis_mas == pytest.approx(axis, rel=1e-10)
    assert recovered.inclination_deg == pytest.approx(inclination, abs=1e-8)
    # Astrometry alone fixes the node only modulo 180 degrees.
    assert _angle_gap(recovered.node_deg, node, 180.0) < 1e-8


def test_omega_and_node_share_a_180_degree_degeneracy():
    """Shifting both angles by 180 degrees leaves the constants unchanged."""
    first = campbell_to_thiele_innes(7.0, 63.0, 40.0, 100.0)
    second = campbell_to_thiele_innes(7.0, 63.0, 220.0, 280.0)
    assert np.allclose(first, second, atol=1e-12)


def test_conversion_matches_the_projects_own_orbit_model():
    """The published constants must describe the same ellipse this project draws."""
    params = BinaryParams(
        period_yr=2.0, t_peri_yr=0.15, eccentricity=0.25, inclination_deg=72.0,
        omega_deg=55.0, node_deg=120.0, m1_msun=1.25, m2_msun=0.85,
        parallax_mas=20.0,
    )
    axis_mas = params.a_rel_au * params.parallax_mas
    A, B, F, G = campbell_to_thiele_innes(
        axis_mas, params.inclination_deg, params.omega_relative_deg, params.node_deg
    )

    times = np.array([0.37, 0.91, 1.44])
    _, E, _ = anomalies(times, params)
    X = np.cos(E) - params.eccentricity
    Y = np.sqrt(1.0 - params.eccentricity**2) * np.sin(E)
    predicted = relative_astrometry_mas(times, params)

    assert np.allclose(predicted[:, 1], A * X + F * Y, atol=1e-9)   # North
    assert np.allclose(predicted[:, 0], B * X + G * Y, atol=1e-9)   # East

    recovered = thiele_innes_to_campbell(A, B, F, G)
    assert recovered.semi_major_axis_mas == pytest.approx(axis_mas, rel=1e-10)
    assert recovered.inclination_deg == pytest.approx(params.inclination_deg, abs=1e-8)


def test_degenerate_constants_are_rejected():
    with pytest.raises(ValueError, match="non-positive axis"):
        thiele_innes_to_campbell(0.0, 0.0, 0.0, 0.0)


def _write_export(path, rows):
    header = (
        "source_id,period,eccentricity,a_thiele_innes,b_thiele_innes,"
        "f_thiele_innes,g_thiele_innes,parallax,ruwe,ipd_frac_multi_peak,"
        "ipd_gof_harmonic_amplitude\n"
    )
    path.write_text(header + "".join(rows))


def test_loading_an_export_and_building_the_campbell_table(tmp_path):
    A, B, F, G = campbell_to_thiele_innes(6.0, 55.0, 30.0, 200.0)
    export = tmp_path / "nss.csv"
    _write_export(export, [
        f"1234567890,410.5,0.31,{A},{B},{F},{G},12.4,3.8,7.5,0.09\n",
        f"2234567890,120.0,0.05,{A},{B},{F},{G},9.1,1.02,0.0,0.001\n",
        "3234567890,300.0,0.2,,,,,8.0,2.2,5.0,0.05\n",   # unusable constants
    ])

    rows = load_nss_csv(export)
    assert len(rows) == 3
    assert isinstance(rows[0]["source_id"], str)
    assert rows[0]["period"] == pytest.approx(410.5)

    table = campbell_table(rows)
    assert len(table) == 2                       # the blank row is dropped
    assert table[0]["inclination_deg"] == pytest.approx(55.0, abs=1e-8)
    assert table[0]["photocentre_semi_major_mas"] == pytest.approx(6.0, rel=1e-10)


def test_candidate_selection_uses_the_published_duplicity_diagnostics(tmp_path):
    A, B, F, G = campbell_to_thiele_innes(6.0, 55.0, 30.0, 200.0)
    export = tmp_path / "nss.csv"
    _write_export(export, [
        f"1,410.5,0.31,{A},{B},{F},{G},12.4,3.8,7.5,0.09\n",     # resolved-ish
        f"2,120.0,0.05,{A},{B},{F},{G},9.1,1.02,0.0,0.001\n",    # clean single
        f"3,300.0,0.20,{A},{B},{F},{G},8.0,2.20,4.0,0.04\n",     # borderline
    ])
    table = campbell_table(load_nss_csv(export))
    selected = marginally_resolved_candidates(table, min_multi_peak_fraction=2.0)
    assert {row["source_id"] for row in selected} == {"1", "3"}

    stricter = marginally_resolved_candidates(
        table, min_multi_peak_fraction=2.0, min_harmonic_amplitude=0.05
    )
    assert {row["source_id"] for row in stricter} == {"1"}


def test_missing_columns_are_reported(tmp_path):
    export = tmp_path / "bad.csv"
    export.write_text("source_id,period\n1,410.5\n")
    with pytest.raises(ValueError, match="missing columns"):
        load_nss_csv(export)


def test_empty_export_is_rejected(tmp_path):
    export = tmp_path / "empty.csv"
    export.write_text("source_id,period\n")
    with pytest.raises(ValueError, match="empty"):
        load_nss_csv(export)


def test_query_targets_astrometric_orbits_with_duplicity_diagnostics():
    query = NSS_ASTROMETRIC_ORBIT_QUERY
    assert "gaiadr3.nss_two_body_orbit" in query
    assert "nss_solution_type = 'Orbital'" in query
    assert "ipd_frac_multi_peak" in query
    assert "ipd_gof_harmonic_amplitude" in query
    for column in ("a_thiele_innes", "b_thiele_innes", "f_thiele_innes", "g_thiele_innes"):
        assert column in query
