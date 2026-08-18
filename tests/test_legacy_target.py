"""Harness for the real-binary V6a fit.

The real GJ 765.2 measurements are not in this repository, so these tests
exercise the loader and the fit path with a synthetic file written in the same
format from a known orbit. That validates the machinery without inventing
measurements: recovering an injected orbit through the legacy reader proves the
polar-to-tangent-plane conversion, the uncertainty handling and the node-branch
logic are right.
"""

from pathlib import Path

import numpy as np
import pytest
from dataclasses import replace

from raa_orbit_model.fit import fit_joint
from raa_orbit_model.kepler import BinaryParams, radial_velocities_kms, relative_astrometry_mas
from raa_orbit_model.legacy_target import (
    initial_guess_from_header,
    legacy_joint_data,
    node_branches,
    parse_legacy_file,
    summarise_fit,
    tangent_plane_offsets_mas,
    total_mass_msun,
    total_mass_uncertainty,
)
from raa_orbit_model.model import GaiaResponseConfig

FREE = (
    "period_yr", "t_peri_yr", "eccentricity", "inclination_deg",
    "omega_deg", "node_deg", "m1_msun", "m2_msun", "parallax_mas", "gamma_kms",
)

TRUTH = BinaryParams(
    period_yr=11.919, t_peri_yr=1.5, eccentricity=0.240, inclination_deg=80.2,
    omega_deg=70.0, node_deg=293.0, m1_msun=0.831, m2_msun=0.763,
    parallax_mas=31.0, gamma_kms=-12.0,
)
EPOCH0 = 1980.0


def _write_synthetic_legacy(path, n_ast=11, n_rv=44, sigma_mas=5.0, sigma_rv=0.35, seed=0):
    """Write a noiseless file in the real legacy CSV layout, built from TRUTH.

    Velocities carry JD-2400000 epochs and visual rows decimal years, matching
    the two time systems the real file uses.
    """
    rng = np.random.default_rng(seed)
    t_ast = np.sort(rng.uniform(0.5, 14.0, n_ast))
    t_rv = np.sort(rng.uniform(0.5, 14.0, n_rv))

    rel = relative_astrometry_mas(t_ast, TRUTH)
    east, north = rel[:, 0], rel[:, 1]
    rho_mas = np.hypot(east, north)
    theta_deg = np.degrees(np.arctan2(east, north)) % 360.0
    rv = radial_velocities_kms(t_rv, TRUTH)
    to_jd = lambda y: (y - 2000.0) * 365.25 + 51545.0

    lines = ["# synthetic legacy file for testing only", "#Objectinfo",
             "Object, TEST", "par,54.27", "", "# RV1 Measurements"]
    for t, (v1, _) in zip(t_rv, rv):
        lines.append(f"{to_jd(t + EPOCH0):.4f},{v1:.5f},{sigma_rv:.3f},Va,COR")
    lines += ["", "# RV2 Measurements"]
    for t, (_, v2) in zip(t_rv, rv):
        lines.append(f"{to_jd(t + EPOCH0):.4f},{v2:.5f},{sigma_rv:.3f},Vb,COR")
    lines += ["", "# Visual Measurements"]
    for t, th, r in zip(t_ast, theta_deg, rho_mas):
        lines.append(f"{t + EPOCH0:.4f},{th:.4f},{r / 1000.0:.6f},{sigma_mas / 1000.0:.6f},I1,M1")
    path.write_text("\n".join(lines) + "\n")
    return t_ast, t_rv


def test_tangent_plane_conversion_uses_north_through_east():
    east, north = tangent_plane_offsets_mas([0.0, 90.0, 180.0], [100.0, 100.0, 100.0])
    assert east == pytest.approx([0.0, 100.0, 0.0], abs=1e-9)
    assert north == pytest.approx([100.0, 0.0, -100.0], abs=1e-9)


def test_parses_blocks_header_and_units(tmp_path):
    path = tmp_path / "legacy.dat"
    _write_synthetic_legacy(path)
    data = parse_legacy_file(path, separation_unit="arcsec")

    assert data.n_astrometry == 11
    assert data.n_rv == 44
    assert data.n_constraints == 2 * 11 + 2 * 44 == 110
    assert data.header_float("par") == pytest.approx(54.27)
    # arcsec input must be scaled to mas
    assert data.separation_mas.max() > 50.0
    assert data.astrometry_sigma_mas == pytest.approx(np.full(11, 5.0))


def test_separation_unit_is_never_guessed(tmp_path):
    path = tmp_path / "legacy.dat"
    _write_synthetic_legacy(path)
    arcsec = parse_legacy_file(path, separation_unit="arcsec")
    mas = parse_legacy_file(path, separation_unit="mas")
    assert arcsec.separation_mas == pytest.approx(1000.0 * mas.separation_mas)
    with pytest.raises(ValueError, match="separation_unit"):
        parse_legacy_file(path, separation_unit="degrees")


def test_gaia_channel_is_empty_so_the_response_is_irrelevant(tmp_path):
    path = tmp_path / "legacy.dat"
    _write_synthetic_legacy(path)
    joint = legacy_joint_data(parse_legacy_file(path))
    assert len(joint.gaia_al.times_yr) == 0
    assert len(joint.relative_astrometry.times_yr) == 11
    assert len(joint.sb2_rv.times_yr) == 44


def test_isotropic_uncertainty_becomes_a_diagonal_covariance(tmp_path):
    path = tmp_path / "legacy.dat"
    _write_synthetic_legacy(path, sigma_mas=5.0)
    joint = legacy_joint_data(parse_legacy_file(path))
    cov = joint.relative_astrometry.covariance_mas2
    assert cov.shape == (11, 2, 2)
    assert cov[:, 0, 0] == pytest.approx(25.0)
    assert cov[:, 1, 1] == pytest.approx(25.0)
    assert cov[:, 0, 1] == pytest.approx(0.0)


def test_node_branches_differ_by_180_in_both_angles():
    first, second = node_branches(TRUTH)
    assert (second.node_deg - first.node_deg) % 360.0 == pytest.approx(180.0)
    assert (second.omega_deg - first.omega_deg) % 360.0 == pytest.approx(180.0)


def test_injected_orbit_is_recovered_through_the_legacy_reader(tmp_path):
    """End-to-end: the harness must reproduce a known orbit from a legacy file."""
    path = tmp_path / "legacy.dat"
    _write_synthetic_legacy(path)
    data = parse_legacy_file(path)
    joint = legacy_joint_data(data)

    start = replace(
        TRUTH, period_yr=TRUTH.period_yr * 1.02, eccentricity=0.27,
        inclination_deg=78.0, m1_msun=0.87, m2_msun=0.72,
        parallax_mas=35.0, gamma_kms=-11.0, t_peri_yr=1.4,
    )
    result = fit_joint(joint, start, GaiaResponseConfig("photocentre"), free_names=FREE)

    assert result.success
    assert result.n_free == 10
    assert result.n_residuals == 110
    assert result.dof == 100
    assert result.params.period_yr == pytest.approx(TRUTH.period_yr, rel=1e-3)
    assert result.params.eccentricity == pytest.approx(TRUTH.eccentricity, abs=5e-3)
    assert result.params.inclination_deg == pytest.approx(TRUTH.inclination_deg, abs=0.5)
    assert total_mass_msun(result.params) == pytest.approx(
        total_mass_msun(TRUTH), rel=5e-3
    )
    assert result.params.parallax_mas == pytest.approx(TRUTH.parallax_mas, rel=0.02)


def test_summary_record_matches_the_manuscript_table_fields(tmp_path):
    path = tmp_path / "legacy.dat"
    _write_synthetic_legacy(path)
    data = parse_legacy_file(path)
    joint = legacy_joint_data(data)
    result = fit_joint(joint, TRUTH, GaiaResponseConfig("photocentre"), free_names=FREE)
    summary = summarise_fit(result, data)
    for field in ("chi2", "dof", "reduced_chi2", "period_yr", "eccentricity",
                  "inclination_deg", "node_deg", "omega_relative_deg",
                  "total_mass_msun", "parallax_mas", "n_constraints"):
        assert field in summary
    assert summary["n_constraints"] == 110
    assert summary["dof"] == 100


def test_total_mass_error_includes_the_component_covariance():
    """The masses are strongly correlated; quadrature alone misstates the total."""
    names = ("m1_msun", "m2_msun")
    # Perfectly correlated: total-mass error must be the linear sum, not sqrt(2)*s.
    covariance = np.array([[0.01, 0.01], [0.01, 0.01]])
    assert total_mass_uncertainty(covariance, names) == pytest.approx(0.2)
    # Anti-correlated: the total is far better determined than either mass.
    covariance = np.array([[0.01, -0.01], [-0.01, 0.01]])
    assert total_mass_uncertainty(covariance, names) == pytest.approx(0.0, abs=1e-12)
    # Uncorrelated reduces to quadrature.
    covariance = np.array([[0.01, 0.0], [0.0, 0.01]])
    assert total_mass_uncertainty(covariance, names) == pytest.approx(np.sqrt(0.02))


def test_total_mass_error_is_undefined_when_masses_are_not_free():
    assert np.isnan(total_mass_uncertainty(np.eye(2), ("period_yr", "eccentricity")))


@pytest.mark.parametrize("body,message", [
    ("# Visual Measurements\n1980.0,10.0,0.2\n", "visual row needs 4"),
    ("# RV1 Measurements\n45000.0,1.0,0.3\n", "no relative astrometry"),
    ("# Visual Measurements\n1980.0,10.0,0.2,0.01\n", "no radial velocities"),
])
def test_malformed_files_are_reported_clearly(tmp_path, body, message):
    path = tmp_path / "bad.dat"
    path.write_text(body)
    with pytest.raises(ValueError, match=message):
        parse_legacy_file(path)


# --- Real GJ 765.2 data (V6a) -------------------------------------------------
# The committed legacy file reproduces the manuscript table tab:gl765. These
# tests pin that result so the real-data section cannot drift silently.

LEGACY_CSV = Path(__file__).resolve().parents[1] / "data" / "gl765_legacy.csv"


def _real_fit():
    data = parse_legacy_file(LEGACY_CSV)
    joint = legacy_joint_data(data)
    start = initial_guess_from_header(data)
    best = None
    for branch in node_branches(start):
        result = fit_joint(joint, branch, GaiaResponseConfig("photocentre"), free_names=FREE)
        if best is None or result.chi2 < best.chi2:
            best = result
    return data, best


def test_real_file_has_the_dimensions_quoted_in_the_manuscript():
    data = parse_legacy_file(LEGACY_CSV)
    assert data.n_astrometry == 11
    assert data.n_rv == 44
    assert data.n_constraints == 110
    # The two 1993 speckle positions carry 5 mas uncertainties.
    assert sorted(data.astrometry_sigma_mas)[:2] == pytest.approx([5.0, 5.0])


def test_header_distinguishes_node_from_argument_of_periastron():
    """W and w differ only by case; conflating them would corrupt the orbit."""
    data = parse_legacy_file(LEGACY_CSV)
    assert data.header_float("W") == pytest.approx(106.34)
    assert data.header_float("w") == pytest.approx(89.4)


def test_velocity_and_visual_epochs_are_brought_to_one_time_system():
    data = parse_legacy_file(LEGACY_CSV)
    assert 1983.0 < data.rv_epoch_yr.min() < 1984.0
    assert 1994.0 < data.rv_epoch_yr.max() < 1995.0
    assert 1971.0 < data.astrometry_epoch_yr.min() < 1972.0


def test_unconverted_velocity_epochs_are_rejected():
    with pytest.raises(ValueError, match="do not overlap"):
        parse_legacy_file(LEGACY_CSV, rv_time_system="decimalyear")


def test_real_fit_reproduces_the_manuscript_table():
    """tab:gl765, regenerated from the committed data."""
    data, result = _real_fit()
    assert result.success
    assert result.n_residuals == 110
    assert result.dof == 100
    assert result.chi2 == pytest.approx(104.6, abs=0.3)
    assert result.reduced_chi2 == pytest.approx(1.046, abs=0.005)
    assert result.params.period_yr == pytest.approx(11.7284, abs=1e-3)
    assert result.params.eccentricity == pytest.approx(0.24888, abs=1e-4)
    assert result.params.inclination_deg == pytest.approx(81.834, abs=0.01)
    assert result.params.node_deg % 360.0 == pytest.approx(289.07, abs=0.05)
    assert result.params.omega_relative_deg == pytest.approx(251.79, abs=0.05)
    assert result.params.parallax_mas == pytest.approx(35.44, abs=0.01)
    assert total_mass_msun(result.params) == pytest.approx(1.5897, abs=1e-3)


def test_total_mass_is_close_to_the_published_solution():
    """Balega et al. (2007) give 1.594 Msun; the manuscript quotes 0.27% apart."""
    _, result = _real_fit()
    assert abs(total_mass_msun(result.params) / 1.594 - 1.0) < 0.005


def test_component_labels_do_not_map_onto_the_published_convention():
    """Va has the larger amplitude, so it is the LESS massive component."""
    data = parse_legacy_file(LEGACY_CSV)
    assert data.header_float("K1") > data.header_float("K2")
    _, result = _real_fit()
    assert result.params.m2_msun > result.params.m1_msun


def test_legacy_header_parallax_is_strongly_disfavoured():
    """The joint data force the scale away from the inconsistent 54.27 mas."""
    data, result = _real_fit()
    joint = legacy_joint_data(data)
    free_no_parallax = tuple(n for n in FREE if n != "parallax_mas")
    held = {}
    for value in (54.27, 31.0):
        fit = fit_joint(joint, replace(result.params, parallax_mas=value),
                        GaiaResponseConfig("photocentre"), free_names=free_no_parallax)
        held[value] = fit.chi2
    assert held[54.27] > result.chi2 + 50.0
    assert held[31.0] == pytest.approx(108.6, abs=0.5)
    assert held[54.27] == pytest.approx(168.8, abs=0.5)
