import numpy as np
import pytest

from raa_orbit_model.kepler import (
    AU_PER_YR_TO_KMS,
    BinaryParams,
    anomalies,
    radial_velocities_kms,
    relative_astrometry_mas,
    relative_position_au,
)


def params(**kw):
    base = dict(
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
    base.update(kw)
    return BinaryParams(**base)


def test_relative_periastron_argument_is_primary_omega_plus_180():
    p = params(omega_deg=55.0)
    assert p.omega_relative_deg == pytest.approx(235.0)


def test_primary_periastron_relative_position_is_opposite_primary_node():
    p = params(omega_deg=0.0, node_deg=0.0)
    east, north = relative_astrometry_mas(np.array([p.t_peri_yr]), p)[0]
    expected = p.a_rel_au * (1.0 - p.eccentricity) * p.parallax_mas
    assert east == pytest.approx(0.0, abs=1e-10)
    assert north == pytest.approx(-expected, rel=1e-12)


def test_node_position_angle_is_measured_north_through_east():
    p = params(omega_deg=0.0, node_deg=90.0)
    east, north = relative_astrometry_mas(np.array([p.t_peri_yr]), p)[0]
    expected = p.a_rel_au * (1.0 - p.eccentricity) * p.parallax_mas
    assert east == pytest.approx(-expected, rel=1e-12)
    assert north == pytest.approx(0.0, abs=1e-10)


def test_relative_astrometry_matches_thiele_innes_north_east_equations():
    p = params()
    t = np.array([0.37, 0.91, 1.44])
    _, E, _ = anomalies(t, p)

    e = p.eccentricity
    X = np.cos(E) - e
    Y = np.sqrt(1.0 - e**2) * np.sin(E)
    a = p.a_rel_au * p.parallax_mas
    Om = np.deg2rad(p.node_deg)
    omega_rel = np.deg2rad(p.omega_relative_deg)
    inc = np.deg2rad(p.inclination_deg)

    A = a * (np.cos(Om) * np.cos(omega_rel) - np.sin(Om) * np.sin(omega_rel) * np.cos(inc))
    B = a * (np.sin(Om) * np.cos(omega_rel) + np.cos(Om) * np.sin(omega_rel) * np.cos(inc))
    F = a * (-np.cos(Om) * np.sin(omega_rel) - np.sin(Om) * np.cos(omega_rel) * np.cos(inc))
    G = a * (-np.sin(Om) * np.sin(omega_rel) + np.cos(Om) * np.cos(omega_rel) * np.cos(inc))

    expected_north = A * X + F * Y
    expected_east = B * X + G * Y
    got = relative_astrometry_mas(t, p)

    assert np.allclose(got[:, 0], expected_east, atol=1e-11)
    assert np.allclose(got[:, 1], expected_north, atol=1e-11)


def test_relative_los_derivative_equals_secondary_minus_primary_rv():
    p = params()
    t0 = 0.73
    dt = 1e-6
    z_minus = relative_position_au(np.array([t0 - dt]), p)[0, 2]
    z_plus = relative_position_au(np.array([t0 + dt]), p)[0, 2]
    dzdt_kms = (z_plus - z_minus) / (2.0 * dt) * AU_PER_YR_TO_KMS

    rv1, rv2 = radial_velocities_kms(np.array([t0]), p)[0]
    assert dzdt_kms == pytest.approx(rv2 - rv1, rel=2e-9, abs=2e-9)


def test_primary_rv_uses_primary_argument_of_periastron():
    p = params(omega_deg=0.0)
    rv1, rv2 = radial_velocities_kms(np.array([p.t_peri_yr]), p)[0]
    assert rv1 > p.gamma_kms
    assert rv2 < p.gamma_kms
