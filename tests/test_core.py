import numpy as np
import pytest

from raa_orbit_model.kepler import (
    BinaryParams, solve_kepler, relative_astrometry_mas,
    radial_velocities_kms, rv_semiamplitudes_kms,
)
from raa_orbit_model.gaia import (
    project_along_scan, photocentre_along_scan, component_al_positions,
    blended_gaussian_peak,
)


def params(**kw):
    base = dict(
        period_yr=2.0,
        t_peri_yr=0.0,
        eccentricity=0.2,
        inclination_deg=70.0,
        omega_deg=40.0,
        node_deg=110.0,
        m1_msun=1.2,
        m2_msun=0.8,
        parallax_mas=20.0,
        gamma_kms=5.0,
        beta_g=0.25,
    )
    base.update(kw)
    return BinaryParams(**base)


def test_kepler_equation():
    e = 0.73
    M = np.linspace(-9 * np.pi, 9 * np.pi, 101)
    E = solve_kepler(M, e)
    assert np.max(np.abs(E - e * np.sin(E) - M)) < 1e-11


def test_invalid_eccentricity_rejected():
    with pytest.raises(ValueError):
        params(eccentricity=1.0).validate()


def test_keplers_third_law_solar_units():
    p = params(period_yr=3.0, m1_msun=1.0, m2_msun=2.0)
    assert p.a_rel_au**3 == pytest.approx(3.0 * 3.0**2)


def test_relative_astrometry_parallax_scaling():
    p1 = params(parallax_mas=10.0)
    p2 = params(parallax_mas=20.0)
    t = np.array([0.3, 0.7])
    assert np.allclose(relative_astrometry_mas(t, p2), 2 * relative_astrometry_mas(t, p1))


def test_sb2_amplitude_ratio_is_inverse_mass_ratio():
    p = params(m1_msun=1.5, m2_msun=0.5)
    K1, K2 = rv_semiamplitudes_kms(p)
    assert K1 / K2 == pytest.approx(p.m2_msun / p.m1_msun)


def test_barycentric_rv_momentum_balance():
    p = params()
    v = radial_velocities_kms(np.linspace(0, p.period_yr, 30), p)
    lhs = p.m1_msun * (v[:, 0] - p.gamma_kms) + p.m2_msun * (v[:, 1] - p.gamma_kms)
    assert np.max(np.abs(lhs)) < 1e-10


def test_scan_angle_convention():
    east = np.array([3.0])
    north = np.array([4.0])
    assert project_along_scan(east, north, 0.0)[0] == pytest.approx(4.0)
    assert project_along_scan(east, north, 90.0)[0] == pytest.approx(3.0)


def test_photocentre_identity():
    d = np.array([-100.0, 20.0, 50.0])
    B, beta = 0.4, 0.1
    x1, x2 = component_al_positions(d, B)
    direct = (1 - beta) * x1 + beta * x2
    assert np.allclose(direct, photocentre_along_scan(d, B, beta))


def test_blended_peak_tends_to_photocentre_when_unresolved():
    B, beta, sigma = 0.4, 0.2, 100.0
    for d in [1e-3, 1e-2, 0.1]:
        peak = blended_gaussian_peak(d, B, beta, sigma)
        ph = photocentre_along_scan(d, B, beta)
        assert peak == pytest.approx(ph, abs=2e-5)


def test_blended_peak_tends_to_primary_when_well_separated():
    B, beta, sigma = 0.4, 0.2, 10.0
    d = 200.0
    x1, _ = component_al_positions(d, B)
    peak = blended_gaussian_peak(d, B, beta, sigma)
    assert peak == pytest.approx(x1, abs=1e-3)
