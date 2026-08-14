import numpy as np
import pytest

from raa_orbit_model.experiments import (
    compare_models_once,
    single_peak_schedule_for_response,
)
from raa_orbit_model.fit import fit_joint
from raa_orbit_model.kepler import BinaryParams
from raa_orbit_model.model import GaiaResponseConfig, predict_relative_astrometry
from raa_orbit_model.scanning import schedule_from_arrays
from raa_orbit_model.synthetic import simulate_joint_data


def truth(**kw):
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


def make_schedule(n=80, span_yr=5.0):
    u = (np.arange(n, dtype=float) + 0.5) / n
    times = span_yr * u**1.35
    angles = np.mod(31.0 + 211.0 * u + 53.0 * np.sin(9.0 * u), 360.0)
    return schedule_from_arrays(
        times,
        angles,
        ra_deg=20.0,
        dec_deg=45.0,
        release="custom",
        source="deterministic-test-fixture",
    )


def test_noise_free_joint_recovery_of_key_physical_parameters():
    p = truth()
    schedule = make_schedule(80)
    data = simulate_joint_data(
        p,
        GaiaResponseConfig("photocentre"),
        schedule,
        seed=4,
        n_ast=30,
        n_rv=60,
        ast_sigma_mas=1e-5,
        rv_sigma_kms=1e-5,
        gaia_sigma_mas=1e-5,
    )
    from raa_orbit_model.model import predict_sb2_rv, predict_gaia_orbital_al
    from raa_orbit_model.synthetic import RelativeAstrometryData, SB2RVData, GaiaALData, JointData
    ast = data.relative_astrometry
    rv = data.sb2_rv
    gaia = data.gaia_al
    data = JointData(
        RelativeAstrometryData(ast.times_yr, predict_relative_astrometry(ast.times_yr, p), ast.covariance_mas2),
        SB2RVData(rv.times_yr, predict_sb2_rv(rv.times_yr, p), rv.sigma_kms),
        GaiaALData(
            gaia.times_yr,
            gaia.scan_angle_deg,
            predict_gaia_orbital_al(
                gaia.times_yr,
                gaia.scan_angle_deg,
                p,
                GaiaResponseConfig("photocentre"),
            ),
            gaia.sigma_mas,
            schedule_source=gaia.schedule_source,
            ra_deg=gaia.ra_deg,
            dec_deg=gaia.dec_deg,
            release=gaia.release,
        ),
    )
    initial = BinaryParams(
        period_yr=2.03,
        t_peri_yr=0.19,
        eccentricity=0.27,
        inclination_deg=74.0,
        omega_deg=59.0,
        node_deg=117.0,
        m1_msun=1.30,
        m2_msun=0.81,
        parallax_mas=20.8,
        gamma_kms=7.2,
        beta_g=0.23,
    )
    result = fit_joint(data, initial, GaiaResponseConfig("photocentre"), max_nfev=5000)
    assert result.success
    assert result.params.period_yr == pytest.approx(p.period_yr, rel=1e-7)
    assert result.params.eccentricity == pytest.approx(p.eccentricity, abs=1e-7)
    assert result.params.m1_msun == pytest.approx(p.m1_msun, rel=1e-6)
    assert result.params.m2_msun == pytest.approx(p.m2_msun, rel=1e-6)
    assert result.params.parallax_mas == pytest.approx(p.parallax_mas, rel=1e-6)
    assert result.params.beta_g == pytest.approx(p.beta_g, abs=1e-6)


def test_single_peak_schedule_filters_aligned_equal_light_epoch():
    p = truth(beta_g=0.5, parallax_mas=100.0)
    t = 0.5
    east, north = predict_relative_astrometry(np.array([t]), p)[0]
    psi_aligned = np.degrees(np.arctan2(east, north))
    psi_perpendicular = psi_aligned + 90.0
    schedule = schedule_from_arrays(
        [t, t],
        [psi_aligned, psi_perpendicular],
        ra_deg=20.0,
        dec_deg=45.0,
        release="custom",
        source="mode-filter-test",
    )
    selection = single_peak_schedule_for_response(
        p,
        schedule,
        GaiaResponseConfig("blended_gaussian_peak", 30.0),
    )
    assert selection.n_total == 2
    assert selection.n_multi_peak == 1
    assert selection.n_single_peak == 1
    assert selection.schedule.n_transits == 1
    assert selection.multi_peak_fraction == pytest.approx(0.5)


def test_resolution_aware_fit_beats_misspecified_photocentre_in_resolved_regime():
    p = truth(parallax_mas=80.0, beta_g=0.20)
    _, photo, raa = compare_models_once(
        p,
        make_schedule(60),
        sigma_response_mas=25.0,
        seed=2,
        n_ast=20,
        n_rv=40,
        ast_sigma_mas=0.15,
        rv_sigma_kms=0.08,
        gaia_sigma_mas=0.05,
    )
    assert photo.success and raa.success
    assert raa.chi2 < photo.chi2
    assert abs(raa.params.beta_g - p.beta_g) < abs(photo.params.beta_g - p.beta_g)
