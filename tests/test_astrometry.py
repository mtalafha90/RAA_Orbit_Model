"""Absolute astrometric motion: parallax factors and joint recovery.

The parallax factor tests are geometric identities that hold independently of
the ephemeris precision, so they remain valid checks even though the Earth
position here is a low-precision analytic series.
"""

import numpy as np
import pytest
from dataclasses import replace

from raa_orbit_model.astrometry import (
    OBLIQUITY_J2000_DEG,
    earth_barycentric_position_au,
    parallax_factors,
)
from raa_orbit_model.experiments import perturbed_start, single_peak_schedule_for_response
from raa_orbit_model.fit import (
    ALL_PARAMETER_NAMES,
    ASTROMETRIC_PARAMETER_NAMES,
    fit_joint,
)
from raa_orbit_model.kepler import BinaryParams
from raa_orbit_model.model import AbsoluteAstrometryConfig, GaiaResponseConfig
from raa_orbit_model.scanning import schedule_from_arrays
from raa_orbit_model.synthetic import simulate_joint_data

YEAR = 2016.0 + np.linspace(0.0, 1.0, 2000)
SIGMA = 50.0


def test_earth_stays_near_one_astronomical_unit():
    r = np.linalg.norm(earth_barycentric_position_au(YEAR), axis=1)
    assert 0.98 < r.min() < 0.99
    assert 1.01 < r.max() < 1.02


def test_parallax_ellipse_is_circular_at_the_ecliptic_pole():
    p_alpha, p_delta = parallax_factors(YEAR, 270.0, 90.0 - OBLIQUITY_J2000_DEG)
    amplitude = np.hypot(p_alpha, p_delta)
    assert amplitude.min() == pytest.approx(1.0, abs=0.02)
    assert amplitude.max() == pytest.approx(1.0, abs=0.02)


def test_parallax_ellipse_degenerates_to_a_line_on_the_ecliptic():
    """On the ecliptic the ellipse collapses to a line tilted by the obliquity."""
    p_alpha, p_delta = parallax_factors(YEAR, 0.0, 0.0)
    ratio = p_delta / p_alpha
    assert ratio.std() < 1e-12
    assert ratio.mean() == pytest.approx(np.tan(np.deg2rad(OBLIQUITY_J2000_DEG)), rel=1e-6)


@pytest.mark.parametrize(
    "ra_deg,dec_deg,ecliptic_latitude_deg",
    [(270.0, 66.5607, 90.0), (0.0, 0.0, 0.0), (90.0, 23.4393, 0.0), (90.0, 66.5607, 43.0)],
)
def test_axis_ratio_tracks_sine_of_ecliptic_latitude(ra_deg, dec_deg, ecliptic_latitude_deg):
    p_alpha, p_delta = parallax_factors(YEAR, ra_deg, dec_deg)
    centred = np.vstack([p_alpha, p_delta])
    centred = centred - centred.mean(axis=1, keepdims=True)
    singular = np.linalg.svd(centred, compute_uv=False)
    measured = singular[1] / singular[0]
    expected = abs(np.sin(np.deg2rad(ecliptic_latitude_deg)))
    assert measured == pytest.approx(expected, abs=0.01)


def test_parallax_factor_amplitude_never_exceeds_unity_by_more_than_eccentricity():
    worst = 0.0
    for ra in range(0, 360, 45):
        for dec in (-60.0, -30.0, 0.0, 30.0, 60.0):
            p_alpha, p_delta = parallax_factors(YEAR, float(ra), dec)
            worst = max(worst, float(np.hypot(p_alpha, p_delta).max()))
    assert worst < 1.03


def _setup(period_yr=2.0):
    rng = np.random.default_rng(5)
    n = 180
    ra, dec, start = 120.0, 30.0, 2014.5
    schedule = schedule_from_arrays(
        np.sort(rng.uniform(0.0, 5.0, n)), rng.uniform(0.0, 360.0, n),
        ra_deg=ra, dec_deg=dec, mission_start_decimalyear=start,
    )
    astrometry = AbsoluteAstrometryConfig(
        ra_deg=ra, dec_deg=dec, mission_start_decimalyear=start
    )
    truth = BinaryParams(
        period_yr=period_yr, t_peri_yr=0.15, eccentricity=0.25, inclination_deg=72.0,
        omega_deg=55.0, node_deg=120.0, m1_msun=1.25, m2_msun=0.85,
        parallax_mas=20.0, gamma_kms=7.0, beta_g=0.25,
        pmra_mas_yr=12.0, pmdec_mas_yr=-8.0,
        delta_alpha_star_mas=3.0, delta_delta_mas=-2.0,
    )
    truth = replace(truth, parallax_mas=SIGMA / truth.a_rel_au)
    return truth, schedule, astrometry


def test_absolute_astrometry_is_off_unless_requested():
    """Frozen results depend on the Gaia channel carrying the orbit alone."""
    truth, schedule, astrometry = _setup()
    injected = GaiaResponseConfig("blended_gaussian_peak", SIGMA)
    selection = single_peak_schedule_for_response(truth, schedule, injected)
    without = simulate_joint_data(truth, injected, selection.schedule, seed=0,
                                  baseline_yr=schedule.mission_span_yr)
    with_astro = simulate_joint_data(truth, injected, selection.schedule, seed=0,
                                     baseline_yr=schedule.mission_span_yr,
                                     astrometry=astrometry)
    assert without.gaia_al.astrometry is None
    assert with_astro.gaia_al.astrometry is astrometry
    assert not np.allclose(without.gaia_al.values_mas, with_astro.gaia_al.values_mas)


def test_orbit_parallax_and_proper_motion_are_recovered_together():
    """The 15-parameter joint fit must recover all of them at once."""
    truth, schedule, astrometry = _setup()
    injected = GaiaResponseConfig("blended_gaussian_peak", SIGMA)
    selection = single_peak_schedule_for_response(truth, schedule, injected)
    data = simulate_joint_data(truth, injected, selection.schedule, seed=0,
                               baseline_yr=schedule.mission_span_yr, astrometry=astrometry)
    start = replace(perturbed_start(truth), pmra_mas_yr=0.0, pmdec_mas_yr=0.0,
                    delta_alpha_star_mas=0.0, delta_delta_mas=0.0)
    result = fit_joint(
        data, start,
        GaiaResponseConfig("blended_gaussian_peak", SIGMA, allow_multi_peak_continuation=True),
        free_names=tuple(ALL_PARAMETER_NAMES) + tuple(ASTROMETRIC_PARAMETER_NAMES),
    )
    assert result.success
    assert result.n_free == 15
    assert result.reduced_chi2 == pytest.approx(1.0, abs=0.25)
    assert result.params.parallax_mas == pytest.approx(truth.parallax_mas, rel=1e-3)
    assert result.params.pmra_mas_yr == pytest.approx(truth.pmra_mas_yr, abs=0.05)
    assert result.params.pmdec_mas_yr == pytest.approx(truth.pmdec_mas_yr, abs=0.05)


def test_photocentre_misspecification_biases_the_fitted_parallax():
    """The mechanism the orbit-only model could not represent."""
    truth, schedule, astrometry = _setup(period_yr=2.0)
    injected = GaiaResponseConfig("blended_gaussian_peak", SIGMA)
    selection = single_peak_schedule_for_response(truth, schedule, injected)
    data = simulate_joint_data(truth, injected, selection.schedule, seed=0,
                               baseline_yr=schedule.mission_span_yr, astrometry=astrometry)
    start = replace(perturbed_start(truth), pmra_mas_yr=0.0, pmdec_mas_yr=0.0,
                    delta_alpha_star_mas=0.0, delta_delta_mas=0.0)
    free = tuple(ALL_PARAMETER_NAMES) + tuple(ASTROMETRIC_PARAMETER_NAMES)
    photocentre = fit_joint(data, start, GaiaResponseConfig("photocentre"), free_names=free)
    fractional_bias = (
        photocentre.params.parallax_mas - truth.parallax_mas
    ) / truth.parallax_mas
    assert abs(fractional_bias) > 1e-3
