import math

import numpy as np
import pytest

from raa_orbit_model.gaia import blended_gaussian_response, project_along_scan
from raa_orbit_model.kepler import BinaryParams
from raa_orbit_model.model import GaiaResponseConfig, predict_gaia_orbital_response
from raa_orbit_model.penoyre import (
    penoyre_effective_width,
    penoyre_oriented_gaussian_peak,
    penoyre_oriented_gaussian_response,
)


def test_effective_width_has_the_penoyre_axis_limits():
    alpha, beta = 20.0, 60.0
    assert penoyre_effective_width(alpha, beta, 0.0) == pytest.approx(alpha)
    assert penoyre_effective_width(alpha, beta, math.pi) == pytest.approx(alpha)
    assert penoyre_effective_width(alpha, beta, math.pi / 2.0) == pytest.approx(beta)
    assert penoyre_effective_width(alpha, beta, -math.pi / 2.0) == pytest.approx(beta)


def test_effective_width_is_pi_periodic():
    phi = np.linspace(-1.2, 1.2, 17)
    g1 = penoyre_effective_width(20.0, 60.0, phi)
    g2 = penoyre_effective_width(20.0, 60.0, phi + math.pi)
    np.testing.assert_allclose(g1, g2, rtol=0.0, atol=1e-12)


def test_invalid_axis_widths_are_rejected():
    with pytest.raises(ValueError):
        penoyre_effective_width(0.0, 60.0, 0.0)
    with pytest.raises(ValueError):
        penoyre_effective_width(20.0, 10.0, 0.0)


def test_infinite_across_scan_width_exactly_recovers_the_1d_response():
    east = np.array([5.0, 20.0, -12.0, 30.0])
    north = np.array([20.0, -8.0, 25.0, 5.0])
    scan = np.array([0.0, 37.0, 121.0, 250.0])
    B, beta_g, alpha = 0.4, 0.25, 50.0

    d_al = project_along_scan(east, north, scan)
    reference = blended_gaussian_response(d_al, B, beta_g, alpha)
    oriented = penoyre_oriented_gaussian_response(
        east, north, scan, B, beta_g, alpha, math.inf
    )

    np.testing.assert_allclose(oriented.al_mas, reference.al_mas, rtol=0.0, atol=1e-12)
    np.testing.assert_array_equal(oriented.n_peaks, reference.n_peaks)
    np.testing.assert_allclose(
        oriented.separation_sigma,
        np.abs(d_al) / alpha,
        rtol=0.0,
        atol=1e-12,
    )


def test_finite_across_scan_width_changes_the_response_at_oblique_orientation():
    # Pair points North.  A 60-degree scan is oblique enough that finite AC
    # width changes the blend while retaining a non-zero AL projection.
    east, north, scan = 0.0, 40.0, 60.0
    B, beta_g, alpha = 0.4, 0.25, 50.0
    one_d = penoyre_oriented_gaussian_response(
        east, north, scan, B, beta_g, alpha, math.inf
    )
    finite = penoyre_oriented_gaussian_response(
        east, north, scan, B, beta_g, alpha, 3.0 * alpha
    )
    assert finite.n_peaks == 1
    assert finite.effective_width_mas > alpha
    assert finite.al_mas != pytest.approx(one_d.al_mas, abs=1e-8)


def test_mode_splitting_depends_on_full_separation_and_orientation():
    # Equal light has the exact critical separation rho=2.  The same 25 mas
    # pair is multi-peak when aligned with alpha=10 mas, but single-peak when
    # perpendicular and therefore seen with beta=30 mas.
    aligned = penoyre_oriented_gaussian_response(
        0.0, 25.0, 0.0, 0.4, 0.5, 10.0, 30.0
    )
    perpendicular = penoyre_oriented_gaussian_response(
        0.0, 25.0, 90.0, 0.4, 0.5, 10.0, 30.0
    )
    assert aligned.separation_sigma == pytest.approx(2.5)
    assert aligned.n_peaks == 2
    assert np.isnan(aligned.al_mas)
    assert perpendicular.separation_sigma == pytest.approx(25.0 / 30.0)
    assert perpendicular.n_peaks == 1
    assert np.isfinite(perpendicular.al_mas)


def test_optimizer_continuation_is_finite_but_scientific_response_is_flagged():
    scientific = penoyre_oriented_gaussian_response(
        0.0, 25.0, 0.0, 0.4, 0.5, 10.0, 30.0
    )
    continuation = penoyre_oriented_gaussian_peak(
        0.0, 25.0, 0.0, 0.4, 0.5, 10.0, 30.0,
        allow_multi_peak_continuation=True,
    )
    assert scientific.n_peaks == 2
    assert np.isnan(scientific.al_mas)
    assert np.isfinite(continuation)


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
        beta_g=0.25,
    )


def test_joint_model_penoyre_mode_reduces_to_existing_m1_in_infinite_ac_limit():
    times = np.linspace(0.0, 4.8, 41)
    scan = np.linspace(3.0, 347.0, len(times))
    params = _truth()
    alpha = 50.0

    m1 = predict_gaia_orbital_response(
        times,
        scan,
        params,
        GaiaResponseConfig("blended_gaussian_peak", alpha),
    )
    m2_limit = predict_gaia_orbital_response(
        times,
        scan,
        params,
        GaiaResponseConfig(
            "penoyre_gaussian_peak",
            alpha,
            sigma_ac_mas=math.inf,
        ),
    )
    np.testing.assert_allclose(m2_limit.values_mas, m1.values_mas, rtol=0.0, atol=1e-12,
                               equal_nan=True)
    np.testing.assert_array_equal(m2_limit.n_peaks, m1.n_peaks)


def test_joint_model_requires_both_penoyre_axis_widths():
    times = np.array([0.0, 0.5])
    scan = np.array([0.0, 90.0])
    with pytest.raises(ValueError, match="sigma_ac_mas"):
        predict_gaia_orbital_response(
            times,
            scan,
            _truth(),
            GaiaResponseConfig("penoyre_gaussian_peak", 50.0),
        )
