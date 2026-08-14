"""Score this project's surrogate against the published gaiamock response.

El-Badry et al. (2024) released ``al_bias_binary``, which places the Gaia
along-scan coordinate at the peak of the combined along-scan flux profile.
These tests establish that this project's surrogate is the same model, so it
can be presented as reproducing a published response rather than as an
invention of this project.
"""

import numpy as np
import pytest

from raa_orbit_model.gaia import blended_gaussian_response, critical_blended_separation_sigma
from raa_orbit_model.gaiamock_reference import (
    al_bias_binary,
    light_ratio_from_beta,
    mass_ratio_from_fraction,
)

U = 90.0  # gaiamock's effective angular resolution, in mas
B = 0.4   # secondary mass fraction


def _both(separation_mas, beta_g):
    raa = blended_gaussian_response(separation_mas, B, beta_g, U)
    published = al_bias_binary(
        separation_mas, mass_ratio_from_fraction(B), light_ratio_from_beta(beta_g), U
    )
    return raa, published


def test_ratio_conversions_round_trip():
    assert light_ratio_from_beta(0.25) == pytest.approx(1.0 / 3.0)
    assert mass_ratio_from_fraction(0.5) == pytest.approx(1.0)


@pytest.mark.parametrize("beta_g", [0.10, 0.25, 0.45])
@pytest.mark.parametrize("ratio", [0.2, 0.5, 1.0, 1.5, 2.0])
def test_surrogate_reproduces_published_response_in_the_single_peak_regime(beta_g, ratio):
    """The decisive check: same model, to the published solver's own tolerance."""
    raa, published = _both(ratio * U, beta_g)
    assert raa.n_peaks == 1
    assert raa.al_mas == pytest.approx(published, rel=1e-4, abs=1e-3)


@pytest.mark.parametrize("beta_g", [0.10, 0.25, 0.45])
def test_published_linearisation_below_a_tenth_of_a_resolution_element(beta_g):
    """Below 0.1u gaiamock deliberately linearises; we solve exactly and agree closely."""
    raa, published = _both(0.05 * U, beta_g)
    assert raa.al_mas == pytest.approx(published, rel=2e-3)


@pytest.mark.parametrize("beta_g", [0.10, 0.25, 0.45])
def test_exact_boundary_is_wider_than_the_published_piecewise_cut(beta_g):
    """gaiamock switches to the primary at (3-f)u, before the profile truly splits.

    The exact equal-width criterion places mode splitting later, so the
    published cut fires while the profile is still single-peaked.
    """
    f = light_ratio_from_beta(beta_g)
    assert critical_blended_separation_sigma(beta_g) > 3.0 - f


def test_published_cut_costs_accuracy_where_the_profile_is_still_single_peaked():
    """At beta=0.10, d=2.9u is past gaiamock's cut but still genuinely single-peaked."""
    beta_g = 0.10
    raa, published = _both(2.9 * U, beta_g)
    assert raa.n_peaks == 1                      # exact criterion: still one mode
    assert 2.9 > 3.0 - light_ratio_from_beta(beta_g)  # but past the published cut
    assert abs(raa.al_mas - published) > 0.1     # so the published value is displaced


def test_surrogate_refuses_where_the_published_model_asserts_the_primary():
    """Beyond mode splitting we return NaN; gaiamock returns the primary position.

    This is a deliberate difference in policy, not a disagreement about physics,
    and it must be stated whenever the two are compared.
    """
    raa, published = _both(2.9 * U, 0.45)
    assert raa.n_peaks == 2
    assert np.isnan(raa.al_mas)
    assert np.isfinite(published)


def test_photocentre_limit_agrees_with_the_published_small_separation_branch():
    beta_g, d = 0.25, 1e-3 * U
    raa, published = _both(d, beta_g)
    photocentre = (beta_g - B) * d
    assert raa.al_mas == pytest.approx(photocentre, rel=1e-6)
    assert published == pytest.approx(photocentre, rel=1e-12)
