"""The orbit and noise axes that methodology.md section 12 names but the
published experiment never varied."""

import argparse

import pytest

from raa_orbit_model.experiment_config import (
    DEFAULT_NOISE,
    DEFAULT_ORBIT,
    add_noise_arguments,
    add_orbit_arguments,
    noise_kwargs_from_args,
    truth_from_args,
)


def _parse(argv):
    parser = argparse.ArgumentParser()
    add_orbit_arguments(parser)
    add_noise_arguments(parser)
    return parser.parse_args(argv)


def test_defaults_reproduce_the_frozen_orbit():
    """Adding these switches must not change any existing experiment."""
    truth = truth_from_args(_parse([]))
    for name, expected in DEFAULT_ORBIT.items():
        assert getattr(truth, name) == pytest.approx(expected)


def test_defaults_reproduce_the_frozen_noise_prescription():
    assert noise_kwargs_from_args(_parse([])) == DEFAULT_NOISE


def test_mass_ratio_sets_the_secondary_mass():
    truth = truth_from_args(_parse(["--m1-msun", "1.25", "--mass-ratio", "0.4"]))
    assert truth.m2_msun == pytest.approx(0.5)
    assert truth.q == pytest.approx(0.4)


def test_mass_ratio_overrides_an_explicit_secondary_mass():
    truth = truth_from_args(_parse(["--m2-msun", "0.85", "--mass-ratio", "0.2"]))
    assert truth.m2_msun == pytest.approx(0.25)


def test_period_and_eccentricity_are_configurable():
    truth = truth_from_args(_parse(["--period-yr", "1.0", "--eccentricity", "0.6"]))
    assert truth.period_yr == pytest.approx(1.0)
    assert truth.eccentricity == pytest.approx(0.6)


@pytest.mark.parametrize("argv", [
    ["--mass-ratio", "1.8"],       # secondary heavier than primary
    ["--mass-ratio", "0.0"],
    ["--eccentricity", "1.2"],
    ["--period-yr", "-1.0"],
    ["--m1-msun", "0.0"],
    ["--m2-msun", "2.0"],          # heavier than the default primary
])
def test_unphysical_orbits_are_rejected(argv):
    with pytest.raises(ValueError):
        truth_from_args(_parse(argv))


@pytest.mark.parametrize("argv", [
    ["--gaia-sigma-mas", "0.0"],
    ["--ast-sigma-mas", "-1.0"],
    ["--rv-sigma-kms", "0.0"],
    ["--n-ast", "-1"],
])
def test_unphysical_noise_settings_are_rejected(argv):
    with pytest.raises(ValueError):
        noise_kwargs_from_args(_parse(argv))


def test_signal_to_noise_axis_is_reachable():
    kwargs = noise_kwargs_from_args(_parse(["--gaia-sigma-mas", "0.02", "--n-rv", "12"]))
    assert kwargs["gaia_sigma_mas"] == pytest.approx(0.02)
    assert kwargs["n_rv"] == 12
