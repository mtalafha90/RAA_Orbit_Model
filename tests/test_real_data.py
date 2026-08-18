from pathlib import Path
import numpy as np

from raa_orbit_model.real_data import fit_visual_sb2, parse_legacy_binary_csv

DATA = Path(__file__).parents[1] / "data" / "real" / "gj7652" / "GL765_Test1.csv"


def test_gl765_parser_and_parallax_alias():
    d = parse_legacy_binary_csv(DATA)
    assert d.object_name == "GL765.2"
    assert d.metadata["parallax"] == "54.27"
    assert len(d.visual_time_yr) == 11
    assert len(d.rv1_time_yr) == 44
    assert len(d.rv2_time_yr) == 44
    assert d.n_scalar_constraints == 110


def test_gl765_v6a_regression():
    fit = fit_visual_sb2(parse_legacy_binary_csv(DATA))
    assert fit.success
    assert fit.dof == 100
    assert np.isclose(fit.chi2, 104.628951586, rtol=0, atol=2e-6)
    assert np.isclose(fit.reduced_chi2, 1.046289516, rtol=0, atol=2e-8)
    assert np.isclose(fit.params.m1_msun + fit.params.m2_msun, 1.5896874326, atol=2e-7)
    assert np.isclose(fit.params.parallax_mas, 35.44258197, atol=2e-6)
    assert "beta_g" not in fit.free_names


def test_gl765_fixed_parallax_controls():
    d = parse_legacy_binary_csv(DATA)
    header = fit_visual_sb2(d, fixed_parallax_mas=54.27)
    balega = fit_visual_sb2(d, fixed_parallax_mas=31.0)
    assert np.isclose(header.chi2, 168.94353, atol=2e-3)
    assert np.isclose(balega.chi2, 108.70277, atol=2e-3)
    assert header.chi2 > balega.chi2
