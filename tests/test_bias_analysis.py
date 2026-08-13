import pandas as pd
import pytest

from raa_orbit_model.bias_analysis import (
    build_paired_results,
    load_bias_results,
    summarize_paired_results,
)


def _row(model, seed, chi2, bias_beta, bias_plx, bias_m1, bias_m2, bias_i):
    return {
        "true_beta_g": 0.20,
        "a_over_sigma": 0.50,
        "seed": seed,
        "model": model,
        "chi2": chi2,
        "reduced_chi2": chi2 / 100.0,
        "bias_beta_g": bias_beta,
        "bias_parallax_mas": bias_plx,
        "bias_m1_msun": bias_m1,
        "bias_m2_msun": bias_m2,
        "bias_inclination_deg": bias_i,
        "true_parallax_mas": 20.0,
        "true_m1_msun": 1.25,
        "true_m2_msun": 0.85,
    }


def test_pairing_uses_same_seed_and_photo_minus_raa_sign(tmp_path):
    rows = [
        _row("photocentre", 0, 130.0, -0.02, 1.0, -0.05, -0.03, 1.0),
        _row("resolution_aware", 0, 100.0, 0.001, 0.1, 0.001, 0.001, 0.1),
        _row("photocentre", 1, 150.0, -0.04, 2.0, -0.10, -0.06, 2.0),
        _row("resolution_aware", 1, 105.0, -0.001, -0.1, -0.001, -0.001, -0.1),
    ]
    path = tmp_path / "bias.csv"
    pd.DataFrame(rows).to_csv(path, index=False)

    df = load_bias_results(path)
    paired = build_paired_results(df)

    assert list(paired["seed"]) == [0, 1]
    assert list(paired["delta_chi2"]) == pytest.approx([30.0, 45.0])
    assert list(paired["pct_bias_beta_g_photocentre"]) == pytest.approx([-10.0, -20.0])
    assert list(paired["pct_bias_parallax_mas_photocentre"]) == pytest.approx([5.0, 10.0])


def test_summary_reports_seed_distribution(tmp_path):
    rows = [
        _row("photocentre", 0, 130.0, -0.02, 1.0, -0.05, -0.03, 1.0),
        _row("resolution_aware", 0, 100.0, 0.001, 0.1, 0.001, 0.001, 0.1),
        _row("photocentre", 1, 150.0, -0.04, 2.0, -0.10, -0.06, 2.0),
        _row("resolution_aware", 1, 105.0, -0.001, -0.1, -0.001, -0.001, -0.1),
    ]
    path = tmp_path / "bias.csv"
    pd.DataFrame(rows).to_csv(path, index=False)

    summary = summarize_paired_results(build_paired_results(load_bias_results(path)))
    row = summary.iloc[0]

    assert row["n_pairs"] == 2
    assert row["delta_chi2_median"] == pytest.approx(37.5)
    assert row["delta_chi2_min"] == pytest.approx(30.0)
    assert row["delta_chi2_max"] == pytest.approx(45.0)


def test_duplicate_model_record_is_rejected(tmp_path):
    row = _row("photocentre", 0, 130.0, -0.02, 1.0, -0.05, -0.03, 1.0)
    rows = [
        row,
        dict(row),
        _row("resolution_aware", 0, 100.0, 0.001, 0.1, 0.001, 0.001, 0.1),
    ]
    path = tmp_path / "bias.csv"
    pd.DataFrame(rows).to_csv(path, index=False)

    with pytest.raises(ValueError, match="duplicate model records"):
        load_bias_results(path)
