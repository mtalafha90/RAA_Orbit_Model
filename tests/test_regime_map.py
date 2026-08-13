import numpy as np
import pandas as pd
import pytest

from raa_orbit_model.bias_analysis import (
    build_regime_boundary_table,
    make_regime_map,
    write_regime_boundary_table,
)


def _summary_grid():
    rows = []
    values = {
        0.10: {
            0.4: (-0.2, -1.0, -2.0),
            0.5: (4.0, 0.5, -0.2),
            0.6: (12.0, 3.0, 0.4),
        },
        0.20: {
            0.4: (1.0, -0.2, -1.0),
            0.5: (15.0, 2.0, 0.1),
            0.6: (30.0, 10.0, 4.0),
        },
    }
    for beta, ratios in values.items():
        for ratio, (median, q16, minimum) in ratios.items():
            rows.append(
                {
                    "true_beta_g": beta,
                    "a_over_sigma": ratio,
                    "n_pairs": 20,
                    "delta_chi2_median": median,
                    "delta_chi2_mean": median,
                    "delta_chi2_q16": q16,
                    "delta_chi2_q84": median + 1.0,
                    "delta_chi2_min": minimum,
                    "delta_chi2_max": median + 2.0,
                }
            )
    return pd.DataFrame(rows)


def test_regime_boundary_table_uses_first_sampled_crossing():
    summary = _summary_grid()
    raw = pd.DataFrame(
        {
            "true_beta_g": [0.10, 0.10, 0.20, 0.20],
            "a_over_sigma": [0.5, 0.6, 0.5, 0.6],
            "gaia_multi_peak_fraction": [0.0, 0.1, 0.0, 0.0],
        }
    )
    boundary = build_regime_boundary_table(summary, raw=raw).set_index("true_beta_g")
    assert boundary.loc[0.10, "q16_positive_a_over_sigma"] == pytest.approx(0.5)
    assert boundary.loc[0.10, "all_positive_a_over_sigma"] == pytest.approx(0.6)
    assert boundary.loc[0.10, "median_gt10_a_over_sigma"] == pytest.approx(0.6)
    assert boundary.loc[0.10, "multi_peak_onset_a_over_sigma"] == pytest.approx(0.6)
    assert boundary.loc[0.20, "q16_positive_a_over_sigma"] == pytest.approx(0.5)
    assert boundary.loc[0.20, "all_positive_a_over_sigma"] == pytest.approx(0.5)
    assert boundary.loc[0.20, "median_gt10_a_over_sigma"] == pytest.approx(0.5)
    assert np.isnan(boundary.loc[0.20, "multi_peak_onset_a_over_sigma"])


def test_regime_map_rejects_incomplete_grid(tmp_path):
    with pytest.raises(ValueError, match="complete rectangular"):
        make_regime_map(_summary_grid().iloc[:-1].copy(), tmp_path)


def test_regime_outputs_are_written(tmp_path):
    import matplotlib
    matplotlib.use("Agg")
    summary = _summary_grid()
    boundary = build_regime_boundary_table(summary)
    csv_path = write_regime_boundary_table(boundary, tmp_path, prefix="smoke")
    outputs = make_regime_map(summary, tmp_path, prefix="smoke", dpi=80)
    assert csv_path.exists() and csv_path.stat().st_size > 0
    assert len(outputs) == 2
    assert {p.suffix for p in outputs} == {".png", ".pdf"}
    assert all(p.exists() and p.stat().st_size > 0 for p in outputs)
