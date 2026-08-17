import math

import numpy as np
import pytest

from raa_orbit_model.kepler import BinaryParams
from raa_orbit_model.response_fidelity import (
    MODEL_ORDER,
    response_fidelity_once,
    response_fidelity_scan,
    summarise_response_fidelity,
)
from raa_orbit_model.scanning import schedule_from_arrays


ALPHA = 50.0


def _truth(a_over_alpha=0.8, beta_g=0.25):
    base = BinaryParams(
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
        beta_g=beta_g,
    )
    return BinaryParams(
        **{
            **base.__dict__,
            "parallax_mas": a_over_alpha * ALPHA / base.a_rel_au,
        }
    )


def _schedule(n=48):
    times = np.linspace(0.03, 5.20, n)
    # Deterministic directional coverage with no 180-degree folding.
    angles = np.mod(np.linspace(7.0, 997.0, n), 360.0)
    return schedule_from_arrays(
        times,
        angles,
        ra_deg=120.0,
        dec_deg=30.0,
        release="dr4",
        source="response-fidelity regression",
    )


def test_response_fidelity_once_fits_one_paired_realization_with_all_three_models():
    rows = response_fidelity_once(
        _truth(),
        _schedule(),
        alpha_mas=ALPHA,
        beta_over_alpha=3.0,
        seed=2,
        n_ast=16,
        n_rv=28,
    )
    assert [row["model"] for row in rows] == list(MODEL_ORDER)
    assert all(row["seed"] == 2 for row in rows)
    assert all(row["gaia_n_multi_peak_injection"] == 0 for row in rows)
    assert all(row["success"] for row in rows)
    assert rows[2]["chi2"] < rows[0]["chi2"]


def test_infinite_ac_limit_makes_m1_and_m2_identical_inside_the_hierarchy():
    rows = response_fidelity_once(
        _truth(a_over_alpha=0.6),
        _schedule(),
        alpha_mas=ALPHA,
        beta_over_alpha=math.inf,
        seed=1,
        n_ast=12,
        n_rv=20,
    )
    by_model = {row["model"]: row for row in rows}
    m1 = by_model["equal_width"]
    m2 = by_model["penoyre_oriented"]
    assert m1["chi2"] == pytest.approx(m2["chi2"], rel=0.0, abs=1e-8)
    for name in ("m1_msun", "m2_msun", "parallax_mas", "beta_g", "inclination_deg"):
        assert m1[f"fit_{name}"] == pytest.approx(m2[f"fit_{name}"], rel=0.0, abs=1e-8)


def test_response_fidelity_scan_and_summary_have_expected_shape():
    rows = response_fidelity_scan(
        _truth(),
        _schedule(32),
        a_over_alpha_values=(0.5, 0.8),
        beta_values=(0.25,),
        beta_over_alpha_values=(1.5, 3.0),
        seeds=(0,),
        alpha_mas=ALPHA,
        n_ast=8,
        n_rv=12,
    )
    assert len(rows) == 2 * 1 * 2 * 1 * 3
    summary = summarise_response_fidelity(rows)
    assert len(summary) == 2 * 1 * 2 * 3
    assert all(item["n"] == 1 for item in summary)
    assert all(0.0 <= item["success_fraction"] <= 1.0 for item in summary)


def test_invalid_elongation_is_rejected():
    with pytest.raises(ValueError, match="beta_over_alpha"):
        response_fidelity_once(
            _truth(), _schedule(), alpha_mas=ALPHA, beta_over_alpha=0.8
        )
