"""Robustness of the photocentre-versus-RAA comparison to a wrong measurement model.

The headline full-sky experiment fits the resolution-aware (RAA) model using
exactly the surrogate that generated the data, with the true response width
supplied. That makes the RAA hypothesis correct by construction and the
photocentre hypothesis the only misspecified one, so the measured
``Delta chi2`` is an upper bound on the advantage available in practice.

This module removes that advantage in two controlled ways.

1. **Width misspecification.** The response width used by the fit differs from
   the injected width, or is fitted from the data rather than asserted.
2. **Shape misspecification.** The data are injected with two Gaussians of
   *different* widths, which is outside the equal-width family the fit uses.
   This is a generic controlled proxy for unmodelled component-dependent or
   profile-shape differences. It is not a calibrated Gaia colour-to-PLSF model.

Under (2) both hypotheses are wrong, which is the situation that will hold
against real Gaia data. A resolution-aware model that only wins when it is
exactly correct has not demonstrated anything useful, so this is an important
control. The Penoyre-style orientation-dependent hierarchy in
``response_fidelity.py`` is now the more literature-grounded response-fidelity
test.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from .experiments import perturbed_start, single_peak_schedule_for_response
from .fit import ALL_PARAMETER_NAMES, fit_joint
from .kepler import BinaryParams
from .model import GaiaResponseConfig
from .scanning import GaiaScanSchedule
from .synthetic import simulate_joint_data


def compare_under_misspecification(
    truth: BinaryParams,
    gaia_schedule: GaiaScanSchedule,
    *,
    sigma_response_mas: float,
    secondary_width_ratio: float = 1.0,
    fit_width_ratio: float = 1.0,
    fit_sigma: bool = False,
    retain_multi_peak: bool = False,
    seed: int = 0,
    n_ast: int = 24,
    n_rv: int = 48,
    ast_sigma_mas: float = 0.20,
    rv_sigma_kms: float = 0.10,
    gaia_sigma_mas: float = 0.10,
) -> dict:
    """Inject with one response and fit with another, then compare hypotheses.

    ``secondary_width_ratio`` scales the injected secondary width relative to
    the primary; 1.0 reproduces the equal-width surrogate used for the frozen
    results. ``fit_width_ratio`` scales the width the RAA fit assumes relative
    to the injected primary width. ``fit_sigma=True`` frees that width instead
    of asserting it.
    """
    if sigma_response_mas <= 0:
        raise ValueError("sigma_response_mas must be > 0")
    if secondary_width_ratio <= 0 or fit_width_ratio <= 0:
        raise ValueError("width ratios must be > 0")

    sigma_secondary = (
        None if secondary_width_ratio == 1.0
        else float(secondary_width_ratio) * sigma_response_mas
    )
    injected = GaiaResponseConfig(
        "blended_gaussian_peak",
        sigma_response_mas,
        allow_multi_peak_continuation=retain_multi_peak,
        sigma_secondary_mas=sigma_secondary,
    )
    selection = single_peak_schedule_for_response(
        truth, gaia_schedule, injected, retain_multi_peak=retain_multi_peak
    )
    data = simulate_joint_data(
        truth,
        injected,
        selection.schedule,
        seed=seed,
        n_ast=n_ast,
        n_rv=n_rv,
        baseline_yr=gaia_schedule.mission_span_yr,
        ast_sigma_mas=ast_sigma_mas,
        rv_sigma_kms=rv_sigma_kms,
        gaia_sigma_mas=gaia_sigma_mas,
    )
    initial = perturbed_start(truth)

    photo = fit_joint(data, initial, GaiaResponseConfig("photocentre"))

    # The RAA hypothesis here uses the equal-width family. Only the injection
    # leaves that family; the newer orientation-dependent hierarchy is separate.
    raa_response = GaiaResponseConfig(
        "blended_gaussian_peak",
        float(fit_width_ratio) * sigma_response_mas,
        allow_multi_peak_continuation=True,
    )
    free = tuple(ALL_PARAMETER_NAMES) + (("sigma_response_mas",) if fit_sigma else ())
    raa = fit_joint(data, initial, raa_response, free_names=free)

    row = {
        "seed": int(seed),
        "true_beta_g": float(truth.beta_g),
        "a_over_sigma": float(truth.a_rel_au * truth.parallax_mas / sigma_response_mas),
        "sigma_response_mas": float(sigma_response_mas),
        "secondary_width_ratio": float(secondary_width_ratio),
        "fit_width_ratio": float(fit_width_ratio),
        "fit_sigma": bool(fit_sigma),
        "retain_multi_peak": bool(retain_multi_peak),
        "gaia_n_transits": selection.n_total,
        "gaia_n_multi_peak_flagged": selection.n_multi_peak,
        "chi2_photocentre": photo.chi2,
        "chi2_resolution_aware": raa.chi2,
        "delta_chi2": photo.chi2 - raa.chi2,
        "raa_favoured": bool(photo.chi2 - raa.chi2 > 0.0),
        "photocentre_success": bool(photo.success),
        "raa_success": bool(raa.success),
        "fitted_sigma_response_mas": raa.fitted_sigma_response_mas,
    }
    for name in ("beta_g", "parallax_mas", "m1_msun", "m2_msun", "inclination_deg"):
        true_value = float(getattr(truth, name))
        for label, result in (("photocentre", photo), ("resolution_aware", raa)):
            fitted = float(getattr(result.params, name))
            row[f"bias_{name}_{label}"] = fitted - true_value
            if true_value != 0:
                row[f"frac_bias_{name}_{label}"] = (fitted - true_value) / true_value
    return row


def shape_misspecification_scan(
    base_truth: BinaryParams,
    gaia_schedule: GaiaScanSchedule,
    *,
    secondary_width_ratios=(1.0, 1.1, 1.2, 1.3, 1.4, 1.6),
    a_over_sigma_values=(1.0,),
    beta_values=(0.25,),
    seeds=(0, 1, 2, 3, 4),
    sigma_response_mas: float = 50.0,
    fit_sigma: bool = False,
    retain_multi_peak: bool = False,
    **kwargs,
) -> list[dict]:
    """Scan the equal-width RAA fit against injected profile-shape error."""
    rows: list[dict] = []
    for beta in beta_values:
        if not (0.0 <= beta <= 0.5):
            raise ValueError("beta_g must satisfy 0 <= beta_g <= 0.5")
        for ratio in a_over_sigma_values:
            if ratio <= 0:
                raise ValueError("a_over_sigma_values must be > 0")
            parallax = float(ratio) * sigma_response_mas / base_truth.a_rel_au
            truth = replace(base_truth, parallax_mas=parallax, beta_g=float(beta))
            for width_ratio in secondary_width_ratios:
                for seed in seeds:
                    rows.append(compare_under_misspecification(
                        truth,
                        gaia_schedule,
                        sigma_response_mas=sigma_response_mas,
                        secondary_width_ratio=float(width_ratio),
                        fit_sigma=fit_sigma,
                        retain_multi_peak=retain_multi_peak,
                        seed=int(seed),
                        **kwargs,
                    ))
    return rows


def summarise_by(rows: list[dict], key: str) -> list[dict]:
    """Median summary of the paired discrepancy grouped by one scan variable."""
    if not rows:
        raise ValueError("rows is empty")
    out = []
    for value in sorted({row[key] for row in rows}):
        group = [row for row in rows if row[key] == value]
        delta = np.array([row["delta_chi2"] for row in group], dtype=float)
        out.append({
            key: value,
            "n": len(group),
            "delta_chi2_median": float(np.median(delta)),
            "delta_chi2_q16": float(np.percentile(delta, 16)),
            "delta_chi2_q84": float(np.percentile(delta, 84)),
            "fraction_raa_favoured": float(np.mean(delta > 0.0)),
            "median_multi_peak_flagged": float(np.median(
                [row["gaia_n_multi_peak_flagged"] for row in group]
            )),
        })
    return out
