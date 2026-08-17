"""Response-fidelity experiments for marginal-resolution orbit inference.

The frozen baseline experiments injected and fitted the same equal-width 1-D
response.  This module deliberately breaks that matched-model assumption.
Synthetic Gaia-like measurements are generated with the finite-elongation
Penoyre-style surrogate (M2) and the *same realization* is then fitted with:

M0  ordinary unresolved photocentre;
M1  equal-width 1-D Lindegren/gaiamock-family peak response;
M2  finite-elongation Penoyre-style response used for injection.

The primary outputs are physical-parameter biases.  Delta-chi2 remains a
secondary diagnostic.  Posterior coverage is a later stage once the M2 forward
model and deterministic hierarchy have passed regression tests.
"""

from __future__ import annotations

from dataclasses import asdict, replace
import math

import numpy as np

from .experiments import perturbed_start, single_peak_schedule_for_response
from .fit import ALL_PARAMETER_NAMES, fit_joint
from .kepler import BinaryParams
from .model import GaiaResponseConfig, predict_gaia_orbital_response
from .scanning import GaiaScanSchedule
from .synthetic import simulate_joint_data


MODEL_ORDER = ("photocentre", "equal_width", "penoyre_oriented")


def _validate_elongation(beta_over_alpha: float) -> float:
    value = float(beta_over_alpha)
    if math.isnan(value) or value < 1.0:
        raise ValueError("beta_over_alpha must be >= 1 or +inf")
    return value


def _science_response(model: str, alpha_mas: float, beta_mas: float) -> GaiaResponseConfig:
    if model == "photocentre":
        return GaiaResponseConfig("photocentre")
    if model == "equal_width":
        return GaiaResponseConfig("blended_gaussian_peak", alpha_mas)
    if model == "penoyre_oriented":
        return GaiaResponseConfig(
            "penoyre_gaussian_peak",
            alpha_mas,
            sigma_ac_mas=beta_mas,
        )
    raise KeyError(model)


def _fit_response(model: str, alpha_mas: float, beta_mas: float) -> GaiaResponseConfig:
    response = _science_response(model, alpha_mas, beta_mas)
    if model == "photocentre":
        return response
    return replace(response, allow_multi_peak_continuation=True)


def response_fidelity_once(
    truth: BinaryParams,
    gaia_schedule: GaiaScanSchedule,
    *,
    alpha_mas: float,
    beta_over_alpha: float = 3.0,
    seed: int = 0,
    free_names=ALL_PARAMETER_NAMES,
    n_ast: int = 24,
    n_rv: int = 48,
    ast_sigma_mas: float = 0.20,
    rv_sigma_kms: float = 0.10,
    gaia_sigma_mas: float = 0.10,
) -> list[dict]:
    """Inject M2 and fit M0, M1, and M2 to one paired realization."""
    if alpha_mas <= 0 or not math.isfinite(float(alpha_mas)):
        raise ValueError("alpha_mas must be finite and > 0")
    elongation = _validate_elongation(beta_over_alpha)
    beta_mas = math.inf if math.isinf(elongation) else float(alpha_mas) * elongation

    injected = GaiaResponseConfig(
        "penoyre_gaussian_peak",
        float(alpha_mas),
        sigma_ac_mas=beta_mas,
    )
    selection = single_peak_schedule_for_response(truth, gaia_schedule, injected)
    data = simulate_joint_data(
        truth,
        injected,
        selection.schedule,
        seed=int(seed),
        n_ast=n_ast,
        n_rv=n_rv,
        baseline_yr=gaia_schedule.mission_span_yr,
        ast_sigma_mas=ast_sigma_mas,
        rv_sigma_kms=rv_sigma_kms,
        gaia_sigma_mas=gaia_sigma_mas,
    )
    initial = perturbed_start(truth)

    rows: list[dict] = []
    for model in MODEL_ORDER:
        result = fit_joint(
            data,
            initial,
            _fit_response(model, float(alpha_mas), beta_mas),
            free_names=free_names,
        )
        final_multi = 0
        scientific_valid = bool(result.success)
        if model != "photocentre":
            final_prediction = predict_gaia_orbital_response(
                data.gaia_al.times_yr,
                data.gaia_al.scan_angle_deg,
                result.params,
                _science_response(model, float(alpha_mas), beta_mas),
            )
            final_multi = int(final_prediction.n_multi_peak)
            scientific_valid = bool(result.success and final_multi == 0)

        row = {
            "model": model,
            "seed": int(seed),
            "injection_model": "penoyre_oriented",
            "alpha_mas": float(alpha_mas),
            "beta_mas": float(beta_mas),
            "beta_over_alpha": float(elongation),
            "a_rel_mas": float(truth.a_rel_au * truth.parallax_mas),
            "a_over_alpha": float(truth.a_rel_au * truth.parallax_mas / alpha_mas),
            "true_beta_g": float(truth.beta_g),
            "chi2": float(result.chi2),
            "reduced_chi2": float(result.reduced_chi2),
            "success": bool(result.success),
            "scientific_valid": scientific_valid,
            "nfev": int(result.nfev),
            "gaia_n_transits_native": int(selection.n_total),
            "gaia_n_single_peak_used": int(selection.n_single_peak),
            "gaia_n_multi_peak_injection": int(selection.n_multi_peak),
            "gaia_final_multi_peak_predicted": int(final_multi),
            "gaia_ra_deg": gaia_schedule.ra_deg,
            "gaia_dec_deg": gaia_schedule.dec_deg,
            "gaia_release": gaia_schedule.release,
            "gaia_schedule_source": gaia_schedule.source,
        }
        for name, true_value in asdict(truth).items():
            fit_value = float(getattr(result.params, name))
            row[f"true_{name}"] = float(true_value)
            row[f"fit_{name}"] = fit_value
            row[f"bias_{name}"] = fit_value - float(true_value)
            if true_value != 0:
                row[f"frac_bias_{name}"] = (fit_value - float(true_value)) / float(true_value)
        rows.append(row)

    # Add paired discrepancy diagnostics without changing the one-row-per-model
    # representation used elsewhere in the repository.
    chi2 = {row["model"]: row["chi2"] for row in rows}
    for row in rows:
        row["delta_chi2_photo_minus_m2"] = chi2["photocentre"] - chi2["penoyre_oriented"]
        row["delta_chi2_m1_minus_m2"] = chi2["equal_width"] - chi2["penoyre_oriented"]
    return rows


def response_fidelity_scan(
    base_truth: BinaryParams,
    gaia_schedule: GaiaScanSchedule,
    *,
    a_over_alpha_values=(0.4, 0.6, 0.8, 1.0),
    beta_values=(0.05, 0.25, 0.45),
    beta_over_alpha_values=(3.0,),
    seeds=(0, 1, 2),
    alpha_mas: float = 50.0,
    free_names=ALL_PARAMETER_NAMES,
    **kwargs,
) -> list[dict]:
    """Run the paired M0/M1/M2 hierarchy over a controlled physical grid."""
    gaia_schedule.validate()
    if alpha_mas <= 0 or not math.isfinite(float(alpha_mas)):
        raise ValueError("alpha_mas must be finite and > 0")

    rows: list[dict] = []
    for beta_g in beta_values:
        if not (0.0 <= float(beta_g) <= 0.5):
            raise ValueError("beta_values must satisfy 0 <= beta_g <= 0.5")
        for ratio in a_over_alpha_values:
            if float(ratio) <= 0:
                raise ValueError("a_over_alpha_values must be > 0")
            parallax = float(ratio) * float(alpha_mas) / base_truth.a_rel_au
            truth = replace(
                base_truth,
                parallax_mas=parallax,
                beta_g=float(beta_g),
            )
            for elongation in beta_over_alpha_values:
                _validate_elongation(elongation)
                for seed in seeds:
                    rows.extend(response_fidelity_once(
                        truth,
                        gaia_schedule,
                        alpha_mas=float(alpha_mas),
                        beta_over_alpha=float(elongation),
                        seed=int(seed),
                        free_names=free_names,
                        **kwargs,
                    ))
    return rows


def summarise_response_fidelity(rows: list[dict]) -> list[dict]:
    """Median physical-bias summary for each response/grid combination."""
    if not rows:
        raise ValueError("rows is empty")
    keys = ("beta_over_alpha", "true_beta_g", "a_over_alpha", "model")
    groups: dict[tuple, list[dict]] = {}
    for row in rows:
        key = tuple(row[name] for name in keys)
        groups.setdefault(key, []).append(row)

    out: list[dict] = []
    for key in sorted(groups, key=lambda item: tuple(str(x) for x in item)):
        group = groups[key]
        summary = dict(zip(keys, key))
        summary["n"] = len(group)
        summary["success_fraction"] = float(np.mean([r["success"] for r in group]))
        summary["scientific_valid_fraction"] = float(
            np.mean([r["scientific_valid"] for r in group])
        )
        summary["chi2_median"] = float(np.median([r["chi2"] for r in group]))
        for name in ("m1_msun", "m2_msun", "parallax_mas", "beta_g", "inclination_deg"):
            values = np.asarray([r[f"frac_bias_{name}"] for r in group], dtype=float)
            summary[f"frac_bias_{name}_median"] = float(np.median(values))
            summary[f"frac_bias_{name}_q16"] = float(np.percentile(values, 16))
            summary[f"frac_bias_{name}_q84"] = float(np.percentile(values, 84))
        out.append(summary)
    return out
