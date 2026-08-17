"""Targeted controls following the 720-fit response-fidelity experiment.

Two questions are isolated here rather than expanding the full response grid.

1. External-information strength: at a fixed difficult M2 injection, how much
   does the protection supplied by resolved relative astrometry and SB2 RVs
   determine the physical mass/parallax bias under M0, M1, and M2?
2. Matched-M2 estimator control: when the injection and fitted response are
   identical, does the small positive parallax median seen in the ten-seed
   response grid persist over many independent noise/time realizations?

These are deterministic-fit controls.  They do not establish posterior
coverage; that remains a later sampling experiment.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import math

import numpy as np

from .experiments import perturbed_start, single_peak_schedule_for_response
from .fit import ALL_PARAMETER_NAMES, fit_joint
from .kepler import BinaryParams
from .model import GaiaResponseConfig, predict_gaia_orbital_response
from .response_fidelity import response_fidelity_once
from .scanning import GaiaScanSchedule
from .synthetic import simulate_joint_data


@dataclass(frozen=True)
class ExternalInformationLevel:
    """Precision/count prescription for the non-Gaia orbit constraints."""

    name: str
    ast_sigma_mas: float
    rv_sigma_kms: float
    n_ast: int = 24
    n_rv: int = 48

    def validate(self) -> None:
        if not self.name:
            raise ValueError("external-information level needs a non-empty name")
        if self.ast_sigma_mas <= 0 or self.rv_sigma_kms <= 0:
            raise ValueError("external-information uncertainties must be > 0")
        if self.n_ast < 0 or self.n_rv < 0:
            raise ValueError("external-information epoch counts must be >= 0")


DEFAULT_EXTERNAL_LEVELS = (
    ExternalInformationLevel("strong", ast_sigma_mas=0.20, rv_sigma_kms=0.10),
    ExternalInformationLevel("medium", ast_sigma_mas=1.00, rv_sigma_kms=0.50),
    ExternalInformationLevel("weak", ast_sigma_mas=2.00, rv_sigma_kms=1.00),
)


def truth_at_response_point(
    base_truth: BinaryParams,
    *,
    alpha_mas: float,
    a_over_alpha: float,
    beta_g: float,
) -> BinaryParams:
    """Move the base orbit to one angular-separation/light-fraction point."""
    if alpha_mas <= 0 or not math.isfinite(float(alpha_mas)):
        raise ValueError("alpha_mas must be finite and > 0")
    if a_over_alpha <= 0:
        raise ValueError("a_over_alpha must be > 0")
    if not (0.0 <= beta_g <= 0.5):
        raise ValueError("beta_g must satisfy 0 <= beta_g <= 0.5")
    parallax = float(a_over_alpha) * float(alpha_mas) / base_truth.a_rel_au
    return replace(base_truth, parallax_mas=parallax, beta_g=float(beta_g))


def external_information_strength_scan(
    base_truth: BinaryParams,
    gaia_schedule: GaiaScanSchedule,
    *,
    levels=DEFAULT_EXTERNAL_LEVELS,
    seeds=tuple(range(30)),
    alpha_mas: float = 50.0,
    a_over_alpha: float = 1.0,
    beta_g: float = 0.25,
    beta_over_alpha: float = 1.5,
    gaia_sigma_mas: float = 0.10,
    free_names=ALL_PARAMETER_NAMES,
) -> list[dict]:
    """Run M0/M1/M2 while varying only the external orbit information.

    For a given level and seed all three models see exactly the same synthetic
    realization.  Across levels, the seed is reused deliberately while the
    stated external uncertainties change.  The Gaia uncertainty and schedule
    stay fixed, so differences between levels isolate the strength of the
    resolved-astrometry/SB2 constraints within this synthetic design.
    """
    gaia_schedule.validate()
    if gaia_sigma_mas <= 0:
        raise ValueError("gaia_sigma_mas must be > 0")
    if beta_over_alpha < 1:
        raise ValueError("beta_over_alpha must be >= 1")

    truth = truth_at_response_point(
        base_truth,
        alpha_mas=alpha_mas,
        a_over_alpha=a_over_alpha,
        beta_g=beta_g,
    )
    rows: list[dict] = []
    for level in tuple(levels):
        level.validate()
        for seed in seeds:
            group = response_fidelity_once(
                truth,
                gaia_schedule,
                alpha_mas=float(alpha_mas),
                beta_over_alpha=float(beta_over_alpha),
                seed=int(seed),
                free_names=free_names,
                n_ast=int(level.n_ast),
                n_rv=int(level.n_rv),
                ast_sigma_mas=float(level.ast_sigma_mas),
                rv_sigma_kms=float(level.rv_sigma_kms),
                gaia_sigma_mas=float(gaia_sigma_mas),
            )
            for row in group:
                row["external_level"] = level.name
                row["n_ast"] = int(level.n_ast)
                row["n_rv"] = int(level.n_rv)
                row["ast_sigma_mas"] = float(level.ast_sigma_mas)
                row["rv_sigma_kms"] = float(level.rv_sigma_kms)
                row["gaia_sigma_mas"] = float(gaia_sigma_mas)
            rows.extend(group)
    return rows


def summarise_external_information(rows: list[dict]) -> list[dict]:
    """Summarise physical bias by external-information level and fitted model."""
    if not rows:
        raise ValueError("rows is empty")
    groups: dict[tuple[str, str], list[dict]] = {}
    for row in rows:
        groups.setdefault((row["external_level"], row["model"]), []).append(row)

    out: list[dict] = []
    order = {level.name: i for i, level in enumerate(DEFAULT_EXTERNAL_LEVELS)}
    for (level, model), group in sorted(
        groups.items(), key=lambda item: (order.get(item[0][0], 999), item[0][1])
    ):
        first = group[0]
        summary = {
            "external_level": level,
            "model": model,
            "n": len(group),
            "n_ast": first["n_ast"],
            "n_rv": first["n_rv"],
            "ast_sigma_mas": first["ast_sigma_mas"],
            "rv_sigma_kms": first["rv_sigma_kms"],
            "gaia_sigma_mas": first["gaia_sigma_mas"],
            "success_fraction": float(np.mean([r["success"] for r in group])),
            "scientific_valid_fraction": float(
                np.mean([r["scientific_valid"] for r in group])
            ),
            "chi2_median": float(np.median([r["chi2"] for r in group])),
            "delta_chi2_photo_minus_m2_median": float(
                np.median([r["delta_chi2_photo_minus_m2"] for r in group])
            ),
            "delta_chi2_m1_minus_m2_median": float(
                np.median([r["delta_chi2_m1_minus_m2"] for r in group])
            ),
        }
        for name in ("m1_msun", "m2_msun", "parallax_mas", "beta_g"):
            values = np.asarray([r[f"frac_bias_{name}"] for r in group], dtype=float)
            summary[f"frac_bias_{name}_median"] = float(np.median(values))
            summary[f"frac_bias_{name}_q16"] = float(np.percentile(values, 16))
            summary[f"frac_bias_{name}_q84"] = float(np.percentile(values, 84))
            summary[f"frac_abs_bias_{name}_median"] = float(np.median(np.abs(values)))
        out.append(summary)
    return out


def matched_m2_once(
    truth: BinaryParams,
    gaia_schedule: GaiaScanSchedule,
    *,
    alpha_mas: float,
    beta_over_alpha: float,
    seed: int,
    n_ast: int = 24,
    n_rv: int = 48,
    ast_sigma_mas: float = 0.20,
    rv_sigma_kms: float = 0.10,
    gaia_sigma_mas: float = 0.10,
    free_names=ALL_PARAMETER_NAMES,
) -> dict:
    """One matched M2 injection/recovery realization with no M0/M1 fits."""
    if alpha_mas <= 0 or not math.isfinite(float(alpha_mas)):
        raise ValueError("alpha_mas must be finite and > 0")
    if beta_over_alpha < 1:
        raise ValueError("beta_over_alpha must be >= 1")
    beta_mas = float(alpha_mas) * float(beta_over_alpha)
    response = GaiaResponseConfig(
        "penoyre_gaussian_peak",
        float(alpha_mas),
        sigma_ac_mas=beta_mas,
    )
    selection = single_peak_schedule_for_response(truth, gaia_schedule, response)
    data = simulate_joint_data(
        truth,
        response,
        selection.schedule,
        seed=int(seed),
        n_ast=int(n_ast),
        n_rv=int(n_rv),
        baseline_yr=gaia_schedule.mission_span_yr,
        ast_sigma_mas=float(ast_sigma_mas),
        rv_sigma_kms=float(rv_sigma_kms),
        gaia_sigma_mas=float(gaia_sigma_mas),
    )
    fit_response = replace(response, allow_multi_peak_continuation=True)
    result = fit_joint(data, perturbed_start(truth), fit_response, free_names=free_names)
    final = predict_gaia_orbital_response(
        data.gaia_al.times_yr,
        data.gaia_al.scan_angle_deg,
        result.params,
        response,
    )
    row = {
        "model": "penoyre_oriented",
        "injection_model": "penoyre_oriented",
        "seed": int(seed),
        "alpha_mas": float(alpha_mas),
        "beta_mas": beta_mas,
        "beta_over_alpha": float(beta_over_alpha),
        "a_rel_mas": float(truth.a_rel_au * truth.parallax_mas),
        "a_over_alpha": float(truth.a_rel_au * truth.parallax_mas / alpha_mas),
        "true_beta_g": float(truth.beta_g),
        "chi2": float(result.chi2),
        "reduced_chi2": float(result.reduced_chi2),
        "success": bool(result.success),
        "scientific_valid": bool(result.success and final.n_multi_peak == 0),
        "nfev": int(result.nfev),
        "gaia_n_transits_native": int(selection.n_total),
        "gaia_n_single_peak_used": int(selection.n_single_peak),
        "gaia_n_multi_peak_injection": int(selection.n_multi_peak),
        "gaia_final_multi_peak_predicted": int(final.n_multi_peak),
        "n_ast": int(n_ast),
        "n_rv": int(n_rv),
        "ast_sigma_mas": float(ast_sigma_mas),
        "rv_sigma_kms": float(rv_sigma_kms),
        "gaia_sigma_mas": float(gaia_sigma_mas),
    }
    for name, true_value in asdict(truth).items():
        fitted = float(getattr(result.params, name))
        row[f"true_{name}"] = float(true_value)
        row[f"fit_{name}"] = fitted
        row[f"bias_{name}"] = fitted - float(true_value)
        if true_value != 0:
            row[f"frac_bias_{name}"] = (fitted - float(true_value)) / float(true_value)
    return row


def matched_m2_scan(
    base_truth: BinaryParams,
    gaia_schedule: GaiaScanSchedule,
    *,
    seeds=tuple(range(100)),
    alpha_mas: float = 50.0,
    a_over_alpha: float = 1.0,
    beta_g: float = 0.25,
    beta_over_alpha: float = 1.5,
    **kwargs,
) -> list[dict]:
    """High-statistics matched M2 control at one physical response point."""
    truth = truth_at_response_point(
        base_truth,
        alpha_mas=alpha_mas,
        a_over_alpha=a_over_alpha,
        beta_g=beta_g,
    )
    return [
        matched_m2_once(
            truth,
            gaia_schedule,
            alpha_mas=float(alpha_mas),
            beta_over_alpha=float(beta_over_alpha),
            seed=int(seed),
            **kwargs,
        )
        for seed in seeds
    ]


def summarise_matched_m2(rows: list[dict]) -> list[dict]:
    """One-row summary focused on whether matched-M2 biases centre on zero."""
    if not rows:
        raise ValueError("rows is empty")
    summary = {
        "n": len(rows),
        "success_fraction": float(np.mean([r["success"] for r in rows])),
        "scientific_valid_fraction": float(np.mean([r["scientific_valid"] for r in rows])),
        "chi2_median": float(np.median([r["chi2"] for r in rows])),
        "reduced_chi2_median": float(np.median([r["reduced_chi2"] for r in rows])),
    }
    for name in ("m1_msun", "m2_msun", "parallax_mas", "beta_g", "inclination_deg"):
        values = np.asarray([r[f"frac_bias_{name}"] for r in rows], dtype=float)
        summary[f"frac_bias_{name}_mean"] = float(np.mean(values))
        summary[f"frac_bias_{name}_median"] = float(np.median(values))
        summary[f"frac_bias_{name}_std"] = float(np.std(values, ddof=1)) if len(values) > 1 else 0.0
        summary[f"frac_bias_{name}_q16"] = float(np.percentile(values, 16))
        summary[f"frac_bias_{name}_q84"] = float(np.percentile(values, 84))
        summary[f"fraction_positive_{name}"] = float(np.mean(values > 0.0))
    return [summary]
