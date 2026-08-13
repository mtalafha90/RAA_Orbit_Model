from __future__ import annotations

from dataclasses import asdict, replace
import csv

from .fit import ALL_PARAMETER_NAMES, fit_joint
from .kepler import BinaryParams
from .model import GaiaResponseConfig
from .scanning import GaiaScanSchedule
from .synthetic import simulate_joint_data


def perturbed_start(truth: BinaryParams) -> BinaryParams:
    """Deterministic near-truth start for measurement-model bias experiments."""
    return replace(
        truth,
        period_yr=truth.period_yr * 1.015,
        t_peri_yr=truth.t_peri_yr + 0.03 * truth.period_yr,
        eccentricity=min(0.95, truth.eccentricity + 0.025),
        inclination_deg=min(175.0, truth.inclination_deg + 2.0),
        omega_deg=truth.omega_deg + 4.0,
        node_deg=truth.node_deg - 3.0,
        m1_msun=truth.m1_msun * 1.05,
        m2_msun=truth.m2_msun * 0.95,
        parallax_mas=truth.parallax_mas * 1.04,
        gamma_kms=truth.gamma_kms + 0.2,
        beta_g=min(0.48, max(0.01, truth.beta_g + 0.03)),
    )


def compare_models_once(
    truth: BinaryParams,
    gaia_schedule: GaiaScanSchedule,
    *,
    sigma_response_mas: float,
    seed: int = 0,
    free_names=ALL_PARAMETER_NAMES,
    n_ast: int = 24,
    n_rv: int = 48,
    ast_sigma_mas: float = 0.20,
    rv_sigma_kms: float = 0.10,
    gaia_sigma_mas: float = 0.10,
):
    """Inject with the resolution-aware surrogate and fit both hypotheses."""
    injected_response = GaiaResponseConfig("blended_gaussian_peak", sigma_response_mas)
    data = simulate_joint_data(
        truth,
        injected_response,
        gaia_schedule,
        seed=seed,
        n_ast=n_ast,
        n_rv=n_rv,
        ast_sigma_mas=ast_sigma_mas,
        rv_sigma_kms=rv_sigma_kms,
        gaia_sigma_mas=gaia_sigma_mas,
    )
    initial = perturbed_start(truth)
    photo = fit_joint(data, initial, GaiaResponseConfig("photocentre"), free_names=free_names)
    raa = fit_joint(data, initial, injected_response, free_names=free_names)
    return data, photo, raa


def _record(
    model_name: str,
    truth: BinaryParams,
    result,
    sigma_response_mas: float,
    seed: int,
    schedule: GaiaScanSchedule,
) -> dict:
    row = {
        "model": model_name,
        "seed": int(seed),
        "sigma_response_mas": float(sigma_response_mas),
        "a_rel_mas": float(truth.a_rel_au * truth.parallax_mas),
        "a_over_sigma": float(truth.a_rel_au * truth.parallax_mas / sigma_response_mas),
        "chi2": result.chi2,
        "reduced_chi2": result.reduced_chi2,
        "success": result.success,
        "nfev": result.nfev,
        "gaia_ra_deg": schedule.ra_deg,
        "gaia_dec_deg": schedule.dec_deg,
        "gaia_release": schedule.release,
        "gaia_schedule_source": schedule.source,
        "gaia_n_transits": schedule.n_transits,
        "gaia_schedule_duration_yr": schedule.duration_yr,
    }
    for name, true_value in asdict(truth).items():
        fit_value = getattr(result.params, name)
        row[f"true_{name}"] = true_value
        row[f"fit_{name}"] = fit_value
        row[f"bias_{name}"] = fit_value - true_value
        if true_value != 0:
            row[f"frac_bias_{name}"] = (fit_value - true_value) / true_value
    return row


def resolution_bias_scan(
    base_truth: BinaryParams,
    gaia_schedule: GaiaScanSchedule,
    *,
    a_over_sigma_values=(0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 4.0),
    beta_values=None,
    seeds=(0, 1, 2),
    sigma_response_mas: float = 50.0,
    free_names=ALL_PARAMETER_NAMES,
    n_ast: int = 24,
    n_rv: int = 48,
    ast_sigma_mas: float = 0.20,
    rv_sigma_kms: float = 0.10,
    gaia_sigma_mas: float = 0.10,
) -> list[dict]:
    """Map bias versus angular resolution using one explicit Gaia schedule."""
    gaia_schedule.validate()
    if sigma_response_mas <= 0:
        raise ValueError("sigma_response_mas must be > 0")
    beta_values = (base_truth.beta_g,) if beta_values is None else tuple(beta_values)
    rows: list[dict] = []

    for beta in beta_values:
        if not (0.0 <= beta <= 0.5):
            raise ValueError("prototype bias scan requires 0 <= beta_g <= 0.5")
        for ratio in a_over_sigma_values:
            if ratio <= 0:
                raise ValueError("a_over_sigma_values must be > 0")
            target_a_mas = float(ratio) * sigma_response_mas
            parallax = target_a_mas / base_truth.a_rel_au
            truth = replace(base_truth, parallax_mas=parallax, beta_g=float(beta))
            for seed in seeds:
                _, photo, raa = compare_models_once(
                    truth,
                    gaia_schedule,
                    sigma_response_mas=sigma_response_mas,
                    seed=int(seed),
                    free_names=free_names,
                    n_ast=n_ast,
                    n_rv=n_rv,
                    ast_sigma_mas=ast_sigma_mas,
                    rv_sigma_kms=rv_sigma_kms,
                    gaia_sigma_mas=gaia_sigma_mas,
                )
                rows.append(_record("photocentre", truth, photo, sigma_response_mas, int(seed), gaia_schedule))
                rows.append(_record("resolution_aware", truth, raa, sigma_response_mas, int(seed), gaia_schedule))
    return rows


def write_rows_csv(rows: list[dict], path: str) -> None:
    if not rows:
        raise ValueError("rows is empty")
    fields = sorted({key for row in rows for key in row})
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
