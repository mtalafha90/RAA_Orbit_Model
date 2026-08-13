from __future__ import annotations

from dataclasses import replace

from .experiments import ComparisonDiagnostics, _record, perturbed_start, single_peak_schedule_for_response
from .fit import ALL_PARAMETER_NAMES, fit_joint
from .kepler import BinaryParams
from .model import GaiaResponseConfig, predict_gaia_orbital_response
from .scanning import GaiaScanSchedule
from .synthetic import simulate_joint_data


def compare_models_fixed_external_baseline(
    truth: BinaryParams,
    gaia_schedule: GaiaScanSchedule,
    *,
    sigma_response_mas: float,
    external_baseline_yr: float,
    seed: int = 0,
    free_names=ALL_PARAMETER_NAMES,
    n_ast: int = 24,
    n_rv: int = 48,
    ast_sigma_mas: float = 0.20,
    rv_sigma_kms: float = 0.10,
    gaia_sigma_mas: float = 0.10,
):
    """Fit both hypotheses while holding external-data time coverage fixed."""
    external_baseline_yr = float(external_baseline_yr)
    if external_baseline_yr <= 0:
        raise ValueError("external_baseline_yr must be > 0")

    injected_response = GaiaResponseConfig("blended_gaussian_peak", sigma_response_mas)
    selection = single_peak_schedule_for_response(truth, gaia_schedule, injected_response)
    data = simulate_joint_data(
        truth,
        injected_response,
        selection.schedule,
        seed=seed,
        n_ast=n_ast,
        n_rv=n_rv,
        baseline_yr=external_baseline_yr,
        ast_sigma_mas=ast_sigma_mas,
        rv_sigma_kms=rv_sigma_kms,
        gaia_sigma_mas=gaia_sigma_mas,
    )
    initial = perturbed_start(truth)
    photo = fit_joint(data, initial, GaiaResponseConfig("photocentre"), free_names=free_names)
    fit_response = GaiaResponseConfig(
        "blended_gaussian_peak",
        sigma_response_mas,
        allow_multi_peak_continuation=True,
    )
    raa = fit_joint(data, initial, fit_response, free_names=free_names)
    final_prediction = predict_gaia_orbital_response(
        data.gaia_al.times_yr,
        data.gaia_al.scan_angle_deg,
        raa.params,
        injected_response,
    )
    diagnostics = ComparisonDiagnostics(
        selection=selection,
        raa_scientific_valid=(raa.success and final_prediction.n_multi_peak == 0),
        raa_final_multi_peak_predicted=final_prediction.n_multi_peak,
    )
    return photo, raa, diagnostics


def sky_resolution_bias_scan(
    base_truth: BinaryParams,
    gaia_schedule: GaiaScanSchedule,
    *,
    external_baseline_yr: float,
    a_over_sigma_values=(0.40, 0.50, 0.60, 0.80, 1.00),
    beta_values=(0.05, 0.15, 0.25, 0.35, 0.45),
    seeds=tuple(range(10)),
    sigma_response_mas: float = 50.0,
    free_names=ALL_PARAMETER_NAMES,
    n_ast: int = 24,
    n_rv: int = 48,
    ast_sigma_mas: float = 0.20,
    rv_sigma_kms: float = 0.10,
    gaia_sigma_mas: float = 0.10,
) -> list[dict]:
    """Sky-position version of the paired bias scan with fixed external data."""
    gaia_schedule.validate()
    external_baseline_yr = float(external_baseline_yr)
    if external_baseline_yr <= 0:
        raise ValueError("external_baseline_yr must be > 0")
    if sigma_response_mas <= 0:
        raise ValueError("sigma_response_mas must be > 0")

    rows: list[dict] = []
    for beta in tuple(beta_values):
        if not (0.0 <= beta <= 0.5):
            raise ValueError("sky bias scan requires 0 <= beta_g <= 0.5")
        for ratio in tuple(a_over_sigma_values):
            if ratio <= 0:
                raise ValueError("a_over_sigma_values must be > 0")
            target_a_mas = float(ratio) * sigma_response_mas
            parallax = target_a_mas / base_truth.a_rel_au
            truth = replace(base_truth, parallax_mas=parallax, beta_g=float(beta))
            for seed in tuple(seeds):
                photo, raa, diagnostics = compare_models_fixed_external_baseline(
                    truth,
                    gaia_schedule,
                    sigma_response_mas=sigma_response_mas,
                    external_baseline_yr=external_baseline_yr,
                    seed=int(seed),
                    free_names=free_names,
                    n_ast=n_ast,
                    n_rv=n_rv,
                    ast_sigma_mas=ast_sigma_mas,
                    rv_sigma_kms=rv_sigma_kms,
                    gaia_sigma_mas=gaia_sigma_mas,
                )
                for model_name, result in (("photocentre", photo), ("resolution_aware", raa)):
                    row = _record(
                        model_name,
                        truth,
                        result,
                        sigma_response_mas,
                        int(seed),
                        gaia_schedule,
                        diagnostics,
                    )
                    row["external_data_baseline_yr"] = external_baseline_yr
                    rows.append(row)
    return rows
