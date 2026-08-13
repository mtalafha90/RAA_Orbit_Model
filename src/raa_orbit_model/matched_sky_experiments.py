from __future__ import annotations

from dataclasses import replace

import numpy as np

from .experiments import ComparisonDiagnostics, _record, perturbed_start, single_peak_schedule_for_response
from .fit import ALL_PARAMETER_NAMES, fit_joint
from .kepler import BinaryParams
from .matched_control import matched_n_schedule
from .model import GaiaResponseConfig, predict_gaia_orbital_response
from .scanning import GaiaScanSchedule
from .synthetic import GaiaALData, JointData, simulate_joint_data


def _subset_gaia_data(data: JointData, indices: np.ndarray, schedule: GaiaScanSchedule) -> JointData:
    gaia = data.gaia_al
    idx = np.asarray(indices, dtype=int)
    subset = GaiaALData(
        np.asarray(gaia.times_yr)[idx],
        np.asarray(gaia.scan_angle_deg)[idx],
        np.asarray(gaia.values_mas)[idx],
        np.asarray(gaia.sigma_mas)[idx],
        schedule_source=schedule.source,
        ra_deg=schedule.ra_deg,
        dec_deg=schedule.dec_deg,
        release=schedule.release,
    )
    return JointData(data.relative_astrometry, data.sb2_rv, subset)


def compare_models_matched_n(
    truth: BinaryParams,
    native_schedule: GaiaScanSchedule,
    *,
    matched_n_transits: int,
    sigma_response_mas: float,
    external_baseline_yr: float,
    seed: int = 0,
    free_names=ALL_PARAMETER_NAMES,
):
    """Fit both hypotheses after deterministic matched-N Gaia subsampling.

    The full native synthetic realization is generated first. The matched Gaia
    epochs are then selected from that realization, so the selected epochs keep
    the same noise draws as the corresponding native experiment for the same
    seed. The external astrometry and SB2 realization is also unchanged.
    """
    external_baseline_yr = float(external_baseline_yr)
    if external_baseline_yr <= 0:
        raise ValueError("external_baseline_yr must be > 0")
    native_schedule.validate()
    matched_schedule, indices = matched_n_schedule(native_schedule, matched_n_transits)

    injected_response = GaiaResponseConfig("blended_gaussian_peak", sigma_response_mas)
    native_selection = single_peak_schedule_for_response(truth, native_schedule, injected_response)
    if native_selection.n_multi_peak != 0:
        raise ValueError(
            "matched-N control is restricted to a native single-peak grid; "
            f"found {native_selection.n_multi_peak} multi-peak epoch(s)"
        )
    matched_selection = single_peak_schedule_for_response(truth, matched_schedule, injected_response)
    if matched_selection.n_multi_peak != 0:
        raise RuntimeError("matched-N subschedule unexpectedly contains multi-peak epochs")

    full_data = simulate_joint_data(
        truth,
        injected_response,
        native_schedule,
        seed=seed,
        baseline_yr=external_baseline_yr,
    )
    data = _subset_gaia_data(full_data, indices, matched_schedule)
    initial = perturbed_start(truth)
    photo = fit_joint(data, initial, GaiaResponseConfig("photocentre"), free_names=free_names)
    raa = fit_joint(
        data,
        initial,
        GaiaResponseConfig("blended_gaussian_peak", sigma_response_mas, allow_multi_peak_continuation=True),
        free_names=free_names,
    )
    final_prediction = predict_gaia_orbital_response(
        data.gaia_al.times_yr,
        data.gaia_al.scan_angle_deg,
        raa.params,
        injected_response,
    )
    diagnostics = ComparisonDiagnostics(
        selection=matched_selection,
        raa_scientific_valid=(raa.success and final_prediction.n_multi_peak == 0),
        raa_final_multi_peak_predicted=final_prediction.n_multi_peak,
    )
    return photo, raa, diagnostics, matched_schedule


def matched_n_sky_resolution_bias_scan(
    base_truth: BinaryParams,
    native_schedule: GaiaScanSchedule,
    *,
    matched_n_transits: int,
    external_baseline_yr: float,
    a_over_sigma_values=(0.40, 0.50, 0.60, 1.00),
    beta_values=(0.05, 0.25, 0.45),
    seeds=tuple(range(10)),
    sigma_response_mas: float = 50.0,
    free_names=ALL_PARAMETER_NAMES,
) -> list[dict]:
    """Reduced matched-N full-sky control."""
    rows: list[dict] = []
    for beta in tuple(beta_values):
        if not (0.0 <= beta <= 0.5):
            raise ValueError("matched-N sky control requires 0 <= beta_g <= 0.5")
        for ratio in tuple(a_over_sigma_values):
            if ratio <= 0:
                raise ValueError("a_over_sigma_values must be > 0")
            target_a_mas = float(ratio) * sigma_response_mas
            truth = replace(
                base_truth,
                parallax_mas=target_a_mas / base_truth.a_rel_au,
                beta_g=float(beta),
            )
            for seed in tuple(seeds):
                photo, raa, diagnostics, matched_schedule = compare_models_matched_n(
                    truth,
                    native_schedule,
                    matched_n_transits=matched_n_transits,
                    sigma_response_mas=sigma_response_mas,
                    external_baseline_yr=external_baseline_yr,
                    seed=int(seed),
                    free_names=free_names,
                )
                for model_name, result in (("photocentre", photo), ("resolution_aware", raa)):
                    row = _record(
                        model_name,
                        truth,
                        result,
                        sigma_response_mas,
                        int(seed),
                        matched_schedule,
                        diagnostics,
                    )
                    row["external_data_baseline_yr"] = external_baseline_yr
                    row["native_gaia_n_transits"] = native_schedule.n_transits
                    row["matched_n_transits"] = matched_n_transits
                    row["control_subsample_strategy"] = "rank_stratified_equal_count"
                    rows.append(row)
    return rows
