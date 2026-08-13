#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

PAIR = ["sky_id", "true_beta_g", "a_over_sigma", "seed"]
META = [
    "sky_id", "ecliptic_lon_deg", "ecliptic_lat_deg", "scan_n_transits",
    "scan_directional_max_gap_deg", "scan_directional_entropy_normalized",
    "scan_directional_resultant_length", "scan_directional_effective_bins",
    "scan_time_angle_r2",
]


def first_ratio(group, mask):
    selected = group.loc[mask].sort_values("a_over_sigma")
    return np.nan if selected.empty else float(selected.iloc[0]["a_over_sigma"])


def sky_plot(frame, value, label, title, path, *, log1p=False):
    values = frame[value].to_numpy(float)
    if log1p:
        values = np.log10(1.0 + np.maximum(values, 0.0))
    lon = ((frame["ecliptic_lon_deg"].to_numpy(float) + 180.0) % 360.0) - 180.0
    fig, ax = plt.subplots(figsize=(8.0, 4.6))
    sc = ax.scatter(lon, frame["ecliptic_lat_deg"], c=values, s=145, marker="s")
    ax.set(xlim=(-195, 195), ylim=(-97, 97), xlabel="J2000 ecliptic longitude [deg]", ylabel="J2000 ecliptic latitude [deg]", title=title)
    ax.set_xticks([-180, -90, 0, 90, 180]); ax.set_yticks([-90, -60, -30, 0, 30, 60, 90])
    fig.colorbar(sc, ax=ax).set_label(label)
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description="Analyze the full-sky RAA bias experiment")
    p.add_argument("input_csv")
    p.add_argument("--output-dir", default="results/sky_position_dr4_full/analysis")
    p.add_argument("--prefix", default="full_sky_dr4")
    p.add_argument("--beta-map", type=float, default=0.25)
    p.add_argument("--a-over-sigma-map", type=float, default=1.0)
    args = p.parse_args()

    raw = pd.read_csv(args.input_csv)
    required = {*PAIR, *META, "model", "chi2", "success", "scientific_valid"}
    missing = required - set(raw.columns)
    if missing:
        raise ValueError(f"missing columns: {sorted(missing)}")
    if not raw["success"].astype(bool).all() or not raw["scientific_valid"].astype(bool).all():
        raise ValueError("input contains failed or scientifically invalid fits")
    if raw.duplicated([*PAIR, "model"]).any():
        raise ValueError("input contains duplicate model records")

    wide = raw.pivot(index=PAIR, columns="model", values="chi2").reset_index()
    wide["delta_chi2"] = wide["photocentre"] - wide["resolution_aware"]
    meta = raw[META].drop_duplicates("sky_id")
    paired = wide.merge(meta, on="sky_id", validate="many_to_one")
    summary = paired.groupby(["sky_id", "true_beta_g", "a_over_sigma"], as_index=False).agg(
        n_pairs=("delta_chi2", "size"),
        delta_chi2_median=("delta_chi2", "median"),
        delta_chi2_mean=("delta_chi2", "mean"),
        delta_chi2_q16=("delta_chi2", lambda x: x.quantile(0.16)),
        delta_chi2_q84=("delta_chi2", lambda x: x.quantile(0.84)),
        delta_chi2_min=("delta_chi2", "min"), delta_chi2_max=("delta_chi2", "max"),
    ).merge(meta, on="sky_id", validate="many_to_one")

    rows = []
    for (sky, beta), g in summary.groupby(["sky_id", "true_beta_g"]):
        g = g.sort_values("a_over_sigma"); m = meta.set_index("sky_id").loc[sky]
        row = {"sky_id": sky, "true_beta_g": beta, **m.to_dict()}
        row["q16_positive_a_over_sigma"] = first_ratio(g, g["delta_chi2_q16"] > 0)
        row["all_positive_a_over_sigma"] = first_ratio(g, g["delta_chi2_min"] > 0)
        row["median_gt10_a_over_sigma"] = first_ratio(g, g["delta_chi2_median"] > 10)
        rows.append(row)
    boundaries = pd.DataFrame(rows)

    selected = summary[np.isclose(summary.true_beta_g, args.beta_map) & np.isclose(summary.a_over_sigma, args.a_over_sigma_map)]
    corr = []
    for metric in META[3:]:
        r = spearmanr(selected[metric], selected["delta_chi2_median"], nan_policy="omit")
        corr.append({"metric": metric, "spearman_rho": r.statistic, "spearman_pvalue": r.pvalue, "n_sky": len(selected)})
    corr = pd.DataFrame(corr)

    out = Path(args.output_dir); out.mkdir(parents=True, exist_ok=True); pre = args.prefix
    paired.to_csv(out / f"{pre}_paired.csv", index=False)
    summary.to_csv(out / f"{pre}_summary.csv", index=False)
    boundaries.to_csv(out / f"{pre}_boundaries.csv", index=False)
    corr.to_csv(out / f"{pre}_correlations.csv", index=False)

    bsel = boundaries[np.isclose(boundaries.true_beta_g, args.beta_map)]
    sky_plot(selected, "delta_chi2_median", r"$\log_{10}[1+\max(\mathrm{median}\,\Delta\chi^2,0)]$", "Photocentre--RAA mismatch across the Gaia sky", out / f"{pre}_delta_chi2_sky", log1p=True)
    sky_plot(bsel, "q16_positive_a_over_sigma", r"first sampled $a/\sigma$ with $Q_{16}(\Delta\chi^2)>0$", "Sampled robust-positive transition boundary", out / f"{pre}_q16_boundary_sky")
    sky_plot(selected, "scan_n_transits", "nominal DR4 transit count", "Native Gaia nominal transit count", out / f"{pre}_native_transit_count_sky")
    sky_plot(selected, "scan_directional_max_gap_deg", "maximum directional scan-angle gap [deg]", "Native directional scan-angle coverage", out / f"{pre}_native_max_gap_sky")

    print(f"fit records: {len(raw)}; sky positions: {raw.sky_id.nunique()}; paired realizations: {len(paired)}")
    print(corr.to_string(index=False))


if __name__ == "__main__":
    main()
