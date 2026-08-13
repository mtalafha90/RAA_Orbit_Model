#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


PAIR_KEYS = [
    "sky_id",
    "ecliptic_lon_deg",
    "ecliptic_lat_deg",
    "true_beta_g",
    "a_over_sigma",
    "seed",
]


def pair_models(df: pd.DataFrame) -> pd.DataFrame:
    required = set(PAIR_KEYS) | {"model", "chi2"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"input is missing required columns: {sorted(missing)}")
    wide = df.pivot(index=PAIR_KEYS, columns="model", values="chi2").reset_index()
    for model in ("photocentre", "resolution_aware"):
        if model not in wide.columns:
            raise ValueError(f"input is missing model {model!r}")
    wide["delta_chi2"] = wide["photocentre"] - wide["resolution_aware"]
    return wide


def sky_summary(paired: pd.DataFrame, beta: float, ratio: float) -> pd.DataFrame:
    sub = paired[
        np.isclose(paired["true_beta_g"], beta)
        & np.isclose(paired["a_over_sigma"], ratio)
    ]
    if sub.empty:
        raise ValueError(f"no paired rows for beta={beta:g}, a/sigma={ratio:g}")
    return (
        sub.groupby(["sky_id", "ecliptic_lon_deg", "ecliptic_lat_deg"])["delta_chi2"]
        .agg(
            delta_chi2_median="median",
            delta_chi2_q16=lambda x: x.quantile(0.16),
            delta_chi2_q84=lambda x: x.quantile(0.84),
            delta_chi2_min="min",
            delta_chi2_max="max",
        )
        .reset_index()
    )


def boundaries(paired: pd.DataFrame, beta: float) -> pd.DataFrame:
    sub = (
        paired[np.isclose(paired["true_beta_g"], beta)]
        .groupby(["sky_id", "ecliptic_lon_deg", "ecliptic_lat_deg", "a_over_sigma"])["delta_chi2"]
        .agg(
            delta_chi2_median="median",
            delta_chi2_q16=lambda x: x.quantile(0.16),
            delta_chi2_min="min",
        )
        .reset_index()
    )
    rows = []
    for keys, group in sub.groupby(["sky_id", "ecliptic_lon_deg", "ecliptic_lat_deg"]):
        group = group.sort_values("a_over_sigma")

        def first(mask: pd.Series) -> float:
            x = group.loc[mask, "a_over_sigma"]
            return float(x.iloc[0]) if len(x) else np.nan

        rows.append(
            {
                "sky_id": keys[0],
                "ecliptic_lon_deg": keys[1],
                "ecliptic_lat_deg": keys[2],
                "q16_positive_a_over_sigma": first(group["delta_chi2_q16"] > 0),
                "all_positive_a_over_sigma": first(group["delta_chi2_min"] > 0),
                "median_gt10_a_over_sigma": first(group["delta_chi2_median"] > 10),
            }
        )
    return pd.DataFrame(rows)


def save_fig(fig, output_dir: Path, stem: str, dpi: int) -> None:
    fig.savefig(output_dir / f"{stem}.png", dpi=dpi, bbox_inches="tight")
    fig.savefig(output_dir / f"{stem}.pdf", bbox_inches="tight")


def make_plots(comp: pd.DataFrame, output_dir: Path, dpi: int) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    sc = ax.scatter(
        comp["delta_chi2_median_native"],
        comp["delta_chi2_median_matched_n53"],
        c=comp["native_n_transits"],
        s=65,
        edgecolor="k",
        linewidth=0.35,
    )
    lims = [
        0.85 * min(comp["delta_chi2_median_native"].min(), comp["delta_chi2_median_matched_n53"].min()),
        1.08 * max(comp["delta_chi2_median_native"].max(), comp["delta_chi2_median_matched_n53"].max()),
    ]
    ax.plot(lims, lims, "--", linewidth=1, label=r"$y=x$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel(r"Native median $\Delta\chi^2$")
    ax.set_ylabel(r"Matched-$N=53$ median $\Delta\chi^2$")
    ax.set_title("Native versus matched transit-count control")
    fig.colorbar(sc, ax=ax, pad=0.02, label=r"Native $N_{\rm transit}$")
    ax.legend(frameon=False)
    ax.grid(alpha=0.2, which="both")
    fig.tight_layout()
    save_fig(fig, output_dir, "native_vs_matched_n53", dpi)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    sc = ax.scatter(
        comp["retained_fraction"],
        comp["delta_chi2_ratio"],
        c=comp["matched_entropy_normalized"],
        s=65,
        edgecolor="k",
        linewidth=0.35,
    )
    lo = 0.9 * min(comp["retained_fraction"].min(), comp["delta_chi2_ratio"].min())
    hi = 1.08 * max(comp["retained_fraction"].max(), comp["delta_chi2_ratio"].max())
    ax.plot([lo, hi], [lo, hi], "--", linewidth=1, label=r"$y=x$")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel(r"Retained-transit fraction $53/N_{\rm native}$")
    ax.set_ylabel(r"$\Delta\chi^2_{N=53}/\Delta\chi^2_{\rm native}$")
    ax.set_title("Approximate information scaling with transit count")
    fig.colorbar(sc, ax=ax, pad=0.02, label="Matched scan-angle entropy")
    ax.legend(frameon=False)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    save_fig(fig, output_dir, "matched_n53_scaling", dpi)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare native full-sky and matched-N=53 controls")
    parser.add_argument("native_bias")
    parser.add_argument("matched_bias")
    parser.add_argument("native_geometry")
    parser.add_argument("matched_geometry")
    parser.add_argument("--output-dir", default="results/native_vs_matched_n53")
    parser.add_argument("--beta", type=float, default=0.25)
    parser.add_argument("--a-over-sigma", type=float, default=1.0)
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    native = pd.read_csv(args.native_bias)
    matched = pd.read_csv(args.matched_bias)
    native_geom = pd.read_csv(args.native_geometry)
    matched_geom = pd.read_csv(args.matched_geometry)

    if not native["success"].all() or not native["scientific_valid"].all():
        raise ValueError("native dataset contains failed or scientifically invalid fits")
    if not matched["success"].all() or not matched["scientific_valid"].all():
        raise ValueError("matched dataset contains failed or scientifically invalid fits")
    if native["gaia_n_multi_peak_flagged"].max() != 0 or matched["gaia_n_multi_peak_flagged"].max() != 0:
        raise ValueError("comparison is restricted to the single-peak domain")

    native_paired = pair_models(native)
    matched_paired = pair_models(matched)
    native_sky = sky_summary(native_paired, args.beta, args.a_over_sigma)
    matched_sky = sky_summary(matched_paired, args.beta, args.a_over_sigma)

    native_geometry_columns = [
        "sky_id",
        "scan_n_transits",
        "scan_directional_max_gap_deg",
        "scan_directional_entropy_normalized",
    ]
    matched_geometry_columns = [
        "sky_id",
        "matched_n_transits",
        "matched_max_gap_deg",
        "matched_entropy_normalized",
    ]
    comp = native_sky.merge(
        matched_sky,
        on=["sky_id", "ecliptic_lon_deg", "ecliptic_lat_deg"],
        suffixes=("_native", "_matched_n53"),
    )
    comp = comp.merge(native_geom[native_geometry_columns], on="sky_id")
    comp = comp.merge(matched_geom[matched_geometry_columns], on="sky_id")
    comp = comp.rename(columns={"scan_n_transits": "native_n_transits"})
    comp["delta_chi2_ratio"] = comp["delta_chi2_median_matched_n53"] / comp["delta_chi2_median_native"]
    comp["retained_fraction"] = comp["matched_n_transits"] / comp["native_n_transits"]
    comp["ratio_over_retained_fraction"] = comp["delta_chi2_ratio"] / comp["retained_fraction"]
    comp = comp.sort_values(["ecliptic_lat_deg", "ecliptic_lon_deg"])
    comp.to_csv(output_dir / "native_vs_matched_n53_by_sky.csv", index=False)

    native_boundary = boundaries(native_paired, args.beta)
    matched_boundary = boundaries(matched_paired, args.beta)
    boundary = native_boundary.merge(
        matched_boundary,
        on=["sky_id", "ecliptic_lon_deg", "ecliptic_lat_deg"],
        suffixes=("_native", "_matched_n53"),
    )
    boundary.to_csv(output_dir / "native_vs_matched_n53_boundaries_by_sky.csv", index=False)

    summary = pd.DataFrame(
        [
            {
                "beta_g": args.beta,
                "a_over_sigma": args.a_over_sigma,
                "sky_positions": len(comp),
                "native_delta_chi2_sky_min": comp["delta_chi2_median_native"].min(),
                "native_delta_chi2_sky_median": comp["delta_chi2_median_native"].median(),
                "native_delta_chi2_sky_max": comp["delta_chi2_median_native"].max(),
                "matched_delta_chi2_sky_min": comp["delta_chi2_median_matched_n53"].min(),
                "matched_delta_chi2_sky_median": comp["delta_chi2_median_matched_n53"].median(),
                "matched_delta_chi2_sky_max": comp["delta_chi2_median_matched_n53"].max(),
                "matched_to_native_ratio_of_sky_medians": comp["delta_chi2_median_matched_n53"].median() / comp["delta_chi2_median_native"].median(),
                "spearman_native_delta_vs_native_n": spearmanr(comp["delta_chi2_median_native"], comp["native_n_transits"]).statistic,
                "spearman_matched_delta_vs_native_n": spearmanr(comp["delta_chi2_median_matched_n53"], comp["native_n_transits"]).statistic,
                "spearman_matched_delta_vs_matched_entropy": spearmanr(comp["delta_chi2_median_matched_n53"], comp["matched_entropy_normalized"]).statistic,
                "spearman_matched_delta_vs_matched_max_gap": spearmanr(comp["delta_chi2_median_matched_n53"], comp["matched_max_gap_deg"]).statistic,
                "spearman_delta_ratio_vs_retained_fraction": spearmanr(comp["delta_chi2_ratio"], comp["retained_fraction"]).statistic,
                "median_ratio_over_retained_fraction": comp["ratio_over_retained_fraction"].median(),
                "median_native_delta_chi2_per_transit": (comp["delta_chi2_median_native"] / comp["native_n_transits"]).median(),
                "median_matched_delta_chi2_per_transit": (comp["delta_chi2_median_matched_n53"] / comp["matched_n_transits"]).median(),
            }
        ]
    )
    summary.to_csv(output_dir / "native_vs_matched_n53_summary.csv", index=False)
    make_plots(comp, output_dir, args.dpi)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
