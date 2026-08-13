from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


PAIR_KEYS = ("true_beta_g", "a_over_sigma", "seed")
MODELS = ("photocentre", "resolution_aware")
FRACTIONAL_PARAMETERS = ("beta_g", "parallax_mas", "m1_msun", "m2_msun")


def load_bias_results(path: str | Path) -> pd.DataFrame:
    """Load and validate a bias-scan CSV produced by ``run_bias_scan.py``."""
    df = pd.read_csv(path)
    required = {
        *PAIR_KEYS,
        "model",
        "chi2",
        "reduced_chi2",
        "bias_beta_g",
        "bias_parallax_mas",
        "bias_m1_msun",
        "bias_m2_msun",
        "bias_inclination_deg",
        "true_parallax_mas",
        "true_m1_msun",
        "true_m2_msun",
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"bias-scan CSV is missing required columns: {sorted(missing)}")

    present_models = set(df["model"].dropna().astype(str))
    missing_models = set(MODELS) - present_models
    if missing_models:
        raise ValueError(f"bias-scan CSV is missing model(s): {sorted(missing_models)}")

    duplicated = df.duplicated([*PAIR_KEYS, "model"], keep=False)
    if duplicated.any():
        sample = df.loc[duplicated, [*PAIR_KEYS, "model"]].head().to_dict("records")
        raise ValueError(f"duplicate model records for paired keys; examples: {sample}")

    out = df.copy()
    true_lookup = {
        "beta_g": "true_beta_g",
        "parallax_mas": "true_parallax_mas",
        "m1_msun": "true_m1_msun",
        "m2_msun": "true_m2_msun",
    }
    for name, true_col in true_lookup.items():
        denominator = out[true_col].replace(0.0, np.nan)
        out[f"pct_bias_{name}"] = 100.0 * out[f"bias_{name}"] / denominator
    return out


def build_paired_results(df: pd.DataFrame) -> pd.DataFrame:
    """Pair photocentre and RAA fits from the same injected realization."""
    metrics = [
        "chi2",
        "reduced_chi2",
        "pct_bias_beta_g",
        "pct_bias_parallax_mas",
        "pct_bias_m1_msun",
        "pct_bias_m2_msun",
        "bias_inclination_deg",
    ]
    missing = set(metrics) - set(df.columns)
    if missing:
        raise ValueError(
            "input dataframe must first be passed through load_bias_results; "
            f"missing derived columns: {sorted(missing)}"
        )

    wide = df.pivot(index=list(PAIR_KEYS), columns="model", values=metrics)
    for model in MODELS:
        if model not in set(wide.columns.get_level_values("model")):
            raise ValueError(f"no paired records found for model {model!r}")

    wide.columns = [f"{metric}_{model}" for metric, model in wide.columns]
    wide = wide.reset_index()
    required_paired = [f"chi2_{model}" for model in MODELS]
    incomplete = wide[required_paired].isna().any(axis=1)
    if incomplete.any():
        bad = wide.loc[incomplete, list(PAIR_KEYS)].head().to_dict("records")
        raise ValueError(f"unpaired model records; examples: {bad}")

    wide["delta_chi2"] = wide["chi2_photocentre"] - wide["chi2_resolution_aware"]
    wide["delta_reduced_chi2"] = (
        wide["reduced_chi2_photocentre"] - wide["reduced_chi2_resolution_aware"]
    )
    return wide.sort_values(list(PAIR_KEYS)).reset_index(drop=True)


def _distribution(values: pd.Series) -> dict[str, float]:
    values = values.dropna()
    if len(values) == 0:
        return {
            "median": np.nan,
            "mean": np.nan,
            "q16": np.nan,
            "q84": np.nan,
            "min": np.nan,
            "max": np.nan,
        }
    return {
        "median": float(values.median()),
        "mean": float(values.mean()),
        "q16": float(values.quantile(0.16)),
        "q84": float(values.quantile(0.84)),
        "min": float(values.min()),
        "max": float(values.max()),
    }


def summarize_paired_results(paired: pd.DataFrame) -> pd.DataFrame:
    """Summarize paired model differences and parameter biases over seeds."""
    bias_metrics = [
        *(f"pct_bias_{name}_{model}" for name in FRACTIONAL_PARAMETERS for model in MODELS),
        *(f"bias_inclination_deg_{model}" for model in MODELS),
    ]
    metrics = ["delta_chi2", "delta_reduced_chi2", *bias_metrics]
    missing = set(metrics) - set(paired.columns)
    if missing:
        raise ValueError(f"paired dataframe is missing columns: {sorted(missing)}")

    rows: list[dict] = []
    for (beta, ratio), group in paired.groupby(["true_beta_g", "a_over_sigma"], sort=True):
        row: dict[str, float | int] = {
            "true_beta_g": float(beta),
            "a_over_sigma": float(ratio),
            "n_pairs": int(len(group)),
        }
        for metric in metrics:
            stats = _distribution(group[metric])
            for stat, value in stats.items():
                row[f"{metric}_{stat}"] = value
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["true_beta_g", "a_over_sigma"]).reset_index(drop=True)


def write_analysis_tables(
    paired: pd.DataFrame,
    summary: pd.DataFrame,
    output_dir: str | Path,
    *,
    prefix: str = "transition",
) -> tuple[Path, Path]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paired_path = output_dir / f"{prefix}_paired.csv"
    summary_path = output_dir / f"{prefix}_summary.csv"
    paired.to_csv(paired_path, index=False)
    summary.to_csv(summary_path, index=False)
    return paired_path, summary_path


def _save_figure(fig, output_dir: Path, stem: str, dpi: int = 300) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = [output_dir / f"{stem}.png", output_dir / f"{stem}.pdf"]
    fig.savefig(paths[0], dpi=dpi, bbox_inches="tight")
    fig.savefig(paths[1], bbox_inches="tight")
    return paths


def _plot_delta_chi2(summary: pd.DataFrame, output_dir: Path, prefix: str, dpi: int) -> list[Path]:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    for beta, group in summary.groupby("true_beta_g", sort=True):
        group = group.sort_values("a_over_sigma")
        x = group["a_over_sigma"].to_numpy()
        y = group["delta_chi2_median"].to_numpy()
        lo = group["delta_chi2_q16"].to_numpy()
        hi = group["delta_chi2_q84"].to_numpy()
        line, = ax.plot(x, y, marker="o", label=rf"$\beta_G={beta:g}$")
        ax.fill_between(x, lo, hi, alpha=0.18, color=line.get_color())
    ax.axhline(0.0, linewidth=1.0)
    ax.set_yscale("symlog", linthresh=1.0)
    ax.set_xlabel(r"$a_{\rm rel,ang}/\sigma$")
    ax.set_ylabel(r"$\Delta\chi^2=\chi^2_{\rm photo}-\chi^2_{\rm RAA}$")
    ax.legend(frameon=False)
    fig.tight_layout()
    paths = _save_figure(fig, output_dir, f"{prefix}_delta_chi2", dpi=dpi)
    plt.close(fig)
    return paths


def _plot_parameter_bias(
    summary: pd.DataFrame,
    output_dir: Path,
    *,
    prefix: str,
    metric: str,
    ylabel: str,
    stem: str,
    dpi: int,
) -> list[Path]:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    for beta, group in summary.groupby("true_beta_g", sort=True):
        group = group.sort_values("a_over_sigma")
        x = group["a_over_sigma"].to_numpy()
        color = None
        for model, linestyle, label_name in (
            ("photocentre", "--", "photocentre"),
            ("resolution_aware", "-", "RAA"),
        ):
            base = f"{metric}_{model}"
            y = group[f"{base}_median"].to_numpy()
            lo = group[f"{base}_q16"].to_numpy()
            hi = group[f"{base}_q84"].to_numpy()
            if color is None:
                line, = ax.plot(
                    x,
                    y,
                    linestyle=linestyle,
                    marker="o",
                    label=rf"$\beta_G={beta:g}$ {label_name}",
                )
                color = line.get_color()
            else:
                line, = ax.plot(
                    x,
                    y,
                    linestyle=linestyle,
                    marker="o",
                    color=color,
                    label=rf"$\beta_G={beta:g}$ {label_name}",
                )
            ax.fill_between(x, lo, hi, alpha=0.12, color=color)
    ax.axhline(0.0, linewidth=1.0)
    ax.set_xlabel(r"$a_{\rm rel,ang}/\sigma$")
    ax.set_ylabel(ylabel)
    ax.legend(frameon=False, fontsize="small", ncol=2)
    fig.tight_layout()
    paths = _save_figure(fig, output_dir, f"{prefix}_{stem}", dpi=dpi)
    plt.close(fig)
    return paths


def make_transition_plots(
    summary: pd.DataFrame,
    output_dir: str | Path,
    *,
    prefix: str = "transition",
    dpi: int = 300,
) -> list[Path]:
    """Create paired transition diagnostics as PNG and vector PDF files."""
    output_dir = Path(output_dir)
    paths: list[Path] = []
    paths.extend(_plot_delta_chi2(summary, output_dir, prefix, dpi))

    specs = (
        ("pct_bias_beta_g", r"$100(\hat\beta_G-\beta_G)/\beta_G$ [\%]", "beta_fractional_bias"),
        ("pct_bias_parallax_mas", r"$100(\hat\varpi-\varpi)/\varpi$ [\%]", "parallax_fractional_bias"),
        ("pct_bias_m1_msun", r"$100(\hat M_1-M_1)/M_1$ [\%]", "m1_fractional_bias"),
        ("pct_bias_m2_msun", r"$100(\hat M_2-M_2)/M_2$ [\%]", "m2_fractional_bias"),
        ("bias_inclination_deg", r"$\hat i-i$ [deg]", "inclination_bias"),
    )
    for metric, ylabel, stem in specs:
        paths.extend(
            _plot_parameter_bias(
                summary,
                output_dir,
                prefix=prefix,
                metric=metric,
                ylabel=ylabel,
                stem=stem,
                dpi=dpi,
            )
        )
    return paths


def compact_delta_chi2_table(summary: pd.DataFrame) -> pd.DataFrame:
    """Return the main paired model-comparison columns for terminal display."""
    columns = [
        "true_beta_g",
        "a_over_sigma",
        "n_pairs",
        "delta_chi2_median",
        "delta_chi2_mean",
        "delta_chi2_q16",
        "delta_chi2_q84",
        "delta_chi2_min",
        "delta_chi2_max",
    ]
    return summary[columns].copy()


def _validate_regime_grid(summary: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    required = {
        "true_beta_g",
        "a_over_sigma",
        "delta_chi2_median",
        "delta_chi2_q16",
        "delta_chi2_min",
    }
    missing = required - set(summary.columns)
    if missing:
        raise ValueError(f"summary is missing regime-map columns: {sorted(missing)}")
    if summary.duplicated(["true_beta_g", "a_over_sigma"]).any():
        raise ValueError("regime map summary contains duplicate grid points")
    betas = np.sort(summary["true_beta_g"].unique().astype(float))
    ratios = np.sort(summary["a_over_sigma"].unique().astype(float))
    expected = len(betas) * len(ratios)
    if len(summary) != expected:
        raise ValueError(
            "regime map requires a complete rectangular (beta_G, a/sigma) grid; "
            f"expected {expected} rows, found {len(summary)}"
        )
    return betas, ratios


def _first_ratio(group: pd.DataFrame, condition: pd.Series) -> float:
    selected = group.loc[condition].sort_values("a_over_sigma")
    if selected.empty:
        return np.nan
    return float(selected.iloc[0]["a_over_sigma"])


def build_regime_boundary_table(
    summary: pd.DataFrame,
    raw: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Build sampled descriptive boundaries for the 2-D response regime."""
    betas, _ = _validate_regime_grid(summary)
    multi = None
    if raw is not None and "gaia_multi_peak_fraction" in raw.columns:
        multi = (
            raw.groupby(["true_beta_g", "a_over_sigma"], as_index=False)
            .agg(gaia_multi_peak_fraction=("gaia_multi_peak_fraction", "max"))
        )

    rows = []
    for beta in betas:
        group = summary[summary["true_beta_g"] == beta].sort_values("a_over_sigma")
        row = {
            "true_beta_g": float(beta),
            "q16_positive_a_over_sigma": _first_ratio(
                group, group["delta_chi2_q16"] > 0.0
            ),
            "all_positive_a_over_sigma": _first_ratio(
                group, group["delta_chi2_min"] > 0.0
            ),
            "median_gt10_a_over_sigma": _first_ratio(
                group, group["delta_chi2_median"] > 10.0
            ),
            "multi_peak_onset_a_over_sigma": np.nan,
        }
        if multi is not None:
            mgroup = multi[multi["true_beta_g"] == beta].sort_values("a_over_sigma")
            if not mgroup.empty:
                row["multi_peak_onset_a_over_sigma"] = _first_ratio(
                    mgroup, mgroup["gaia_multi_peak_fraction"] > 0.0
                )
        rows.append(row)
    return pd.DataFrame(rows)


def write_regime_boundary_table(
    boundary: pd.DataFrame,
    output_dir: str | Path,
    *,
    prefix: str = "transition",
) -> Path:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f"{prefix}_regime_boundaries.csv"
    boundary.to_csv(path, index=False)
    return path


def _grid_edges(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if len(values) == 1:
        return np.array([values[0] - 0.5, values[0] + 0.5])
    mids = 0.5 * (values[:-1] + values[1:])
    return np.concatenate((
        [values[0] - 0.5 * (values[1] - values[0])],
        mids,
        [values[-1] + 0.5 * (values[-1] - values[-2])],
    ))


def make_regime_map(
    summary: pd.DataFrame,
    output_dir: str | Path,
    *,
    raw: pd.DataFrame | None = None,
    prefix: str = "transition",
    dpi: int = 300,
) -> list[Path]:
    """Create the 2-D single-peak photocentre-mismatch regime map.

    The heatmap uses log10(1 + max(median Delta chi2, 0)); negative medians are
    plotted at zero rather than passed to a logarithm. Curves mark the first
    sampled a/sigma with q16>0, all seeds >0, and, when present in the raw scan,
    any multi-peak surrogate epoch. These are descriptive grid boundaries, not
    formal significance levels or Gaia instrumental resolution thresholds.
    """
    import matplotlib.pyplot as plt

    betas, ratios = _validate_regime_grid(summary)
    pivot = summary.pivot(
        index="true_beta_g", columns="a_over_sigma", values="delta_chi2_median"
    ).reindex(index=betas, columns=ratios)
    if pivot.isna().any().any():
        raise ValueError("regime map contains missing median Delta chi2 values")

    values = np.log10(1.0 + np.clip(pivot.to_numpy(dtype=float), 0.0, None))
    boundary = build_regime_boundary_table(summary, raw=raw)

    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    mesh = ax.pcolormesh(
        _grid_edges(ratios),
        _grid_edges(betas),
        values,
        shading="flat",
    )
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label(r"$\log_{10}[1+\max(\mathrm{median}\ \Delta\chi^2,0)]$")

    for column, label, linestyle in (
        ("q16_positive_a_over_sigma", r"$Q_{16}(\Delta\chi^2)>0$", "-"),
        ("all_positive_a_over_sigma", r"all seeds: $\Delta\chi^2>0$", "--"),
        ("multi_peak_onset_a_over_sigma", "multi-peak onset", ":"),
    ):
        points = boundary[["true_beta_g", column]].dropna()
        if len(points) >= 2:
            ax.plot(
                points[column].to_numpy(),
                points["true_beta_g"].to_numpy(),
                marker="o",
                linestyle=linestyle,
                label=label,
            )
        elif len(points) == 1:
            ax.scatter(
                points[column].to_numpy(),
                points["true_beta_g"].to_numpy(),
                label=label,
            )

    ax.set_xlabel(r"$a_{\rm rel,ang}/\sigma$")
    ax.set_ylabel(r"secondary light fraction $\beta_G$")
    ax.set_title("Single-peak photocentre-mismatch regime")
    handles, _ = ax.get_legend_handles_labels()
    if handles:
        ax.legend(frameon=False, fontsize="small")
    fig.tight_layout()

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = [
        output_dir / f"{prefix}_regime_map.png",
        output_dir / f"{prefix}_regime_map.pdf",
    ]
    fig.savefig(paths[0], dpi=dpi, bbox_inches="tight")
    fig.savefig(paths[1], bbox_inches="tight")
    plt.close(fig)
    return paths
