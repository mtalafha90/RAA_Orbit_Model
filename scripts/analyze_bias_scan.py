#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from raa_orbit_model.bias_analysis import (
    build_paired_results,
    compact_delta_chi2_table,
    load_bias_results,
    make_transition_plots,
    summarize_paired_results,
    write_analysis_tables,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Pair photocentre and resolution-aware fits by injection seed, summarize "
            "the transition, and generate publication-ready diagnostics."
        )
    )
    parser.add_argument("input", help="bias-scan CSV produced by scripts/run_bias_scan.py")
    parser.add_argument(
        "--output-dir",
        default=None,
        help="analysis directory; default is <input stem>_analysis beside the input CSV",
    )
    parser.add_argument("--prefix", default="transition", help="output filename prefix")
    parser.add_argument("--dpi", type=int, default=300, help="PNG resolution")
    parser.add_argument("--no-plots", action="store_true", help="write paired/summary CSVs only")
    args = parser.parse_args()

    input_path = Path(args.input)
    if args.dpi <= 0:
        parser.error("--dpi must be > 0")
    output_dir = (
        Path(args.output_dir)
        if args.output_dir is not None
        else input_path.with_name(f"{input_path.stem}_analysis")
    )

    df = load_bias_results(input_path)
    paired = build_paired_results(df)
    summary = summarize_paired_results(paired)
    paired_path, summary_path = write_analysis_tables(
        paired,
        summary,
        output_dir,
        prefix=args.prefix,
    )

    print(compact_delta_chi2_table(summary).to_string(index=False))
    print()
    print(f"paired records: {paired_path}")
    print(f"summary table:  {summary_path}")

    if not args.no_plots:
        plot_paths = make_transition_plots(
            summary,
            output_dir,
            prefix=args.prefix,
            dpi=args.dpi,
        )
        for path in plot_paths:
            print(f"plot:           {path}")


if __name__ == "__main__":
    main()
