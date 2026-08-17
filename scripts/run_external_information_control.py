#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from raa_orbit_model.experiment_config import add_orbit_arguments, truth_from_args
from raa_orbit_model.experiments import write_rows_csv
from raa_orbit_model.response_controls import (
    DEFAULT_EXTERNAL_LEVELS,
    external_information_strength_scan,
    summarise_external_information,
)
from raa_orbit_model.scanning import schedule_from_csv, schedule_from_gaiascanlaw


def main():
    parser = argparse.ArgumentParser(
        description=(
            "At one difficult M2 injection, vary the precision of resolved astrometry and "
            "SB2 RVs while fitting M0/M1/M2. The default strong/medium/weak levels are "
            "(0.2 mas,0.1 km/s), (1 mas,0.5 km/s), and (2 mas,1 km/s)."
        )
    )
    parser.add_argument("--schedule-file", help="archived Gaia schedule CSV")
    parser.add_argument("--ra-deg", type=float)
    parser.add_argument("--dec-deg", type=float)
    parser.add_argument("--release", choices=("dr1", "dr2", "dr3", "dr4", "dr5"), default="dr4")
    parser.add_argument("--output", default="results/external_information_control.csv")
    parser.add_argument("--summary-output", default=None)
    parser.add_argument("--seeds", type=int, default=30)
    parser.add_argument("--alpha-mas", type=float, default=50.0)
    parser.add_argument("--a-over-alpha", type=float, default=1.0)
    parser.add_argument("--beta-g", type=float, default=0.25)
    parser.add_argument("--beta-over-alpha", type=float, default=1.5)
    parser.add_argument("--gaia-sigma-mas", type=float, default=0.10)
    add_orbit_arguments(parser)
    args = parser.parse_args()

    if args.seeds <= 0:
        parser.error("--seeds must be > 0")
    if args.alpha_mas <= 0 or args.a_over_alpha <= 0:
        parser.error("--alpha-mas and --a-over-alpha must be > 0")
    if not (0 <= args.beta_g <= 0.5):
        parser.error("--beta-g must be in [0,0.5]")
    if args.beta_over_alpha < 1:
        parser.error("--beta-over-alpha must be >= 1")
    if args.gaia_sigma_mas <= 0:
        parser.error("--gaia-sigma-mas must be > 0")

    if args.schedule_file:
        schedule = schedule_from_csv(args.schedule_file)
    else:
        if args.ra_deg is None or args.dec_deg is None:
            parser.error("--schedule-file or both --ra-deg/--dec-deg are required")
        schedule = schedule_from_gaiascanlaw(args.ra_deg, args.dec_deg, release=args.release)

    truth = truth_from_args(args, parser)
    expected = len(DEFAULT_EXTERNAL_LEVELS) * args.seeds * 3
    print(
        f"schedule N={schedule.n_transits}; central point: beta_G={args.beta_g:g}, "
        f"a/alpha={args.a_over_alpha:g}, beta_PSF/alpha={args.beta_over_alpha:g}"
    )
    for level in DEFAULT_EXTERNAL_LEVELS:
        print(
            f"{level.name}: N_ast={level.n_ast}, sigma_ast={level.ast_sigma_mas:g} mas; "
            f"N_RV={level.n_rv}, sigma_RV={level.rv_sigma_kms:g} km/s"
        )
    print(f"seeds 0..{args.seeds-1}; expected fit records: {expected}")

    rows = external_information_strength_scan(
        truth,
        schedule,
        seeds=tuple(range(args.seeds)),
        alpha_mas=float(args.alpha_mas),
        a_over_alpha=float(args.a_over_alpha),
        beta_g=float(args.beta_g),
        beta_over_alpha=float(args.beta_over_alpha),
        gaia_sigma_mas=float(args.gaia_sigma_mas),
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    write_rows_csv(rows, str(output))
    print(f"wrote {len(rows)} fit records to {output}")

    summary_path = Path(args.summary_output) if args.summary_output else output.with_name(output.stem + "_summary.csv")
    summary = summarise_external_information(rows)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    write_rows_csv(summary, str(summary_path))
    print(f"wrote {len(summary)} summary rows to {summary_path}")


if __name__ == "__main__":
    main()
