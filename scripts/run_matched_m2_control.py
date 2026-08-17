#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from raa_orbit_model.experiment_config import add_orbit_arguments, truth_from_args
from raa_orbit_model.experiments import write_rows_csv
from raa_orbit_model.response_controls import matched_m2_scan, summarise_matched_m2
from raa_orbit_model.scanning import schedule_from_csv, schedule_from_gaiascanlaw


def main():
    parser = argparse.ArgumentParser(
        description=(
            "High-statistics matched-M2 injection/recovery control. Fits only the "
            "correct Penoyre-style response to test whether small residual parameter "
            "offsets persist over many independent noise/time realizations."
        )
    )
    parser.add_argument("--schedule-file", help="archived Gaia schedule CSV")
    parser.add_argument("--ra-deg", type=float)
    parser.add_argument("--dec-deg", type=float)
    parser.add_argument("--release", choices=("dr1", "dr2", "dr3", "dr4", "dr5"), default="dr4")
    parser.add_argument("--output", default="results/matched_m2_100seed.csv")
    parser.add_argument("--summary-output", default=None)
    parser.add_argument("--seeds", type=int, default=100)
    parser.add_argument("--alpha-mas", type=float, default=50.0)
    parser.add_argument("--a-over-alpha", type=float, default=1.0)
    parser.add_argument("--beta-g", type=float, default=0.25)
    parser.add_argument("--beta-over-alpha", type=float, default=1.5)
    parser.add_argument("--n-ast", type=int, default=24)
    parser.add_argument("--n-rv", type=int, default=48)
    parser.add_argument("--ast-sigma-mas", type=float, default=0.20)
    parser.add_argument("--rv-sigma-kms", type=float, default=0.10)
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
    if args.n_ast < 0 or args.n_rv < 0:
        parser.error("--n-ast and --n-rv must be >= 0")
    if min(args.ast_sigma_mas, args.rv_sigma_kms, args.gaia_sigma_mas) <= 0:
        parser.error("all uncertainties must be > 0")

    if args.schedule_file:
        schedule = schedule_from_csv(args.schedule_file)
    else:
        if args.ra_deg is None or args.dec_deg is None:
            parser.error("--schedule-file or both --ra-deg/--dec-deg are required")
        schedule = schedule_from_gaiascanlaw(args.ra_deg, args.dec_deg, release=args.release)

    truth = truth_from_args(args, parser)
    print(
        f"schedule N={schedule.n_transits}; matched M2 point: beta_G={args.beta_g:g}, "
        f"a/alpha={args.a_over_alpha:g}, beta_PSF/alpha={args.beta_over_alpha:g}"
    )
    print(
        f"external data: N_ast={args.n_ast}, sigma_ast={args.ast_sigma_mas:g} mas; "
        f"N_RV={args.n_rv}, sigma_RV={args.rv_sigma_kms:g} km/s; "
        f"sigma_Gaia={args.gaia_sigma_mas:g} mas"
    )
    print(f"seeds 0..{args.seeds-1}; expected fit records: {args.seeds}")

    rows = matched_m2_scan(
        truth,
        schedule,
        seeds=tuple(range(args.seeds)),
        alpha_mas=float(args.alpha_mas),
        a_over_alpha=float(args.a_over_alpha),
        beta_g=float(args.beta_g),
        beta_over_alpha=float(args.beta_over_alpha),
        n_ast=int(args.n_ast),
        n_rv=int(args.n_rv),
        ast_sigma_mas=float(args.ast_sigma_mas),
        rv_sigma_kms=float(args.rv_sigma_kms),
        gaia_sigma_mas=float(args.gaia_sigma_mas),
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    write_rows_csv(rows, str(output))
    print(f"wrote {len(rows)} fit records to {output}")

    summary_path = Path(args.summary_output) if args.summary_output else output.with_name(output.stem + "_summary.csv")
    summary = summarise_matched_m2(rows)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    write_rows_csv(summary, str(summary_path))
    print(f"wrote {len(summary)} summary row to {summary_path}")


if __name__ == "__main__":
    main()
