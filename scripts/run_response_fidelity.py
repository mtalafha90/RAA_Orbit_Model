#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from raa_orbit_model.experiment_config import (
    add_noise_arguments,
    add_orbit_arguments,
    noise_kwargs_from_args,
    truth_from_args,
)
from raa_orbit_model.experiments import write_rows_csv
from raa_orbit_model.response_fidelity import (
    response_fidelity_scan,
    summarise_response_fidelity,
)
from raa_orbit_model.scanning import (
    schedule_from_csv,
    schedule_from_gaiascanlaw,
    write_schedule_csv,
)


DEFAULT_A_OVER_ALPHA = (0.4, 0.6, 0.8, 1.0)
DEFAULT_BETA_VALUES = (0.05, 0.25, 0.45)
DEFAULT_ELONGATION = (3.0,)


def _positive(parser, values, name):
    out = tuple(float(v) for v in values)
    if not out or any(v <= 0 for v in out):
        parser.error(f"{name} must contain one or more values > 0")
    return out


def _light_fractions(parser, values):
    out = tuple(float(v) for v in values)
    if not out or any(v < 0 or v > 0.5 for v in out):
        parser.error("--beta-values must contain one or more values in [0, 0.5]")
    return out


def _elongations(parser, values):
    out = tuple(float(v) for v in values)
    if not out or any(v < 1 for v in out):
        parser.error("--beta-over-alpha-values must contain values >= 1")
    return out


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Inject the finite-elongation Penoyre-style response (M2) and fit the same "
            "synthetic realization with M0 photocentre, M1 equal-width 1-D response, "
            "and M2. Widths are research-surrogate parameters, not Gaia calibration values."
        )
    )
    parser.add_argument("--output", default="results/response_fidelity.csv")
    parser.add_argument("--summary-output", default=None)
    parser.add_argument("--seeds", type=int, default=3, help="number of seeds starting at zero")
    parser.add_argument(
        "--alpha-mas", type=float, default=50.0,
        help="idealised AL/narrow Gaussian width alpha; research scale only",
    )
    parser.add_argument(
        "--beta-over-alpha-values", nargs="+", type=float, default=DEFAULT_ELONGATION,
        metavar="E", help="idealised PSF elongation beta/alpha >=1; default: 3",
    )
    parser.add_argument(
        "--a-over-alpha-values", nargs="+", type=float, default=DEFAULT_A_OVER_ALPHA,
        metavar="R", help="angular relative semimajor axis divided by alpha",
    )
    parser.add_argument(
        "--beta-values", nargs="+", type=float, default=DEFAULT_BETA_VALUES,
        metavar="BETA", help="secondary Gaia-band light fractions in [0, 0.5]",
    )
    parser.add_argument("--ra-deg", type=float)
    parser.add_argument("--dec-deg", type=float)
    parser.add_argument("--release", choices=("dr1", "dr2", "dr3", "dr4", "dr5"), default="dr4")
    parser.add_argument("--schedule-file", help="reuse an exact schedule CSV archived by this project")
    parser.add_argument("--write-schedule", help="archive the resolved Gaia schedule to this CSV")
    parser.add_argument(
        "--apply-astrometry-gaps", action="store_true",
        help="apply gaiascanlaw's published astrometric gap mask when generating a schedule",
    )
    add_orbit_arguments(parser)
    add_noise_arguments(parser)
    args = parser.parse_args()

    if args.seeds <= 0:
        parser.error("--seeds must be > 0")
    if args.alpha_mas <= 0:
        parser.error("--alpha-mas must be > 0")
    a_values = _positive(parser, args.a_over_alpha_values, "--a-over-alpha-values")
    beta_values = _light_fractions(parser, args.beta_values)
    elongations = _elongations(parser, args.beta_over_alpha_values)

    if args.schedule_file:
        schedule = schedule_from_csv(args.schedule_file)
    else:
        if args.ra_deg is None or args.dec_deg is None:
            parser.error("--ra-deg and --dec-deg are required unless --schedule-file is supplied")
        schedule = schedule_from_gaiascanlaw(
            args.ra_deg,
            args.dec_deg,
            release=args.release,
            obstype="astrometry" if args.apply_astrometry_gaps else None,
        )
    if args.write_schedule:
        write_schedule_csv(schedule, args.write_schedule)

    truth = truth_from_args(args, parser)
    noise = noise_kwargs_from_args(args, parser)
    n_expected = len(a_values) * len(beta_values) * len(elongations) * args.seeds * 3
    print(
        f"schedule: N={schedule.n_transits} RA={schedule.ra_deg} Dec={schedule.dec_deg} "
        f"release={schedule.release}; {schedule.source}"
    )
    print("M2 injection: finite-elongation Penoyre-style Gaussian surrogate")
    print(f"alpha={args.alpha_mas:g} mas (research scale only)")
    print("beta/alpha:", " ".join(f"{x:g}" for x in elongations))
    print("a/alpha:", " ".join(f"{x:g}" for x in a_values))
    print("beta_G:", " ".join(f"{x:g}" for x in beta_values))
    print(f"seeds: 0..{args.seeds - 1}; expected fit records: {n_expected}")

    rows = response_fidelity_scan(
        truth,
        schedule,
        a_over_alpha_values=a_values,
        beta_values=beta_values,
        beta_over_alpha_values=elongations,
        seeds=tuple(range(args.seeds)),
        alpha_mas=float(args.alpha_mas),
        **noise,
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    write_rows_csv(rows, str(output))
    print(f"wrote {len(rows)} fit records to {output}")

    summary_path = (
        Path(args.summary_output)
        if args.summary_output
        else output.with_name(output.stem + "_summary.csv")
    )
    summary = summarise_response_fidelity(rows)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    write_rows_csv(summary, str(summary_path))
    print(f"wrote {len(summary)} summary rows to {summary_path}")


if __name__ == "__main__":
    main()
