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
from raa_orbit_model.experiments import resolution_bias_scan, write_rows_csv
from raa_orbit_model.scanning import (
    schedule_from_csv,
    schedule_from_gaiascanlaw,
    write_schedule_csv,
)


DEFAULT_A_OVER_SIGMA = (0.05, 0.10, 0.20, 0.50, 1.00, 2.00, 4.00)
DEFAULT_BETA_VALUES = (0.05, 0.20, 0.40)


def _positive_values(parser: argparse.ArgumentParser, values, name: str) -> tuple[float, ...]:
    values = tuple(float(v) for v in values)
    if not values or any(v <= 0 for v in values):
        parser.error(f"{name} must contain one or more values > 0")
    return values


def _beta_values(parser: argparse.ArgumentParser, values) -> tuple[float, ...]:
    values = tuple(float(v) for v in values)
    if not values or any(v < 0 or v > 0.5 for v in values):
        parser.error("--beta-values must contain one or more values in [0, 0.5]")
    return values


def main():
    parser = argparse.ArgumentParser(
        description="Run the RAA photocentre-vs-resolution experiment on a Gaia scanning-law schedule"
    )
    parser.add_argument("--output", default="bias_scan.csv")
    parser.add_argument("--seeds", type=int, default=3, help="number of seeds, starting at zero")
    parser.add_argument("--sigma-response-mas", type=float, default=50.0)
    parser.add_argument(
        "--a-over-sigma-values",
        nargs="+",
        type=float,
        default=DEFAULT_A_OVER_SIGMA,
        metavar="R",
        help=(
            "dimensionless angular-separation grid a_rel,ang/sigma. "
            "Default: 0.05 0.10 0.20 0.50 1 2 4"
        ),
    )
    parser.add_argument(
        "--beta-values",
        nargs="+",
        type=float,
        default=DEFAULT_BETA_VALUES,
        metavar="BETA",
        help="secondary Gaia-band light fractions in [0, 0.5]. Default: 0.05 0.20 0.40",
    )
    parser.add_argument("--ra-deg", type=float, help="ICRS right ascension; required unless --schedule-file is used")
    parser.add_argument("--dec-deg", type=float, help="ICRS declination; required unless --schedule-file is used")
    parser.add_argument("--release", choices=("dr1", "dr2", "dr3", "dr4", "dr5"), default="dr4")
    parser.add_argument(
        "--schedule-file",
        help="reuse an exact schedule CSV previously written by this project instead of querying gaiascanlaw",
    )
    parser.add_argument(
        "--write-schedule",
        help="write the exact resolved Gaia schedule to this CSV for reproducibility",
    )
    parser.add_argument(
        "--apply-astrometry-gaps",
        action="store_true",
        help=(
            "apply gaiascanlaw's published astrometric gap mask; the default is the "
            "uninterrupted nominal scanning law used for mission simulations"
        ),
    )
    add_orbit_arguments(parser)
    add_noise_arguments(parser)
    args = parser.parse_args()

    if args.seeds <= 0:
        parser.error("--seeds must be > 0")
    if args.sigma_response_mas <= 0:
        parser.error("--sigma-response-mas must be > 0")
    a_over_sigma_values = _positive_values(parser, args.a_over_sigma_values, "--a-over-sigma-values")
    beta_values = _beta_values(parser, args.beta_values)

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

    print(
        f"Gaia schedule: {schedule.n_transits} transits; "
        f"RA={schedule.ra_deg:.6f} deg Dec={schedule.dec_deg:.6f} deg; "
        f"release={schedule.release}; source={schedule.source}"
    )
    print("a/sigma grid:", " ".join(f"{v:g}" for v in a_over_sigma_values))
    print("beta_G grid:", " ".join(f"{v:g}" for v in beta_values))
    print(f"seeds: 0..{args.seeds - 1}")

    truth = truth_from_args(args, parser)
    print(
        f"orbit: P={truth.period_yr} yr e={truth.eccentricity} i={truth.inclination_deg} deg "
        f"M1={truth.m1_msun} M2={truth.m2_msun} Msun (q={truth.q:.4f})"
    )
    rows = resolution_bias_scan(
        truth,
        schedule,
        a_over_sigma_values=a_over_sigma_values,
        beta_values=beta_values,
        seeds=tuple(range(args.seeds)),
        sigma_response_mas=args.sigma_response_mas,
        **noise_kwargs_from_args(args, parser),
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    write_rows_csv(rows, str(output))
    print(f"wrote {len(rows)} fit records to {output}")


if __name__ == "__main__":
    main()
