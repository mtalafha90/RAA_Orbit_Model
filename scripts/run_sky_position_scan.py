#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from raa_orbit_model.experiments import write_rows_csv
from raa_orbit_model.kepler import BinaryParams
from raa_orbit_model.scanning import schedule_from_csv, schedule_from_gaiascanlaw, write_schedule_csv
from raa_orbit_model.experiment_config import (
    add_noise_arguments,
    add_orbit_arguments,
    noise_kwargs_from_args,
    truth_from_args,
)
from raa_orbit_model.sky_experiments import sky_resolution_bias_scan
from raa_orbit_model.sky_study import (
    DEFAULT_ECLIPTIC_LATITUDES_DEG,
    DEFAULT_ECLIPTIC_LONGITUDES_DEG,
    FULL_SKY_ECLIPTIC_LATITUDES_DEG,
    FULL_SKY_POSITION_COUNT,
    augment_bias_rows,
    build_ecliptic_grid,
    position_record,
    scan_geometry_diagnostics,
)

DEFAULT_A_OVER_SIGMA = (0.40, 0.50, 0.60, 0.80, 1.00)
DEFAULT_BETA_VALUES = (0.05, 0.15, 0.25, 0.35, 0.45)


def _write_dict_rows(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _read_rows(path: Path) -> list[dict]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def _release_window_yr(release: str) -> float:
    try:
        import gaiascanlaw as gsl
    except ImportError as exc:
        raise ImportError("sky-position runs require the gaiascanlaw package") from exc
    end_name = f"t{release}"
    if not hasattr(gsl, "tstart") or not hasattr(gsl, end_name):
        raise RuntimeError("installed gaiascanlaw does not expose the expected release time constants")
    return float(getattr(gsl, end_name) - gsl.tstart)


def _validate_existing_schedule(schedule, position, release: str) -> None:
    if schedule.release != release:
        raise ValueError(f"existing schedule {position.sky_id} has release={schedule.release}, expected {release}")
    if not np.isclose(schedule.ra_deg, position.ra_deg, atol=1e-8) or not np.isclose(schedule.dec_deg, position.dec_deg, atol=1e-8):
        raise ValueError(f"existing schedule {position.sky_id} does not match the requested ICRS position")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the RAA transition experiment over a controlled Gaia sky-position grid")
    parser.add_argument("--output-dir", default="results/sky_position_dr4")
    parser.add_argument("--release", choices=("dr1", "dr2", "dr3", "dr4", "dr5"), default="dr4")
    parser.add_argument("--ecliptic-lat-values", nargs="+", type=float, default=DEFAULT_ECLIPTIC_LATITUDES_DEG,
                        help="default is the superseded 25-position northern pilot; see --full-sky")
    parser.add_argument("--ecliptic-lon-values", nargs="+", type=float, default=DEFAULT_ECLIPTIC_LONGITUDES_DEG)
    parser.add_argument("--full-sky", action="store_true",
                        help="use the published pole-to-pole latitude grid, giving the "
                             "46 sky positions reported in the manuscript; overrides "
                             "--ecliptic-lat-values")
    parser.add_argument("--a-over-sigma-values", nargs="+", type=float, default=DEFAULT_A_OVER_SIGMA)
    parser.add_argument("--beta-values", nargs="+", type=float, default=DEFAULT_BETA_VALUES)
    parser.add_argument("--seeds", type=int, default=10)
    parser.add_argument("--sigma-response-mas", type=float, default=50.0)
    parser.add_argument("--entropy-bins", type=int, default=12)
    parser.add_argument("--external-baseline-yr", type=float, default=None, help="fixed resolved-astrometry/RV baseline; default is the gaiascanlaw release window")
    parser.add_argument("--apply-astrometry-gaps", action="store_true")
    parser.add_argument("--schedule-only", action="store_true", help="generate/archive schedules and scan diagnostics without fitting")
    parser.add_argument("--overwrite", action="store_true", help="regenerate schedules and rerun completed per-position fit files")
    add_orbit_arguments(parser)
    add_noise_arguments(parser)
    args = parser.parse_args()

    if args.seeds <= 0:
        parser.error("--seeds must be > 0")
    if args.sigma_response_mas <= 0:
        parser.error("--sigma-response-mas must be > 0")
    if args.entropy_bins < 2:
        parser.error("--entropy-bins must be >= 2")
    if any(v <= 0 for v in args.a_over_sigma_values):
        parser.error("--a-over-sigma-values must all be > 0")
    if any(v < 0 or v > 0.5 for v in args.beta_values):
        parser.error("--beta-values must all lie in [0, 0.5]")
    if args.external_baseline_yr is not None and args.external_baseline_yr <= 0:
        parser.error("--external-baseline-yr must be > 0")

    output_dir = Path(args.output_dir)
    schedule_dir = output_dir / "schedules"
    position_dir = output_dir / "per_position"
    schedule_dir.mkdir(parents=True, exist_ok=True)
    position_dir.mkdir(parents=True, exist_ok=True)

    latitudes = FULL_SKY_ECLIPTIC_LATITUDES_DEG if args.full_sky else args.ecliptic_lat_values
    positions = build_ecliptic_grid(latitudes, args.ecliptic_lon_values)
    print(f"sky positions: {len(positions)}")
    if args.full_sky and len(positions) != FULL_SKY_POSITION_COUNT:
        parser.error(
            f"--full-sky expects {FULL_SKY_POSITION_COUNT} positions but the requested "
            f"longitude grid gives {len(positions)}; the published run uses the default longitudes"
        )

    resolved = {}
    scan_rows = []
    for index, position in enumerate(positions, start=1):
        schedule_path = schedule_dir / f"{position.sky_id}_{args.release}.csv"
        if schedule_path.exists() and not args.overwrite:
            schedule = schedule_from_csv(schedule_path)
            _validate_existing_schedule(schedule, position, args.release)
        else:
            schedule = schedule_from_gaiascanlaw(
                position.ra_deg,
                position.dec_deg,
                release=args.release,
                obstype="astrometry" if args.apply_astrometry_gaps else None,
            )
            write_schedule_csv(schedule, schedule_path)
        diagnostics = scan_geometry_diagnostics(schedule, entropy_bins=args.entropy_bins)
        resolved[position.sky_id] = (position, schedule, diagnostics)
        scan_rows.append(position_record(position, diagnostics))
        print(
            f"[{index:02d}/{len(positions):02d}] {position.sky_id}: "
            f"RA={position.ra_deg:.6f} Dec={position.dec_deg:.6f}; "
            f"N={schedule.n_transits}; max_gap={diagnostics.directional_max_gap_deg:.2f} deg; "
            f"Hnorm={diagnostics.directional_entropy_normalized:.3f}"
        )

    _write_dict_rows(output_dir / "sky_grid.csv", [position_record(p, resolved[p.sky_id][2]) for p in positions])
    _write_dict_rows(output_dir / "scan_geometry.csv", scan_rows)
    if args.schedule_only:
        print("schedule-only mode: no orbit fits were run")
        return

    external_baseline_yr = args.external_baseline_yr
    if external_baseline_yr is None:
        external_baseline_yr = _release_window_yr(args.release)
    print(f"fixed external astrometry/RV baseline: {external_baseline_yr:.9f} yr")

    truth = truth_from_args(args, parser)
    noise_kwargs = noise_kwargs_from_args(args, parser)
    print(
        f"orbit: P={truth.period_yr} yr e={truth.eccentricity} i={truth.inclination_deg} deg "
        f"M1={truth.m1_msun} M2={truth.m2_msun} Msun (q={truth.q:.4f})"
    )
    expected_per_position = len(args.a_over_sigma_values) * len(args.beta_values) * args.seeds * 2
    expected_total = expected_per_position * len(positions)
    print(f"expected fit records: {expected_total} ({expected_per_position} per sky position)")

    for index, position in enumerate(positions, start=1):
        _, schedule, diagnostics = resolved[position.sky_id]
        result_path = position_dir / f"{position.sky_id}_{args.release}_bias.csv"
        if result_path.exists() and not args.overwrite:
            n_existing = len(_read_rows(result_path))
            if n_existing != expected_per_position:
                raise RuntimeError(
                    f"incomplete existing result file {result_path}: {n_existing} rows; "
                    f"expected {expected_per_position}. Use --overwrite to rerun it."
                )
            print(f"[{index:02d}/{len(positions):02d}] resume: {position.sky_id} already complete")
            continue

        rows = sky_resolution_bias_scan(
            truth,
            schedule,
            external_baseline_yr=external_baseline_yr,
            a_over_sigma_values=tuple(args.a_over_sigma_values),
            beta_values=tuple(args.beta_values),
            seeds=tuple(range(args.seeds)),
            sigma_response_mas=args.sigma_response_mas,
            **noise_kwargs,
        )
        augment_bias_rows(rows, position, diagnostics)
        write_rows_csv(rows, str(result_path))
        print(f"[{index:02d}/{len(positions):02d}] wrote {len(rows)} rows: {result_path}")

    combined_rows = []
    for position in positions:
        result_path = position_dir / f"{position.sky_id}_{args.release}_bias.csv"
        rows = _read_rows(result_path)
        if len(rows) != expected_per_position:
            raise RuntimeError(f"cannot combine incomplete file: {result_path}")
        combined_rows.extend(rows)
    combined_path = output_dir / "sky_position_bias.csv"
    _write_dict_rows(combined_path, combined_rows)
    print(f"combined fit records: {len(combined_rows)}")
    print(f"combined output: {combined_path}")
    print(f"scan diagnostics: {output_dir / 'scan_geometry.csv'}")


if __name__ == "__main__":
    main()
