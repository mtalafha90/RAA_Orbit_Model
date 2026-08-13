#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from raa_orbit_model.experiments import resolution_bias_scan, write_rows_csv
from raa_orbit_model.kepler import BinaryParams
from raa_orbit_model.scanning import (
    schedule_from_csv,
    schedule_from_gaiascanlaw,
    write_schedule_csv,
)


def main():
    parser = argparse.ArgumentParser(
        description="Run the RAA photocentre-vs-resolution experiment on a Gaia scanning-law schedule"
    )
    parser.add_argument("--output", default="bias_scan.csv")
    parser.add_argument("--seeds", type=int, default=3, help="number of seeds, starting at zero")
    parser.add_argument("--sigma-response-mas", type=float, default=50.0)
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
        "--nominal-no-gaps",
        action="store_true",
        help="use the uninterrupted nominal law; default applies the gaiascanlaw astrometric gap mask",
    )
    args = parser.parse_args()

    if args.schedule_file:
        schedule = schedule_from_csv(args.schedule_file)
    else:
        if args.ra_deg is None or args.dec_deg is None:
            parser.error("--ra-deg and --dec-deg are required unless --schedule-file is supplied")
        schedule = schedule_from_gaiascanlaw(
            args.ra_deg,
            args.dec_deg,
            release=args.release,
            obstype=None if args.nominal_no_gaps else "astrometry",
        )

    if args.write_schedule:
        write_schedule_csv(schedule, args.write_schedule)

    print(
        f"Gaia schedule: {schedule.n_transits} transits; "
        f"RA={schedule.ra_deg:.6f} deg Dec={schedule.dec_deg:.6f} deg; "
        f"release={schedule.release}; source={schedule.source}"
    )

    truth = BinaryParams(
        period_yr=2.0,
        t_peri_yr=0.15,
        eccentricity=0.25,
        inclination_deg=72.0,
        omega_deg=55.0,
        node_deg=120.0,
        m1_msun=1.25,
        m2_msun=0.85,
        parallax_mas=20.0,  # replaced internally for each a/sigma ratio
        gamma_kms=7.0,
        beta_g=0.20,
    )
    rows = resolution_bias_scan(
        truth,
        schedule,
        a_over_sigma_values=(0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 4.0),
        beta_values=(0.05, 0.20, 0.40),
        seeds=tuple(range(args.seeds)),
        sigma_response_mas=args.sigma_response_mas,
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    write_rows_csv(rows, str(output))
    print(f"wrote {len(rows)} fit records to {output}")


if __name__ == "__main__":
    main()
