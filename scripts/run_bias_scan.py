#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from raa_orbit_model.experiments import resolution_bias_scan, write_rows_csv
from raa_orbit_model.kepler import BinaryParams


def main():
    parser = argparse.ArgumentParser(description="Run the first RAA photocentre-vs-resolution bias experiment")
    parser.add_argument("--output", default="bias_scan.csv")
    parser.add_argument("--seeds", type=int, default=3, help="number of seeds, starting at zero")
    parser.add_argument("--sigma-response-mas", type=float, default=50.0)
    args = parser.parse_args()

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
