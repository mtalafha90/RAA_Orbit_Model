#!/usr/bin/env python
"""Validation step V6: consistency against published Gaia DR3 NSS solutions.

The Gaia archive is not reachable from every environment, and it is not
reachable from the one this script was written in, so downloading is not
automated. Export the catalogue once with the query printed by ``--show-query``
and pass the resulting CSV here.

    python scripts/validate_against_dr3.py --show-query
    python scripts/validate_against_dr3.py nss_export.csv --output-dir results/dr3

The output converts the published Thiele-Innes constants to Campbell elements
in this project's convention and ranks the sample by the catalogue duplicity
diagnostics that the resolution-aware surrogate predicts.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from raa_orbit_model.dr3_validation import (  # noqa: E402
    NSS_ASTROMETRIC_ORBIT_QUERY,
    campbell_table,
    load_nss_csv,
    marginally_resolved_candidates,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare published Gaia DR3 NSS astrometric orbits with this model"
    )
    parser.add_argument("nss_csv", nargs="?", help="CSV exported from the Gaia archive")
    parser.add_argument("--output-dir", default="results/dr3_validation")
    parser.add_argument("--show-query", action="store_true",
                        help="print the ADQL to run at the Gaia archive, then exit")
    parser.add_argument("--min-multi-peak-fraction", type=float, default=2.0,
                        help="minimum ipd_frac_multi_peak, in per cent")
    parser.add_argument("--min-harmonic-amplitude", type=float, default=0.0,
                        help="minimum ipd_gof_harmonic_amplitude")
    args = parser.parse_args()

    if args.show_query:
        print("Run this at https://gea.esac.esa.int/archive/ and export as CSV:")
        print(NSS_ASTROMETRIC_ORBIT_QUERY)
        return

    if not args.nss_csv:
        parser.error("provide an exported NSS CSV, or use --show-query")

    rows = load_nss_csv(args.nss_csv)
    print(f"loaded {len(rows)} DR3 NSS rows from {args.nss_csv}")

    table = campbell_table(rows)
    print(f"converted {len(table)} rows to Campbell elements")

    candidates = marginally_resolved_candidates(
        table,
        min_multi_peak_fraction=args.min_multi_peak_fraction,
        min_harmonic_amplitude=args.min_harmonic_amplitude,
    )
    print(f"marginally resolved candidates: {len(candidates)}")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, data in (("dr3_campbell.csv", table), ("dr3_candidates.csv", candidates)):
        if not data:
            print(f"skipped {name}: nothing to write")
            continue
        path = output_dir / name
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=sorted(data[0]))
            writer.writeheader()
            writer.writerows(data)
        print(f"wrote {len(data)} rows to {path}")

    print(
        "\nNote: DR3 publishes catalogue-level orbits only. Epoch astrometry needed "
        "for a measurement-level test is a DR4 product (validation step V7)."
    )


if __name__ == "__main__":
    main()
