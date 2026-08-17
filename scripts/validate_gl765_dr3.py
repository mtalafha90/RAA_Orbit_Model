#!/usr/bin/env python
"""Gaia DR3 catalogue-level validation for GJ 765.2 / HIP 96656.

Usage
-----
Print the target-specific ADQL query::

    python scripts/validate_gl765_dr3.py --show-query

Run the Gaia Archive query, export its result as CSV, then analyse it::

    python scripts/validate_gl765_dr3.py gj765_dr3.csv \
        --output results/dr3_validation/gj765_dr3_summary.json

This is validation step V6b: a catalogue/IPD consistency check. It is not the
DR4 epoch-level M0/M1/M2 response test.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from raa_orbit_model.dr3_target import (  # noqa: E402
    GJ765_ORBITAL_PARALLAX_MAS,
    GJ765_PERIOD_YR,
    gj765_target_query,
    load_target_export,
    summarize_target_rows,
    write_summary_json,
)


def fmt(value: float, digits: int = 3) -> str:
    return "null" if not math.isfinite(value) else f"{value:.{digits}f}"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Catalogue-level Gaia DR3 validation for GJ 765.2 / HIP 96656"
    )
    parser.add_argument("csv", nargs="?", help="CSV exported from the Gaia Archive target query")
    parser.add_argument("--show-query", action="store_true", help="print exact target ADQL and exit")
    parser.add_argument(
        "--output",
        default="results/dr3_validation/gj765_dr3_summary.json",
        help="JSON summary path",
    )
    parser.add_argument(
        "--radius-arcsec", type=float, default=5.0,
        help="target cone radius used by --show-query",
    )
    args = parser.parse_args()

    if args.show_query:
        print(gj765_target_query(args.radius_arcsec))
        return
    if not args.csv:
        parser.error("provide a Gaia Archive CSV or use --show-query")

    rows = load_target_export(args.csv)
    summary = summarize_target_rows(rows)
    write_summary_json(summary, args.output)

    print("GJ 765.2 / HIP 96656 — Gaia DR3 catalogue-level validation")
    print(f"source_id: {summary.source_id}")
    print(f"coordinate separation: {fmt(summary.separation_arcsec, 4)} arcsec")
    print(f"DR3 parallax: {fmt(summary.parallax_mas)} ± {fmt(summary.parallax_error_mas)} mas")
    if math.isfinite(summary.parallax_mas) and math.isfinite(summary.parallax_error_mas):
        delta = summary.parallax_mas - GJ765_ORBITAL_PARALLAX_MAS
        sigma = delta / summary.parallax_error_mas if summary.parallax_error_mas > 0 else math.nan
        print(
            f"vs Balega orbital parallax {GJ765_ORBITAL_PARALLAX_MAS:.1f} mas: "
            f"Δ={delta:+.3f} mas ({fmt(sigma, 2)} catalogue σ; external-orbit uncertainty not folded in)"
        )
    print(f"RUWE: {fmt(summary.ruwe)}")
    print(f"ipd_frac_multi_peak: {fmt(summary.ipd_frac_multi_peak_percent)} %")
    print(f"ipd_gof_harmonic_amplitude: {fmt(summary.ipd_gof_harmonic_amplitude)}")
    print(f"ipd_gof_harmonic_phase: {fmt(summary.ipd_gof_harmonic_phase_deg)} deg")
    print(f"non_single_star flag: {summary.non_single_star_flag}")
    print(f"NSS solution types: {summary.nss_solution_types or 'none'}")

    print("\nExternal M0 benchmark (not an RAA response fit):")
    print(f"mass fraction B = {summary.mass_fraction_secondary:.5f}")
    print(f"beta_V proxy = {summary.beta_V_proxy:.5f}")
    print(
        "predicted ordinary-photocentre axis = "
        f"{summary.predicted_M0_photocentre_axis_mas:.3f} mas"
    )

    if math.isfinite(summary.nearest_nss_photocentre_axis_mas):
        print("\nAstrometric NSS orbit found:")
        print(f"a0 = {summary.nearest_nss_photocentre_axis_mas:.3f} mas")
        print(f"period = {summary.nearest_nss_period_days:.3f} d")
        print(f"published visual/SB2 period = {GJ765_PERIOD_YR * 365.25:.3f} d")
        print(f"eccentricity = {summary.nearest_nss_eccentricity:.4f}")
        print(f"inclination = {summary.nearest_nss_inclination_deg:.3f} deg")
        print(f"omega_relative = {summary.nearest_nss_omega_relative_deg:.3f} deg")
        print(f"node (astrometric branch) = {summary.nearest_nss_node_deg:.3f} deg")
        print(
            "Reminder: the Thiele-Innes astrometric orbit has the usual 180-degree "
            "node ambiguity; SB2 velocities select the physical ascending-node branch."
        )
    else:
        print(
            "\nNo usable astrometric Thiele-Innes NSS orbit is present in this export. "
            "That is not a falsification of RAA because DR3 NSS processing is a selected "
            "sample and partially resolved doubles were filtered upstream."
        )

    print(f"\nwrote {args.output}")
    print(
        "Scope: DR3 tests catalogue/IPD consistency only. The direct M0/M1/M2 "
        "measurement-response likelihood requires released epoch astrometry/images."
    )


if __name__ == "__main__":
    main()
