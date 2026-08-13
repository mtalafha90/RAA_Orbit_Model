#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from raa_orbit_model.experiments import write_rows_csv
from raa_orbit_model.kepler import BinaryParams
from raa_orbit_model.matched_control import matched_n_schedule
from raa_orbit_model.matched_sky_experiments import matched_n_sky_resolution_bias_scan
from raa_orbit_model.scanning import schedule_from_csv, write_schedule_csv
from raa_orbit_model.sky_study import SkyPosition, augment_bias_rows, scan_geometry_diagnostics


def main():
    p = argparse.ArgumentParser(description="Run the matched-N full-sky control")
    p.add_argument("--native-dir", default="results/sky_position_dr4_full")
    p.add_argument("--output-dir", default="results/sky_position_dr4_full_matched_n53")
    p.add_argument("--release", default="dr4"); p.add_argument("--matched-n", type=int, default=53)
    p.add_argument("--a-over-sigma-values", nargs="+", type=float, default=(0.40, 0.50, 0.60, 1.00))
    p.add_argument("--beta-values", nargs="+", type=float, default=(0.05, 0.25, 0.45))
    p.add_argument("--seeds", type=int, default=10); p.add_argument("--sigma-response-mas", type=float, default=50.0)
    p.add_argument("--overwrite", action="store_true"); args = p.parse_args()
    if args.matched_n <= 0 or args.seeds <= 0 or args.sigma_response_mas <= 0: p.error("matched-n, seeds and sigma must be > 0")
    if any(v <= 0 for v in args.a_over_sigma_values): p.error("a/sigma values must be > 0")
    if any(v < 0 or v > 0.5 for v in args.beta_values): p.error("beta values must lie in [0,0.5]")

    native = Path(args.native_dir); out = Path(args.output_dir)
    grid_path = native / "sky_grid.csv"; bias_path = native / "sky_position_bias.csv"; schedules = native / "schedules"
    for path in (grid_path, bias_path, schedules):
        if not path.exists(): raise FileNotFoundError(path)
    grid = pd.read_csv(grid_path)
    baselines = pd.read_csv(bias_path, usecols=["external_data_baseline_yr"])["external_data_baseline_yr"].unique()
    if len(baselines) != 1: raise ValueError("native experiment must have one fixed external-data baseline")
    baseline = float(baselines[0])
    positions = [SkyPosition(r.sky_id, r.ecliptic_lon_deg, r.ecliptic_lat_deg, r.ra_deg, r.dec_deg) for r in grid.itertuples()]
    per = out / "per_position"; sched_out = out / "schedules"; per.mkdir(parents=True, exist_ok=True); sched_out.mkdir(parents=True, exist_ok=True)
    nper = len(args.a_over_sigma_values) * len(args.beta_values) * args.seeds * 2
    print(f"sky positions: {len(positions)}; matched N: {args.matched_n}; expected records: {len(positions)*nper}")
    print(f"fixed external astrometry/RV baseline: {baseline:.9f} yr")

    truth = BinaryParams(period_yr=2.0, t_peri_yr=0.15, eccentricity=0.25, inclination_deg=72.0,
        omega_deg=55.0, node_deg=120.0, m1_msun=1.25, m2_msun=0.85, parallax_mas=20.0, gamma_kms=7.0, beta_g=0.20)
    geometry = []
    for k, pos in enumerate(positions, 1):
        full = schedule_from_csv(schedules / f"{pos.sky_id}_{args.release}.csv")
        if full.n_transits < args.matched_n: raise ValueError(f"{pos.sky_id}: only {full.n_transits} native transits")
        matched, _ = matched_n_schedule(full, args.matched_n)
        write_schedule_csv(matched, sched_out / f"{pos.sky_id}_{args.release}_matched_n{args.matched_n}.csv")
        d0, d1 = scan_geometry_diagnostics(full), scan_geometry_diagnostics(matched)
        geometry.append({"sky_id":pos.sky_id,"ecliptic_lon_deg":pos.ecliptic_lon_deg,"ecliptic_lat_deg":pos.ecliptic_lat_deg,
            "native_n_transits":full.n_transits,"matched_n_transits":matched.n_transits,
            "native_max_gap_deg":d0.directional_max_gap_deg,"matched_max_gap_deg":d1.directional_max_gap_deg,
            "native_entropy_normalized":d0.directional_entropy_normalized,"matched_entropy_normalized":d1.directional_entropy_normalized})
        result = per / f"{pos.sky_id}_{args.release}_matched_n{args.matched_n}_bias.csv"
        if result.exists() and not args.overwrite:
            if len(pd.read_csv(result)) != nper: raise RuntimeError(f"incomplete result: {result}; use --overwrite")
            print(f"[{k:02d}/{len(positions):02d}] resume: {pos.sky_id}"); continue
        rows = matched_n_sky_resolution_bias_scan(truth, full, matched_n_transits=args.matched_n,
            external_baseline_yr=baseline, a_over_sigma_values=tuple(args.a_over_sigma_values),
            beta_values=tuple(args.beta_values), seeds=tuple(range(args.seeds)), sigma_response_mas=args.sigma_response_mas)
        augment_bias_rows(rows, pos, d1)
        for row in rows:
            row["native_scan_n_transits"] = full.n_transits
            row["native_scan_directional_max_gap_deg"] = d0.directional_max_gap_deg
            row["native_scan_directional_entropy_normalized"] = d0.directional_entropy_normalized
        write_rows_csv(rows, str(result)); print(f"[{k:02d}/{len(positions):02d}] wrote {len(rows)} rows: {result}")

    pd.DataFrame(geometry).to_csv(out / f"matched_n{args.matched_n}_scan_geometry.csv", index=False)
    combined = []
    for pos in positions:
        f = per / f"{pos.sky_id}_{args.release}_matched_n{args.matched_n}_bias.csv"; d = pd.read_csv(f)
        if len(d) != nper: raise RuntimeError(f"cannot combine incomplete file: {f}")
        combined.append(d)
    final = pd.concat(combined, ignore_index=True); final.to_csv(out / f"matched_n{args.matched_n}_bias.csv", index=False)
    print(f"combined fit records: {len(final)}")

if __name__ == "__main__": main()
