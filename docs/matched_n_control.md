# Matched-transit-count control

The native full-sky experiment includes both real variation in Gaia transit count and scan-angle/time geometry. This control fixes the Gaia observation count to N=53, the smallest count in the 46-position full-sky grid.

For each native ordered schedule, the control divides the transit sequence into 53 equal-count strata and selects the central-ranked transit from each stratum. This deterministic rule spans the full schedule. It is a control-device subsampling rule, not a Gaia source-selection model.

For a given seed and injected system, the code first regenerates the full native synthetic realization and then selects the matched Gaia epochs from it. The resolved astrometry, SB2 data, and selected Gaia noise draws therefore remain identical to the native realization for the same seed; only the retained Gaia epoch/scan geometry changes.

Default reduced grid:

- beta_G = 0.05, 0.25, 0.45
- a_rel,ang/sigma = 0.40, 0.50, 0.60, 1.00
- 10 seeds
- 46 sky positions
- 2 measurement models

Total: 11,040 fit records.

Run:

```bash
python scripts/run_matched_n_control.py \
  --native-dir results/sky_position_dr4_full \
  --output-dir results/sky_position_dr4_full_matched_n53 \
  --release dr4 \
  --matched-n 53
```

The run is restart-safe at the per-position level and writes `matched_n53_bias.csv` after all 46 positions are complete.
