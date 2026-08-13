# Full-sky nominal-DR4 experiment

This note freezes the corrected 46-position full-sky development result after the orbit-convention and multi-peak-validity fixes.

## Integrity

- 46 distinct J2000 barycentric-mean-ecliptic sky positions.
- 23,000 fit records = 46 positions x 5 light fractions x 5 separation ratios x 10 seeds x 2 models.
- All fits succeeded and all rows are scientifically valid.
- The tested range a_rel,ang/sigma <= 1 stayed single-peak in the surrogate.
- The external resolved-astrometry/SB2 baseline is fixed at 5.4915748 yr.
- Nominal transit counts range from 53 to 310.

## Main result

At all 46 sky positions and all 5 tested separation ratios, the largest median paired model discrepancy occurs at beta_G = 0.25: 230/230 sky/separation combinations.

At beta_G = 0.25 and a_rel,ang/sigma = 1, the sky distribution of median Delta chi2 = chi2_photocentre - chi2_RAA has minimum 834.141, median 2195.391, and maximum 8766.639. The minimum sampled location is (+30 deg, 180 deg) in ecliptic latitude/longitude and the maximum is (+45 deg, 90 deg).

For beta_G = 0.25, the descriptive Q16(Delta chi2) > 0 boundary occurs at the smallest sampled ratio, a/sigma = 0.40, for 41/46 sky positions; the other 5 first cross at 0.50. Requiring all 10 realizations to have Delta chi2 > 0 gives 18/46 first crossing at 0.40, 23/46 at 0.50, and 5/46 at 0.60.

These are sampled surrogate-model transition descriptors, not formal significance thresholds and not Gaia instrumental resolution limits.

At beta_G = 0.25 and a/sigma = 1, descriptive Spearman correlations across the 46 sky positions are approximately +0.771 with transit count, -0.640 with maximum directional scan-angle gap, and +0.358 with normalized directional entropy. Transit count matters strongly but does not fully determine the mismatch strength.

## Reproduction

```bash
python scripts/analyze_sky_position_scan.py \
  results/sky_position_dr4_full/sky_position_bias.csv \
  --output-dir results/sky_position_dr4_full/analysis \
  --prefix full_sky_dr4 \
  --beta-map 0.25 \
  --a-over-sigma-map 1.0
```

The CLI writes paired, summary, boundary, and correlation CSVs plus PNG/PDF sampled-sky plots for median Delta chi2, the Q16 transition boundary, native transit count, and native maximum directional gap.
