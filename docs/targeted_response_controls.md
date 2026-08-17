# Targeted controls after the 720-fit response-fidelity grid

The 720-fit M0/M1/M2 grid showed two results that require targeted follow-up before the manuscript Results section is rewritten:

1. the ordinary photocentre can produce percent-level component-mass and parallax bias at the difficult single-peak point, while the simpler resolution-aware M1 response removes most of the mass bias even when M1 is measurably misspecified relative to M2;
2. the matched M2 fits showed a small positive parallax median across the ten reused seed realizations, which cannot yet be interpreted as a systematic estimator bias.

This document freezes the two controls used to distinguish those effects.

## Control A: external-information strength

The physical response point is fixed at

- `beta_G = 0.25`,
- `a/alpha = 1.0`,
- `beta_PSF/alpha = 1.5`,
- `alpha = 50 mas` as a research-surrogate scale,
- Gaia AL uncertainty `0.10 mas`,
- the same archived nominal Gaia schedule used for the response-fidelity experiment.

The same M2 injection is fit with M0, M1, and M2. Only the precision of the independent resolved-astrometry and SB2 constraints changes:

| level | resolved astrometry | SB2 RV |
|---|---:|---:|
| strong | 0.20 mas | 0.10 km/s |
| medium | 1.00 mas | 0.50 km/s |
| weak | 2.00 mas | 1.00 km/s |

The epoch counts remain `N_ast=24` and `N_RV=48`. Thirty seeds give

`3 levels x 30 seeds x 3 fitted response models = 270 fits`.

Run:

```bash
python -u scripts/run_external_information_control.py \
  --schedule-file schedules/ra120_dec30_dr4.csv \
  --seeds 30 \
  --alpha-mas 50 \
  --a-over-alpha 1.0 \
  --beta-g 0.25 \
  --beta-over-alpha 1.5 \
  --output results/external_information_control.csv
```

The companion summary file groups by external-information level and fitted model and reports median/16th/84th-percentile fractional bias in `M1`, `M2`, parallax, and `beta_G`, plus paired Delta-chi2 diagnostics.

### Interpretation target

This control asks when Gaia response fidelity begins to control the dynamical masses rather than primarily the light fraction. It is not a sky-position experiment and does not establish a universal precision threshold.

## Control B: 100-seed matched-M2 estimator check

The second control uses exactly the matched M2 injection and fit at the same difficult response point, with the original strong external constraints. It fits **M2 only** for 100 independent seeds, avoiding the unnecessary M0/M1 fits.

Run:

```bash
python -u scripts/run_matched_m2_control.py \
  --schedule-file schedules/ra120_dec30_dr4.csv \
  --seeds 100 \
  --alpha-mas 50 \
  --a-over-alpha 1.0 \
  --beta-g 0.25 \
  --beta-over-alpha 1.5 \
  --output results/matched_m2_100seed.csv
```

The one-row summary reports mean, median, standard deviation, 16th/84th percentiles, and fraction positive for the fractional biases in the two component masses, parallax, light fraction, and inclination.

### Interpretation target

If the matched-M2 parallax distribution recentres around zero over 100 seeds, the small ten-seed positive median was sampling noise. If a comparable positive offset persists, the next task is to separate nonlinear-estimator bias from optimization or model-parameterization effects before interpreting the response-fidelity grid.

## Scope

Both controls use the deterministic validation fitter. They do not test posterior calibration or credible-interval coverage. Posterior coverage should be attempted only after these deterministic controls establish which measurement-response hierarchy and external-information regime are scientifically informative.
