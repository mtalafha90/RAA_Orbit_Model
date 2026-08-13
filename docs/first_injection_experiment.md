# First injection/recovery experiment

## Status

This is a **code-validation experiment**, not an empirical Gaia result and not a claimed astrophysical threshold.

The purpose is to verify that the proposed comparison behaves correctly in controlled limiting regimes before adding a Gaia scanning law, calibrated PSF/LSF information, posterior sampling, or real target data.

## Injection model

The injected binary is a Newtonian SB2 with

- `P = 2.0 yr`
- `T_peri = 0.15 yr`
- `e = 0.25`
- `i = 72 deg`
- `omega = 55 deg`
- `Omega = 120 deg`
- `M1 = 1.25 Msun`
- `M2 = 0.85 Msun`
- `gamma = 7 km/s`

The angular scale is varied through the dimensionless ratio

\[
R = a_{\rm rel,ang}/\sigma,
\]

where \(\sigma=50\) mas is the width of the **Gaussian research surrogate**, not a calibrated Gaia LSF width. The injected parallax is changed to achieve each requested \(R\).

Three secondary Gaia-band light fractions are tested:

\[
\beta_G \in \{0.05, 0.20, 0.40\}.
\]

Each synthetic realization contains

- 24 two-dimensional relative-astrometry epochs with 0.20 mas isotropic noise;
- 48 SB2 epochs with 0.10 km/s noise on each component;
- 72 synthetic along-scan epochs with 0.10 mas noise;
- scan angles drawn uniformly from 0 to 180 degrees.

The uniform scan-angle distribution is intentionally **not** a Gaia scanning-law simulation.

All 11 baseline physical parameters are free in each fit. Three random seeds (0, 1, 2) are used.

## Competing models

Every dataset is injected with the two-profile resolution-aware surrogate and then fitted independently with:

1. `photocentre` — assumes the unresolved photocentre at every epoch;
2. `resolution_aware` — uses the same two-profile response used for injection.

Thus this experiment tests whether a misspecified photocentre likelihood becomes detectably biased as the pair becomes more resolved along scan.

## Pilot results

The table reports the median over the three seeds. Define

\[
\Delta\chi^2 = \chi^2_{\rm photocentre}-\chi^2_{\rm resolution-aware}.
\]

Positive values favour the correctly specified resolution-aware surrogate. `|bias beta|` is the absolute fitted-minus-injected error in \(\beta_G\).

| beta_G | a/sigma | median Delta chi2 | median |bias beta| photocentre | median |bias beta| resolution-aware |
|---:|---:|---:|---:|---:|
| 0.05 | 0.05 | -0.002 | 0.004856 | 0.004892 |
| 0.05 | 0.10 | -0.021 | 0.002482 | 0.002632 |
| 0.05 | 0.20 | -0.146 | 0.000948 | 0.001558 |
| 0.05 | 0.50 | 2.157 | 0.003225 | 0.000482 |
| 0.05 | 1.00 | 165.604 | 0.012466 | 0.000502 |
| 0.05 | 2.00 | 2794.442 | 0.032501 | 0.000794 |
| 0.05 | 4.00 | 5119.561 | 0.046160 | 0.000637 |
| 0.20 | 0.05 | -0.006 | 0.002648 | 0.002736 |
| 0.20 | 0.10 | -0.053 | 0.001194 | 0.001546 |
| 0.20 | 0.20 | -0.329 | 0.001552 | 0.001025 |
| 0.20 | 0.50 | 17.894 | 0.008327 | 0.000411 |
| 0.20 | 1.00 | 1793.453 | 0.034064 | 0.000170 |
| 0.20 | 2.00 | 55148.048 | 0.116202 | 0.000481 |
| 0.20 | 4.00 | 118099.835 | 0.182545 | 0.000512 |
| 0.40 | 0.05 | -0.003 | 0.006684 | 0.006641 |
| 0.40 | 0.10 | -0.028 | 0.003444 | 0.003276 |
| 0.40 | 0.20 | -0.208 | 0.002275 | 0.001610 |
| 0.40 | 0.50 | 2.843 | 0.004524 | 0.000470 |
| 0.40 | 1.00 | 838.870 | 0.020677 | 0.000251 |
| 0.40 | 2.00 | 264184.634 | 0.161441 | 0.000065 |
| 0.40 | 4.00 | 833704.351 | 0.351917 | 0.000305 |

## Interpretation limited to the surrogate experiment

The results show the required limiting behaviour:

- at `a/sigma <= 0.2`, the two models are effectively indistinguishable at the level of this three-seed pilot;
- by `a/sigma = 0.5`, the misspecified photocentre response begins to produce detectable degradation for these simulated conditions;
- at larger `a/sigma`, forcing a photocentre response produces increasingly poor fits and substantial bias in the recovered light fraction.

These numbers **must not be converted into a Gaia angular-resolution threshold**. The response width is a free Gaussian surrogate and the scan angles are not generated from the Gaia scanning law.

The experiment only establishes that the code can expose a resolution-dependent model-misspecification regime, which is the prerequisite for the next, more realistic stages.

## Reproduction

After installing the package:

```bash
python scripts/run_bias_scan.py --seeds 3 --output bias_scan.csv
```

## Next required upgrades

1. replace uniform scan angles with a documented Gaia scanning-law calculation for selected sky positions;
2. increase Monte Carlo seeds and map uncertainty on the bias boundary;
3. add heteroscedastic and correlated relative-astrometry errors;
4. vary period, eccentricity, mass ratio, phase coverage, and noise levels;
5. separate optimization failure from forward-model failure using multiple starts;
6. add posterior inference only after the deterministic recovery surface is understood;
7. do not fit real Gaia epoch data until appropriate epoch products and calibration information are available.
