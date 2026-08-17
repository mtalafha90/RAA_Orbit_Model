# Response-fidelity manuscript results

This note records the quantitative results used in the revised manuscript after the M0/M1/M2 response-fidelity programme.

## Model hierarchy

- **M0**: ordinary unresolved photocentre.
- **M1**: equal-width 1-D blended-profile peak response, matched to the published `gaiamock` response family.
- **M2**: finite-elongation orientation-dependent Gaussian response using the exact Penoyre (2026) effective width.

The response widths are research-surrogate parameters, not calibrated Gaia PLSF widths.

## 720-fit response-fidelity grid

Grid:
- beta_G = 0.05, 0.25, 0.45
- a/alpha = 0.4, 0.6, 0.8, 1.0
- beta_PSF/alpha = 1.5, 3.0
- 10 seeds
- 3 fitted models
- 87 Gaia-like transits on the archived RA=120 deg, Dec=30 deg nominal DR4 schedule.

All 720 fits succeeded and were scientifically valid; no multi-peak injections or final predictions occurred.

At beta_G=0.25, a/alpha=1, beta_PSF/alpha=1.5, the median biases are:

| model | M1 mass | M2 mass | parallax | beta_G | median chi2 |
|---|---:|---:|---:|---:|---:|
| M0 photocentre | -1.255% | -1.277% | +0.935% | -15.343% | 1494.5 |
| M1 equal-width | -0.026% | -0.095% | +0.132% | -3.303% | 352.4 |
| M2 oriented | +0.070% | -0.014% | +0.083% | +0.053% | 207.2 |

The paired medians at this point are Delta-chi2(M0-M2)=1272.7 and Delta-chi2(M1-M2)=129.3. At beta_PSF/alpha=3, Delta-chi2(M1-M2) falls to 5.0, consistent with convergence toward the M1 limit.

Among the three sampled light fractions, beta_G=0.25 maximizes the ten-seed median M0-M2 discrepancy at all eight sampled a/alpha and elongation combinations.

## 270-fit external-information control

Fixed response point:
- beta_G=0.25
- a/alpha=1
- beta_PSF/alpha=1.5
- Gaia AL sigma=0.10 mas
- 24 resolved-astrometry epochs and 48 SB2 epochs.

Only the external precision changes:
- strong: 0.20 mas, 0.10 km/s
- medium: 1.00 mas, 0.50 km/s
- weak: 2.00 mas, 1.00 km/s

Photocentre median mass biases grow from about 1.2% under strong external constraints to about 5--6% under weak constraints. M1 mass biases grow from below 0.1% to about 1--1.4%. The correctly specified M2 medians remain close to the injected values relative to the increasing random scatter.

Median Delta-chi2(M0-M2) decreases from 1289.8 (strong) to 887.4 (medium) and 807.6 (weak) even while the physical photocentre mass bias increases. This demonstrates that residual mismatch and physical bias need not move together.

## 100-seed matched-M2 control

At the same fixed response point with strong external constraints, 100 independent M2->M2 realizations give:
- median parallax fractional bias: +0.0167%
- q16/q84: -0.1274% / +0.1779%
- mean: +0.0214%
- standard deviation: 0.1509%
- fraction positive: 0.53
- median reduced chi2: 0.982

Median component-mass biases are -0.0348% (M1) and -0.0234% (M2), and the median beta_G bias is -0.0189%. The small positive parallax median seen in the earlier ten-seed grid is therefore not evidence of a coherent matched-model offset.

## Interpretation

The main manuscript result is that the response fidelity required for unbiased dynamical masses depends on the strength of independent orbit information. Strong visual-SB2 constraints can protect the masses against moderate Gaia response misspecification, while the light ratio remains substantially more sensitive. As the independent orbit weakens, Gaia response misspecification propagates into percent-level mass and parallax biases.
