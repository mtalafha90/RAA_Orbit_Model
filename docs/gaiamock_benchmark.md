# Benchmark against the published gaiamock along-scan response

## Why this benchmark exists

El-Badry et al. (2024), *Open Journal of Astrophysics* **7**, DOI 10.33232/001c.125461, released `gaiamock`, which contains a function `al_bias_binary`. Its docstring states that it predicts epoch astrometry for a binary by placing the one-dimensional centroid at the peak of the combined along-scan flux profile.

That is the same physical statement this project makes about its own baseline along-scan surrogate. A scan-angle- and separation-dependent Gaia response is therefore **already published and openly implemented**, and this project must not claim it as new. See `docs/literature_gap.md` for the surviving candidate inference gap.

The purpose of this benchmark is narrower: establish that the baseline implementation used in this repository is a correct reimplementation of the published equal-width response rather than an unrelated in-house construction.

## The two implementations are the same equal-width model

The published routine solves

```text
x = f * xi / (f + exp(xi^2 / 2 - xi * x))
```

by fixed-point iteration, where `xi` is the projected along-scan separation in units of an effective angular resolution `u`, and `f = F2 / F1`.

That is the stationary-point condition for the maximum of

```text
I(x) = exp(-x^2 / 2) + f * exp(-(x - xi)^2 / 2),
```

an equal-width two-Gaussian blend whose measured coordinate is the profile peak.

This project builds the same profile with flux fractions `1 - beta_G` and `beta_G` and width `sigma`, and locates its maximum. The parameter correspondence is

| this project | gaiamock |
|---|---|
| `sigma` (surrogate width) | `u` (effective angular-resolution scale) |
| `beta_G = F2 / (F1 + F2)` | `f = F2 / F1 = beta_G / (1 - beta_G)` |
| `B = M2 / (M1 + M2)` | `q = M2 / M1 = B / (1 - B)` |

`src/raa_orbit_model/gaiamock_reference.py` reproduces the published function for comparison. It is never used in the inference path.

## Numerical agreement

With `sigma = u = 90 mas` and `B = 0.4`, offsets in mas:

| `beta_G` | `d / u` | this project | gaiamock | relative difference |
|---|---|---|---|---|
| 0.10 | 0.20 | −5.42585 | −5.42587 | < 0.01% |
| 0.10 | 1.00 | −29.94903 | −29.94912 | < 0.01% |
| 0.10 | 2.50 | −88.87230 | −88.87234 | < 0.01% |
| 0.25 | 0.50 | −7.29373 | −7.29377 | < 0.01% |
| 0.25 | 2.00 | −62.50065 | −62.50076 | < 0.01% |
| 0.45 | 1.00 | +3.02715 | +3.02712 | < 0.01% |
| 0.45 | 2.00 | −37.11342 | −37.11358 | < 0.01% |

Agreement is at the level of the published solver's own convergence tolerance (`tol = 1e-6`) everywhere both models are valid. **The baseline surrogate reproduces the published equal-width response.**

## Differences from gaiamock's piecewise implementation

Two `gaiamock` branches are approximations that the numerical peak solver here does not require.

**1. Linearisation below `0.1 u`.** `gaiamock` returns the flux-weighted photocentre below a tenth of a resolution element. This project solves the profile peak directly. The small difference in that range is an approximation difference, not new measurement physics.

**2. A hard switch to the primary above `(3 - f) u`.** Above that separation `gaiamock` assigns the primary's position. For the restricted equal-width model, this project instead identifies the actual mode-splitting point and refuses to interpret a genuinely bimodal profile as one astrometric coordinate.

For the tested light fractions, the exact equal-width mode-splitting point lies beyond the `gaiamock` hard cut:

| `beta_G` | published cut `(3 - f)` | equal-width mode splitting | gap |
|---|---|---|---|
| 0.10 | 2.889 | 3.315 | 0.426 |
| 0.25 | 2.667 | 2.845 | 0.178 |
| 0.45 | 2.182 | 2.279 | 0.097 |

This comparison is useful for defining the domain of the baseline experiment. It is **not a claim that this repository introduced general blended-source resolvability theory**. Penoyre (2026), *RAS Techniques and Instruments* 5, rzaf062 (corrected by rzag016), gives a broader analytical treatment of blended Gaussian-source position and resolvability, including an elongated PSF with orientation-dependent effective width.

**3. Behaviour beyond mode splitting.** Once the restricted profile is genuinely bimodal this project flags the transit rather than returning one coordinate; `gaiamock` returns the primary's position. This is a policy difference. A real Gaia analysis requires the actual detection/window/IPD behaviour rather than either simplified policy.

## What this benchmark does and does not support

It supports the following statements:

1. the repository correctly reproduces the published `gaiamock` equal-width blended response over the shared single-peak domain;
2. its baseline numerical solver avoids `gaiamock`'s small-separation linearisation and uses the actual equal-width mode boundary instead of the `(3-f)` implementation cut;
3. the frozen experiments are therefore grounded in a published response family rather than in an arbitrary coordinate formula.

It does **not** support novelty of the response or resolvability calculation. It also does not show that an equal-width Gaussian is an accurate Gaia PLSF. Penoyre (2026) already motivates orientation-dependent effective width even in an idealized Gaussian treatment, and Rowell et al. (2026) describe the substantially richer DR4 PLSF.

## Reproducing this

```bash
pytest -q tests/test_gaiamock_benchmark.py
```

## Caveats

`u` in `gaiamock` and `sigma` here are effective model scales, not calibrated Gaia line-spread widths. Numerical agreement establishes implementation consistency only. The `gaiamock` default `u = 90 mas` can be used to reproduce that published simulation convention, but it must not be relabelled as a physical Gaia resolution threshold.

## References

- El-Badry, K. et al. (2024), *Open Journal of Astrophysics* **7**, DOI 10.33232/001c.125461, arXiv:2411.00088. Source: <https://github.com/kareemelbadry/gaiamock>.
- Penoyre, Z. (2026), "The position and resolvability of blended point sources", *RAS Techniques and Instruments* **5**, rzaf062, DOI 10.1093/rasti/rzaf062; correction rzag016, DOI 10.1093/rasti/rzag016.
- Rowell, N. et al. (2026), *A&A* **708**, A174, DOI 10.1051/0004-6361/202558618.
