# Benchmark against the published gaiamock along-scan response

## Why this benchmark exists

El-Badry et al. (2024), *Open Journal of Astrophysics* **7**, DOI 10.33232/001c.125461, released `gaiamock`, which contains a function `al_bias_binary`. Its docstring states that it

> predicts the epoch astrometry for a binary assuming that the 1D centroid is at the peak of the combined AL flux profile, following the model from Lindegren+2022

That is the same physical statement this project makes about its own along-scan surrogate. A scan-angle- and separation-dependent Gaia response is therefore **already published and openly implemented**, and this project must not claim it as new. See `docs/literature_gap.md` for the revised claim.

The remaining question is whether this project's surrogate *agrees* with the published one. If it does, the surrogate stops being an in-house invention and becomes a reimplementation of a published response. This note records that comparison.

## The two implementations are the same model

The published routine solves

```text
x = f * xi / (f + exp(xi^2 / 2 - xi * x))
```

by fixed-point iteration, where `xi` is the projected along-scan separation in units of an effective angular resolution `u`, and `f = F2 / F1`.

That is exactly the stationary-point condition for the maximum of

```text
I(x) = exp(-x^2 / 2) + f * exp(-(x - xi)^2 / 2),
```

an **equal-width two-Gaussian blend** whose measured coordinate is the profile peak. Setting `d(ln I)/dx = 0` and rearranging gives the published fixed-point form directly.

This project builds the same profile with flux fractions `1 - beta_G` and `beta_G` and width `sigma`, and locates its maximum. The parameter correspondence is

| this project | gaiamock |
|---|---|
| `sigma` (surrogate width) | `u` (effective angular resolution) |
| `beta_G = F2 / (F1 + F2)` | `f = F2 / F1 = beta_G / (1 - beta_G)` |
| `B = M2 / (M1 + M2)` | `q = M2 / M1 = B / (1 - B)` |

`src/raa_orbit_model/gaiamock_reference.py` reproduces the published function verbatim for comparison. It is never used in the inference path.

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

Agreement is at the level of the published solver's own convergence tolerance (`tol = 1e-6`) everywhere both models are valid. **The surrogate reproduces the published response.**

## Where they differ, and why this project is the more exact of the two

The published implementation is piecewise. Two of its three branches are approximations that this project does not need.

**1. Linearisation below `0.1 u`.** `gaiamock` returns the flux-weighted photocentre below a tenth of a resolution element, for numerical stability. This project solves the peak everywhere. The two differ by 0.03–0.08% at `d = 0.05 u`, which is the published approximation error, not a disagreement.

**2. A hard switch to the primary above `(3 - f) u`.** Above that separation `gaiamock` assigns the primary's position. This project instead applies the exact mode-splitting criterion derived in `src/raa_orbit_model/gaia.py`, which solves `sinh(2u) - 2u = ln((1 - beta_G) / beta_G)` for `d_crit / sigma = 2 cosh(u)`.

The exact boundary is **always wider** than the published cut:

| `beta_G` | published cut `(3 - f)` | exact mode splitting | gap |
|---|---|---|---|
| 0.10 | 2.889 | 3.315 | 0.426 |
| 0.25 | 2.667 | 2.845 | 0.178 |
| 0.45 | 2.182 | 2.279 | 0.097 |

So the published cut fires while the profile is **still genuinely single-peaked**, and in that gap it returns a displaced value. At `beta_G = 0.10` and `d = 2.9 u` the true peak is at −103.96 mas while `gaiamock` returns −104.40 mas, the primary's position: a 0.42% displacement in a regime where an exact single-peak solution exists.

**3. Behaviour beyond mode splitting.** Once the profile is genuinely bimodal this project returns `NaN` and flags the transit, because a single-coordinate response is undefined; `gaiamock` returns the primary's position. This is a difference in **policy**, not in physics, and it must be stated whenever the two are compared. The published choice is closer to what Gaia's own pipeline does — the documented behaviour is to exclude samples around a detected secondary peak and still fit a single line-spread function, which yields a finite, biased position rather than no position at all.

## What this project can and cannot claim

**Cannot claim:** that a scan-angle- and separation-dependent Gaia response is novel. It is published (El-Badry et al. 2024) and follows Lindegren (2022).

**Can claim, and now demonstrate:**

1. this project's surrogate reproduces the published response to the published solver's own tolerance;
2. it removes both of the published piecewise approximations, solving the peak exactly at all separations below mode splitting;
3. it supplies an **exact** mode-splitting boundary where the published implementation uses the approximation `(3 - f)`, and that approximation is systematically early.

## Reproducing this

```bash
pytest -q tests/test_gaiamock_benchmark.py
```

## Caveats

`u` in `gaiamock` and `sigma` here are both effective research scales, not calibrated Gaia line-spread widths. Agreement between the two models says nothing about whether either matches the real Gaia instrument; it establishes only that this project has implemented the published response correctly. The `gaiamock` default `u = 90 mas` is the closest thing to a published anchor for that scale and is worth adopting in place of an arbitrary value.

## Reference

- El-Badry, K., Lam, C., Holl, B., Halbwachs, J.-L., Rix, H.-W., Mazeh, T., & Shahaf, S. (2024), *Open Journal of Astrophysics* **7**, DOI 10.33232/001c.125461, arXiv:2411.00088. Source: <https://github.com/kareemelbadry/gaiamock>.
- The `Lindegren (2022)` reference credited in the `al_bias_binary` docstring has **not** been resolved to a specific document. Confirm it from the reference list of Holl et al. (2023), A&A **674**, A25, before citing it.
