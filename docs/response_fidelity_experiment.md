# Response-fidelity experiment

## Scientific question

The baseline full-sky experiment injected and fitted the same equal-width 1-D blended response. That is useful for code validation and for demonstrating that the photocentre approximation can be misspecified, but the resolution-aware model is correct by construction.

The next experiment asks a stronger question:

> **What biases in individual stellar masses and parallax arise when a finite-elongation marginal-resolution response is fitted with progressively simpler measurement models?**

The primary endpoint is physical-parameter bias and, after posterior sampling is enabled for this hierarchy, 68/95% posterior coverage. `Delta chi2` is retained only as a model-mismatch diagnostic.

## Measurement-model hierarchy

The same synthetic realization is fitted three times.

- **M0 — photocentre:** ordinary unresolved flux-weighted photocentre.
- **M1 — equal-width 1-D response:** the Lindegren/`gaiamock`-family blended-profile peak used by the frozen experiments.
- **M2 — finite-elongation response:** an idealised elongated Gaussian following the reduction in Penoyre (2026), using its orientation-dependent effective width.

M2 is the injection model for this stage.

## Penoyre-style effective width

For an idealised elongated Gaussian with along-scan/narrow width `alpha` and across-scan/long width `beta`, Penoyre (2026) gives

```text
gamma(phi) = alpha / sqrt(1 + [(alpha/beta)^2 - 1] sin(phi)^2)
```

where `phi` is the angle between the source-pair separation vector and the scan axis. The extrema lie on the line joining the two sources. Once distances are normalized by `gamma`, the blend has the same one-dimensional functional form as the un-elongated case.

The implementation therefore does **not** use Penoyre's low-separation series expansion. It computes `gamma(phi)` exactly and feeds `r/gamma` into the repository's numerical equal-width peak solver. This avoids dependence on the sign errors corrected in Penoyre's published correction rzag016.

### Relation between M1 and M2

M1 is **not** the `beta=alpha` limit. Penoyre notes that the 1-D Lindegren-style Gaia calculation corresponds to an effectively unconstrained across-scan direction. The code therefore regression-tests the correct relation:

```text
M2(beta -> infinity) == M1
```

for both the predicted AL coordinate and the single/multi-peak classification.

For finite `beta`, mode splitting depends on the full on-sky pair separation through

```text
rho = r / gamma(phi),
```

not merely on the projected AL separation.

## Instrument caveat

`alpha`, `beta`, and their ratio are **research-surrogate parameters**. They are not calibrated Gaia PLSF widths. Penoyre uses `beta/alpha=3` as an approximate Gaia-like elongation example, but this repository treats elongation as a configurable experiment axis rather than a fixed Gaia calibration. The final instrument-facing model must follow the released Gaia DR4 PLSF/epoch-image calibration products.

## First regression pilot

Use the already archived pilot DR4 nominal schedule:

```bash
python scripts/run_response_fidelity.py \
  --schedule-file schedules/ra120_dec30_dr4.csv \
  --a-over-alpha-values 0.6 1.0 \
  --beta-values 0.25 \
  --beta-over-alpha-values 1.5 3.0 \
  --seeds 3 \
  --alpha-mas 50 \
  --output results/response_fidelity_pilot.csv
```

This produces

```text
2 separations x 1 light fraction x 2 elongations x 3 seeds x 3 models = 36 fit records.
```

The pilot is intentionally small. Before expanding the grid, verify:

1. every M2-injected retained epoch is single-peaked;
2. all three deterministic fits converge;
3. final M1/M2 solutions remain scientifically valid;
4. the `beta/alpha -> infinity` regression test recovers M1 exactly;
5. M2 gives the best self-consistent fit to M2 injections without pathological parameter excursions;
6. the mass/parallax bias ordering is stable across more than one seed.

## Next grid after the pilot

If the pilot passes, a useful deterministic map is

```text
beta_G = 0.05, 0.25, 0.45
a/alpha = 0.4, 0.6, 0.8, 1.0
beta/alpha = 1.5, 3.0
10 seeds
3 fitted models
```

for `720` fit records on one exact Gaia schedule. The same grid should then be repeated on a small set of deliberately different scan geometries before posterior coverage is attempted.

The objective is **not** to declare `beta/alpha=3` the Gaia PSF. It is to determine whether the qualitative mass-bias correction survives when the data-generating response is more general than M1.

## Posterior stage

Only after the deterministic hierarchy is validated should the response-fidelity experiment be moved to the posterior sampler. The quantities to record are at minimum

```text
(M1_fit - M1_true) / M1_true
(M2_fit - M2_true) / M2_true
(parallax_fit - parallax_true) / parallax_true
(beta_G_fit - beta_G_true) / beta_G_true
```

plus empirical inclusion of the truth in nominal 68% and 95% credible intervals. The strongest scientific result would be a map of where M0 and M1 lose calibration and whether M2 restores it.

## References

- Penoyre, Z. (2026), *RAS Techniques and Instruments* 5, rzaf062, DOI 10.1093/rasti/rzaf062; correction rzag016.
- El-Badry, K. et al. (2024), *Open Journal of Astrophysics* 7, DOI 10.33232/001c.125461.
- Rowell, N. et al. (2026), *A&A* 708, A174, DOI 10.1051/0004-6361/202558618.
