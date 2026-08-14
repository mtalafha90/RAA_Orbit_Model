# Literature gap investigation

## Scope

This note records the literature status of the proposed RAA Orbit Model. It deliberately separates established results from the **candidate** research gap. It does not claim priority.

## Established capabilities

The following are already present in the literature and therefore are **not** claimed as novel contributions here.

### Joint astrometry and spectroscopy

BINARYS (Leclerc et al. 2023) combines Hipparcos/Gaia absolute astrometry, relative astrometry, and radial velocities. Its RV interface can identify measurements of either component, allowing treatment of SB1 and SB2 systems. Therefore, “relative astrometry + SB2 + Gaia” is not itself a gap.

### Gaia DR3 astrometric binaries

The Gaia DR3 non-single-star astrometric processing (Halbwachs et al. 2023) fits Keplerian orbital solutions, but partially resolved double stars were removed upstream from the astrometric-binary processing. This is direct evidence that the difficult transition regime was not part of the published DR3 orbital pipeline.

### Scan-angle-dependent close-pair systematics

Holl et al. (2023) showed that close source structure can produce scan-angle-dependent biases in Gaia image parameter determination and derived time series. They provide an analytical approximation, based on Lindegren (2022), for the along-scan bias of a binary as a function of projected along-scan separation and flux ratio.

### Scan-angle- and resolution-aware measurement response

**This is published, openly implemented, and must not be claimed as novel by this project.**

El-Badry et al. (2024) released `gaiamock`, which contains `al_bias_binary`. Its docstring states that it "predicts the epoch astrometry for a binary assuming that the 1D centroid is at the peak of the combined AL flux profile, following the model from Lindegren+2022". The measured along-scan coordinate is an explicitly nonlinear function of projected along-scan separation and light ratio, reducing to the flux-weighted photocentre only in the small-separation limit. That is precisely the physical statement this project set out to make.

The published routine solves `x = f xi / (f + exp(xi^2/2 - xi x))`, which is the stationary-point condition for the maximum of an equal-width two-Gaussian blend. It is therefore the *same model* as this project's surrogate. `docs/gaiamock_benchmark.md` records the numerical comparison: the two agree to the published solver's own convergence tolerance everywhere both are valid.

Their stated purpose was catalogue selection-function modelling rather than a target-level joint posterior using independent resolved astrometry and both SB2 velocity curves.

### The published response is used for simulation, not for inference

Verified by reading the `gaiamock` source: `al_bias_binary` is called only by the prediction functions `predict_astrometry_luminous_binary` and `predict_astrometry_and_rvs_simultaneously`. **No fitting routine uses it.** `fit_5par_solution_only`, `check_7par`, `check_9par`, `fit_orbital_solution_nonlinear`, `fit_full_astrometric_cascade` and `mcmc_fit_with_thiele_innes_elements` all fit a plain photocentre model on a Thiele-Innes design matrix.

That is a deliberate design choice, because `gaiamock` emulates Gaia's own pipeline. The published state of the art is therefore: **generate with a resolution-aware response, fit with a photocentre model.**

### Gaia PSF/LSF calibration

Gaia image parameter determination rests on calibrated PSF/LSF models; Rowell et al. (2021) document the EDR3 PSF/LSF calibration. Therefore a final image-level RAA model must not treat a Gaussian profile as an exact Gaia instrument model.

## Candidate gap, as revised

The original wording of this section listed the resolution-aware measurement response itself as one of the missing ingredients. **That was wrong**, and the revision below is the direct result of finding `al_bias_binary` in `gaiamock`. The claim is now considerably narrower.

What is **not** a gap:

- combining Gaia or Hipparcos astrometry with relative astrometry and radial velocities (BINARYS);
- a scan-angle- and separation-dependent Gaia along-scan response (El-Badry et al. 2024, following Lindegren 2022);
- demonstrating that fitting a photocentre model to resolution-aware data produces biased or spurious astrometric solutions — El-Badry et al. (2024) show this at population level.

What we did **not find** published, after the searches carried out for this project:

1. the resolution-aware response used as the **inference** model rather than only as a generative one;
2. that inference performed at **target level**, jointly with independent resolved relative astrometry and **both** SB2 velocity curves;
3. propagation of the resulting measurement-model choice to **component masses** specifically, as opposed to parallax and proper motion.

This is a **methodological** claim about where an existing measurement model is applied, not a claim to new measurement physics. It remains a candidate: no priority claim is made, and the manuscript must not use "first", "novel", or "unprecedented".

## Verification status of this note

The `gaiamock` findings above were verified by reading the source directly at <https://github.com/kareemelbadry/gaiamock>, and are reproduced for comparison in `src/raa_orbit_model/gaiamock_reference.py`.

The following leads surfaced during the literature review but **could not be read**, because scholarly domains were unreachable from the working environment. They must be checked against the published papers before citation:

- **Penoyre, Z. (2025/26), "The position and resolvability of blended point sources", RAS Techniques and Instruments**, DOI 10.1093/rasti/rzaf062, arXiv:2512.05047. Reportedly derives analytic blended-source positions for an *elongated* Gaia-like PSF, including a resolvability criterion and an orientation-dependent effective width. If accurate, this project's `critical_blended_separation_sigma` is the equal-width special case of a published general result, and the constant surrogate width `sigma` should become an orientation-dependent `sigma_eff(psi)`. **Highest-priority check.**
- **Baycroft, Faria & Delisle (2026)**, `kima` with Gaia epoch astrometry and RV joint Bayesian fitting, reportedly photocentre-based. Closest competitor on the inference side.
- **Harrison et al. (2023), "SEAPipe: the source environment analysis pipeline", A&A 679, A158.** If SEAPipe feeds DR4, the statement that partially resolved doubles are excluded upstream may no longer hold for DR4 and must be hedged.
- **Rowell et al. (2026)**, DR4 PSF drift-scan modelling, arXiv:2602.20906, which would supersede Rowell et al. (2021) as the PSF anchor for DR4.
- The **`Lindegren (2022)`** reference credited in the `al_bias_binary` docstring is **unresolved**. Two independent searches returned contradictory candidates. Settle it from the reference list of Holl et al. (2023), A&A 674, A25.
- **Gaia DR4 is not released.** Expected 2 December 2026. Product names in the DR4 section below (`epoch_astrometry`, `epoch_image`, `rvs_epoch_data_double`) should be checked against the released DR4 datamodel; searches returned conflicting evidence on whether `epoch_image` exists under that name, and the RVS double-lined product may be named `EPOCH_PARAMETERS_RVS_DOUBLE`.

## Why the transition regime is distinct

For a fully unresolved binary, the photocentre relative to the barycentre is

\[
\mathbf r_{\rm ph}=(\beta-B)\,\mathbf r,
\]

where

\[
B=\frac{M_2}{M_1+M_2},\qquad
\beta=\frac{F_2}{F_1+F_2},
\]

and \(\mathbf r=\mathbf r_2-\mathbf r_1\).

For Gaia observation angle \(\psi\), the projected along-scan component of the true relative separation is

\[
\Delta\eta = \Delta\alpha^*\sin\psi + \Delta\delta\cos\psi,
\]

up to the adopted scan-angle convention.

In the marginal-resolution regime, the effective measured location can depend nonlinearly on \(\Delta\eta\) and the flux ratio. Thus a single photocentre orbit is no longer a complete observation model.

## DR4 motivation

The official Gaia DR4 expected-content documentation lists planned epoch-level products including `epoch_astrometry`, `epoch_image`, and double-lined epoch spectroscopy (`rvs_epoch_data_double`). These products are planned and explicitly subject to processing/validation changes. If released as expected, they would allow the candidate model to be tested much closer to the measurement level than is possible with DR3 catalogue products.

Until then, the resolution-aware component of this project must be validated with synthetic/injection tests and with catalogue-level consistency tests only. It must not be presented as empirically validated against unavailable Gaia epoch images.

## References

- Leclerc, A. et al. (2023), “Combining Hipparcos and Gaia data for the study of binaries: The BINARYS tool”, *A&A*, 672, A82. DOI: 10.1051/0004-6361/202244144.
- Halbwachs, J.-L. et al. (2023), “Gaia Data Release 3: Astrometric binary star processing”, *A&A*, 674, A9. arXiv:2206.05726.
- Holl, B. et al. (2023), “Gaia Data Release 3: Gaia scan-angle-dependent signals and spurious periods”, *A&A*, 674, A25. DOI: 10.1051/0004-6361/202245353.
- El-Badry, K. et al. (2024), Gaia astrometric-orbit catalogue selection-function generative modelling, *Open Journal of Astrophysics*, 7. DOI: 10.33232/001c.125461. arXiv:2411.00088. Source: <https://github.com/kareemelbadry/gaiamock>. **Contains `al_bias_binary`, the published scan-angle- and separation-dependent along-scan response benchmarked in `docs/gaiamock_benchmark.md`.**
- Rowell, N. et al. (2021), “Gaia Early Data Release 3: Point and line spread function modelling and calibration”, *A&A*, 649, A11. DOI: 10.1051/0004-6361/202039448.
