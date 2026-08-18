# 01 — Literature gap and prior-art audit

## Scope

This document records the defensible scientific gap for the RAA Orbit Model. It is intentionally conservative: related capabilities are treated as prior art whenever they are demonstrated in published literature or public software. The manuscript must not use priority language such as “first”, “novel”, or “unprecedented” for the combined method unless a later exhaustive review establishes it.

## Established prior art

The following ingredients are already established and therefore are **not** the research gap:

- heterogeneous joint orbit inference from astrometry, resolved relative astrometry and radial velocities (e.g. BINARYS, orvara);
- Gaia astrometry + SB2 spectroscopy for individual stellar masses and luminosities;
- nonlinear separation-dependent Gaia-like blended-source response (`gaiamock`);
- response-aware orbital inference in a broader astronomical context, including the Gaia analysis of binary asteroid (4337) Arecibo;
- analytical blended-Gaussian resolvability and orientation-dependent elongated-profile response (Penoyre 2026);
- Bayesian Gaia epoch astrometry + radial-velocity inference (`kima`);
- general joint-fitting infrastructure containing Gaia epoch astrometry, relative astrometry and RV likelihoods (e.g. current Octofitter infrastructure).

BINARYS is particularly important because it already has the broad architecture

`binary dynamics -> instrument response -> astrometric measurement -> joint orbit fit`

for Hipparcos and explicitly identifies luminous partially resolved Gaia systems as a regime requiring more detailed Gaia response information. Gaia DR3 NSS processing also treated partially resolved doubles as problematic and filtered them upstream of the published astrometric-binary solutions. These facts motivate the project but do not establish priority.

## Response theory is not the contribution

The equal-width M1 model used here is the same response family implemented in `gaiamock`; it is an implementation/benchmark, not new response theory. The finite-elongation M2 model is an idealized application of the orientation-dependent effective-width framework of Penoyre (2026), not a calibrated Gaia PSF/LSF model. Gaia DR4 processing uses a substantially richer calibrated PLSF with colour, focal-plane, time and drift-scan dependences.

Therefore:

- no surrogate width may be described as a Gaia angular-resolution threshold;
- the equal-width single/multi-peak criterion is a validity boundary for the surrogate only;
- M2 is a response-fidelity experiment, not a hardware calibration model.

## Surviving candidate gap

Based on the published literature and publicly available implementations examined for this project, we did not identify a demonstrated **stellar target-level** analysis that simultaneously contains all of the following:

1. epoch-level Gaia astrometric inference with a physical marginal-resolution blended-source observation operator **inside the likelihood**;
2. the same stellar dynamical model constrained by **independent resolved relative astrometry**;
3. **both radial-velocity curves of an SB2** in the same inference;
4. explicit propagation of response misspecification into individual masses, parallax, light ratio and posterior calibration.

This remains a **candidate gap**, not a proof of absence or priority claim.

## Stronger and more durable scientific question

The paper should be framed around response fidelity rather than a checklist novelty claim:

> What biases in individual stellar dynamical masses, parallax and light ratio arise when marginal-resolution Gaia measurements are treated as an ordinary photocentre, and how accurate must the astrometric response model be before those biases disappear for a given strength of external orbit information?

That question remains scientifically meaningful even if another package later supports the same data types.

## Capability summary

| Capability | Established elsewhere? | Role here |
|---|---:|---|
| joint astrometry + RV inference | yes | infrastructure, not novelty |
| relative astrometry + SB2 | yes | external dynamical constraint |
| Gaia-like nonlinear blend response | yes | M1 benchmark |
| orientation-dependent blend response | yes analytically | M2 fidelity surrogate |
| response-aware orbital inference | yes in other contexts | prior art |
| stellar Gaia response + both SB2 curves + resolved astrometry in one target likelihood | not identified in audited sources | candidate intersection |
| systematic mass/parallax bias versus response fidelity and external-information strength | central question here | principal contribution |

## Claims that should not appear

Do not claim:

- “first resolution-aware Gaia binary model”;
- “new Gaia blended-source response”;
- “first resolution-aware Gaia+RV model”;
- “first Gaia astrometry + SB2 mass inference”;
- “first orbital MCMC with a partially resolved Gaia response”;
- “new resolvability criterion” without explicitly restricting it to the equal-width validation surrogate;
- any physical Gaia resolution threshold derived from the research width parameters.

## DR3/DR4 relevance

DR3 provides catalogue-level IPD and NSS diagnostics but not the general stellar epoch AL measurements required for a direct M0/M1/M2 likelihood comparison. DR4 expected-content documentation lists `epoch_astrometry`, `epoch_image`, and `rvs_epoch_data_double`; these products are the natural future measurement-level test. The final DPAC DR4 binary-processing response must not be guessed before the relevant public processing documentation is available.

## Submission-time verification

Immediately before submission, repeat the literature/source audit for BINARYS, `gaiamock`, `kima`, Octofitter, Gaia DR4 processing papers and recent stellar dynamical-mass work. If another study occupies the full candidate intersection, the manuscript remains viable because its primary contribution is the quantified response-fidelity/mass-bias problem rather than priority.
