# Literature gap investigation

## Scope and standard of evidence

This note records the literature status of the proposed RAA Orbit Model as of 2026-08-17. Its purpose is to identify a defensible research gap, not to protect an earlier novelty claim. Claims are separated into **established prior art**, a **surviving candidate gap**, and **open verification tasks**. No priority claim is made, and the manuscript should not use "first", "novel", or "unprecedented" unless a later exhaustive review justifies it.

## Established prior art

### Combined astrometry, relative astrometry, and spectroscopy are established

BINARYS (Leclerc et al. 2023) combines Hipparcos/Gaia absolute astrometry, relative astrometry, and radial velocities, and its RV interface can distinguish measurements of the two components. Therefore, "Gaia + relative astrometry + SB2" is not itself a gap.

orvara (Brandt et al. 2021) likewise fits combinations of radial velocity, absolute astrometry, and relative astrometry and performs posterior inference. Its Gaia information is based on Hipparcos--Gaia absolute astrometry / scan information rather than a marginal-resolution image-response likelihood. It is nevertheless important prior art for dynamical-mass inference from heterogeneous orbit data.

### Gaia DR3 orbital processing did not target partially resolved doubles

The Gaia DR3 non-single-star astrometric processing (Halbwachs et al. 2023) fits Keplerian orbital solutions, but partially resolved double stars were discarded in a preliminary treatment. This supports the relevance of the transition regime, but it is not evidence that no external method has treated that regime.

### Scan-angle-dependent close-pair systematics are established

Holl et al. (2023) showed that close source structure can generate scan-angle-dependent image-parameter biases and spurious time-series signals in Gaia. The existence of scan-angle-dependent astrometric bias is therefore not a new result of this project.

### A nonlinear, resolution-aware along-scan response is already published

El-Badry et al. (2024) released `gaiamock`, whose `al_bias_binary` predicts the one-dimensional measured position from the peak of the blended along-scan flux profile. The response is nonlinear in projected separation and light ratio and tends to the ordinary flux-weighted photocentre at small separation.

The equal-width Gaussian response used by this project is mathematically the same response family. `docs/gaiamock_benchmark.md` records a direct source-level and numerical comparison. The two implementations agree to the published solver tolerance over their common single-peak domain. The measurement-response idea itself must therefore not be claimed as new.

### The blended-source/resolvability mathematics is now broader in the literature

Penoyre (2026), *RAS Techniques and Instruments* 5, rzaf062, derives analytical expressions and numerical recipes for the effective position and resolvability of blended Gaussian point sources. Importantly for Gaia, the paper treats an elongated PSF and shows that the problem reduces to a one-dimensional blend with an **orientation-dependent effective width**. A correction to several signs in the published equations appeared as rzag016 on 2026-03-19 and the corrected article should be used.

Consequently, the constant-width, equal-profile surrogate in this repository should be presented as a **controlled restricted case** of a more general published blended-source framework. The exact equal-width mode-splitting calculation remains useful for code validation and for defining the domain of the frozen synthetic experiments, but it is not a claim to general resolvability theory.

### Bayesian Gaia epoch astrometry + radial-velocity inference is established

Baycroft, Faria & Delisle (2026), arXiv:2606.24132, present open-source Bayesian analysis of Gaia epoch astrometry, both alone and simultaneously with radial velocities, using `kima` and diffusive nested sampling. Thus neither "Bayesian Gaia epoch astrometry" nor "Gaia epoch astrometry + RV joint inference" is a remaining gap.

The source reviewed for this project establishes those capabilities. It does **not by itself establish** the more specific combination sought here: a marginal-resolution blended-source measurement response fitted inside the likelihood together with independent resolved relative astrometry and both component RV curves. That narrower combination remains a search target rather than a confirmed priority claim.

### Gaia DR4 PSF modelling is substantially more realistic than the present surrogate

Rowell et al. (2026), *A&A* 708, A174, describe the PSF model deployed in the processing of Gaia DR4. It includes drift-scan effects and calibrated dependences on source colour and focal-plane position. This supersedes Rowell et al. (2021) as the primary instrument-model reference for DR4-era motivation.

The present fixed Gaussian width therefore cannot be interpreted as a physical Gaia LSF/PSF width. The next instrument-facing model must be driven by the released DR4 measurement/calibration products rather than by tuning the surrogate transition scale.

## Literature capability matrix

The table below states only capabilities supported by the sources reviewed; a blank or qualified entry is not evidence of absence.

| Work | Gaia/Hipparcos astrometry in orbit inference | Independent relative astrometry | RVs / SB2 | Resolution-aware blended response | Response fitted inside target likelihood? | Main role |
|---|---|---|---|---|---|---|
| Brandt et al. 2021, orvara | yes | yes | yes | not the focus | no evidence in reviewed source | combined dynamical orbit inference |
| Leclerc et al. 2023, BINARYS | yes | yes | yes, including component-labelled RVs | not the focus | no evidence in reviewed source | combined binary orbit inference |
| Holl et al. 2023 | Gaia image/scan diagnostics | no | no | analytical scan-angle bias | not a target orbit likelihood | close-pair systematics |
| El-Badry et al. 2024, gaiamock | simulated Gaia epochs and Gaia-like fitting | no independent resolved orbit data in the target fit | RV simulation capability | yes | **no: response is generative; fit path is photocentre-based** | catalogue selection-function simulation |
| Penoyre 2026 | repeated astrometric observations analytically | no orbit combination | no | **yes, including elongated/orientation-dependent Gaussian PSF** | no joint binary-orbit likelihood | blended-source position/resolvability theory |
| Baycroft et al. 2026, kima | **yes, epoch astrometry** | not established by the source reviewed here | **yes, joint RV** | not established by the source reviewed here | Bayesian orbit likelihood, but resolution-aware blending not established | Gaia epoch + RV posterior inference |
| This project, candidate target | Gaia-like epoch AL | **yes** | **both SB2 curves** | **yes** | **yes** | measurement-aware target-level combined orbit inference |

## Surviving candidate gap

After the 2026 literature update, the scientifically defensible question is no longer "can a scan-angle-dependent blended response be modelled?" and no longer "can Gaia epoch astrometry be fitted jointly with RVs?" Both are established.

The **surviving candidate gap** is the intersection of capabilities:

1. use a resolution-aware blended-source response as the **fitted measurement model**, not merely as the generator of synthetic Gaia data;
2. perform that inference at **target level** with Gaia epoch along-scan data;
3. constrain the same physical orbit simultaneously with **independent resolved relative astrometry**;
4. include **both component RV curves of an SB2** rather than a single host-star RV time series;
5. propagate measurement-model uncertainty/misspecification into **individual component masses, parallax, and light ratio**.

The project should be presented as testing whether this particular intersection is scientifically useful. It should not claim that any one ingredient is individually new.

### Strongest current wording for the paper

> We investigate a candidate inference gap at the intersection of existing methods: fitting a published resolution-aware blended-source astrometric response at target level while simultaneously constraining the orbit with independent resolved relative astrometry and both SB2 velocity curves, and quantifying how measurement-model assumptions propagate into component masses.

### Claims that should not appear

- "first resolution-aware Gaia binary model";
- "new Gaia blended-source response";
- "first joint Gaia epoch astrometry and RV fit";
- "new resolvability criterion" without explicitly restricting it to the equal-width validation case and acknowledging Penoyre (2026);
- any translation of the surrogate `sigma` or `a/sigma` values into Gaia angular-resolution thresholds.

## Why the present experiments are still useful

The 23,000-fit native full-sky experiment and the 11,040-fit matched-transit control remain valid as **baseline equal-width synthetic experiments**. They demonstrate how a deliberately controlled measurement-model mismatch propagates through the joint orbit fit. They do not validate the response against real Gaia observations and should not be rewritten as if they used the newer absolute-astrometry, posterior-sampling, or unequal-width extensions added later to the code.

The strongest next experiment is therefore not another ordinary sky grid. It is a **response-misspecification study**: inject a more realistic profile family (at minimum unequal and/or orientation-dependent widths motivated by Penoyre 2026 and Rowell et al. 2026), fit with successively simpler response models, and determine when the inference-side advantage survives model error. This directly tests whether the candidate gap has practical value rather than only an advantage under matched injection/recovery.

## Remaining verification tasks

The gap is still a candidate. Before submission, carry out citation-chain searches around the following topics and record any counterexamples:

- target-level Gaia epoch astrometry + **resolved relative astrometry** + spectroscopy after 2024;
- SB2-specific Gaia epoch orbit inference;
- marginally resolved / blended Gaia epoch likelihoods rather than generative bias models;
- methods using Penoyre's blended-source response inside an orbital posterior;
- DR4-oriented methods that fit image/window-level data for binaries.

If a paper is found that already combines all five capabilities above, the project should pivot to a narrower technical gap, such as robust response misspecification, colour-dependent component profiles, or component-mass bias under an instrument-calibrated likelihood.

## DR4 motivation and real-data status

Gaia DR4 is expected in **December 2026**. The official DR4 expected-content page states that lower-level individual observations will be released, subject to processing and validation. Exact product names should be taken from the released datamodel rather than guessed in advance.

Until suitable epoch/image products are public, the resolution-aware inference must be described as synthetic/methodological. DR3 catalogue diagnostics such as `ipd_frac_multi_peak` and `ipd_gof_harmonic_amplitude` can provide qualitative consistency checks, but they are not direct measurements of this surrogate's mode count or likelihood residuals.

## References used for the gap assessment

- Brandt, T. D. et al. (2021), "orvara: An Efficient Code to Fit Orbits using Radial Velocity, Absolute, and/or Relative Astrometry", *AJ* / arXiv:2105.11671.
- Leclerc, A. et al. (2023), "Combining Hipparcos and Gaia data for the study of binaries: the BINARYS tool", *A&A* 672, A82, DOI 10.1051/0004-6361/202244144.
- Halbwachs, J.-L. et al. (2023), "Gaia Data Release 3: Astrometric binary star processing", *A&A* 674, A9, DOI 10.1051/0004-6361/202243969.
- Holl, B. et al. (2023), "Gaia Data Release 3: Gaia scan-angle-dependent signals and spurious periods", *A&A* 674, A25, DOI 10.1051/0004-6361/202245353.
- El-Badry, K. et al. (2024), "A generative model for Gaia astrometric orbit catalogs...", *Open Journal of Astrophysics* 7, DOI 10.33232/001c.125461, arXiv:2411.00088.
- Penoyre, Z. (2026), "The position and resolvability of blended point sources", *RAS Techniques and Instruments* 5, rzaf062, DOI 10.1093/rasti/rzaf062; corrected by rzag016, DOI 10.1093/rasti/rzag016.
- Baycroft, T. A., Faria, J. P. & Delisle, J.-B. (2026), "Bayesian analysis of Gaia epoch astrometry and radial velocities with kima", arXiv:2606.24132.
- Rowell, N. et al. (2026), "Gaia Data Release 4: Modelling of drift-scan related effects in Gaia's point spread function", *A&A* 708, A174, DOI 10.1051/0004-6361/202558618.
