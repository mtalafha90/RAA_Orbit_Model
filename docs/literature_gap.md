# Literature gap investigation

## Scope and standard of evidence

This note records the literature and public-software audit for the proposed RAA Orbit Model as of 2026-08-17. Its purpose is to identify a defensible scientific gap rather than to protect an earlier novelty claim. Statements are separated into established prior art, source-level software observations, a surviving **candidate gap**, and remaining blind spots. No priority claim is made. The manuscript should not use "first", "novel", or "unprecedented" for the combined method unless a later exhaustive review establishes that wording.

The central distinction throughout this audit is between (i) a resolution-aware response used to generate or diagnose data and (ii) that response being fitted **inside the inverse orbit likelihood**. A second distinction is between a single host-star RV series and the two component velocity curves of an SB2.

## Established prior art

### Combined astrometry, relative astrometry, and SB2 spectroscopy are established

BINARYS (Leclerc et al. 2023) combines Hipparcos/Gaia absolute astrometry, independent relative astrometry, and radial velocities. Its RV interface can label the measured component and therefore supports SB1/SB2 information. Consequently, neither heterogeneous joint orbit fitting nor "Gaia + relative astrometry + SB2" is a research gap.

This point is reinforced by Chevalier et al. (2023), who used BINARYS to combine Gaia DR3 non-single-star astrometric solutions with SB2 information from SB9 and APOGEE, obtaining component masses and luminosities for dozens of systems. Gaia astrometry + SB2 -> individual stellar masses is therefore established.

orvara (Brandt et al. 2021) independently establishes posterior inference from radial velocity, absolute astrometry, and relative astrometry. Its Gaia use is not the marginal-resolution epoch-image response considered here, but it is essential prior art for dynamical-mass inference from heterogeneous data.

Other Gaia DR3 work also infers component properties from unresolved astrometric orbits. ESMORGA (Campo et al. 2024) derives possible relative orbits and individual masses using Gaia DR3 astrometric, photometric, and spectroscopic information. Bailer-Jones & Kreidberg (2026) infer component masses for unresolved binaries from Gaia astrometric orbits plus three-band photometry and a mass--flux relation. Thus component-mass inference from a Gaia photocentre orbit and flux information is not itself new.

### BINARYS provides the closest conceptual bridge--and explicitly leaves the Gaia transition regime open

BINARYS is especially relevant because it already embeds non-trivial measurement response physics for Hipparcos transit data when the secondary contributes significant light, while also combining relative astrometry and component-labelled RVs. The broad architecture

`binary dynamics -> instrument response -> astrometric measurement -> joint orbit fit`

is therefore established.

For Gaia, however, Leclerc et al. explicitly limit the treatment of luminous systems that are not cleanly resolved. Their discussion notes that partially resolved systems with non-negligible secondary light require modelling the flux contamination/line-spread-function fitting and await more detailed Gaia information. This is a much stronger motivation for the present project than a generic absence claim: a mature combined-orbit framework identifies the same luminous partially resolved Gaia regime as untreated.

### Gaia DR3 NSS did not model this regime in the orbital pipeline

Halbwachs et al. (2023) state that a preliminary treatment was used to **discard partially resolved double stars** before the DR3 astrometric binary models were applied. Therefore the population of interest was explicitly problematic for the published DR3 NSS processing. This establishes relevance, not priority: external methods can still exist.

Holl et al. (2023) further established that close source structure produces scan-angle-dependent Gaia image-parameter biases and spurious time-series signals. Scan-angle-dependent close-pair bias is therefore prior art, not a result unique to this project.

### A nonlinear Gaia-like blended-source response is established

El-Badry et al. (2024) released `gaiamock`. Its public `al_bias_binary` routine predicts the one-dimensional measured position from the peak of a blended along-scan light profile as a nonlinear function of projected separation, mass ratio, light ratio, and an effective angular scale. The equal-width response used in the baseline experiments of this repository is the same response family; `docs/gaiamock_benchmark.md` records the numerical agreement over the common single-peak domain.

The deeper source audit also found `predict_astrometry_and_rvs_simultaneously(...)` in the public `gaiamock` source. That routine uses the nonlinear blended astrometric response and predicts a **primary-star RV curve** for the same stellar binary. Therefore it is not safe to claim that resolution-aware stellar Gaia modelling has never been combined with radial velocities. In the audited source this is a forward-prediction path. We did **not identify** a corresponding target-level posterior that simultaneously fits the response-aware Gaia epochs, external resolved relative astrometry, and two SB2 velocity curves. This is a source-audit observation, not proof that no such implementation exists elsewhere.

### Response-aware orbital inference itself is not unique to stars or to this project

Liu et al. (2024) modelled the binary asteroid (4337) Arecibo using Gaia Focused Product Release astrometry. Their analysis fits the heliocentric centre-of-mass orbit and then the relative binary orbit, uses a Lindegren-style unresolved/partially-resolved/resolved astrometric response, includes external stellar-occultation relative positions, and infers orbital/physical quantities including flux and component-mass information. The analysis is sequential rather than a stellar SB2 joint posterior, and there are no two stellar RV curves, but it removes any defensible claim that "putting a partially resolved Gaia response inside an orbital MCMC" is new in general.

### The blended-source/resolvability mathematics is broader than our baseline surrogate

Penoyre (2026), *RAS Techniques and Instruments* 5, rzaf062, derives analytical expressions and numerical recipes for blended Gaussian sources and treats an elongated Gaia-like PSF. The two-dimensional problem reduces to a one-dimensional blend with an **orientation-dependent effective width**. The paper's correction rzag016 (2026-03-19) changes signs in Eq. 18 and Appendix equations B7, B8, and B13; the corrected article must be the source for any implementation.

Consequently, the constant-width equal-profile response in this repository is a controlled restricted case, not new general resolvability theory and not a calibrated Gaia PSF/LSF. Penoyre also points out that changing orientation can add information about separation and light ratio, motivating the next response-fidelity experiment.

### Bayesian Gaia epoch astrometry + RV inference is established

Baycroft, Faria & Delisle (2026) present Bayesian Gaia epoch-astrometry inference, both alone and simultaneous with RVs, in `kima` using diffusive nested sampling. Thus Bayesian Gaia epoch astrometry and Gaia epoch astrometry + RV are established capabilities.

A source-level audit of the current `RVGAIAmodel` found a Gaia AL model with position, proper motion, parallax and Keplerian terms, plus an optional phenomenological scan-angle-dependent harmonic signal. In that **specific audited Gaia+RV path**, we identified one Gaia data object and one RV data vector and did not identify double-lined handling. This must not be generalized into a statement that `kima` as a whole cannot model SB2 data. More importantly, the audited Gaia likelihood uses an orbital AL signal/phenomenological scan-angle terms rather than the physical blended-image operator sought here.

### Octofitter closes most of the infrastructure gap, but not the audited response gap

The current public Octofitter source contains Gaia DR4 epoch-astrometry likelihood machinery as well as relative astrometry, interferometry, radial-velocity, sampling, and joint-fitting infrastructure. Its audited `GaiaDR4Astrom` likelihood predicts AL residuals from sky-plane offsets, proper motion, parallax, and orbital perturbations,

`eta = dRA* sin(psi) + dDec cos(psi) + parallax_factor_AL * parallax`,

with no marginal-resolution blended-profile operator identified in that likelihood path. This means multi-dataset Bayesian infrastructure is not a contribution of this project; the potentially missing element is the nonlinear measurement operator for a luminous close pair.

### Gaia DR4 PSF modelling is substantially more realistic than our surrogate

Rowell et al. (2026), *A&A* 708, A174, describe the PSF model deployed in Gaia DR4 processing, including drift-scan effects and calibrated dependences on source colour and focal-plane position. The present fixed Gaussian width therefore cannot be interpreted as a physical Gaia LSF/PSF width.

## Capability matrix after the deeper audit

The table states only capabilities supported by the sources or public code paths reviewed. "Not identified" means exactly that--it is not proof of absence.

| Work / implementation | Gaia epoch AL | Physical marginal-resolution response | RV | SB2 | Independent resolved relative astrometry | Response inside target inverse orbit fit? |
|---|---:|---:|---:|---:|---:|---:|
| Gaia DR3 NSS (Halbwachs et al. 2023) | internal processing | problematic partially resolved doubles were filtered | combined NSS products exist | some spectroscopy products | catalogue context | not for the filtered transition population |
| BINARYS (Leclerc et al. 2023) | Gaia catalogue/absolute astrometry; Hipparcos transit path | non-trivial Hipparcos response; luminous partially resolved Gaia case explicitly limited | yes | **yes** | **yes** | not established for marginally resolved Gaia |
| Chevalier et al. 2023 | Gaia DR3 NSS solutions | standard Gaia DR3 astrometric solution | yes | **yes** | some resolved systems separately | no marginal-response fit |
| `gaiamock` (El-Badry et al. 2024) | simulated epoch AL | **yes** | **primary RV forward predictor identified** | two fitted RV curves not identified in audited path | not identified in target fit | response-aware forward model identified; full joint inverse not identified |
| Liu et al. 2024, Arecibo | **yes, Gaia FPR** | **yes** | no stellar RV | no | **yes, occultation positions** | **yes, in relative-orbit inference; asteroid and sequential analysis** |
| Penoyre 2026 | repeated astrometric observations analytically | **yes, including orientation-dependent elongated Gaussian** | no | no | no | no joint binary-orbit posterior |
| `kima` Gaia+RV (Baycroft et al. 2026) | **yes** | physical blend operator not identified in audited path | **yes** | not identified in the audited Gaia+RV path | not established there | Bayesian joint Gaia+RV likelihood, but physical blending not identified |
| Octofitter current Gaia DR4 likelihood | **yes** | not identified in audited Gaia likelihood | yes via package modules | broad RV machinery | **yes** | joint framework exists; physical Gaia blend operator not identified |
| Bailer-Jones & Kreidberg 2026 | Gaia DR3 astrometric orbit | ordinary unresolved photocentre premise | no requirement | no | no | component-mass posterior from astrometry+photometry, not marginal-response epochs |
| This project, candidate target | Gaia-like epoch AL | **yes** | **yes** | **both curves** | **yes** | **yes, intended shared stellar posterior** |

## Surviving candidate gap

After this audit, the project must not be sold as a new response, a new Gaia+RV fitter, a new combined-orbit framework, a new route to component masses, or the first orbital use of a partially resolved response.

The **surviving candidate gap** is the complete stellar intersection:

1. epoch-level Gaia AL inference with a physical marginal-resolution blended-source observation operator **inside the likelihood**;
2. one shared stellar dynamical model constrained simultaneously by that Gaia likelihood and **independent resolved relative astrometry**;
3. **both component velocity curves of an SB2** in the same inference;
4. propagation of response assumptions and response misspecification into `M1`, `M2`, parallax, light ratio, and posterior coverage.

Based on the published literature and public implementations audited here, we did not identify a demonstrated stellar target-level analysis occupying all four items simultaneously. This is still a **candidate gap**, not a claim of priority.

### Strongest current manuscript wording

> Based on the published literature and publicly available implementations examined here, we did not identify a stellar target-level analysis in which a physical marginal-resolution Gaia blended-source response is placed directly inside an epoch-astrometric likelihood and inferred simultaneously with independent resolved relative astrometry and both radial-velocity curves of an SB2. Related ingredients exist separately: BINARYS combines astrometry, relative astrometry, and SB2 velocities but explicitly leaves luminous partially resolved Gaia systems untreated; `gaiamock` implements a resolution-aware stellar Gaia response and an astrometry-plus-primary-RV forward predictor; Bayesian Gaia epoch-astrometry-plus-RV inference exists in `kima`; and response-aware orbital inference has been demonstrated for the binary asteroid Arecibo. The scientific question is therefore not whether these ingredients exist separately, but how marginal-resolution measurement physics propagates into stellar dynamical-mass inference when all component-resolving constraints are combined.

### Claims that should not appear

- "first resolution-aware Gaia binary model";
- "new Gaia blended-source response";
- "first resolution-aware Gaia+RV model";
- "first joint Gaia epoch astrometry and RV fit";
- "first Gaia astrometry + SB2 mass inference";
- "first orbital MCMC with a partially resolved Gaia response";
- "new resolvability criterion" without restricting it to the equal-width code-validation case and citing Penoyre (2026);
- any conversion of surrogate `sigma` or `a/sigma` into a Gaia angular-resolution threshold.

## The stronger scientific gap: response fidelity and mass bias

A simple checklist combination is vulnerable to becoming incremental as DR4-era tools mature. The stronger scientific question is instead:

> **What biases in individual stellar dynamical masses and parallax arise when marginal-resolution Gaia measurements are treated as an ordinary photocentre, and how accurate must the astrometric response model be before those biases and coverage failures disappear?**

This question remains meaningful even if another package later gains the same data types. It also turns the current equal-width work into a baseline rather than the endpoint.

The decisive experiment is a response hierarchy on the **same synthetic observations**:

- `M0`: ordinary photocentre;
- `M1`: Lindegren/`gaiamock`-family equal-width peak response;
- `M2`: corrected Penoyre-style orientation-dependent effective-width response;
- later `M3`: response parameters treated as uncertain nuisance parameters rather than asserted.

The primary outputs should be fractional bias in `M1`, `M2`, parallax and light fraction, together with 68%/95% posterior coverage. `Delta chi2` remains useful diagnostically but is not the principal scientific endpoint.

A particularly interesting identifiability question follows from the external data. SB2 amplitudes strongly constrain the mass ratio, while resolved relative astrometry constrains the relative orbit. A well-measured visual SB2 can therefore act as a **calibrator of the mapping from known projected component separation to Gaia's measured blended coordinate**. A possible later question is whether an ensemble of such systems can empirically constrain response nuisance parameters rather than assuming them.

## Why the frozen experiments remain useful

The 23,000-fit native full-sky experiment and 11,040-fit matched-transit control remain valid as baseline equal-width synthetic experiments. They establish that, in a controlled single-peak surrogate, an incorrect photocentre likelihood can produce coherent joint-fit mismatch and that nominal sky dependence is strongly affected by transit count. They do **not** establish a physical Gaia resolution threshold, do not validate the response against real Gaia images, and do not prove that the advantage survives response misspecification.

They should remain frozen and labelled as the `M1` matched-surrogate baseline. The next experiment must deliberately generate data from a more general response than the simple fitted model.

## Gaia DR4: opportunity and principal overlap risk

The official Gaia DR4 expected-content page already lists provisional product names including:

- `epoch_astrometry`: individual astrometric measurements;
- `epoch_image`: preprocessed, sky-projected individual CCD sample values;
- `rvs_epoch_data_double`: FoV-transit-level information for double-lined RVS transits;
- non-single-star and mass/multiplicity products in the planned DR4 archive.

The same official page states that the expected content remains under development and subject to processing/validation changes. Detailed DR4 processing documentation and the final data model are forthcoming. Therefore it is incorrect to say that these product names are unconfirmed guesses; they are **official expected names**, but they are not yet final released interfaces.

This is also the largest remaining novelty blind spot. The public material reviewed here does **not establish** whether the eventual DPAC NSS orbit likelihood uses a physical marginal-resolution image/PLSF response internally. It would be unsafe to claim that DPAC does or does not implement such an inversion until the DR4 processing papers and documentation are public.

`epoch_image` also motivates a future image-domain extension. Once the blended profile becomes genuinely multi-peaked, reducing a transit to one astrometric coordinate is intrinsically fragile. A later likelihood can work directly with CCD samples rather than inventing a single-coordinate continuation.

## Remaining verification tasks before submission

- Re-run citation-chain searches immediately before submission for stellar Gaia epoch AL + resolved relative astrometry + SB2.
- Check all DR4 NSS/epoch-astrometry processing papers as soon as they appear.
- Search specifically for papers that cite Penoyre (2026) and use its response inside an orbital posterior.
- Re-audit current `gaiamock`, `kima`, Octofitter, BINARYS and other public orbit software shortly before submission because these packages are active.
- Search for image/window-level binary inference once DR4 previews/data-model documentation are released.

If a direct stellar precedent appears, the paper should pivot without resistance to the response-fidelity / posterior-coverage question, which remains scientifically useful independent of software priority.

## References used for this gap assessment

- Brandt, T. D. et al. (2021), `orvara`, *AJ* 162, 186, DOI 10.3847/1538-3881/ac042e.
- Leclerc, A. et al. (2023), BINARYS, *A&A* 672, A82, DOI 10.1051/0004-6361/202244144.
- Halbwachs, J.-L. et al. (2023), Gaia DR3 astrometric binary processing, *A&A* 674, A9, DOI 10.1051/0004-6361/202243969.
- Holl, B. et al. (2023), scan-angle-dependent Gaia signals, *A&A* 674, A25, DOI 10.1051/0004-6361/202245353.
- Chevalier, S. et al. (2023), "Binary masses and luminosities with Gaia DR3", *A&A* 678, A19, DOI 10.1051/0004-6361/202347111.
- Liu, Z. et al. (2024), "Asteroid (4337) Arecibo: Two ice-rich bodies forming a binary -- Based on Gaia astrometric data", *A&A* 688, L23, DOI 10.1051/0004-6361/202450586.
- El-Badry, K. et al. (2024), `gaiamock` / generative Gaia orbit-catalogue model, *Open Journal of Astrophysics* 7, DOI 10.33232/001c.125461.
- Campo, P. P. et al. (2024), ESMORGA methodology, *A&A* 682, A12.
- Penoyre, Z. (2026), "The position and resolvability of blended point sources", *RAS Techniques and Instruments* 5, rzaf062, DOI 10.1093/rasti/rzaf062; correction rzag016, DOI 10.1093/rasti/rzag016.
- Baycroft, T. A., Faria, J. P. & Delisle, J.-B. (2026), Gaia epoch astrometry + RV with `kima`, arXiv:2606.24132.
- Bailer-Jones, C. A. L. & Kreidberg, L. (2026), "Component masses in stellar and sub-stellar binaries from Gaia astrometry and photometry", *A&A* 708, A249, DOI 10.1051/0004-6361/202659004.
- Rowell, N. et al. (2026), Gaia DR4 PSF modelling, *A&A* 708, A174, DOI 10.1051/0004-6361/202558618.
- ESA Gaia COSMOS, "Gaia DR4 content", official expected-content table (content explicitly marked under development / subject to change).
