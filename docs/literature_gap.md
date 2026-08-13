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

### Marginal-resolution generative modelling

El-Badry et al. (2024) modelled Gaia DR3 astrometric-orbit selection using a forward generative calculation that includes marginally resolved binaries. Their purpose was primarily catalogue selection-function modelling rather than a target-level joint posterior using independent resolved astrometry and both SB2 velocity curves.

### Gaia PSF/LSF calibration

Gaia image parameter determination rests on calibrated PSF/LSF models; Rowell et al. (2021) document the EDR3 PSF/LSF calibration. Therefore a final image-level RAA model must not treat a Gaussian profile as an exact Gaia instrument model.

## Candidate gap

After the targeted literature search carried out for this project, we did **not find** a published target-level framework explicitly documented as simultaneously fitting all of the following:

1. independent resolved relative astrometry of a luminous binary;
2. both SB2 radial-velocity curves;
3. Gaia epoch-level astrometric or image information;
4. a scan-angle- and resolution-aware Gaia measurement response that permits the measured location to differ from the simple unresolved photocentre;
5. full joint uncertainty propagation to orbital elements, component masses, distance/parallax, and passband flux ratio.

This is a **candidate methodological gap**, not a proven first-in-literature claim. A broader ADS/citation-chain review is required before any manuscript uses words such as “first”, “novel”, or “unprecedented”.

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
- El-Badry, K. et al. (2024), Gaia astrometric-orbit catalogue selection-function generative modelling, *Open Journal of Astrophysics*, 7. DOI: 10.33232/001c.125461.
- Rowell, N. et al. (2021), “Gaia Early Data Release 3: Point and line spread function modelling and calibration”, *A&A*, 649, A11. DOI: 10.1051/0004-6361/202039448.
