# Real-data validation: status of V6 and V7

`docs/methodology.md` section 11 defines a validation ladder. Steps V0 to V5 are synthetic and are exercised by the test suite. The final two steps involve real Gaia measurements. **Neither has been executed.** This note records exactly why, and what has been built so that each can be run without overstating current validation.

## V6 — DR3 catalogue consistency

**Status: blocked on data access in the development environment, not on implementation.**

The Gaia archive was unreachable from the environment in which the validation harness was developed, so no DR3 row was retrieved there and **no claim of consistency with real Gaia data is made anywhere in this project.**

What exists instead is the machinery to run V6 as soon as a DR3 export is available:

- `src/raa_orbit_model/dr3_validation.py`
  - `NSS_ASTROMETRIC_ORBIT_QUERY` — ADQL selecting DR3 astrometric orbits, joined to `gaia_source` for duplicity diagnostics.
  - `thiele_innes_to_campbell` / `campbell_to_thiele_innes` — DR3 publishes astrometric orbits as Thiele-Innes constants, not Campbell elements, so this conversion is required before any comparison. It is written in the same North/East convention asserted in `tests/test_orbit_conventions.py`.
  - `campbell_table`, `marginally_resolved_candidates` — conversion and selection helpers.
- `scripts/validate_against_dr3.py` — run `--show-query`, export from the archive, then pass the CSV.
- `tests/test_dr3_validation.py` — offline regression tests.

The conversion is verified two ways without network access: it round-trips against itself for a range of geometries, and it reproduces the ellipse drawn by this project's own orbit model to numerical precision.

### Why these DR3 columns

A catalogue-level test cannot see the along-scan measurements, so it cannot test the measurement model directly. It can only test diagnostics that are qualitatively associated with partially resolved structure.

- `ipd_frac_multi_peak` is the percentage of successful-IPD windows in which the Gaia IPD processing identified more than one peak. The surrogate also tracks whether its simplified blended profile is single- or multi-peaked, so the two quantities are conceptually related. They are **not equivalent**: the Gaia statistic additionally depends on the real PLSF, detection thresholds, windowing, gating, background, calibration, and IPD processing.
- `ipd_gof_harmonic_amplitude` measures scan-angle-dependent structure in the IPD goodness of fit. A close pair can generate such angular dependence, but this statistic is not a direct likelihood residual of the RAA surrogate.

These columns therefore provide possible falsification/selection diagnostics, not a direct calibration target.

### Known limitation of this step

Halbwachs et al. (2023) state that partially resolved doubles were filtered upstream of the DR3 astrometric-binary processing. The DR3 orbit catalogue is therefore depleted in exactly the regime this project targets, and a null result in V6 would partly reflect selection. Any V6 analysis must account for that fact and should distinguish the all-source IPD diagnostics from the much more selected NSS orbital-solution sample.

## V7 — DR4 epoch and image validation

**Status: not possible yet. Gaia DR4 has not been released.**

The official Gaia release scenario currently places DR4 in approximately 2026 (the release-scenario page states not before mid-2026). The exact release date should not be hard-coded here until ESA publishes it.

The official Gaia DR4 **expected-content** page already provides provisional archive product names. Of direct relevance to this project are:

- `epoch_astrometry` — individual astrometric measurements for sources in the main astrometric processing;
- `epoch_image` — preprocessed, sky-projected individual CCD sample values;
- `rvs_epoch_data_double` — FoV-transit-level information for double-lined RVS transits;
- planned non-single-star, mass, and multiplicity products.

These names are therefore not guesses. However, the same official page states that the content is under development and changes can be expected, while the detailed DR4 processing documentation and final data model are still forthcoming. Code that ingests DR4 should consequently isolate these interfaces behind adapters and verify them against the released schema rather than assuming today's expected table layout is final.

This creates the largest remaining literature/implementation blind spot. The public DR4 material reviewed for this project does **not establish** whether the eventual DPAC NSS orbit likelihood internally uses a physical marginal-resolution PLSF/image response for close binaries. Until the processing papers and documentation are available, the manuscript should state neither that DPAC has such an inference nor that it does not.

### Why `epoch_image` matters beyond scalar epoch astrometry

The current single-coordinate surrogate is intentionally invalid once the blended profile becomes genuinely multi-peaked. If DR4 `epoch_image` provides the expected CCD sample information at useful fidelity, an image/sample-domain likelihood becomes the natural extension:

`orbit -> component positions + fluxes -> PLSF/profile samples -> CCD data`.

That would avoid assigning one astrometric coordinate to a transit whose brightness profile has multiple competing maxima. It is therefore a scientifically stronger long-term endpoint than extending the present one-coordinate surrogate into the resolved regime.

The DR4 instrument reference has also changed. Rowell et al. (2026, A&A 708, A174) describe the PLSF model deployed in DR4 processing, including drift-scan effects and calibrated dependences on source colour and focal-plane position. Any measurement-level implementation should be checked against released calibration products rather than treating a fixed Gaussian width as the Gaia instrument model.

## What can be claimed today

Every headline quantitative result in the manuscript is synthetic. The equal-width surrogate reproduces the public `gaiamock` blended response over their common validity domain, but the deeper source audit also shows that `gaiamock` already contains a response-aware stellar astrometry-plus-primary-RV **forward predictor**. The measurement response and that forward combination are therefore not novelty claims.

Likewise, Liu et al. (2024) already used a partially resolved Gaia response inside orbital inference for the binary asteroid (4337) Arecibo. What remains under investigation is the narrower stellar inverse problem and, more importantly, the response-fidelity question: whether marginal-resolution model error biases individual stellar masses/parallax and how accurate the response must be for calibrated inference.

Neither the `gaiamock` benchmark nor synthetic recovery is a test against real Gaia measurements. The project should continue to say this plainly until V6 and, ultimately, V7 are actually executed.
