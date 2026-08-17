# Real-data validation: status of V6 and V7

`docs/methodology.md` section 11 defines a validation ladder. Steps V0 to V5 are synthetic and are exercised by the test suite. The final two steps involve real Gaia measurements. **Neither has been executed.** This note records exactly why, and what has been built so that each can be run without further development.

## V6 — DR3 catalogue consistency

**Status: blocked on data access, not on implementation.**

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

Current Gaia/DR4 sources describe DR4 as expected in **December 2026**. This repository does not assume a specific day. The official expected-content page states that lower-level individual observations are planned, but the released datamodel should be treated as the authority for exact table/product names.

Until release, this project should avoid hard-coding scientific claims around provisional names such as `epoch_image` or `rvs_epoch_data_double`. The relevant requirement is functional rather than nominal: target-level validation needs the epoch observation times, scan geometry, astrometric measurements or window/image information, calibration metadata, and (for the SB2 application) component-resolved spectroscopy or external SB2 radial velocities.

The DR4 instrument reference has also changed. Rowell et al. (2026, A&A 708, A174) describe the PLSF model deployed in DR4 processing, including drift-scan effects and calibrated dependences on source colour and focal-plane position. Any measurement-level RAA implementation should be checked against those released calibration products rather than treating a fixed Gaussian width as the Gaia instrument model.

## What can be claimed today

Every headline quantitative result in the manuscript is synthetic. The equal-width surrogate has been shown to reproduce the published `gaiamock` blended response over their common validity domain, and the inference machinery has been shown to recover controlled injections. Penoyre (2026) supplies a broader published treatment of blended-source position/resolvability, so the present constant-width model should be viewed as a restricted baseline rather than a new resolvability theory.

Neither the `gaiamock` benchmark nor synthetic recovery is a test against real Gaia measurements. The project should continue to say this plainly until V6 and, ultimately, V7 are actually executed.
