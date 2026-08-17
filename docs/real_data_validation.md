# Real-data validation status

The validation ladder now has three distinct real-data stages.

## V6a — real resolved-astrometry + SB2 core check: executed

A legacy measurement subset for **GJ 765.2 = HD 186922 = HIP 96656** has been fitted with the physical Newtonian two-body core, resolved relative astrometry, and both SB2 velocity curves. The input contains 11 visual/speckle relative positions and 44 paired CORAVEL epochs for each RV component, but no Gaia epoch astrometry.

The joint fit uses 110 scalar constraints and 10 free physical parameters and gives

- chi2 = 104.629
- dof = 100
- reduced chi2 = 1.046
- P = 11.7284 +/- 0.0732 yr
- e = 0.24888 +/- 0.01005
- i = 81.834 +/- 1.371 deg
- ascending-node branch Omega = 289.07 +/- 3.27 deg
- relative-orbit omega = 251.79 +/- 2.25 deg
- total mass = 1.5897 Msun
- orbital parallax = 35.44 +/- 2.24 mas

The later combined interferometric/spectroscopic solution of Balega et al. (2007, A&A 464, 635; DOI 10.1051/0004-6361:20066224) gives P = 11.919 yr, e = 0.240, i = 80.2 deg, Omega = 293.0 deg, omega = 250.0 deg, masses 0.831 and 0.763 Msun, and orbital parallax 31.0 +/- 0.5 mas. The legacy-subset total mass differs from the later total mass (1.594 Msun) by about 0.27%.

This is a **consistency/implementation check**, not an independent astrophysical determination, because the legacy CORAVEL data overlap with those used in the later combined solution.

The legacy file header contains parallax = 54.27 mas. When that value is held fixed, the fit degrades to chi2 = 168.94 (reduced chi2 = 1.673 for 101 dof). Holding the later orbital parallax 31.0 mas fixed gives chi2 = 108.70 (reduced chi2 = 1.076). The joint resolved-astrometry/SB2 data therefore force the physical scale away from the inconsistent legacy header value.

### Scope

V6a validates the non-Gaia physical orbit engine on a real stellar binary. It does **not** validate the marginal-resolution Gaia response.

## V6b — Gaia DR3 catalogue/IPD consistency: target workflow implemented, target row not yet retrieved here

DR3 can provide a scientifically useful intermediate test even though it cannot supply the general stellar epoch along-scan measurements required to fit the RAA response hierarchy directly.

The target-specific implementation is now in:

- `src/raa_orbit_model/dr3_target.py`
- `scripts/validate_gl765_dr3.py`
- `tests/test_dr3_target.py`

The query uses the SIMBAD position of HD 186922 / HIP 96656, approximately RA = 294.7765558 deg and Dec = +76.4220233 deg. It deliberately does **not** use the incorrect coordinates stored in the legacy GL765 input header.

Run

```bash
python scripts/validate_gl765_dr3.py --show-query
```

and execute the printed ADQL in the Gaia Archive, then export the result as CSV and run

```bash
python scripts/validate_gl765_dr3.py gj765_dr3.csv
```

The target query retrieves, where available:

- source ID and coordinate separation from the external target position
- parallax and proper motion
- RUWE, along-scan observation counts, astrometric GoF/chi2, and excess noise
- `ipd_frac_multi_peak`
- `ipd_gof_harmonic_amplitude` and `ipd_gof_harmonic_phase`
- scan-direction moment strengths and mean directions
- `duplicated_source`, `non_single_star`, and `has_epoch_rv`
- every matching `nss_two_body_orbit` solution and its Thiele-Innes constants

The official DR3 documentation states that `ipd_frac_multi_peak` is the percentage of successful-IPD windows for which the IPD algorithm identified more than one peak. `ipd_gof_harmonic_amplitude` measures the scan-angle-dependent amplitude of the IPD goodness of fit and can indicate non-isotropic spatial structure such as a partially resolved binary. Its phase has a physical but non-trivial relation to binary position angle that depends on the resolution regime. These diagnostics are therefore highly relevant to RAA, but they are **not equivalent** to the mode classification or residuals of the research surrogate.

The ordinary unresolved-photocentre benchmark from the Balega et al. orbit is also encoded. With M1 = 0.831 Msun, M2 = 0.763 Msun, a_rel = 189 mas, and Delta V = 0.65 mag,

- secondary mass fraction B = 0.47867
- beta_V = 0.35465, used only as a temporary optical light-fraction proxy
- predicted M0 photocentre semi-major axis = |beta_V - B| a_rel = 23.44 mas

`beta_V` must not be described as the Gaia G-band light fraction. The catalogue comparison becomes stronger once a component-resolved G-band flux ratio is available.

If an astrometric NSS row exists, its A/B/F/G Thiele-Innes constants are converted to Campbell elements and compared with the external visual-SB2 orbit. If no NSS row exists, that result is **not** by itself a falsification: DR3 NSS is a selected sample, and partially resolved doubles were filtered upstream of the DR3 astrometric-binary processing.

### Current execution status

The exact Gaia DR3 source row for HIP 96656 has not been inserted into this repository. In the current development environment, dynamic Gaia Archive/VizieR catalogue queries are not reachable even though the public documentation is reachable. No RUWE, IPD value, DR3 source ID, or NSS solution is therefore inferred or invented here. V6b should remain labelled **query-ready / data-pull pending** until the exported target row is ingested and archived.

## V7 — Gaia DR4 epoch/image validation: pending public release

The official Gaia DR4 expected-content page lists the products required for direct measurement-level validation:

- `epoch_astrometry`: individual astrometric measurements for all sources treated in the main astrometric processing (DL2)
- `epoch_image`: pre-processed, sky-projected individual CCD sample values for those sources (DL1)
- `rvs_epoch_data_double`: FoV-transit-level information for double-lined RVS transits (DL2)

Official source: https://www.cosmos.esa.int/web/gaia/dr4

The same page states that DR4 is based on 66 months of observations from 2014-07-25 10:30 UTC to 2020-01-20 22:00 UTC, with reference epoch J2017.5, and that the expected release volume is about 500 TB. The content is explicitly described as under development and may change before release.

The DR4 overview still labels the release **Coming up**. The official release-scenario page currently states DR4 is **not before mid 2026**; no firm public release date is stated there.

## GJ 765.2 as a DR4 feasibility target

Using the Balega et al. (2007) orbit over the official 66-month DR4 interval gives a sky-projected component separation ranging from about 24.7 to 199.1 mas. Using the published V-band component magnitude difference (0.65 mag) only as a temporary light-fraction proxy gives beta_V = 0.3547.

For the current equal-width research baseline alpha = 50 mas, the exact surrogate mode-splitting threshold is |d_AL| about 128.5 mas at that light fraction. On a uniform time/scan-angle envelope over the DR4 interval:

- 24.1% of time-angle combinations are multi-peak in the surrogate
- for 58.0% of the interval, at least some scan orientations can be multi-peak
- at the most favorable epoch, 55.3% of scan orientations cross the threshold

These are **not Gaia transit statistics** and must not be presented as such. The real test requires the released DR4 epoch products. The 50 mas width is a research surrogate, not a calibrated Gaia PLSF resolution.

If an observed transit is genuinely multi-peaked, `epoch_image` is the scientifically preferred endpoint because the likelihood can be written in the sample/image domain instead of assigning an artificial unique astrometric coordinate.
