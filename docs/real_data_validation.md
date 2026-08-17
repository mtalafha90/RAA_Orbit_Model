# Real-data validation status

The validation ladder now has two distinct real-data stages.

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

## V6b — DR3 catalogue consistency: still limited

The existing `dr3_validation.py` machinery remains useful for catalogue-level consistency and candidate selection, including Thiele-Innes/Campbell conversion and IPD diagnostics. A catalogue-level DR3 test cannot directly validate the close-pair along-scan response because DR3 does not publish the stellar epoch astrometry needed by the likelihood, and partially resolved doubles were filtered upstream of the DR3 astrometric-binary processing.

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
