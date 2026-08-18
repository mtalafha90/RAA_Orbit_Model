# 15 — Real-data validation

The real-data programme has three distinct stages and strict claim boundaries.

## V6a — GJ 765.2 visual + SB2 physical-core validation: completed

Target: **GJ 765.2 = HD 186922 = HIP 96656**.

The exact legacy input is versioned at:

`data/real/gj7652/GL765_Test1.csv`

It contains 11 visual/speckle relative positions and 44 `Va` + 44 `Vb` CORAVEL radial-velocity measurements. The legacy header is preserved for provenance but is not treated as authoritative catalogue metadata.

Reproducibility components:

- parser/fitter: `src/raa_orbit_model/real_data.py`;
- runner: `scripts/fit_gl765_visual_sb2.py`;
- regression tests: `tests/test_real_data.py`;
- frozen results: `results/real/gj7652/`;
- orbit figure: `figures/05_gl765_visual_orbit.svg`;
- SB2 figure: `figures/06_gl765_sb2_rv.svg`.

The free-parallax fit uses 110 scalar constraints and 10 free physical parameters and gives

- chi2 = 104.62895;
- dof = 100;
- reduced chi2 = 1.04629;
- P = 11.72839 +/- 0.07319 yr;
- e = 0.248881 +/- 0.010052;
- i = 81.8338 +/- 1.3705 deg;
- ascending-node branch = 289.072 +/- 3.270 deg;
- relative-orbit omega = 251.790 +/- 2.246 deg;
- M1 = 0.78232 +/- 0.02817 Msun;
- M2 = 0.80737 +/- 0.02657 Msun;
- total mass = 1.58969 Msun;
- orbital parallax = 35.4426 +/- 2.2445 mas;
- gamma = -4.12478 +/- 0.05694 km/s.

Balega et al. (2007) give a later combined total mass of 1.594 Msun and orbital parallax 31.0 +/- 0.5 mas. The legacy-subset total mass differs by about 0.27%. Because the CORAVEL measurements overlap data used in the later combined solution, this is an implementation/consistency check rather than an independent astrophysical measurement.

Parallax controls:

- fixed legacy-header parallax 54.27 mas: chi2 = 168.94353, reduced chi2 = 1.67271;
- fixed Balega parallax 31.0 mas: chi2 = 108.70277, reduced chi2 = 1.07627.

V6a validates the Newtonian visual/SB2 physical core. It does **not** validate the Gaia marginal-resolution response.

## V6b — Gaia DR3 catalogue/IPD bridge: query-ready

The target workflow is implemented in `src/raa_orbit_model/dr3_target.py` and `scripts/validate_gl765_dr3.py`. Target identification uses the Gaia DR3 Hipparcos-2 best-neighbour cross-match for HIP 96656 rather than the incorrect coordinates in the legacy header.

The archived ADQL query retrieves, where available, parallax, proper motion, RUWE, astrometric diagnostics, `ipd_frac_multi_peak`, `ipd_gof_harmonic_amplitude`, `ipd_gof_harmonic_phase`, source flags and any `nss_two_body_orbit` solution.

Using the published masses and Delta V = 0.65 mag only as a temporary optical proxy gives

- mass fraction B = 0.47867;
- beta_V = 0.35465;
- ordinary unresolved-photocentre benchmark = 23.44 mas.

`beta_V` is not the Gaia G-band light fraction.

No validated target-specific DR3 row is currently frozen in the repository. Therefore V6b remains **query-ready / data-pull pending**. No DR3 source ID, RUWE, IPD value or NSS orbit is inferred or invented.

## V7 — Gaia DR4 measurement-level response validation: pending

Official Gaia DR4 expected-content documentation lists `epoch_astrometry`, `epoch_image`, and `rvs_epoch_data_double` among the planned products. These products will permit a direct measurement-level response test once publicly available. No firm DR4 release date is asserted here; official release status should be checked again immediately before manuscript submission or analysis.

GJ 765.2 is a useful prospective target because the Balega orbit spans approximately 24.7–199.1 mas in projected component separation over the 66-month DR4 observation interval. With the current research surrogate alpha = 50 mas and beta_V proxy, the equal-width mode-splitting threshold is about 128.5 mas. A uniform time/orientation envelope gives 24.1% multi-peak combinations, 58.0% of the interval with at least some multi-peak orientations, and a maximum orientation fraction of 55.3% at one epoch.

These percentages are **not Gaia transit statistics**. They are only a pre-release feasibility envelope. Real V7 inference must use released epoch products, and genuinely multi-peaked observations should be handled in the image/sample domain rather than forced into one astrometric coordinate.
