# Real-data validation: status of V6 and V7

`docs/methodology.md` section 11 defines a validation ladder. Steps V0 to V5 are synthetic and are exercised by the test suite. The final two steps involve real Gaia measurements. **Neither has been executed.** This note records exactly why, and what has been built so that each can be run without further development.

## V6 — DR3 catalogue consistency

**Status: blocked on data access, not on implementation.**

The Gaia archive is unreachable from the environment this work was carried out in. The gateway refuses the connection outright:

```text
gea.esac.esa.int:443   403 to CONNECT (policy denial)
www.cosmos.esa.int:443 403 to CONNECT (policy denial)
vizier.cds.unistra.fr:443 403 to CONNECT (policy denial)
```

So no DR3 row has been retrieved, and **no claim of consistency with real Gaia data is made anywhere in this project.**

What exists instead is the machinery to run V6 as soon as the catalogue is reachable:

- `src/raa_orbit_model/dr3_validation.py`
  - `NSS_ASTROMETRIC_ORBIT_QUERY` — ADQL selecting DR3 astrometric orbits, joined to `gaia_source` for the duplicity diagnostics, ordered by scan-angle-dependent fit structure.
  - `thiele_innes_to_campbell` / `campbell_to_thiele_innes` — DR3 publishes astrometric orbits as Thiele-Innes constants, not Campbell elements, so this conversion is required before any comparison. It is written in the same North/East convention asserted in `tests/test_orbit_conventions.py`.
  - `campbell_table`, `marginally_resolved_candidates` — conversion and selection.
- `scripts/validate_against_dr3.py` — run `--show-query`, export once from the archive, then pass the CSV.
- `tests/test_dr3_validation.py` — 13 tests, all offline.

The conversion is verified two ways without any network access: it round-trips against itself for a range of geometries, and it reproduces the ellipse drawn by this project's own orbit model to 1e-9 mas.

### Why these DR3 columns

A catalogue-level test cannot see the along-scan measurements, so it cannot test the measurement model directly. It can test the *symptoms* the surrogate predicts. Two published columns are the natural targets:

- `ipd_frac_multi_peak` — the percentage of windows in which image parameter determination found a second peak. This project's surrogate already computes a per-transit mode count, so it predicts this quantity almost directly.
- `ipd_gof_harmonic_amplitude` — the amplitude of the scan-angle-dependent component of the fit quality. A fixed pair produces an even-harmonic signature in scan angle, which is a structural prediction of any resolution-aware response.

Testing against these turns the surrogate from something falsifiable only in principle into something falsifiable against roughly 1.8 billion published sources, before DR4 ships.

### Known limitation of this step

Halbwachs et al. (2023) state that partially resolved doubles were filtered upstream of the DR3 astrometric-binary processing. The DR3 orbit catalogue is therefore *depleted* in exactly the regime this project targets, and a null result in V6 would be partly a selection effect rather than evidence against the model. Any V6 analysis must model that selection, or restrict itself to the diagnostics above, which are published for all sources rather than only for those with orbital solutions.

## V7 — DR4 epoch and image validation

**Status: not possible yet. Gaia DR4 has not been released.**

The expected release date is 2 December 2026. Until then there are no general epoch astrometry products to fit, so the measurement-level test this project is ultimately aimed at cannot be attempted.

Two caveats on the product names used elsewhere in this repository, neither of which could be checked against the DR4 datamodel from this environment:

- `epoch_astrometry` is expected to carry, per transit, the observation time, the along-scan position, the scan angle and the parallax factor. That is exactly the observable this project models, and it is the reason the along-scan channel now includes absolute astrometry.
- The names `epoch_image` and `rvs_epoch_data_double` are **unconfirmed**. Searches returned conflicting evidence, including a possible `residual_image` product and a possible `EPOCH_PARAMETERS_RVS_DOUBLE` retrieval type. Check these against the released datamodel before relying on them.

## What can be claimed today

Every quantitative result in this project is synthetic. The surrogate has been shown to reproduce a published measurement response (`docs/gaiamock_benchmark.md`), and the inference machinery has been shown to recover injected parameters. Neither of those is a test against real Gaia measurements, and the manuscript should continue to say so plainly.
