# RAA Orbit Model

**Resolution-Aware Joint Orbit Inference for Marginally Resolved Gaia Binaries**

This repository contains the code, controlled experiments, manuscript, frozen result summaries, and real-data validation work for a study of how marginal-resolution astrometric response affects binary-star dynamical inference.

## Scientific question

The project asks how much bias is introduced into component masses, parallax, and light ratio when a luminous marginally resolved binary is treated as an ordinary flux-weighted photocentre, and how much measurement-response fidelity is required when Gaia-like astrometry is combined with independent resolved relative astrometry and both SB2 radial-velocity curves.

Three response levels are used in the controlled experiments:

- **M0** — ordinary photocentre;
- **M1** — equal-width nonlinear blended response benchmarked against the published `gaiamock` implementation;
- **M2** — an idealised finite-elongation, scan-orientation-dependent Gaussian response motivated by Penoyre (2026).

The response widths used here are research-surrogate parameters, not calibrated Gaia PLSF resolution scales.

## Current scientific status

The completed synthetic programme includes a 720-fit response-fidelity hierarchy, a 270-fit external-information-strength control, a 100-seed matched-M2 control, and supporting full-sky/matched-transit experiments.

The real-data validation ladder is now:

- **V6a — completed:** resolved relative astrometry + both SB2 RV curves for GJ 765.2 / HIP 96656. The legacy-subset fit gives reduced chi-square 1.046, total mass 1.5897 Msun, and orbital parallax 35.44 +/- 2.24 mas. The total mass differs by about 0.27% from the later combined value of 1.594 Msun.
- **V6b — query-ready:** Gaia DR3 catalogue/IPD validation for HIP 96656, using the Hipparcos-2 best-neighbour cross-match, IPD diagnostics, and any available NSS two-body orbit. The external ordinary-photocentre benchmark is 23.44 mas when the published V-band light ratio is used only as a temporary proxy.
- **V7 — pending public DR4 epoch products:** direct measurement-level comparison of M0/M1/M2 using released epoch astrometry and, where needed, epoch-image samples.

DR3 catalogue diagnostics are an intermediate consistency/falsification layer; they do not replace an epoch-level response likelihood.

## Start here

The canonical manuscript entry point is:

```text
main.tex
```

The documentation has a stable numbered reading order in:

```text
docs/00_DOCUMENT_INDEX.md
```

The journal manuscript is assembled from:

```text
manuscript/introduction.tex
manuscript/model.tex
manuscript/experiments.tex
manuscript/results.tex
manuscript/real_data_validation.tex
manuscript/discussion_conclusions.tex
```

## Numbered documentation order

| No. | Document | Purpose |
|---:|---|---|
| 00 | `docs/00_DOCUMENT_INDEX.md` | Canonical document order and repository policy |
| 01 | `docs/literature_gap.md` | Literature/software audit and candidate gap |
| 02 | `docs/methodology.md` | Mathematical and validation methodology |
| 03 | `docs/orbit_conventions.md` | Orbit, node, RV, and coordinate conventions |
| 04 | `docs/gaia_scanning_law.md` | Mission-grounded scan schedule and caveats |
| 05 | `docs/gaiamock_benchmark.md` | Equal-width response benchmark |
| 06 | `docs/first_injection_experiment.md` | Historical uniform-angle validation pilot |
| 07 | `docs/sky_position_study.md` | Full-sky experiment design |
| 08 | `docs/full_sky_results.md` | Frozen full-sky baseline results |
| 09 | `docs/matched_n_control.md` | Matched-transit-count control design |
| 10 | `docs/native_vs_matched_n53_results.md` | Native versus matched-N frozen result |
| 11 | `docs/multi_peak_validity.md` | Single-/multi-peak validity boundary |
| 12 | `docs/response_fidelity_experiment.md` | M0/M1/M2 response-fidelity design |
| 13 | `docs/response_fidelity_results.md` | Frozen response-fidelity results |
| 14 | `docs/targeted_response_controls.md` | External-information and high-statistics controls |
| 15 | `docs/real_data_validation.md` | GJ 765.2, DR3 bridge, and DR4 validation plan |

The numbering is intentionally maintained in the index rather than renaming historical files, so frozen experiment provenance and existing citations/links remain stable.

## Repository structure

```text
RAA_Orbit_Model/
├── main.tex                         # canonical manuscript entry point
├── raa.cls                          # RAA journal class required for compilation
├── raa.bst                          # RAA bibliography style
├── raa_orbit_refs.bib               # main project references
├── raa_realdata_refs.bib            # real-binary reference(s)
├── README.md
├── pyproject.toml
├── docs/                            # numbered scientific documentation
├── manuscript/                      # manuscript section files
├── figures/                         # project figures only
├── results/
│   ├── frozen/                      # compact immutable manuscript-support tables
│   └── dr3_validation/              # target-specific DR3 query/validation products
├── scripts/                         # command-line experiment/analysis runners
├── src/raa_orbit_model/             # reusable Python package
├── tests/                           # regression and scientific-consistency tests
└── .github/workflows/tests.yml      # Python 3.10-3.12 CI matrix
```

The unrelated RAA example manuscript, example quasar bibliography, example figures, and unused alternate class have been removed from the research repository.

## Core implementation

The package includes Newtonian two-body dynamics, resolved tangent-plane astrometry, physically derived SB2 RV amplitudes, Gaia along-scan projection, M0/M1/M2 response dispatch, proper-motion/parallax terms, covariance whitening, deterministic joint fitting, posterior sampling, scanning-law adapters, synthetic data generation, response-fidelity experiments, sky-position studies, matched-transit controls, paired bias analysis, DR3 NSS conversion, and the target-specific GJ 765.2 DR3 workflow.

Important modules include:

```text
src/raa_orbit_model/kepler.py
src/raa_orbit_model/astrometry.py
src/raa_orbit_model/gaia.py
src/raa_orbit_model/model.py
src/raa_orbit_model/fit.py
src/raa_orbit_model/synthetic.py
src/raa_orbit_model/experiments.py
src/raa_orbit_model/response_fidelity.py
src/raa_orbit_model/bias_analysis.py
src/raa_orbit_model/dr3_validation.py
src/raa_orbit_model/dr3_target.py
src/raa_orbit_model/legacy_target.py
```

## Install and test

Core development environment:

```bash
python -m pip install -e ".[test]"
pytest -q
```

For scanning-law generation and analysis/plotting:

```bash
python -m pip install -e ".[test,scanlaw,analysis]"
```

The GitHub Actions workflow tests Python 3.10, 3.11, and 3.12.

## Main experiment runners

Response-fidelity hierarchy:

```bash
python scripts/run_response_fidelity.py \
  --schedule-file schedules/ra120_dec30_dr4.csv \
  --a-over-alpha-values 0.6 1.0 \
  --beta-values 0.25 \
  --beta-over-alpha-values 1.5 3.0 \
  --seeds 3 \
  --alpha-mas 50 \
  --output results/response_fidelity_pilot.csv
```

External-information control:

```bash
python scripts/run_external_information_control.py
```

High-statistics matched-M2 control:

```bash
python scripts/run_matched_m2_control.py
```

Mission-grounded baseline bias experiment:

```bash
python scripts/run_bias_scan.py \
  --ra-deg 120.0 --dec-deg 30.0 --release dr4 \
  --write-schedule schedules/ra120_dec30_dr4.csv \
  --seeds 3 --output results/bias_ra120_dec30_dr4.csv
```

Full-sky and matched-transit controls:

```bash
python scripts/run_sky_position_scan.py
python scripts/run_matched_n_control.py
python scripts/compare_native_matched_n.py
```

## Real-binary validation (V6a)

The GJ 765.2 legacy visual + SB2 fit reported in the manuscript is regenerated
from the committed measurements at `data/gl765_legacy.csv`:

```bash
python scripts/run_legacy_target_fit.py data/gl765_legacy.csv \
  --fixed-parallax-mas 54.27 31.0
```

It writes `results/frozen/legacy_target_fit.csv` and reproduces the manuscript
table to four or five significant figures, uncertainties included. See
`docs/real_data_validation.md`.

## Gaia DR3 validation

The generic DR3 NSS catalogue harness is:

```bash
python scripts/validate_against_dr3.py --show-query
```

The target-specific GJ 765.2 / HIP 96656 workflow is:

```bash
python scripts/validate_gl765_dr3.py --show-query
```

The exact archived ADQL is also stored at:

```text
results/dr3_validation/gj765_dr3_query.adql
```

After exporting the Gaia Archive result as CSV:

```bash
python scripts/validate_gl765_dr3.py gj765_dr3.csv \
  --output results/dr3_validation/gj765_dr3_summary.json
```

No target-specific DR3 numerical values should be quoted until a validated catalogue export has been ingested.

## Frozen results and reproducibility

`results/frozen/` contains compact summaries used directly in the manuscript. These files are retained intentionally and should not be overwritten by exploratory runs. New exploratory/raw outputs should be written elsewhere under `results/` and regenerated from the documented scripts.

The original uniform-angle pilot is retained only as historical code-validation provenance. Scientific interpretation uses the later mission-grounded experiments described in the numbered documentation.

## Manuscript compilation

From the repository root:

```bash
xelatex main.tex
bibtex main
xelatex main.tex
xelatex main.tex
```

The repository retains only the RAA style files required for this manuscript (`raa.cls` and `raa.bst`); unrelated journal example assets are intentionally excluded.

## Scope and caution

The current M1/M2 image-response families are controlled scientific surrogates. They establish inference behaviour under known response misspecification but do not claim to reproduce Gaia's calibrated PLSF or processing chain. DR3 can provide catalogue/IPD diagnostics and NSS orbit comparisons. Direct empirical validation of the response hierarchy requires released Gaia epoch astrometry and, for genuinely multi-peaked cases, sample/image-domain information.
