# RAA Orbit Model

**Resolution-Aware Joint Orbit Inference for Marginally Resolved Gaia Binaries**

This repository contains the physical model, controlled response-fidelity experiments, real-binary validation, manuscript, compact result products and tests for a study of how marginal-resolution astrometric response propagates into binary-star dynamical inference.

## Scientific question

How much bias is introduced into component masses, parallax and light ratio when a luminous marginally resolved binary is treated as an ordinary flux-weighted photocentre, and how much response fidelity is required when Gaia-like astrometry is combined with independent resolved relative astrometry and both SB2 radial-velocity curves?

The controlled response hierarchy is:

- **M0** — ordinary photocentre;
- **M1** — published equal-width nonlinear blended response benchmarked against the `gaiamock` response family;
- **M2** — idealized finite-elongation, scan-orientation-dependent Gaussian response motivated by Penoyre (2026).

The response widths are research-surrogate parameters, not calibrated Gaia PLSF resolution scales.

## Current validation status

- **V0–V5 — completed:** internal/synthetic validation, response-fidelity hierarchy and targeted controls.
- **V6a — completed and reproducible:** real visual astrometry + both SB2 RV curves for GJ 765.2 / HIP 96656. The 110-constraint fit gives reduced chi-square 1.04629, total mass 1.58969 Msun and orbital parallax 35.4426 +/- 2.2445 mas.
- **V6b — query-ready:** Gaia DR3 catalogue/IPD validation for HIP 96656. No target-specific DR3 row is claimed until a validated export is ingested.
- **V7 — pending:** direct measurement-level Gaia epoch/image response validation when the relevant public products are available.

The real V6a test validates the Newtonian resolved-astrometry+SB2 core only; it does not validate the marginal-resolution Gaia response.

## Start here

The manuscript entry point is:

```text
main.tex
```

The numbered scientific documentation begins at:

```text
docs/00_DOCUMENT_INDEX.md
```

Canonical documentation is physically numbered `00` through `15`, from the literature audit and methodology through the real-data validation programme.

## Repository structure

```text
RAA_Orbit_Model/
├── main.tex
├── README.md
├── LICENSE
├── CITATION.cff
├── raa.cls
├── raa.bst
├── raa_orbit_refs.bib
├── raa_realdata_refs.bib
├── pyproject.toml
├── docs/                         # numbered 00–15 scientific documentation
├── data/
│   └── real/gj7652/              # exact V6a input + provenance
├── manuscript/                   # paper sections and figure snippets
├── figures/                      # project figures
├── results/
│   ├── frozen/                   # compact synthetic manuscript-support tables
│   ├── real/gj7652/              # frozen V6a outputs
│   └── dr3_validation/           # DR3 query/validation products
├── scripts/                      # experiment and validation runners
├── src/raa_orbit_model/          # reusable Python package
├── tests/                        # regression/scientific tests
└── .github/workflows/tests.yml   # Python and manuscript CI
```

## Reproduce the GJ 765.2 V6a result

The exact legacy input is committed at `data/real/gj7652/GL765_Test1.csv`. Its provenance and interpretation are documented next to the data. Run:

```bash
python -m pip install -e ".[test,analysis]"
python scripts/fit_gl765_visual_sb2.py
```

The real-data implementation is in `src/raa_orbit_model/real_data.py`, regression tests are in `tests/test_real_data.py`, and the compact reference outputs are committed under `results/real/gj7652/`.

The legacy header is preserved for provenance but is not silently treated as authoritative metadata. In particular, the header parallax 54.27 mas is tested explicitly as a fixed-parallax control rather than imposed on the free solution.

## Run the full test suite

```bash
python -m pip install -e ".[test,analysis]"
pytest -q
```

GitHub Actions tests Python 3.10, 3.11 and 3.12. CI also compiles the RAA manuscript with XeLaTeX/latexmk so missing citations, broken includes and LaTeX errors are caught independently of the Python tests.

## Main experiment runners

Response-fidelity hierarchy:

```bash
python scripts/run_response_fidelity.py --help
```

External-information control:

```bash
python scripts/run_external_information_control.py --help
```

Matched-M2 high-statistics control:

```bash
python scripts/run_matched_m2_control.py --help
```

Full-sky and matched-transit supporting studies:

```bash
python scripts/run_sky_position_scan.py --help
python scripts/run_matched_n_control.py --help
```

## Gaia DR3 target bridge

Print the preferred HIP 96656 cross-match query with:

```bash
python scripts/validate_gl765_dr3.py --show-query
```

After exporting a validated Gaia Archive result as CSV:

```bash
python scripts/validate_gl765_dr3.py gj765_dr3.csv
```

The current repository intentionally does not invent a source ID, RUWE, IPD value or NSS orbit before that target row is retrieved and verified.

## Scientific claim boundaries

The project does not claim a new Gaia blended-source response, a new general astrometry+RV framework, the first Gaia+SB2 mass inference, or a calibrated Gaia resolution threshold. The conservative prior-art/claim audit is maintained in `docs/01_literature_gap.md`.

The strongest contribution is the quantified response-fidelity problem: how measurement-model misspecification propagates into stellar dynamical masses, parallax and light ratio as the independent orbit information changes.

## Citation and license

Citation metadata are provided in `CITATION.cff`. The source code and repository materials are released under the MIT License; third-party journal style files and cited external data/software remain subject to their own terms where applicable.
