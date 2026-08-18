# 00 — Documentation index

This is the canonical reading order for the scientific documentation. The filenames are now physically numbered so GitHub displays them in research order rather than alphabetically by topic.

| No. | Canonical document | Purpose |
|---:|---|---|
| 00 | `00_DOCUMENT_INDEX.md` | Reading order, repository policy, validation map |
| 01 | `01_literature_gap.md` | Prior-art audit, claim boundaries, surviving candidate gap |
| 02 | `02_methodology.md` | Mathematical model, response hierarchy, validation ladder |
| 03 | `03_orbit_conventions.md` | Node, periastron, coordinate and RV conventions |
| 04 | `04_gaia_scanning_law.md` | Nominal scan schedules, angle convention and limitations |
| 05 | `05_gaiamock_benchmark.md` | Equal-width response benchmark and interpretation |
| 06 | `06_first_injection_experiment.md` | Historical development pilot; not manuscript evidence |
| 07 | `07_sky_position_study.md` | Full-sky supporting experiment design |
| 08 | `08_full_sky_results.md` | Frozen full-sky baseline interpretation |
| 09 | `09_matched_n_control.md` | Matched-transit-count control design |
| 10 | `10_native_vs_matched_n53_results.md` | Native versus matched-53 results |
| 11 | `11_multi_peak_validity.md` | Single-coordinate validity boundary |
| 12 | `12_response_fidelity_experiment.md` | M0/M1/M2 response-fidelity experiment |
| 13 | `13_response_fidelity_results.md` | Central response-fidelity results |
| 14 | `14_targeted_response_controls.md` | External-information and 100-seed controls |
| 15 | `15_real_data_validation.md` | GJ 765.2 V6a, DR3 bridge, DR4 plan |

## Validation ladder

- **V0–V5:** synthetic/internal validation and response-fidelity controls.
- **V6a — completed:** real resolved astrometry + both SB2 curves for GJ 765.2/HIP 96656. Exact input, parser, runner, tests and frozen results are committed.
- **V6b — query-ready:** Gaia DR3 catalogue/IPD bridge; target-specific numerical row is not yet frozen and no value is guessed.
- **V7 — pending:** direct public Gaia epoch/image response validation.

## Manuscript assembly

`main.tex` is the canonical paper entry point. It currently assembles:

1. `manuscript/01_introduction.tex`
2. `manuscript/model.tex`
3. `manuscript/experiments.tex`
4. `manuscript/results.tex`
5. `manuscript/04_result_figures.tex`
6. `manuscript/real_data_validation.tex`
7. `manuscript/05a_real_figures.tex`
8. `manuscript/discussion_conclusions.tex`

The scientific section numbering is controlled by LaTeX; auxiliary figure snippets are named by their insertion point.

## Repository policy

- `docs/`: current numbered scientific reasoning and frozen interpretation.
- `data/real/`: exact real inputs with provenance; source metadata are preserved but not silently treated as authoritative priors.
- `src/raa_orbit_model/`: reusable physical/inference code.
- `scripts/`: reproducible command-line runners.
- `tests/`: regression and scientific-consistency tests.
- `results/frozen/`: compact synthetic tables directly supporting manuscript claims.
- `results/real/`: compact real-data validation outputs.
- `results/dr3_validation/`: target-specific DR3 query/validation products.
- `figures/`: project figures only.
- `raa.cls` and `raa.bst`: journal files required to compile `main.tex`.

Superseded unnumbered documentation files are not part of the canonical set and are removed once internal references are migrated.
