# Documentation index

This file defines the canonical reading order for the research documentation in this repository. The numbers are stable document identifiers; the underlying filenames are kept descriptive so existing links and experiment provenance remain valid.

## Foundation and scientific scope

**00 — Documentation index**  
`docs/00_DOCUMENT_INDEX.md`  
This index and the repository reading order.

**01 — Literature gap and prior-art audit**  
`docs/literature_gap.md`  
Deep literature/software audit, prior art, novelty cautions, and the surviving candidate scientific gap.

**02 — Methodology**  
`docs/methodology.md`  
Mathematical formulation, likelihood structure, inference strategy, and validation ladder.

**03 — Orbit conventions**  
`docs/orbit_conventions.md`  
Coordinate, node, argument-of-periastron, RV, and visual-orbit conventions used throughout the code and manuscript.

**04 — Gaia scanning law**  
`docs/gaia_scanning_law.md`  
Mission-grounded scan scheduling, GOST/gaiascanlaw use, conventions, and limitations.

**05 — gaiamock benchmark**  
`docs/gaiamock_benchmark.md`  
Direct benchmark of the equal-width response against the published gaiamock implementation.

## Baseline and control experiments

**06 — First injection experiment [historical validation]**  
`docs/first_injection_experiment.md`  
Original uniform-angle injection/recovery pilot. Preserved for provenance; superseded for scientific interpretation by the mission-grounded experiments.

**07 — Sky-position study design**  
`docs/sky_position_study.md`  
Definition of the sky-grid experiment and scan-geometry diagnostics.

**08 — Full-sky results**  
`docs/full_sky_results.md`  
Frozen 46-position full-sky equal-width result summary.

**09 — Matched-transit-count control**  
`docs/matched_n_control.md`  
Design of the matched-N control used to separate transit-count effects from scan geometry.

**10 — Native versus matched-N results**  
`docs/native_vs_matched_n53_results.md`  
Frozen comparison between native schedules and the matched 53-transit schedules.

**11 — Multi-peak validity**  
`docs/multi_peak_validity.md`  
Definition and validation of the single-peak/multi-peak boundary for the restricted response surrogate.

## Response-fidelity programme

**12 — Response-fidelity experiment**  
`docs/response_fidelity_experiment.md`  
M0/M1/M2 experiment design, parameter grid, and scientific interpretation rules.

**13 — Response-fidelity results**  
`docs/response_fidelity_results.md`  
Frozen summary of the central response-fidelity results used in the manuscript.

**14 — Targeted response controls**  
`docs/targeted_response_controls.md`  
External-information-strength and high-statistics matched-M2 controls.

## Real-data validation

**15 — Real-data validation**  
`docs/real_data_validation.md`  
Validation status for GJ 765.2, the Gaia DR3 catalogue/IPD bridge, and the planned DR4 epoch/image test.

## Canonical manuscript order

The journal manuscript is assembled by `main.tex` in the following order:

1. `manuscript/introduction.tex`
2. `manuscript/model.tex`
3. `manuscript/experiments.tex`
4. `manuscript/results.tex`
5. `manuscript/real_data_validation.tex`
6. `manuscript/discussion_conclusions.tex`

## Repository policy

- `docs/` contains scientific reasoning, experiment design, frozen interpretation, and validation notes.
- `src/raa_orbit_model/` contains reusable library code.
- `scripts/` contains command-line experiment and validation runners.
- `tests/` contains regression and scientific-consistency tests.
- `results/frozen/` contains compact tables that support manuscript claims and must remain reproducible.
- `results/dr3_validation/` contains target-specific DR3 query/validation products.
- `figures/` contains project figures only; journal-template example figures are not retained.
- RAA journal style files `raa.cls` and `raa.bst` are kept at the repository root so `main.tex` compiles directly.
