# RAA Orbit Model

**Resolution-Aware Astrometric Orbit Model for luminous binaries**

This repository is the working research repository for a candidate methodology aimed at the difficult transition between unresolved and resolved binary-star astrometry in Gaia.

## Scientific target

The project asks whether it is useful to fit a **resolution-aware blended-source response inside the orbit likelihood** for a luminous binary while simultaneously using:

1. Gaia-like along-scan epoch astrometry,
2. independent resolved relative astrometry, and
3. both radial-velocity curves of a double-lined spectroscopic binary (SB2),

with the measurement-model choice propagated to parallax, light ratio, and the two component masses.

This is deliberately an **intersection-of-capabilities** question. The individual ingredients are not claimed as new:

- BINARYS and orvara already combine astrometric, relative-astrometric, and radial-velocity information for dynamical orbit inference; BINARYS supports component-labelled RVs.
- Chevalier et al. (2023) already combine Gaia DR3 non-single-star astrometry with SB2 information to obtain individual stellar masses and luminosities.
- El-Badry et al. (2024) `gaiamock` already implements a scan-angle- and separation-dependent blended along-scan response. Its public source also contains a response-aware stellar astrometry + primary-RV forward predictor. This repository reproduces the equal-width response over the common single-peak domain; see `docs/gaiamock_benchmark.md`.
- Liu et al. (2024) already use a partially resolved Gaia response inside orbital inference for the binary asteroid (4337) Arecibo, together with external occultation constraints.
- Penoyre (2026) gives a more general treatment of blended Gaussian source positions and resolvability, including an elongated Gaia-like PSF with orientation-dependent effective width. The constant-width baseline used here is therefore a restricted case, not new resolvability theory.
- Baycroft, Faria & Delisle (2026) already provide Bayesian Gaia epoch-astrometry inference, including simultaneous radial velocities, in `kima`.
- Rowell et al. (2026) describe the substantially more realistic PLSF model deployed in Gaia DR4 processing.

The **surviving candidate gap** is narrow: based on the literature and public implementations audited, we did not identify a demonstrated *stellar target-level* inference that places a physical marginal-resolution Gaia response directly inside an epoch-astrometric likelihood while simultaneously constraining the same orbit with independent resolved relative astrometry and both SB2 velocity curves. **No novelty or priority claim is made.** The evidence matrix and remaining searches are maintained in `docs/literature_gap.md`.

The more durable scientific question is stronger than that checklist:

> **What biases in individual stellar dynamical masses and parallax arise when marginal-resolution Gaia measurements are treated as an ordinary photocentre, and how accurate must the response model be before those biases and posterior-coverage failures disappear?**

The current development therefore pivots from matched-surrogate `Delta chi2` alone to a **response-fidelity hierarchy** judged primarily by physical-parameter bias and, after deterministic validation, posterior coverage.

## Current implementation status

The prototype now contains:

- a physically constrained Newtonian two-body SB2 forward model;
- resolved relative astrometry in the tangent plane;
- physical RV semi-amplitudes derived from component masses;
- Gaia along-scan projection geometry;
- the published equal-width 1-D blended-response family used by `gaiamock`;
- an idealised finite-elongation Penoyre-style response with scan-orientation-dependent effective width;
- explicit recovery of the existing 1-D response in the Penoyre `beta_PSF -> infinity` limit;
- optional unequal component widths for an independent profile-shape misspecification control;
- a paired M0/M1/M2 response-fidelity experiment: photocentre, equal-width response, and finite-elongation response;
- absolute astrometric position, proper motion, and approximate parallax terms for controlled tests;
- full 2-D covariance whitening for resolved astrometry;
- a bounded deterministic joint fitter and normalized joint likelihood;
- a validation-grade affine-invariant posterior sampler;
- direct injection/recovery comparisons on identical synthetic realizations;
- a sky-position-dependent **GOST-derived Gaia Nominal Scanning Law** adapter through `gaiascanlaw`;
- portable CSV archiving/reloading of exact transit times and directional scan angles;
- configurable orbit, noise, resolution, light-fraction, and response-elongation axes;
- paired seed-by-seed model comparison and physical-bias summaries;
- a DR3 catalogue-validation harness for orbit conversion and IPD diagnostic selection;
- regression tests for orbital identities, scan-angle convention, response limits, schedule conversion, exact synthetic recovery, response misspecification, posterior sampling, and paired-analysis sign conventions.

The blended-image responses remain **prototype surrogates**, not Gaia's calibrated PSF/LSF. Their widths are research parameters. Penoyre (2026) shows how an idealised elongated Gaussian introduces an orientation-dependent effective width, while Rowell et al. (2026) document DR4 PLSF dependences on drift-scan motion, colour, and focal-plane position. The final instrument-facing model must follow released Gaia epoch/image and calibration products.

The original uniform-angle pilot is preserved in `docs/first_injection_experiment.md` as a code-validation result only. It is superseded for scientific bias experiments by the mission-grounded schedule implementation documented in `docs/gaia_scanning_law.md`. The 23,000-fit native full-sky and 11,040-fit matched-N result sets are frozen baseline equal-width experiments and are not retroactively reinterpreted as using later code extensions.

## Repository layout

```text
docs/
  literature_gap.md             Deep literature/software audit and surviving candidate gap
  gaiamock_benchmark.md         Benchmark against the published equal-width response
  response_fidelity_experiment.md  M0/M1/M2 response-fidelity design and pilot
  methodology.md                Mathematical and validation methodology
  real_data_validation.md       Status of validation steps V6 and V7
  orbit_conventions.md          Orientation and angle conventions
  multi_peak_validity.md        Mode-splitting boundary of the restricted surrogate
  gaia_scanning_law.md          Mission-grounded schedule choice and caveats
  sky_position_study.md         Sky-grid experiment design
  full_sky_results.md           Frozen 46-position full-sky result
  matched_n_control.md          Matched-transit-count control design
  native_vs_matched_n53_results.md  Frozen native versus matched-N comparison
  first_injection_experiment.md Original controlled uniform-angle pilot (superseded)
src/raa_orbit_model/
  kepler.py                     Newtonian binary dynamics, relative astrometry, SB2 RVs
  gaia.py                       1-D AL projection and equal-width blended surrogate
  penoyre.py                    Idealised orientation-dependent elongated-Gaussian response
  gaiamock_reference.py         Published gaiamock response, for benchmarking only
  astrometry.py                 Parallax factors, proper motion, position offset
  scanning.py                   GOST/gaiascanlaw schedule adapter and CSV archive
  likelihoods.py                Gaussian likelihood building blocks
  model.py                      Joint forward model and measurement-response dispatch
  synthetic.py                  Synthetic datasets on explicit Gaia schedules
  fit.py                        Bounded deterministic joint fitting and log-likelihood
  sampling.py                   Affine-invariant ensemble posterior sampler
  experiments.py                Baseline injection/recovery and bias-grid experiments
  response_fidelity.py          M2 injection fitted with M0/M1/M2 on paired realizations
  experiment_config.py          Orbit and noise axes shared by the runners
  robustness.py                 Unequal-component-width misspecification controls
  sky_study.py                  Ecliptic sky grid and scan-geometry diagnostics
  sky_experiments.py            Sky-position bias scan
  matched_control.py            Stratified matched-transit-count subsetting
  matched_sky_experiments.py    Matched-N sky experiment
  bias_analysis.py              Paired model comparison and transition diagnostics
  dr3_validation.py             DR3 NSS query, Thiele-Innes conversion, selection
scripts/
  run_bias_scan.py              Baseline bias scan on a Gaia schedule
  run_response_fidelity.py      M0/M1/M2 finite-elongation response experiment
  analyze_bias_scan.py          Paired summaries and publication-ready transition plots
  run_sky_position_scan.py      Full-sky scan over an ecliptic grid
  analyze_sky_position_scan.py  Sky maps, boundaries and correlations
  run_matched_n_control.py      Matched-transit-count control run
  compare_native_matched_n.py   Native versus matched-N comparison and figures
  validate_against_dr3.py       Validation step V6 against an exported DR3 table
```

Every experiment runner exposes the orbit (`--period-yr`, `--eccentricity`, `--mass-ratio`, …) and the per-channel noise levels (`--gaia-sigma-mas`, …). The defaults for the original runners reproduce the frozen baseline results exactly.

## Install and test

Core tests and the analysis-table tests:

```bash
python -m pip install -e ".[test]"
pytest -q
```

For the full mission-grounded workflow including scanning-law generation and plots:

```bash
python -m pip install -e ".[test,scanlaw,analysis]"
```

## Run the response-fidelity pilot

The first deterministic pilot deliberately injects the finite-elongation M2 response and fits the same realization with M0, M1, and M2:

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

This is a 36-fit regression pilot. `alpha`, `beta`, and `beta/alpha` are idealised surrogate parameters, not calibrated Gaia PLSF values. See `docs/response_fidelity_experiment.md` before expanding the grid.

## Run a mission-grounded baseline bias experiment

RA and Dec are deliberately explicit; the code does not hide an arbitrary representative sky position.

```bash
python scripts/run_bias_scan.py \
  --ra-deg 120.0 --dec-deg 30.0 --release dr4 \
  --write-schedule schedules/ra120_dec30_dr4.csv \
  --seeds 3 --output results/bias_ra120_dec30_dr4.csv
```

Reuse exactly the archived scan geometry with:

```bash
python scripts/run_bias_scan.py \
  --schedule-file schedules/ra120_dec30_dr4.csv \
  --seeds 3 --output results/bias_ra120_dec30_dr4_repeat.csv
```

The default is the uninterrupted **nominal** scanning law. `--apply-astrometry-gaps` optionally applies the gap information exposed by `gaiascanlaw`; it must not be interpreted as a complete source-level detection model for later releases.

## Dense transition scan

The resolution and light-fraction grids are configurable from the command line. For the currently interesting transition regime:

```bash
python scripts/run_bias_scan.py \
  --schedule-file schedules/ra120_dec30_dr4.csv \
  --a-over-sigma-values \
    0.20 0.25 0.30 0.35 0.40 0.45 0.50 \
    0.55 0.60 0.70 0.80 0.90 1.00 1.20 \
  --beta-values 0.05 0.20 0.40 \
  --seeds 20 \
  --sigma-response-mas 50 \
  --output results/bias_transition_ra120_dec30_dr4_seed20.csv
```

Each model is fit to the same synthetic realization for a given `(beta_G, a/sigma, seed)`. This permits a genuinely paired comparison rather than comparing independent noise realizations.

## Paired analysis and plots

Install the analysis extra if needed:

```bash
python -m pip install -e ".[analysis]"
```

Then analyze a completed baseline scan:

```bash
python scripts/analyze_bias_scan.py \
  results/bias_transition_ra120_dec30_dr4_seed20.csv \
  --output-dir results/transition_ra120_dec30_dr4_seed20 \
  --prefix ra120_dec30_dr4
```

The script writes paired tables and PNG/PDF transition/bias figures. The baseline comparison uses

```text
Delta chi2 = chi2_photocentre - chi2_resolution_aware
```

so positive values favour the resolution-aware measurement model for that injected realization. `Delta chi2` is descriptive; it is not converted into a Gaussian-significance claim.

## Core references

- Brandt et al. (2021), *AJ* **162**, 186, DOI: 10.3847/1538-3881/ac042e — orvara combined RV/absolute/relative astrometry inference.
- Leclerc et al. (2023), *A&A* **672**, A82, DOI: 10.1051/0004-6361/202244144 — BINARYS and its explicit limitation for luminous partially resolved Gaia systems.
- Halbwachs et al. (2023), *A&A* **674**, A9, DOI: 10.1051/0004-6361/202243969 — Gaia DR3 astrometric binary processing and filtering of partially resolved doubles.
- Holl et al. (2023), *A&A* **674**, A25, DOI: 10.1051/0004-6361/202245353 — scan-angle-dependent close-pair biases.
- Chevalier et al. (2023), *A&A* **678**, A19, DOI: 10.1051/0004-6361/202347111 — Gaia DR3 + SB2 component masses/luminosities.
- Liu et al. (2024), *A&A* **688**, L23, DOI: 10.1051/0004-6361/202450586 — partially resolved Gaia response in binary-asteroid orbit inference.
- El-Badry et al. (2024), *Open Journal of Astrophysics* **7**, DOI: 10.33232/001c.125461 — Gaia generative modelling and published blended AL response.
- Penoyre (2026), *RAS Techniques and Instruments* **5**, rzaf062, DOI: 10.1093/rasti/rzaf062; corrected by rzag016 — blended-source position and resolvability with elongated PSFs.
- Baycroft, Faria & Delisle (2026), arXiv:2606.24132 — Bayesian Gaia epoch astrometry and radial velocities with `kima`.
- Bailer-Jones & Kreidberg (2026), *A&A* **708**, A249, DOI: 10.1051/0004-6361/202659004 — component masses from unresolved Gaia astrometry and photometry.
- Rowell et al. (2026), *A&A* **708**, A174, DOI: 10.1051/0004-6361/202558618 — DR4 drift-scan PLSF modelling.
- Gaia DR4 expected-content documentation — official expected tables include epoch astrometry, epoch images, and double-lined RVS epoch products, while the content remains under development and subject to change.
