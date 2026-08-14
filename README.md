# RAA Orbit Model

**Resolution-Aware Astrometric Orbit Model for luminous binaries**

This repository is the working research repository for a candidate methodology aimed at the difficult transition between unresolved and resolved binary-star astrometry in Gaia.

## Scientific target

The project investigates whether a target-level joint inference can combine:

1. resolved relative astrometry (interferometry, speckle, direct imaging),
2. both radial-velocity curves of a double-lined spectroscopic binary (SB2), and
3. Gaia along-scan epoch astrometry or epoch images using a **scan-angle- and resolution-aware measurement model** rather than assuming that Gaia always measures the photocentre.

Two things here are **not novel**, and the project does not claim them.

1. The broad combination “Gaia + RV + relative astrometry”. BINARYS already combines Hipparcos/Gaia absolute astrometry with relative astrometry and/or RV, and can handle SB2 velocities.
2. **The scan-angle- and resolution-aware measurement response itself.** El-Badry et al. (2024) released `gaiamock`, whose `al_bias_binary` already places the measured along-scan coordinate at the peak of the combined flux profile, following Lindegren (2022). This project's surrogate reproduces that published response to better than 0.01% — see `docs/gaiamock_benchmark.md`.

The candidate gap is narrower than originally recorded: in `gaiamock` the resolution-aware response is used **only to generate** data, and every fitting routine uses a plain photocentre model. The unaddressed question is the inference side — *fitting* that response at target level, jointly with resolved relative astrometry and both SB2 velocity curves, and propagating the measurement-model choice through to component masses.

**No novelty/priority claim is made at this stage.** The literature record, including leads that could not be verified, is maintained in `docs/literature_gap.md`.

## Current implementation status

The prototype now contains:

- a physically constrained Newtonian two-body SB2 forward model;
- resolved relative astrometry in the tangent plane;
- physical RV semi-amplitudes derived from component masses;
- Gaia along-scan projection geometry;
- an explicit two-profile blended-image surrogate for the marginal-resolution response;
- full 2-D covariance whitening for resolved astrometry;
- a bounded deterministic joint fitter for all 11 baseline physical parameters;
- a direct injection/recovery comparison between photocentre and resolution-aware hypotheses;
- a sky-position-dependent **GOST-derived Gaia Nominal Scanning Law** adapter through `gaiascanlaw`;
- portable CSV archiving/reloading of the exact transit times and directional scan angles used in an experiment;
- configurable resolution/light-fraction grids for dense transition experiments;
- paired seed-by-seed model comparison, physical-bias summaries, and PNG/PDF transition plots;
- regression tests for orbital identities, scan-angle convention, schedule conversion, exact synthetic recovery, model-misspecification detection, and paired-analysis sign conventions.

The blended-image response remains a **prototype surrogate**, not Gaia's calibrated PSF/LSF. Its width is an explicit research parameter. Gaia PSF/LSF calibration is time-, colour-, and instrument-state dependent (Rowell et al. 2021); the final image-level model must follow released Gaia epoch/image calibration products.

The original uniform-angle pilot is preserved in `docs/first_injection_experiment.md` as a code-validation result only. It is now superseded for scientific bias experiments by the mission-grounded schedule implementation documented in `docs/gaia_scanning_law.md`.

## Repository layout

```text
docs/
  literature_gap.md             Candidate gap, revised, with unverified leads marked
  gaiamock_benchmark.md         Scoring the surrogate against the published response
  methodology.md                Mathematical and validation methodology
  real_data_validation.md       Status of validation steps V6 and V7
  orbit_conventions.md          Orientation and angle conventions
  multi_peak_validity.md        Mode-splitting boundary of the surrogate
  gaia_scanning_law.md          Mission-grounded schedule choice and caveats
  sky_position_study.md         Sky-grid experiment design
  full_sky_results.md           Frozen 46-position full-sky result
  matched_n_control.md          Matched-transit-count control design
  native_vs_matched_n53_results.md  Frozen native versus matched-N comparison
  first_injection_experiment.md Original controlled uniform-angle pilot (superseded)
src/raa_orbit_model/
  kepler.py                     Newtonian binary dynamics, relative astrometry, SB2 RVs
  gaia.py                       Along-scan projection and resolution-aware surrogate
  gaiamock_reference.py         Published gaiamock response, for benchmarking only
  astrometry.py                 Parallax factors, proper motion, position offset
  scanning.py                   GOST/gaiascanlaw schedule adapter and CSV archive
  likelihoods.py                Gaussian likelihood building blocks
  model.py                      Joint forward model
  synthetic.py                  Synthetic datasets on explicit Gaia schedules
  fit.py                        Bounded deterministic joint fitting and log-likelihood
  sampling.py                   Affine-invariant ensemble posterior sampler
  experiments.py                Injection/recovery and bias-grid experiments
  experiment_config.py          Orbit and noise axes shared by the runners
  robustness.py                 Measurement-model misspecification controls
  sky_study.py                  Ecliptic sky grid and scan-geometry diagnostics
  sky_experiments.py            Sky-position bias scan
  matched_control.py            Stratified matched-transit-count subsetting
  matched_sky_experiments.py    Matched-N sky experiment
  bias_analysis.py              Paired model comparison and transition diagnostics
  dr3_validation.py             DR3 NSS query, Thiele-Innes conversion, selection
scripts/
  run_bias_scan.py              Bias scan on a specified Gaia sky position/schedule
  analyze_bias_scan.py          Paired summaries and publication-ready transition plots
  run_sky_position_scan.py      Full-sky scan over an ecliptic grid
  analyze_sky_position_scan.py  Sky maps, boundaries and correlations
  run_matched_n_control.py      Matched-transit-count control run
  compare_native_matched_n.py   Native versus matched-N comparison and figures
  validate_against_dr3.py       Validation step V6 against published DR3 orbits
```

Every experiment runner exposes the orbit (`--period-yr`, `--eccentricity`,
`--mass-ratio`, …) and the per-channel noise levels (`--gaia-sigma-mas`, …).
The defaults reproduce the frozen results exactly.

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

## Run a mission-grounded bias experiment

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

Then analyze a completed scan:

```bash
python scripts/analyze_bias_scan.py \
  results/bias_transition_ra120_dec30_dr4_seed20.csv \
  --output-dir results/transition_ra120_dec30_dr4_seed20 \
  --prefix ra120_dec30_dr4
```

The script writes:

- `<prefix>_paired.csv` — one row per matched `(beta_G, a/sigma, seed)` realization;
- `<prefix>_summary.csv` — median, mean, 16th/84th percentiles, extrema, and pair counts;
- `<prefix>_delta_chi2.{png,pdf}` — paired `chi2_photocentre - chi2_RAA` transition;
- fractional-bias plots for `beta_G`, parallax, `M1`, and `M2`;
- inclination-bias plots in degrees.

The paired analysis uses

```text
Delta chi2 = chi2_photocentre - chi2_resolution_aware
```

so positive values favour the resolution-aware measurement model for that injected realization. The plotting code reports seed distributions; it does not convert `Delta chi2` into a Gaussian-significance claim.

## Core references

- Leclerc et al. (2023), *A&A* **672**, A82, DOI: 10.1051/0004-6361/202244144 — BINARYS.
- Halbwachs et al. (2023), *A&A* **674**, A9, arXiv:2206.05726 — Gaia DR3 astrometric binary processing.
- Holl et al. (2023), *A&A* **674**, A25, DOI: 10.1051/0004-6361/202245353 — scan-angle-dependent close-pair biases.
- El-Badry et al. (2024), *Open Journal of Astrophysics* **7**, DOI: 10.33232/001c.125461 — Gaia astrometric-orbit selection-function generative modelling.
- Guerriero, Penoyre & Brown (2026), *MNRAS* **548**, stag654, DOI: 10.1093/mnras/stag654 — full-mission binary simulations using `gaiascanlaw`.
- Gaia DR3 `commanded_scan_law` documentation — official scan-angle convention and commanded-attitude caveats.
- ESA Gaia Community Tools — lists `gaiascanlaw` as a full nominal mission times/scan-angle interface.
