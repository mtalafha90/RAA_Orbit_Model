# RAA Orbit Model

**Resolution-Aware Astrometric Orbit Model for luminous binaries**

This repository is the working research repository for a candidate methodology aimed at the difficult transition between unresolved and resolved binary-star astrometry in Gaia.

## Scientific target

The project investigates whether a target-level joint inference can combine:

1. resolved relative astrometry (interferometry, speckle, direct imaging),
2. both radial-velocity curves of a double-lined spectroscopic binary (SB2), and
3. Gaia along-scan epoch astrometry or epoch images using a **scan-angle- and resolution-aware measurement model** rather than assuming that Gaia always measures the photocentre.

The broad combination “Gaia + RV + relative astrometry” is **not novel**. BINARYS already combines Hipparcos/Gaia absolute astrometry with relative astrometry and/or RV, and can handle SB2 velocities. The specific candidate gap recorded here is narrower: target-level inference in the luminous, marginally/partially resolved regime where the Gaia measured position depends on projected separation and scan angle and is not generally equal to a simple photocentre.

**No novelty/priority claim is made at this stage.** The literature record is maintained in `docs/literature_gap.md`.

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
- regression tests for orbital identities, scan-angle convention, schedule conversion, exact synthetic recovery, and model-misspecification detection.

The blended-image response remains a **prototype surrogate**, not Gaia's calibrated PSF/LSF. Its width is an explicit research parameter. Gaia PSF/LSF calibration is time-, colour-, and instrument-state dependent (Rowell et al. 2021); the final image-level model must follow released Gaia epoch/image calibration products.

The original uniform-angle pilot is preserved in `docs/first_injection_experiment.md` as a code-validation result only. It is now superseded for scientific bias experiments by the mission-grounded schedule implementation documented in `docs/gaia_scanning_law.md`.

## Repository layout

```text
docs/
  literature_gap.md             Evidence for and limits of the candidate research gap
  methodology.md                Mathematical and validation methodology
  first_injection_experiment.md Original controlled uniform-angle pilot
  gaia_scanning_law.md          Mission-grounded schedule choice and caveats
src/raa_orbit_model/
  kepler.py                     Newtonian binary dynamics, relative astrometry, SB2 RVs
  gaia.py                       Along-scan projection and resolution-aware surrogate
  scanning.py                   GOST/gaiascanlaw schedule adapter and CSV archive
  likelihoods.py                Likelihood building blocks
  model.py                      Joint forward model
  synthetic.py                  Synthetic datasets on explicit Gaia schedules
  fit.py                        Bounded deterministic joint fitting
  experiments.py                Injection/recovery and bias-grid experiments
scripts/
  run_bias_scan.py              Bias scan on a specified Gaia sky position/schedule
```

## Install and test

```bash
python -m pip install -e ".[test]"
pytest -q
```

To generate GOST-derived schedules locally, install the optional scanning-law dependency:

```bash
python -m pip install -e ".[scanlaw]"
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

## Core references

- Leclerc et al. (2023), *A&A* **672**, A82, DOI: 10.1051/0004-6361/202244144 — BINARYS.
- Halbwachs et al. (2023), *A&A* **674**, A9, arXiv:2206.05726 — Gaia DR3 astrometric binary processing.
- Holl et al. (2023), *A&A* **674**, A25, DOI: 10.1051/0004-6361/202245353 — scan-angle-dependent close-pair biases.
- El-Badry et al. (2024), *Open Journal of Astrophysics* **7**, DOI: 10.33232/001c.125461 — Gaia astrometric-orbit selection-function generative modelling.
- Guerriero, Penoyre & Brown (2026), *MNRAS* **548**, stag654, DOI: 10.1093/mnras/stag654 — full-mission binary simulations using `gaiascanlaw`.
- Gaia DR3 `commanded_scan_law` documentation — official scan-angle convention and commanded-attitude caveats.
- ESA Gaia Community Tools — lists `gaiascanlaw` as a full nominal mission times/scan-angle interface.
