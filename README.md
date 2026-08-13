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

The first prototype contains:

- a physically constrained Newtonian two-body SB2 forward model;
- resolved relative astrometry in the tangent plane;
- physical RV semi-amplitudes derived from component masses;
- the published piecewise Gaia along-scan response approximation used by Holl et al. (2023), based on Lindegren (2022), including unresolved, marginally resolved, and primary-resolved regimes;
- full 2-D Gaussian likelihood support for resolved astrometry;
- 1-D Gaussian likelihood support with optional jitter;
- unit tests for orbital, SB2, scan-projection, photocentre, and resolved-limit identities.

The 90-mas resolution response is a **prototype surrogate**, not a final Gaia image model. Gaia PSF/LSF calibration is time-, colour-, and instrument-state dependent; the final DR4 image-level model must follow the released DR4 data model and calibration products.

## Repository layout

```text
docs/
  literature_gap.md      Evidence for and limits of the candidate research gap
  methodology.md         Mathematical and validation methodology
src/raa_orbit_model/
  kepler.py              Newtonian binary dynamics, relative astrometry, SB2 RVs
  gaia.py                Along-scan projection and resolution-aware response
  likelihoods.py         Likelihood building blocks
  model.py               Joint forward model
tests/
  test_core.py           Physics and limiting-case regression tests
```

## Run tests

```bash
python -m pip install -e ".[test]"
pytest -q
```

## Core references

- Leclerc et al. (2023), *A&A* **672**, A82, DOI: 10.1051/0004-6361/202244144 — BINARYS.
- Halbwachs et al. (2023), *A&A* **674**, A9, arXiv:2206.05726 — Gaia DR3 astrometric binary processing and rejection of partially resolved doubles.
- Holl et al. (2023), *A&A* **674**, A25, DOI: 10.1051/0004-6361/202245353 — scan-angle-dependent signals and the analytical along-scan bias model.
- El-Badry et al. (2024), *Open Journal of Astrophysics* **7**, DOI: 10.33232/001c.125461 — Gaia astrometric-orbit selection-function generative model including marginal-resolution effects.
- Rowell et al. (2021), *A&A* **649**, A11, DOI: 10.1051/0004-6361/202039448 — Gaia EDR3 PSF/LSF modelling and calibration.
- ESA Gaia DR4 expected-content page — planned `epoch_astrometry`, `epoch_image`, and `rvs_epoch_data_double` products; content is explicitly subject to processing/validation changes.
