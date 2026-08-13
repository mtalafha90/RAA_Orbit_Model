# Gaia scanning-law schedule used by RAA Orbit Model

## Decision

Scientific injection/recovery runs no longer draw Gaia scan angles uniformly.
They require an explicit sky-position-dependent schedule. The default schedule
source is the **Gaia Nominal Scanning Law derived from the Gaia Observation
Forecast Tool (GOST)**, accessed through Zephyr Penoyre's `gaiascanlaw` package.
The exact schedule used by a run can be archived to CSV and reused without the
external package.

This is a mission-grounded *nominal* schedule. It is not a claim to reproduce
the exact set of CCD observations that survived acquisition, downlink,
calibration, source matching, and filtering for a real Gaia source.

## Why this source was selected

### 1. Official Gaia scanning-law definition

Gaia's routine scanning combines a six-hour spacecraft spin with a roughly
63-day precession of the spin axis around the solar direction at a fixed
45-degree solar-aspect angle. The two astrometric fields of view are separated
by the 106.5-degree basic angle. These facts are part of the Gaia mission and
archive documentation.

The official Gaia DR3 auxiliary table `gaiadr3.commanded_scan_law` provides the
commanded attitude over the 34-month DR3 interval. Its documentation defines
`scan_angle_fov1` and `scan_angle_fov2` as the position angle of the along-scan
direction, with 0 degrees toward local North and 90 degrees toward local East.
That is exactly the convention used by `raa_orbit_model.gaia.project_along_scan`.
The table is sampled every 10 seconds. The documentation also warns that the
actual attitude could deviate from commanded attitude by about 30 arcsec and
that the table does not encode mission data interruptions.

Official documentation:

- Gaia DR3 commanded scanning law: https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_auxiliary_tables/ssec_dm_commanded_scan_law.html
- Gaia scanning-law concepts: https://gea.esac.esa.int/archive/documentation/GEDR3/Introduction/chap_cu0int/cu0int_sec_mission/cu0int_ssec_scanning_law_concepts.html

### 2. GOST

The ESA/DPAC Gaia Observation Forecast Tool predicts field-of-view crossings for
a requested sky position from the Gaia scanning operations. ESA warns that GOST
is a forecast: operational activities and focal-plane gaps can remove predicted
observations, and the predicted time can differ from the actual observation.
The GOST result page recommends the commanded scan-law table when commanded
attitude is needed.

- GOST: https://gaia.esac.esa.int/gost/index.jsp

### 3. `gaiascanlaw`

ESA's Gaia community-tools page lists **Gaiascanlaw** as a Python interface to
Gaia scanning-law times and scan angles for the full nominal mission.
The package stores GOST-derived schedules on a level-6 HEALPix grid and returns
transit times and scan angles for an input RA/Dec. It also contains optional
published data-gap masks.

- ESA community tools: https://www.cosmos.esa.int/web/gaia/community-tools
- Code: https://github.com/zpenoyre/gaiascanlaw

Guerriero, Penoyre & Brown (2026, MNRAS 548, stag654; DOI
10.1093/mnras/stag654) use this same package for Gaia binary simulations. Their
method obtains observation epochs and scan angles as a function of sky position
and uses the along-scan projection

\[
x=\Delta\alpha^*\sin\psi+\Delta\delta\cos\psi+\epsilon.
\]

They use nominal baselines of 34 months for DR3, 66 months for DR4, and 126
months for the full mission simulation. They also note that nominal observation
counts can be overestimates because losses and rejected observations are not
fully modelled.

- Paper: https://doi.org/10.1093/mnras/stag654

### 4. Independent precedent in `gaiamock`

El-Badry et al.'s public `gaiamock` implementation also uses GOST-derived
position-dependent schedules. Its code reads `scanAngle[rad]`,
`parallaxFactorAlongScan`, and barycentric observation time from precomputed
GOST tables. This provides an independent methodological precedent for using
GOST-derived scan geometry in Gaia binary simulations.

- Code: https://github.com/kareemelbadry/gaiamock

## Data-release baselines

The adapter exposes `dr1` through `dr5` using the release endpoints encoded by
`gaiascanlaw`. The main RAA experiment defaults to `dr4`, corresponding to the
66-month DR4 baseline used in current mission simulations. ESA currently states
that DR4 is based on 66 months of data and is expected in December 2026.

Gaia stopped nominal science observations on 15 January 2025. Therefore no
schedule after that date is treated as science observing time in this project.
For the final full-mission experiment, the endpoint is the January 2025 end of
science operations encoded by the current `gaiascanlaw` package.

- ESA mission status / releases: https://www.esa.int/Science_Exploration/Space_Science/Gaia
- End of observations: https://www.cosmos.esa.int/web/gaia/end-of-observations

## Gap treatment

The default RAA schedule is the uninterrupted nominal law. This matches the
usage in Guerriero et al. (2026) and avoids pretending that the currently
available DR3 gap mask is a complete loss model for DR4 or the full mission.

For experiments specifically intended to mimic DR3 data availability,
`--apply-astrometry-gaps` passes `obstype='astrometry'` to `gaiascanlaw`. The
result is labelled in the schedule metadata. Such a schedule is still a model,
not a list of actual source-level detections.

## Reproducibility

Generate and archive one schedule:

```bash
python -m pip install -e ".[scanlaw]"
python scripts/run_bias_scan.py \
  --ra-deg 120.0 --dec-deg 30.0 --release dr4 \
  --write-schedule schedules/ra120_dec30_dr4.csv \
  --output results/bias_ra120_dec30_dr4.csv
```

Reuse exactly the same schedule later:

```bash
python scripts/run_bias_scan.py \
  --schedule-file schedules/ra120_dec30_dr4.csv \
  --output results/bias_ra120_dec30_dr4_repeat.csv
```

The CSV stores mission-relative transit time, directional scan angle, RA, Dec,
release label, provider description, and provider mission-start epoch.

## Required next experiment

One sky position is not sufficient because Gaia cadence and angle coverage vary
strongly over the sky. The next bias map should therefore be repeated at a
controlled set of ecliptic latitudes or on a sky grid. The resulting model bias
must be reported as a distribution over scanning-law geometries rather than as
a single universal threshold.

When real epoch astrometry becomes available, source-level epochs and scan
geometry from the released data should replace the nominal GOST schedule for
empirical target fits.
