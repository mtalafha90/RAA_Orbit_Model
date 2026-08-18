# Sky-position dependence experiment

## Purpose

This experiment tests whether the photocentre-versus-RAA transition measured at the pilot sky position is stable under the Gaia nominal scanning law or depends materially on sky position.

The experiment changes only the Gaia schedule. The physical binary, noise model, reduced `(beta_G, a/sigma)` grid, and random seeds are held fixed. Resolved relative astrometry and SB2 data use one fixed release-wide baseline at every sky position so the external-data realization is identical for a given seed.

## Coordinate grid

The controlled grid is defined in J2000 barycentric mean ecliptic coordinates and transformed to ICRS with Astropy's `BarycentricMeanEcliptic` to ICRS transformation.

Default ecliptic latitudes: `0, 15, 30, 45, 60, 75, 90 deg`.
Default ecliptic longitudes at each non-polar latitude: `0, 90, 180, 270 deg`.

Longitude is degenerate at the ecliptic pole, so the pole is represented once. The default grid contains 25 sky positions. Negative latitudes are supported through the command line but are not part of the first controlled run.

> **Superseded.** This note describes the original 25-position northern pilot. The published experiment is the **46-position full-sky** run documented in `docs/full_sky_results.md`, which extends the latitude grid from pole to pole. The command-line default remains the 25-position pilot so earlier runs stay reproducible, so the published result is **not** obtained from the defaults. Use `--full-sky`, or `FULL_SKY_ECLIPTIC_LATITUDES_DEG` in `sky_study.py`. The physical grid and record counts below likewise refer to the pilot, not to the published run.

## Physical grid

The first sky experiment uses:

- `beta_G = 0.05, 0.15, 0.25, 0.35, 0.45`
- `a/sigma = 0.40, 0.50, 0.60, 0.80, 1.00`
- 10 seeds by default
- photocentre and RAA fits for each realization

For 25 sky positions this gives 12,500 fit records.

The Gaussian response width remains a research-surrogate parameter. No numerical `a/sigma` value from this study is a physical Gaia resolution threshold.

## Directional scan diagnostics

Gaia scan angles remain directional 0--360 degree quantities and are never folded modulo 180 degrees. Each sky position records transit count, time span, circular resultant length, directional circular variance, largest directional gap, Shannon entropy in configurable equal-width directional bins, normalized entropy, effective bin count, occupied bin count, and a descriptive time-angle `R^2` from regression of centered transit time on `sin(psi)` and `cos(psi)`.

The entropy metrics are bin-dependent descriptive diagnostics, not inferential statistics.

## Reproducibility and restart behavior

Each exact Gaia schedule is archived before fitting. Results are also written per sky position. By default, a rerun skips a per-position file only if its row count matches the expected grid size. An incomplete file stops the run and must be replaced with `--overwrite`.

Outputs:

- `sky_grid.csv`
- `scan_geometry.csv`
- `schedules/<sky_id>_<release>.csv`
- `per_position/<sky_id>_<release>_bias.csv`
- `sky_position_bias.csv` after all positions are complete

Use `--schedule-only` to inspect the scan-law geometry before launching the expensive orbit fits.

## Interpretation

The target is to measure `Delta chi2` and parameter bias as functions of `(a/sigma, beta_G, sky position)` while separately recording scan-law geometry. The first analysis should test whether the transition boundary tracks ecliptic latitude, transit count, directional coverage, maximum gap, or time-angle coupling.

Nominal `gaiascanlaw`/GOST epochs are mission-model schedules and must not be described as actual observed Gaia epochs for a real target.
