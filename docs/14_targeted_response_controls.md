# 14 — Targeted response controls

Two controls test whether the apparent success/failure of a response model is driven by measurement physics or by external-orbit information and finite ensemble size.

## External-information-strength control

At the fixed difficult response point

\[
\beta_G=0.25,\quad a_{\rm rel,ang}/\alpha=1,\quad \beta_{\rm PSF}/\alpha=1.5,
\]

the Gaia-like uncertainty remains 0.10 mas while only the independent visual/SB2 precision is weakened:

| Level | resolved astrometry | SB2 RV |
|---|---:|---:|
| strong | 0.20 mas | 0.10 km/s |
| medium | 1.00 mas | 0.50 km/s |
| weak | 2.00 mas | 1.00 km/s |

Thirty paired seeds and three fitted models produce 270 fits.

The photocentre component-mass bias grows from about 1.2% under strong external constraints to roughly 5–6% in the weak case. M1 mass biases grow from below 0.1% to about 1–1.4%. The matched M2 solution remains centered relative to its increasing random scatter.

An important result is that fit mismatch and physical bias do not move monotonically together: as external constraints weaken, the wrong model can distort the physical orbit and reduce residual mismatch while increasing mass/parallax bias.

## 100-seed matched-M2 control

A separate 100-realization M2→M2 control at the central point tests whether the small offsets in the ten-seed grid are coherent estimator biases.

The parallax fractional bias has median +0.0167%, 16th–84th percentiles -0.1274% to +0.1779%, mean +0.0214%, and standard deviation 0.1509%. The component masses and light fraction are likewise centered close to zero.

This supports the interpretation that the large M0 and residual M1 errors are measurement-model effects rather than a coherent offset intrinsic to the fitted M2 parameterization.

Frozen summaries are stored in `results/frozen/external_information_control_summary.csv` and `results/frozen/matched_m2_100seed_summary.csv`.
