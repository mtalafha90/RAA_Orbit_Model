# Methodology: Resolution-Aware Joint Orbit Inference

## 1. Goal

Build a target-level generative model for a luminous binary that predicts, from one physical parameter set:

- resolved relative astrometry;
- primary and secondary radial velocities;
- Gaia along-scan epoch measurements or, later, Gaia epoch image samples;
- derived masses, orbital parallax/distance, and Gaia-band flux ratio.

The two-body dynamics remain Newtonian/Keplerian unless data require additional physics. The methodological focus is the **observation model** in the transition between unresolved and resolved Gaia measurements.

## 2. Physical parameterization

A baseline SB2 parameter vector is

\[
\Theta_{\rm phys}=\{P,T_0,e,i,\Omega,\omega,M_1,M_2,\varpi,\gamma,\beta_G\},
\]

with constraints

\[
P>0,\quad 0\le e<1,\quad M_1>0,\quad M_2>0,\quad \varpi>0,\quad 0\le\beta_G\le1.
\]

For numerical inference, bounded/transformed parameters should be used rather than allowing the optimizer or sampler to enter an unphysical domain. The mass ratio is \(q=M_2/M_1\). The physical relative semi-major axis follows from Kepler's third law,

\[
a_{\rm rel}^3=\frac{G(M_1+M_2)P^2}{4\pi^2},
\]

with barycentric component semi-major axes

\[
a_1=a_{\rm rel}\frac{M_2}{M_1+M_2},\qquad a_2=a_{\rm rel}\frac{M_1}{M_1+M_2}.
\]

The RV semi-amplitudes are derived rather than independently fitted,

\[
K_j=\frac{2\pi a_j\sin i}{P\sqrt{1-e^2}},
\]

so the masses, orbit size, inclination and SB2 amplitudes remain physically consistent.

## 3. Kepler solution

For epoch \(t\),

\[
M(t)=2\pi\frac{t-T_0}{P},\qquad M=E-e\sin E,
\]

and

\[
\nu=2\arctan2\!\left(\sqrt{1+e}\sin(E/2),\sqrt{1-e}\cos(E/2)\right),
\]

with instantaneous radius \(r=a_{\rm rel}(1-e\cos E)\).

## 4. Relative astrometry

The model uses tangent-plane coordinates \((\Delta\alpha^*,\Delta\delta)=(\mathrm{East},\mathrm{North})\), with the primary spectroscopic periastron convention documented in `docs/03_orbit_conventions.md`. Historical \((\rho,\theta)\) measurements can be transformed to the tangent plane, with covariance propagated when the source data permit it. The implementation supports full covariance whitening.

## 5. SB2 radial velocities

Define

\[
F(t)=\cos[\nu(t)+\omega_1]+e\cos\omega_1,
\]

then

\[
v_1(t)=\gamma+K_1F(t),\qquad v_2(t)=\gamma-K_2F(t).
\]

Instrument-dependent offsets or excess-noise terms may be introduced only when a real data set supports them; they are not automatic degrees of freedom.

## 6. Unresolved photocentre limit

Let

\[
B=\frac{M_2}{M_1+M_2},\qquad \beta_G=\frac{F_{2,G}}{F_{1,G}+F_{2,G}}.
\]

For \(\mathbf r=\mathbf r_2-\mathbf r_1\),

\[
\mathbf r_{\rm ph}=(\beta_G-B)\mathbf r.
\]

Any resolution-aware Gaia observation model must recover this relation in the fully unresolved limit.

## 7. Gaia along-scan geometry

For directional scan angle \(\psi\),

\[
\Delta\eta=\Delta\alpha^*\sin\psi+\Delta\delta\cos\psi.
\]

The unresolved orbital along-scan displacement is therefore \((\beta_G-B)\Delta\eta\).

## 8. Response hierarchy

The scientific experiments deliberately separate measurement-response fidelity from binary dynamics:

- **M0:** ordinary unresolved photocentre;
- **M1:** equal-width nonlinear two-profile peak response, benchmarked against the published `gaiamock` response family;
- **M2:** finite-elongation, scan-orientation-dependent Gaussian surrogate based on the exact effective-width construction of Penoyre (2026).

The research widths used by M1/M2 are not calibrated Gaia PLSF resolution scales. The response hierarchy is intended to quantify inference bias caused by measurement-model misspecification, not to claim a new Gaia calibration model.

Once a surrogate profile is genuinely multi-peaked, a single-coordinate astrometric likelihood is no longer the natural representation. The long-term DR4 extension should use the released epoch-image/sample information directly.

## 9. Absolute astrometric motion

The barycentric sky position can include position offset, proper motion and parallax factors in addition to the binary perturbation. The code supports these terms for Gaia-like epoch experiments; real Gaia validation must use released source-level epochs/calibration rather than inferred scan coordinates.

## 10. Joint objective

The target posterior/likelihood factorization is

\[
p(\Theta\mid D)\propto L_{\rm rel}L_{\rm RV1}L_{\rm RV2}L_{\rm Gaia}p(\Theta).
\]

Optional external measurements enter as explicit constraints rather than silently fixed metadata. Deterministic weighted least squares is used for the current response-bias experiments; sampling infrastructure is retained for later posterior-coverage work.

## 11. Validation ladder

### V0--V5 — synthetic/internal validation

These stages cover orbital identities, relative-astrometry/SB2 recovery, unresolved and resolved response limits, scan-angle behaviour, response benchmarks, and injection/recovery under M0/M1/M2. The completed response-fidelity and control experiments are documented in `docs/12_response_fidelity_experiment.md` through `docs/14_targeted_response_controls.md`.

### V6a — real visual + SB2 physical-core validation: completed

The exact legacy GJ 765.2 / HIP 96656 input is versioned under `data/real/gj7652/`. The parser/fitter is in `src/raa_orbit_model/real_data.py`, the reproducibility runner is `scripts/fit_gl765_visual_sb2.py`, and regression tests are in `tests/test_real_data.py`.

The 110-constraint, 10-parameter free-parallax fit gives \(\chi^2=104.62895\), \(\chi^2_\nu=1.04629\), total mass \(1.58969\,M_\odot\), and orbital parallax \(35.4426\) mas. This is a real-data consistency check of the Newtonian resolved-astrometry+SB2 core; it is not a validation of the marginal-resolution Gaia response.

### V6b — Gaia DR3 catalogue/IPD bridge: query-ready

The target-specific HIP 96656 workflow uses the Gaia Hipparcos-2 best-neighbour cross-match, DR3 IPD diagnostics, and any available `nss_two_body_orbit` solution. The query/parser/conversion code is versioned, but no validated target-specific DR3 row is yet frozen in the repository. Therefore V6b remains **query-ready / data-pull pending** and no RUWE, IPD, source ID or NSS orbit is guessed.

### V7 — Gaia DR4 epoch/image validation: pending release

Official Gaia DR4 expected-content documentation lists `epoch_astrometry`, `epoch_image`, and `rvs_epoch_data_double` among the planned products. No firm public release date is asserted in this repository. Once the products are publicly available, V7 will compare the response hierarchy at the measurement level and use image/sample-domain inference where a unique astrometric coordinate is not defensible.

See `docs/15_real_data_validation.md` for the current status and strict claim boundaries.

## 12. Falsifiable scientific criterion

The response-aware methodology is scientifically useful only if there is a measurable regime in which the ordinary photocentre approximation creates meaningful bias and a more faithful response removes that bias without unacceptable loss of identifiability. The current experiments map that question as a function of light fraction, angular scale, response elongation and external-orbit information. Future real-data work must additionally explore period, eccentricity, inclination, magnitude/noise, sky position, scan geometry, cadence and calibrated Gaia response information.
