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
P>0,\quad 0\le e<1,\quad M_1>0,\quad M_2>0,\quad \varpi>0,\quad 0<\beta_G<1.
\]

For numerical inference, bounded/transformed parameters should be used rather than allowing the optimizer or sampler to enter an unphysical domain.

The mass ratio is

\[
q=M_2/M_1.
\]

The physical relative semi-major axis follows from Kepler's third law,

\[
a_{\rm rel}^3=\frac{G(M_1+M_2)P^2}{4\pi^2}.
\]

The barycentric component semi-major axes are

\[
a_1=a_{\rm rel}\frac{M_2}{M_1+M_2},\qquad
 a_2=a_{\rm rel}\frac{M_1}{M_1+M_2}.
\]

The RV semi-amplitudes are derived quantities,

\[
K_j=\frac{2\pi a_j\sin i}{P\sqrt{1-e^2}}.
\]

This parameterization enforces physical consistency between orbital size, masses, inclination, and SB2 velocity amplitudes.

## 3. Kepler solution

For epoch \(t\),

\[
M(t)=2\pi\frac{t-T_0}{P},
\]

and the eccentric anomaly satisfies

\[
M=E-e\sin E.
\]

The true anomaly is

\[
\nu=2\arctan2\!\left(\sqrt{1+e}\sin(E/2),\sqrt{1-e}\cos(E/2)\right).
\]

The instantaneous orbital radius is

\[
r=a_{\rm rel}(1-e\cos E).
\]

## 4. Relative astrometry

Define sky-plane coordinates \((\Delta\alpha^*,\Delta\delta)\), with \(\Delta\alpha^*=\Delta\alpha\cos\delta\). The orbital-plane position is rotated by \(\omega,i,\Omega\) into the sky plane. The model returns the angular relative separation using the parallax/distance conversion.

For a measurement vector

\[
\mathbf y_k=
\begin{pmatrix}
\Delta\alpha_k^*\\
\Delta\delta_k
\end{pmatrix},
\]

and full covariance \(C_k\),

\[
\ln L_{\rm rel,k}=-\frac12\left[
\mathbf r_k^T C_k^{-1}\mathbf r_k+\ln|C_k|+2\ln(2\pi)
\right],
\]

where \(\mathbf r_k=\mathbf y_k-\mathbf y_{k,\rm model}\).

Historical \((\rho,\theta)\) data may be transformed into tangent-plane coordinates, with covariance propagated when the source data permit it.

## 5. SB2 radial velocities

Define

\[
F(t)=\cos[\nu(t)+\omega]+e\cos\omega.
\]

The baseline component velocities are

\[
v_1(t)=\gamma + K_1F(t),
\]

\[
v_2(t)=\gamma - K_2F(t).
\]

For instrument \(j\), optional nuisance terms may be added:

\[
v_{1,j}=v_1+\Delta\gamma_{1,j},\qquad
v_{2,j}=v_2+\Delta\gamma_{2,j}.
\]

An optional excess-noise term \(s_{c,j}\) for component \(c\) gives

\[
\sigma_{\rm eff}^2=\sigma_{\rm quoted}^2+s_{c,j}^2.
\]

These terms must be data-supported and should not be included automatically if the dataset cannot constrain them.

## 6. Unresolved photocentre limit

Let the secondary mass fraction be

\[
B=\frac{M_2}{M_1+M_2},
\]

and the secondary Gaia-band light fraction be

\[
\beta_G=\frac{F_{2,G}}{F_{1,G}+F_{2,G}}.
\]

With relative vector \(\mathbf r=\mathbf r_2-\mathbf r_1\), the photocentre relative to the barycentre is

\[
\mathbf r_{\rm ph}=(\beta_G-B)\mathbf r.
\]

Any resolution-aware Gaia measurement model must recover this expression in the fully unresolved limit.

## 7. Gaia along-scan geometry

For scan angle \(\psi_k\), define the projected along-scan separation

\[
\Delta\eta_k=
\Delta\alpha_k^*\sin\psi_k+
\Delta\delta_k\cos\psi_k,
\]

using one explicit angle convention throughout the code and tests.

The projected component locations relative to the barycentre are

\[
\eta_{1,k}=-B\Delta\eta_k,\qquad
\eta_{2,k}=(1-B)\Delta\eta_k.
\]

In the unresolved limit the photocentre along scan is

\[
\eta_{{\rm ph},k}=\eta_{1,k}+\beta_G\Delta\eta_k=(\beta_G-B)\Delta\eta_k.
\]

## 8. Resolution-aware Gaia response: prototype stage

The first code stage implements the analytical/piecewise approximation used in Holl et al. (2023), based on the Lindegren (2022) close-double response treatment. It provides a deterministic along-scan bias as a function of projected separation and flux ratio and transitions between unresolved, marginally resolved, and resolved-primary regimes.

This is **not** treated as the final Gaia instrument model.

The final DR4 image-level model should instead predict calibrated sample values schematically as

\[
D_{kj}^{\rm model}=F_{1,k}L_k(x_j-x_{1,k})+F_{2,k}L_k(x_j-x_{2,k})+b_{kj},
\]

where \(L_k\) is the appropriate calibrated line-spread/point-spread response and \(b_{kj}\) is background. The exact likelihood must be derived from the released DR4 data model and calibration metadata rather than assumed in advance.

## 9. Absolute astrometric motion

The barycentric sky position may include

\[
\alpha^*(t)=\alpha_0^*+\mu_{\alpha^*}(t-t_{\rm ref})+\varpi P_\alpha(t)+\Delta\alpha_{\rm bin}^*(t),
\]

\[
\delta(t)=\delta_0+\mu_\delta(t-t_{\rm ref})+\varpi P_\delta(t)+\Delta\delta_{\rm bin}(t),
\]

where \(P_\alpha,P_\delta\) are parallax factors. The current prototype does not yet compute spacecraft ephemerides/parallax factors; that belongs in the Gaia-data integration stage.

## 10. Joint likelihood

The target posterior is

\[
p(\Theta\mid D)\propto
L_{\rm rel}
L_{\rm RV1}
L_{\rm RV2}
L_{\rm Gaia}
p(\Theta).
\]

Optional external measurements such as an independent parallax enter as explicit likelihood/prior terms rather than fixed metadata.

The inference backend should remain sampler-agnostic. A deterministic MAP/least-squares solution can initialize MCMC or nested sampling, but no sampler is assumed superior without benchmark tests.

## 11. Identifiability and validation strategy

The model must be validated in increasing complexity.

### V0 — orbital identities

- Kepler equation residuals;
- circular and face-on limiting cases;
- \(K_1/K_2=M_2/M_1\);
- Kepler's third law.

### V1 — relative astrometry + SB2 synthetic recovery

Generate synthetic data from known parameters and test parameter recovery under realistic noise and incomplete phase coverage.

### V2 — unresolved Gaia limit

Verify numerically that the Gaia response approaches

\[
(\beta_G-B)\Delta\eta
\]

as projected separation tends to zero.

### V3 — resolved-primary limit

At sufficiently large projected separation, verify that the model approaches the selected component location rather than an unresolved photocentre.

### V4 — scan-angle experiment

For the same physical binary at fixed orbital phase, vary \(\psi\) and verify the predicted transition between blended and more-resolved along-scan configurations.

### V5 — injection/recovery

Inject synthetic Gaia along-scan observations generated with the resolution-aware response. Fit them with both:

1. a simple photocentre model;
2. the resolution-aware model.

Measure bias in masses, parallax, orbital elements, and flux ratio as a function of angular separation, flux ratio, scan-angle coverage, and orbital phase coverage.

This comparison is central to demonstrating whether the additional measurement physics produces a scientifically significant improvement.

### V6 — DR3 catalogue consistency

Where suitable DR3 NSS systems exist, compare posterior predictions against published catalogue-level solutions while explicitly acknowledging that epoch data are unavailable.

### V7 — DR4 epoch/image validation

Only after DR4 epoch products and calibration documentation are released should the image/epoch model be validated against real Gaia measurement-level data.

## 12. Falsifiable scientific criterion

The proposed methodology is useful only if there exists a measurable regime in which the ordinary photocentre approximation yields statistically significant bias and the resolution-aware model removes that bias without unacceptable loss of identifiability.

The study should map this regime rather than assume it exists everywhere.

Primary axes of the experiment are

\[
\rho,\quad \Delta G\;(\text{or }\beta_G),\quad q,\quad P,\quad e,\quad \text{scan-angle coverage},\quad \text{S/N}.
\]

If injection tests show negligible improvement across realistic Gaia conditions, the proposed complexity is not justified and the methodology should be rejected or narrowed.
