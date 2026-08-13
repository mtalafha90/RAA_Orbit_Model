# Orbit conventions for joint relative astrometry and SB2 radial velocities

## Why this document exists

Joint visual/interferometric astrometry and spectroscopy is vulnerable to a 180-degree convention error in the argument of periastron. The package therefore fixes one explicit convention and tests it numerically.

The convention follows the primary-star Campbell parameterization used by BINARYS (Leclerc et al. 2023):

- `omega_deg` is the primary argument of periastron, `omega_1`;
- `node_deg` is the position angle of the primary ascending node, measured from local North toward East;
- positive spectroscopic radial velocity means recession;
- relative astrometry is component 2 relative to component 1, `r_rel = r2-r1`.

BINARYS states that the secondary periastron argument is related to the primary by

\[
\omega_2=\omega_1+\pi.
\]

Fabry et al. (2021), in a combined interferometric+SB2 analysis of 9 Sgr, likewise note explicitly that the periastron argument of the relative orbit is shifted by 180 degrees from the primary spectroscopic periastron argument.

## Relative orbit used by the code

For the relative vector

\[
\mathbf r_{\rm rel}=\mathbf r_2-\mathbf r_1,
\]

we therefore use

\[
\omega_{\rm rel}=\omega_1+\pi.
\]

The public `BinaryParams.omega_deg` field remains `omega_1`; the derived property `omega_relative_deg` returns the corresponding relative-orbit value modulo 360 degrees.

The tangent-plane coordinate order is

\[
(\Delta\alpha^*,\Delta\delta)=(\mathrm{East},\mathrm{North}).
\]

With `Omega = node_deg`, measured from North through East, the relative orbit is evaluated with `u = nu + omega_rel` as

\[
\Delta N=r[\cos\Omega\cos u-\sin\Omega\sin u\cos i],
\]

\[
\Delta E=r[\sin\Omega\cos u+\cos\Omega\sin u\cos i].
\]

These are equivalent to the standard Thiele-Innes expressions used by Fabry et al. (2021), whose Appendix B writes `Delta N = A X + F Y` and `Delta E = B X + G Y`.

The third coordinate returned by `relative_position_au` is line of sight positive away from the observer. It is tested to satisfy

\[
\frac{d z_{\rm rel}}{dt}=RV_2-RV_1.
\]

## SB2 radial velocities

The primary spectroscopic equation is

\[
F_1(t)=\cos[\nu(t)+\omega_1]+e\cos\omega_1,
\]

\[
RV_1=\gamma+K_1F_1,
\qquad
RV_2=\gamma-K_2F_1.
\]

The second equation is equivalent to using `omega_2 = omega_1 + pi` for the secondary.

This convention also guarantees barycentric momentum balance,

\[
M_1(RV_1-\gamma)+M_2(RV_2-\gamma)=0.
\]

## Regression tests

`tests/test_orbit_conventions.py` locks the convention by checking:

1. `omega_rel = omega_1 + 180 deg`;
2. the North-through-East definition of `Omega`;
3. direct agreement with the published Thiele-Innes North/East equations;
4. the finite-difference identity `d z_rel/dt = RV2-RV1`;
5. the sign of the primary RV curve at primary periastron for `omega_1=0`.

These tests are intentionally redundant. An internally self-consistent synthetic injection can otherwise hide a 180-degree or East/North error.

## Consequence for previous synthetic results

Runs produced before this convention audit used the same numerical `omega_deg` in the relative sky orbit and primary RV curve and also labeled the standard North/East projection terms in reverse order. Those runs were internally self-consistent for injection/recovery, but their orientation relative to a fixed Gaia scanning-law schedule was not tied to the standard astronomical convention.

Therefore, pre-audit transition numbers should be treated as development results. After this change, the mission-grounded transition and boundary experiments must be rerun before their numerical values are used in a manuscript.

## References

- Leclerc, A., et al. 2023, *Astronomy & Astrophysics*, 672, A82, "Combining Hipparcos and Gaia data for the study of binaries: the BINARYS tool", DOI: 10.1051/0004-6361/202244144, arXiv:2209.04210.
- Fabry, M., et al. 2021, *Astronomy & Astrophysics*, 651, A119, "Resolving the dynamical mass tension of the massive binary 9 Sagittarii", arXiv:2105.09968. See Table 1 note and Appendix B.
- Householder, A. & Weiss, L. 2022, arXiv:2212.06966, "The Inconsistent use of omega in the RV Equation". This paper documents the broader 180-degree omega ambiguity when RV and astrometric conventions are mixed.
