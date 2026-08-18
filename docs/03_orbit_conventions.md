# 03 — Orbit conventions

Joint visual astrometry and SB2 spectroscopy is vulnerable to 180-degree convention errors. The project fixes one explicit convention and regression-tests it.

- `omega_deg` is the **primary spectroscopic** argument of periastron, \(\omega_1\).
- `node_deg` is the position angle of the primary ascending node measured from local North through East.
- positive RV means recession.
- relative astrometry is component 2 relative to component 1, \(\mathbf r_{\rm rel}=\mathbf r_2-\mathbf r_1\).
- tangent-plane order is \((\Delta\alpha^*,\Delta\delta)=(\mathrm{East},\mathrm{North})\).

The relative orbit therefore uses

\[
\omega_{\rm rel}=\omega_1+\pi.
\]

With \(u=\nu+\omega_{\rm rel}\),

\[
\Delta N=r(\cos\Omega\cos u-\sin\Omega\sin u\cos i),
\]

\[
\Delta E=r(\sin\Omega\cos u+\cos\Omega\sin u\cos i).
\]

The line-of-sight coordinate returned by the code is positive away from the observer and is tested to satisfy

\[
\frac{dz_{\rm rel}}{dt}=RV_2-RV_1.
\]

The SB2 model is

\[
F(t)=\cos[\nu(t)+\omega_1]+e\cos\omega_1,
\]

\[
RV_1=\gamma+K_1F(t),\qquad RV_2=\gamma-K_2F(t),
\]

which satisfies barycentric momentum balance.

`tests/test_orbit_conventions.py` checks the 180-degree relation, North-through-East node definition, Thiele–Innes agreement, the line-of-sight finite-difference identity, and the primary RV sign convention. Results produced before this convention audit are development provenance only; manuscript results use the audited convention.
