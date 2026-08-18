# 13 — Response-fidelity results

All 720 response-fidelity fits converge and satisfy the final scientific-validity checks.

At the central difficult point

\[
\beta_G=0.25,\qquad a_{\rm rel,ang}/\alpha=1,\qquad \beta_{\rm PSF}/\alpha=1.5,
\]

the median fit statistics are approximately:

| Model | median chi2 | M1 bias | M2 bias | parallax bias | light-fraction bias |
|---|---:|---:|---:|---:|---:|
| M0 photocentre | 1494.5 | -1.255% | -1.277% | +0.935% | -15.343% |
| M1 equal-width | 352.4 | -0.026% | -0.095% | +0.132% | -3.303% |
| M2 oriented | 207.2 | +0.070% | -0.014% | +0.083% | +0.053% |

The paired median mismatches are

\[
\mathrm{median}\,\Delta\chi^2_{0,2}\simeq1272.7,
\qquad
\mathrm{median}\,\Delta\chi^2_{1,2}\simeq129.3.
\]

Thus M1 is much closer to the M2 injection than the ordinary photocentre model but remains measurably misspecified at this elongation. The component masses are much less sensitive than the light fraction when the independent visual-SB2 orbit is strong.

As \(\beta_{\rm PSF}/\alpha\) increases from 1.5 to 3.0, M2 approaches the M1 limiting response: at the same light fraction and angular scale, \(\Delta\chi^2_{1,2}\) falls to about 5 and the M1 component-mass biases remain below about 0.05% in magnitude.

For the photocentre model, the primary-mass bias steepens with angular scale, from about -0.23% at \(a/\alpha=0.4\) to -1.26% at \(a/\alpha=1\) in the central light-fraction/elongation case.

The compact source tables are frozen under `results/frozen/response_fidelity_*`, and `figures/01_response_fidelity_central.svg` visualizes the central bias comparison.
