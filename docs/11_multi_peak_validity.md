# 11 — Single-peak / multi-peak validity

The M1/M2 surrogate likelihood assigns one measured coordinate to a blended profile. That representation is physically ambiguous once the profile has more than one local maximum of comparable relevance. Multi-peak states are therefore a **validity boundary**, not a branch to be arbitrarily continued.

For the equal-width two-Gaussian mixture, the exact critical separation is obtained from the simultaneous first- and second-derivative conditions. With light fraction \(\beta\), define \(u\) through

\[
\sinh(2u)-2u=\ln\left(\frac{1-\beta}{\beta}\right),
\]

then

\[
\frac{s_{\rm crit}}{\alpha}=2\cosh u.
\]

This criterion applies to the restricted equal-width surrogate only. It is not a calibrated Gaia resolution law.

For the orientation-dependent M2 surrogate, the same dimensionless mode criterion is evaluated after scaling the component separation by the effective width \(\gamma(\phi)\). Scientific fits are accepted only if all retained epochs are single-peaked in the final fitted response.

The reported 720 response-fidelity fits, 270 external-information fits and 100 matched-M2 controls retain all 87 Gaia-like epochs and have zero multi-peak injections/final predictions. If real Gaia epoch data are genuinely multi-peaked, the project will move to an image/sample-domain likelihood rather than force a single coordinate.
