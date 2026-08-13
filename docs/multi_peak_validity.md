# Multi-peak validity boundary of the Gaussian response surrogate

## Scope

The current RAA image-response prototype represents a Gaia along-scan profile as the sum of two equal-width one-dimensional Gaussian components. This is a research surrogate only; its width `sigma_response_mas` is not a calibrated Gaia LSF width.

The boundary experiment showed that defining the measurement as an unconstrained global peak becomes ill-defined when the profile develops two modes. The code therefore now treats single- and multi-peak epochs as different regimes.

## Exact mode-splitting criterion

For component flux fractions

\[
F_1 = 1-\beta_G,\qquad F_2=\beta_G,\qquad 0\le\beta_G\le 0.5,
\]

and common Gaussian width \(\sigma\), let \(s\) be the absolute projected along-scan component separation.

At the transition from one to two modes, a stationary point satisfies both the first- and second-derivative equations. If its normalized distances from the two component centres are \(a\) and \(b\), the simultaneous conditions imply

\[
ab=1.
\]

Writing

\[
a=e^u,\qquad b=e^{-u},
\]

gives

\[
\sinh(2u)-2u = \ln\!\left(\frac{1-\beta_G}{\beta_G}\right),
\]

and therefore

\[
\frac{s_{\rm crit}}{\sigma}=2\cosh u.
\]

The implementation solves this scalar equation numerically. The two symmetry limits are exact:

- \(\beta_G=0\): \(s_{\rm crit}/\sigma\to\infty\); a dark secondary cannot create a second light-profile mode.
- \(\beta_G=0.5\): \(s_{\rm crit}/\sigma=2\); the equal-light profile splits into two modes only for \(s>2\sigma\).

This criterion is a mathematical property of the equal-width Gaussian mixture. It is **not** a Gaia angular-resolution threshold.

## Scientific handling of multi-peak epochs

Before an injection/recovery experiment, every epoch in the supplied Gaia scanning-law schedule is classified using the injected physical orbit, the scan angle, the light fraction, and the surrogate width.

Epochs with a single peak remain in the one-dimensional AL likelihood. Epochs with multiple peaks are removed from that likelihood rather than being assigned an arbitrary left/right/global peak. Both the photocentre and resolution-aware fits use the same retained epoch set.

The selection is only a validity operation for the current surrogate. It is not a model of Gaia detection probability, window assignment, source matching, IPD multi-peak processing, or completeness.

## Output diagnostics

Bias-scan CSVs now record:

- `gaia_n_transits`: number of transits in the original mission schedule;
- `gaia_n_single_peak_used`: transits retained in the current one-coordinate likelihood;
- `gaia_n_multi_peak_flagged`: transits outside that likelihood domain;
- `gaia_multi_peak_fraction`: flagged fraction;
- `gaia_critical_separation_sigma`: exact mixture mode-splitting threshold for the injected \(\beta_G\);
- `gaia_max_projected_separation_sigma`: largest projected AL component separation in units of the surrogate width;
- `gaia_final_multi_peak_predicted`: number of retained epochs that the final RAA fit would classify as multi-peak;
- `scientific_valid`: false for a resolution-aware final fit that leaves the single-peak response domain.

The analysis CLI refuses to summarize a new file containing `scientific_valid=False` rows.

## Optimizer continuation

The deterministic least-squares prototype is allowed to use the old continuous peak branch internally while traversing parameter space. This is a numerical continuation device only. A final RAA fit is accepted for scientific analysis only if the validity-aware response predicts one peak at every retained Gaia epoch.

This separation prevents the optimizer from being confused by a discontinuous hard boundary while ensuring that an invalid final solution cannot enter the scientific summary.

## Validation tests

Regression tests cover:

1. the unresolved photocentre limit;
2. the exact \(\beta_G=0\) dark-secondary limit;
3. the exact \(\beta_G=0.5\), \(s_{\rm crit}=2\sigma\) limit;
4. explicit rejection of an equal-light profile above \(2\sigma\);
5. mission-epoch filtering using aligned and perpendicular scan directions.
