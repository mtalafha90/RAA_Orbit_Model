# 12 — Response-fidelity experiment

## Scientific question

How much physical bias is produced when data generated with a more faithful marginal-resolution response are fitted with a simpler measurement model?

## Response hierarchy

- **M0 — photocentre:** ordinary flux-weighted unresolved coordinate.
- **M1 — equal-width blend:** published `gaiamock`-family nonlinear one-dimensional peak response.
- **M2 — oriented finite-elongation surrogate:** exact effective-width reduction of an elongated Gaussian response motivated by Penoyre (2026).

M2 generates the central experiment and the same noisy realization is fitted independently by M0, M1 and M2.

## Grid

The response-fidelity grid samples

- \(\beta_G=0.05,0.25,0.45\);
- \(a_{\rm rel,ang}/\alpha=0.40,0.60,0.80,1.00\);
- \(\beta_{\rm PSF}/\alpha=1.5,3.0\);
- 10 paired noise seeds;
- 3 fitted response models.

This gives

\[
3\times4\times2\times10\times3=720
\]

fit records. The reference schedule contains 87 nominal directional transits, and all reported fits remain in the single-peak validity regime.

## Reference external data

The central synthetic binary uses 24 resolved relative-astrometry epochs with 0.20 mas one-sigma coordinate errors, 48 SB2 epochs with 0.10 km s\(^{-1}\) uncertainty per component, and 0.10 mas Gaia-like along-scan uncertainty.

## Primary outputs

The main diagnostics are fractional bias in \(M_1\), \(M_2\), parallax and \(\beta_G\). Paired \(\Delta\chi^2\) values are used only as descriptive model-mismatch measures, not converted to Gaussian detection significance because the response models are not treated as nested.

The key design principle is that response-model fidelity changes while the physical binary, external data and nominal scan schedule remain fixed.
