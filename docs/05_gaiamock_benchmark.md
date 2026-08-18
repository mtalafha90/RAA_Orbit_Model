# 05 — `gaiamock` equal-width response benchmark

The M1 equal-width response is not introduced as new measurement physics. It is benchmarked against the published `gaiamock` blended-source response family over their common single-peak domain.

Both models represent the measured one-dimensional coordinate through the maximum of a blended two-component profile, with the response depending nonlinearly on projected component separation and light ratio. The repository regression tests compare the internal implementation with the reference behaviour and verify the expected unresolved and brighter-component limits.

The benchmark establishes implementation consistency only. It does **not** imply that the equal-width Gaussian width is a calibrated Gaia PLSF scale, and it does not validate the response against real Gaia epoch images.

The scientific role of M1 is therefore:

1. a published resolution-aware baseline between the ordinary photocentre M0 and the more general orientation-dependent M2 surrogate;
2. a way to quantify how much mass/parallax bias remains when the response is closer to the injection but still misspecified;
3. a limiting model recovered by the M2 implementation as the across-scan width tends to infinity.

Relevant regression tests are in `tests/test_gaiamock_benchmark.py` and `tests/test_penoyre_response.py`.
