# 09 — Matched-transit-count control

The native full-sky experiment mixes scan geometry with a strong variation in the number of nominal transits. To isolate those effects, every sky position is deterministically reduced to 53 mission-spanning transits.

The selection preserves mission coverage while equalizing observation count. The same binary configurations, noise prescription and equal-width response comparison are then repeated, producing 11,040 matched-N fits.

The matched-N experiment is a control, not a replacement for real source-level Gaia observations. Its purpose is to determine how much of the native sky dependence is explained by transit count and how much remains associated with scan-direction/timing geometry.

The comparison and frozen numerical summaries are documented in `docs/10_native_vs_matched_n53_results.md` and `results/frozen/`.
