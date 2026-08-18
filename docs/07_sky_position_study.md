# 07 — Sky-position study design

The supporting sky-position experiment asks how strongly a fixed response mismatch accumulates under different nominal Gaia scan schedules.

A 46-position ecliptic/sky grid is generated with mission-grounded nominal scanning-law schedules. For each sky position, the same physical binary and response configuration are used while transit number, scan-direction distribution and timing change with position.

The original equal-width study samples five light fractions, five angular-scale ratios and ten noise seeds and fits both the ordinary photocentre and equal-width resolution-aware model, producing 23,000 native-schedule fits.

The experiment is not intended to map a physical Gaia angular resolution. Its role is to establish that the *observability* of response mismatch is schedule dependent and therefore that a single sky location cannot define a universal detection threshold.

A separate matched-transit experiment (`09`–`10`) controls the number of retained transits so count effects can be separated from residual scan-geometry effects.
