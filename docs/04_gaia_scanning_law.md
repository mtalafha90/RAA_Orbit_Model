# 04 — Gaia scanning law and schedule policy

## Purpose

The synthetic response experiments require mission-grounded observation times and scan directions. The code therefore separates **nominal scan scheduling** from **actual accepted source observations**.

## Current schedule source

The central response-fidelity and targeted-control experiments use an archived nominal schedule at ICRS RA = 120 deg, Dec = +30 deg generated through the public `gaiascanlaw` interface to GOST-derived nominal scanning-law tables. It contains 87 directional transits over the adopted DR4-like interval.

The supporting full-sky experiment uses the same scheduling machinery over a 46-position sky grid. A matched-transit control deterministically reduces each position to 53 mission-spanning transits so observation-count effects can be separated from scan-direction geometry.

## Angle convention

Directional scan angle \(\psi\) is retained over 0–360 deg, with

- 0 deg toward local North;
- 90 deg toward local East.

The along-scan projection is

\[
\Delta\eta=\Delta\alpha^*\sin\psi+\Delta\delta\cos\psi.
\]

Although an unoriented line has 180-degree symmetry, retaining the directional angle avoids silently mixing sign conventions in the orbital along-scan displacement.

## What nominal schedules are not

GOST/gaiascanlaw schedules are not measured source-level Gaia epoch astrometry. They predict nominal field-of-view crossings and do not encode all effects of acquisition, downlink, source matching, window assignment, gating, calibration or quality filtering. Therefore:

- nominal transit counts must not be presented as observed counts for a real source;
- a GOST scan angle must not be described as an actual accepted DR3/DR4 measurement;
- real response validation requires released epoch-level Gaia products.

## Use in this project

The response-fidelity hierarchy deliberately fixes one nominal schedule so M0/M1/M2 differences are not mixed with sky-position variation. The earlier full-sky and matched-N studies are retained only as supporting evidence that accumulated model mismatch depends strongly on observation count and scan geometry.

For GJ 765.2, the pre-DR4 transition calculation uses a uniform time/orientation envelope only as a feasibility diagnostic. It is explicitly not a Gaia transit statistic.
