# GJ 765.2 / HIP 96656 legacy visual + SB2 subset

`GL765_Test1.csv` is the exact legacy input supplied for the V6a real-binary consistency test. It contains 11 visual/speckle relative-astrometry measurements and 44 `Va` + 44 `Vb` CORAVEL radial-velocity measurements.

The file is preserved for provenance. Its header is **not** treated as authoritative catalogue metadata: in particular `RA=19.404`, `Dec=76.1812`, and `par=54.27` are retained because they are present in the source file, but the Gaia target workflow identifies the system independently as HIP 96656. The V6a fit uses only the measurement rows; the header parallax is tested separately as a fixed-parallax control.

Parser semantics are documented and tested in `src/raa_orbit_model/real_data.py`. The legacy numerical RV epoch is converted using the historical PySVOrbit convention

`year = 1900 + (epoch - 15020.31352) / 365.242198781`.

For the visual rows, position angle is North through East and the stored separation/uncertainty in arcsec are converted to tangent-plane East/North coordinates in mas. The quoted positional uncertainty is represented as an isotropic one-sigma uncertainty in East and North, matching the first-order convention used by the legacy workflow.

Reproduce the V6a fit with:

```bash
python scripts/fit_gl765_visual_sb2.py
```

The later comparison orbit is Balega et al. (2007), A&A 464, 635–640, DOI 10.1051/0004-6361:20066224. Because the CORAVEL measurements overlap the data used in that later combined solution, V6a is an implementation/consistency check rather than a statistically independent astrophysical measurement.
