"""Gaia DR3 catalogue-level consistency helpers.

The real-data validation ladder now distinguishes three stages. V6a has already
checked the Newtonian resolved-astrometry + SB2 core on a legacy GJ 765.2 data
set. This module supports V6b: comparison against published Gaia DR3 NSS
catalogue solutions and DR3 IPD diagnostics. V7 remains the future direct
measurement-level test using released Gaia epoch astrometry/images.

DR3 publishes astrometric orbits as **Thiele-Innes** constants
``(a_thiele_innes, b_thiele_innes, f_thiele_innes, g_thiele_innes)`` rather
than as Campbell elements. This module supplies the conversion in the same
North/East convention the rest of the project uses, together with a broad
archive query and loader. Target-specific GJ 765.2 machinery is in
``dr3_target.py`` and ``scripts/validate_gl765_dr3.py``.

Convention, matching `tests/test_orbit_conventions.py`:

    north = A X + F Y,      east = B X + G Y

with ``X = cos E - e``, ``Y = sqrt(1 - e^2) sin E``, node measured from North
through East, and the *relative* argument of periastron.

**Downloading is deliberately not automated here.** The Gaia archive is not
reachable from every environment, so `scripts/validate_against_dr3.py` reads a
CSV that has already been exported. The query to produce a broad candidate
sample is below.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

#: ADQL selecting DR3 astrometric orbits, favouring marginally resolved candidates.
#: Run at https://gea.esac.esa.int/archive/ and export as CSV.
NSS_ASTROMETRIC_ORBIT_QUERY = """
SELECT TOP 2000
    nss.source_id,
    nss.period,
    nss.eccentricity,
    nss.a_thiele_innes,
    nss.b_thiele_innes,
    nss.f_thiele_innes,
    nss.g_thiele_innes,
    nss.parallax,
    nss.t_periastron,
    nss.goodness_of_fit,
    src.ruwe,
    src.phot_g_mean_mag,
    src.ipd_frac_multi_peak,
    src.ipd_gof_harmonic_amplitude
FROM gaiadr3.nss_two_body_orbit AS nss
JOIN gaiadr3.gaia_source AS src USING (source_id)
WHERE nss.nss_solution_type = 'Orbital'
  AND nss.parallax > 5
  AND src.ipd_frac_multi_peak > 2
ORDER BY src.ipd_gof_harmonic_amplitude DESC
"""

REQUIRED_COLUMNS = (
    "source_id", "period", "eccentricity",
    "a_thiele_innes", "b_thiele_innes", "f_thiele_innes", "g_thiele_innes",
)


@dataclass(frozen=True)
class CampbellElements:
    """Geometric orbit elements recovered from Thiele-Innes constants."""

    semi_major_axis_mas: float
    inclination_deg: float
    omega_relative_deg: float
    node_deg: float


def thiele_innes_to_campbell(A: float, B: float, F: float, G: float) -> CampbellElements:
    """Convert Thiele-Innes constants to Campbell elements.

    Uses the standard decomposition. ``omega + Omega`` and ``omega - Omega``
    follow from the sums and differences of the constants, and the inclination
    from the ratio of the resulting semi-axes.

    The inclination is returned in ``[0, 180)`` degrees and the two angles in
    ``[0, 360)``.

    Thiele-Innes constants are invariant under simultaneously shifting both
    ``omega`` and ``Omega`` by 180 degrees, so astrometry alone cannot
    distinguish the ascending node from the descending one. This routine
    returns one branch consistently; recovering the correct branch requires
    radial velocities, which is precisely what the SB2 channel supplies. A
    round trip therefore reproduces the axis, the inclination and the node
    modulo 180 degrees exactly, and may return the angles shifted as a pair.
    """
    A, B, F, G = float(A), float(B), float(F), float(G)

    u = 0.5 * (A * A + B * B + F * F + G * G)
    v = A * G - B * F
    semi_major = math.sqrt(abs(u) + math.sqrt(max(u * u - v * v, 0.0)))
    if semi_major <= 0.0:
        raise ValueError("degenerate Thiele-Innes constants give a non-positive axis")

    omega_plus_node = math.atan2(B - F, A + G)
    omega_minus_node = math.atan2(-B - F, A - G)
    omega = 0.5 * (omega_plus_node + omega_minus_node)
    node = 0.5 * (omega_plus_node - omega_minus_node)

    # cos i from the ratio of the two half-axes of the apparent ellipse.
    denominator = semi_major * semi_major
    cos_inclination = v / denominator if denominator > 0 else 0.0
    cos_inclination = float(np.clip(cos_inclination, -1.0, 1.0))
    inclination = math.acos(cos_inclination)

    return CampbellElements(
        semi_major_axis_mas=semi_major,
        inclination_deg=math.degrees(inclination) % 180.0,
        omega_relative_deg=math.degrees(omega) % 360.0,
        node_deg=math.degrees(node) % 360.0,
    )


def campbell_to_thiele_innes(
    semi_major_axis_mas: float,
    inclination_deg: float,
    omega_relative_deg: float,
    node_deg: float,
) -> tuple[float, float, float, float]:
    """Inverse of :func:`thiele_innes_to_campbell`, in the project's convention."""
    a = float(semi_major_axis_mas)
    node = math.radians(node_deg)
    omega = math.radians(omega_relative_deg)
    inclination = math.radians(inclination_deg)
    cos_i = math.cos(inclination)
    return (
        a * (math.cos(node) * math.cos(omega) - math.sin(node) * math.sin(omega) * cos_i),
        a * (math.sin(node) * math.cos(omega) + math.cos(node) * math.sin(omega) * cos_i),
        a * (-math.cos(node) * math.sin(omega) - math.sin(node) * math.cos(omega) * cos_i),
        a * (-math.sin(node) * math.sin(omega) + math.cos(node) * math.cos(omega) * cos_i),
    )


def load_nss_csv(path: str | Path) -> list[dict]:
    """Load a DR3 NSS export produced by :data:`NSS_ASTROMETRIC_ORBIT_QUERY`."""
    path = Path(path)
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"NSS export is empty: {path}")
    missing = set(REQUIRED_COLUMNS) - set(rows[0])
    if missing:
        raise ValueError(f"NSS export is missing columns: {sorted(missing)}")

    parsed = []
    for row in rows:
        record = dict(row)
        for key, value in row.items():
            if key == "source_id":
                continue
            try:
                record[key] = float(value) if value not in ("", None) else math.nan
            except (TypeError, ValueError):
                record[key] = math.nan
        parsed.append(record)
    return parsed


def campbell_table(rows: list[dict]) -> list[dict]:
    """Attach Campbell elements to each usable DR3 astrometric orbit."""
    out = []
    for row in rows:
        constants = (
            row["a_thiele_innes"], row["b_thiele_innes"],
            row["f_thiele_innes"], row["g_thiele_innes"],
        )
        if any(not np.isfinite(c) for c in constants):
            continue
        try:
            elements = thiele_innes_to_campbell(*constants)
        except ValueError:
            continue
        out.append({
            "source_id": row["source_id"],
            "period_days": row.get("period", math.nan),
            "eccentricity": row.get("eccentricity", math.nan),
            "parallax_mas": row.get("parallax", math.nan),
            "photocentre_semi_major_mas": elements.semi_major_axis_mas,
            "inclination_deg": elements.inclination_deg,
            "omega_relative_deg": elements.omega_relative_deg,
            "node_deg": elements.node_deg,
            "ruwe": row.get("ruwe", math.nan),
            "ipd_frac_multi_peak": row.get("ipd_frac_multi_peak", math.nan),
            "ipd_gof_harmonic_amplitude": row.get("ipd_gof_harmonic_amplitude", math.nan),
        })
    return out


def marginally_resolved_candidates(
    table: list[dict],
    *,
    min_multi_peak_fraction: float = 2.0,
    min_harmonic_amplitude: float = 0.0,
) -> list[dict]:
    """Select rows whose published duplicity diagnostics suggest partial resolution.

    ``ipd_frac_multi_peak`` is the percentage of windows in which the image
    parameter determination found a second peak, and
    ``ipd_gof_harmonic_amplitude`` measures the scan-angle-dependent component
    of the fit quality. Both are the catalogue-level symptoms this project's
    surrogate predicts, which is what makes them the natural first target for
    a consistency test that does not require epoch data.
    """
    selected = []
    for row in table:
        multi = row.get("ipd_frac_multi_peak", math.nan)
        harmonic = row.get("ipd_gof_harmonic_amplitude", math.nan)
        if not np.isfinite(multi) or not np.isfinite(harmonic):
            continue
        if multi >= min_multi_peak_fraction and harmonic >= min_harmonic_amplitude:
            selected.append(row)
    return selected
