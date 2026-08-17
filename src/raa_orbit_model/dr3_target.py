"""Target-specific Gaia DR3 catalogue validation for GJ 765.2 / HIP 96656.

This module deliberately separates what DR3 can test from what it cannot.
DR3 catalogue products can test the published source-level astrometry, IPD
scan-angle diagnostics, and any NSS orbit. They do not provide the general
stellar epoch along-scan astrometry needed to fit the RAA response hierarchy
M0/M1/M2 directly.

The primary identification route uses Gaia DR3's own Hipparcos-2 best-neighbour
crossmatch for HIP 96656. This is preferable to a nearest-neighbour cone search
for a high-proper-motion binary because Gaia's crossmatch propagates available
astrometric parameters between catalogue epochs. The SIMBAD sky position is
retained only as an independent sanity check / fallback.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from .dr3_validation import thiele_innes_to_campbell


GJ765_HIP = 96656
GJ765_RA_DEG = 294.7765558333333
GJ765_DEC_DEG = 76.42202333333333
GJ765_SEARCH_RADIUS_ARCSEC = 10.0

# Published combined-orbit quantities used only to define an external benchmark.
# Balega et al. (2007), A&A 464, 635--640.
GJ765_PERIOD_YR = 11.919
GJ765_ECCENTRICITY = 0.240
GJ765_INCLINATION_DEG = 80.2
GJ765_NODE_DEG = 293.0
GJ765_OMEGA_REL_DEG = 250.0
GJ765_A_REL_MAS = 189.0
GJ765_M1_MSUN = 0.831
GJ765_M2_MSUN = 0.763
GJ765_DELTA_V_MAG = 0.65
GJ765_ORBITAL_PARALLAX_MAS = 31.0


_SOURCE_AND_NSS_COLUMNS = """
    src.source_id,
    src.designation,
    src.ra,
    src.dec,
    src.parallax,
    src.parallax_error,
    src.pmra,
    src.pmra_error,
    src.pmdec,
    src.pmdec_error,
    src.phot_g_mean_mag,
    src.astrometric_n_obs_al,
    src.astrometric_n_good_obs_al,
    src.astrometric_gof_al,
    src.astrometric_chi2_al,
    src.astrometric_excess_noise,
    src.visibility_periods_used,
    src.ruwe,
    src.ipd_frac_multi_peak,
    src.ipd_gof_harmonic_amplitude,
    src.ipd_gof_harmonic_phase,
    src.scan_direction_strength_k1,
    src.scan_direction_strength_k2,
    src.scan_direction_strength_k3,
    src.scan_direction_strength_k4,
    src.scan_direction_mean_k1,
    src.scan_direction_mean_k2,
    src.scan_direction_mean_k3,
    src.scan_direction_mean_k4,
    src.duplicated_source,
    src.non_single_star,
    src.has_epoch_rv,
    nss.nss_solution_type,
    nss.period AS nss_period_days,
    nss.eccentricity AS nss_eccentricity,
    nss.t_periastron AS nss_t_periastron,
    nss.parallax AS nss_parallax_mas,
    nss.goodness_of_fit AS nss_goodness_of_fit,
    nss.a_thiele_innes,
    nss.b_thiele_innes,
    nss.f_thiele_innes,
    nss.g_thiele_innes
""".strip()


def gj765_target_query() -> str:
    """Return the preferred DR3 ADQL, identified by Gaia's HIP2 crossmatch."""
    return f"""
SELECT
    hip.angular_distance AS separation_arcsec,
    hip.number_of_neighbours AS hipparcos2_number_of_neighbours,
    hip.gaia_astrometric_params AS hipparcos2_gaia_astrometric_params,
    {_SOURCE_AND_NSS_COLUMNS}
FROM gaiadr3.hipparcos2_best_neighbour AS hip
JOIN gaiadr3.gaia_source AS src
    ON hip.source_id = src.source_id
LEFT OUTER JOIN gaiadr3.nss_two_body_orbit AS nss
    ON src.source_id = nss.source_id
WHERE hip.original_ext_source_id = {GJ765_HIP}
ORDER BY nss.nss_solution_type ASC
""".strip()


def gj765_cone_query(radius_arcsec: float = GJ765_SEARCH_RADIUS_ARCSEC) -> str:
    """Fallback coordinate query around the SIMBAD position.

    This should not be used to identify the target by simple nearest distance
    alone unless the catalogue-epoch/proper-motion issue has been checked.
    """
    radius_deg = float(radius_arcsec) / 3600.0
    if radius_deg <= 0:
        raise ValueError("radius_arcsec must be > 0")
    return f"""
SELECT
    DISTANCE(
        POINT('ICRS', src.ra, src.dec),
        POINT('ICRS', {GJ765_RA_DEG:.12f}, {GJ765_DEC_DEG:.12f})
    ) * 3600.0 AS separation_arcsec,
    0 AS hipparcos2_number_of_neighbours,
    0 AS hipparcos2_gaia_astrometric_params,
    {_SOURCE_AND_NSS_COLUMNS}
FROM gaiadr3.gaia_source AS src
LEFT OUTER JOIN gaiadr3.nss_two_body_orbit AS nss
    ON src.source_id = nss.source_id
WHERE 1 = CONTAINS(
    POINT('ICRS', src.ra, src.dec),
    CIRCLE('ICRS', {GJ765_RA_DEG:.12f}, {GJ765_DEC_DEG:.12f}, {radius_deg:.12f})
)
ORDER BY separation_arcsec ASC
""".strip()


def flux_fraction_from_delta_mag(delta_mag: float) -> float:
    """Secondary light fraction for m2-m1 = ``delta_mag``."""
    return 1.0 / (1.0 + 10.0 ** (0.4 * float(delta_mag)))


def secondary_mass_fraction(m1_msun: float, m2_msun: float) -> float:
    total = float(m1_msun) + float(m2_msun)
    if total <= 0:
        raise ValueError("total mass must be positive")
    return float(m2_msun) / total


def photocentre_axis_mas(
    relative_axis_mas: float,
    m1_msun: float,
    m2_msun: float,
    secondary_light_fraction: float,
) -> float:
    """Ordinary unresolved-photocentre semi-major axis |beta-B| a_rel."""
    beta = float(secondary_light_fraction)
    if not 0.0 <= beta <= 1.0:
        raise ValueError("secondary_light_fraction must be in [0, 1]")
    B = secondary_mass_fraction(m1_msun, m2_msun)
    return abs(beta - B) * float(relative_axis_mas)


def gj765_photocentre_benchmark() -> dict[str, float]:
    """Return the M0 benchmark using Delta V only as a pre-G-band proxy.

    ``beta_V`` must not be described as a Gaia G-band light fraction. The value
    is retained because it makes the catalogue-level expectation reproducible
    while a component-resolved G-band flux ratio is unavailable.
    """
    beta_v = flux_fraction_from_delta_mag(GJ765_DELTA_V_MAG)
    B = secondary_mass_fraction(GJ765_M1_MSUN, GJ765_M2_MSUN)
    return {
        "mass_fraction_secondary": B,
        "beta_V_proxy": beta_v,
        "predicted_M0_photocentre_axis_mas": photocentre_axis_mas(
            GJ765_A_REL_MAS,
            GJ765_M1_MSUN,
            GJ765_M2_MSUN,
            beta_v,
        ),
    }


def _as_float(value) -> float:
    if value in (None, "", "null", "NULL", "NaN", "nan"):
        return math.nan
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def _as_int(value) -> int | None:
    if value in (None, "", "null", "NULL"):
        return None
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def _is_astrometric_nss_type(value: str) -> bool:
    value = str(value).strip()
    return (
        value == "Orbital"
        or value.startswith("OrbitalAlternative")
        or value.startswith("OrbitalTargetedSearch")
        or value == "AstroSpectroSB1"
    )


@dataclass(frozen=True)
class DR3TargetSummary:
    source_id: str
    separation_arcsec: float
    hipparcos2_number_of_neighbours: int | None
    hipparcos2_gaia_astrometric_params: int | None
    parallax_mas: float
    parallax_error_mas: float
    ruwe: float
    ipd_frac_multi_peak_percent: float
    ipd_gof_harmonic_amplitude: float
    ipd_gof_harmonic_phase_deg: float
    non_single_star_flag: int | None
    nss_solution_types: tuple[str, ...]
    nss_orbit_count: int
    nearest_nss_photocentre_axis_mas: float
    nearest_nss_period_days: float
    nearest_nss_eccentricity: float
    nearest_nss_inclination_deg: float
    nearest_nss_omega_relative_deg: float
    nearest_nss_node_deg: float
    predicted_M0_photocentre_axis_mas: float
    beta_V_proxy: float
    mass_fraction_secondary: float


def load_target_export(path: str | Path) -> list[dict[str, str]]:
    path = Path(path)
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"DR3 target export is empty: {path}")
    required = {"source_id", "separation_arcsec", "parallax", "ruwe"}
    missing = required - set(rows[0])
    if missing:
        raise ValueError(f"DR3 target export is missing columns: {sorted(missing)}")
    return rows


def summarize_target_rows(rows: list[dict[str, str]]) -> DR3TargetSummary:
    """Summarise the source returned by the HIP2 crossmatch / fallback export."""
    if not rows:
        raise ValueError("rows must not be empty")

    source_ids = {str(r.get("source_id", "")).strip() for r in rows if str(r.get("source_id", "")).strip()}
    if not source_ids:
        raise ValueError("target export contains no source_id")

    # Preferred HIP2 query should return one Gaia source, possibly repeated for
    # several NSS rows. For a fallback cone export, preserve the older nearest
    # selection while requiring the user to verify the identification.
    if len(source_ids) == 1:
        source_id = next(iter(source_ids))
        target_rows = rows
        nearest_row = rows[0]
    else:
        finite_sep = [(_as_float(r.get("separation_arcsec")), r) for r in rows]
        finite_sep = [(s, r) for s, r in finite_sep if np.isfinite(s)]
        if not finite_sep:
            raise ValueError("multiple sources and no finite separation_arcsec values")
        _, nearest_row = min(finite_sep, key=lambda item: item[0])
        source_id = str(nearest_row.get("source_id", "")).strip()
        target_rows = [r for r in rows if str(r.get("source_id", "")).strip() == source_id]

    nss_rows = [r for r in target_rows if str(r.get("nss_solution_type", "")).strip()]
    nss_types = tuple(sorted({str(r["nss_solution_type"]).strip() for r in nss_rows}))
    astrometric_nss = next(
        (r for r in nss_rows if _is_astrometric_nss_type(r.get("nss_solution_type", ""))),
        None,
    )

    campbell = None
    if astrometric_nss is not None:
        constants = [_as_float(astrometric_nss.get(k)) for k in (
            "a_thiele_innes", "b_thiele_innes", "f_thiele_innes", "g_thiele_innes"
        )]
        if all(np.isfinite(constants)):
            try:
                campbell = thiele_innes_to_campbell(*constants)
            except ValueError:
                campbell = None

    bench = gj765_photocentre_benchmark()
    return DR3TargetSummary(
        source_id=source_id,
        separation_arcsec=_as_float(nearest_row.get("separation_arcsec")),
        hipparcos2_number_of_neighbours=_as_int(nearest_row.get("hipparcos2_number_of_neighbours")),
        hipparcos2_gaia_astrometric_params=_as_int(nearest_row.get("hipparcos2_gaia_astrometric_params")),
        parallax_mas=_as_float(nearest_row.get("parallax")),
        parallax_error_mas=_as_float(nearest_row.get("parallax_error")),
        ruwe=_as_float(nearest_row.get("ruwe")),
        ipd_frac_multi_peak_percent=_as_float(nearest_row.get("ipd_frac_multi_peak")),
        ipd_gof_harmonic_amplitude=_as_float(nearest_row.get("ipd_gof_harmonic_amplitude")),
        ipd_gof_harmonic_phase_deg=_as_float(nearest_row.get("ipd_gof_harmonic_phase")),
        non_single_star_flag=_as_int(nearest_row.get("non_single_star")),
        nss_solution_types=nss_types,
        nss_orbit_count=len(nss_rows),
        nearest_nss_photocentre_axis_mas=(campbell.semi_major_axis_mas if campbell else math.nan),
        nearest_nss_period_days=_as_float(astrometric_nss.get("nss_period_days")) if astrometric_nss else math.nan,
        nearest_nss_eccentricity=_as_float(astrometric_nss.get("nss_eccentricity")) if astrometric_nss else math.nan,
        nearest_nss_inclination_deg=(campbell.inclination_deg if campbell else math.nan),
        nearest_nss_omega_relative_deg=(campbell.omega_relative_deg if campbell else math.nan),
        nearest_nss_node_deg=(campbell.node_deg if campbell else math.nan),
        predicted_M0_photocentre_axis_mas=bench["predicted_M0_photocentre_axis_mas"],
        beta_V_proxy=bench["beta_V_proxy"],
        mass_fraction_secondary=bench["mass_fraction_secondary"],
    )


def write_summary_json(summary: DR3TargetSummary, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(asdict(summary), indent=2, allow_nan=True) + "\n")
