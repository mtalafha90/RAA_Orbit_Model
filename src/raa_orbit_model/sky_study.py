from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Callable, Iterable

import numpy as np

from .scanning import GaiaScanSchedule

DEFAULT_ECLIPTIC_LATITUDES_DEG = (0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0)
DEFAULT_ECLIPTIC_LONGITUDES_DEG = (0.0, 90.0, 180.0, 270.0)

@dataclass(frozen=True)
class SkyPosition:
    sky_id: str
    ecliptic_lon_deg: float
    ecliptic_lat_deg: float
    ra_deg: float
    dec_deg: float

@dataclass(frozen=True)
class ScanGeometryDiagnostics:
    n_transits: int
    duration_yr: float
    mission_span_yr: float
    directional_resultant_length: float
    directional_circular_variance: float
    directional_max_gap_deg: float
    directional_entropy_bins: int
    directional_entropy_normalized: float
    directional_effective_bins: float
    directional_occupied_bins: int
    time_angle_r2: float

def _sky_id(lon_deg: float, lat_deg: float) -> str:
    lat_tag = ("p" if lat_deg >= 0 else "m") + f"{abs(lat_deg):05.1f}".replace(".", "p")
    lon_tag = f"{lon_deg % 360.0:06.1f}".replace(".", "p")
    return f"ecl_b{lat_tag}_l{lon_tag}"

def ecliptic_to_icrs(lon_deg: float, lat_deg: float) -> tuple[float, float]:
    try:
        import astropy.units as u
        from astropy.coordinates import BarycentricMeanEcliptic, SkyCoord
        from astropy.time import Time
    except ImportError as exc:
        raise ImportError("Ecliptic sky-grid generation requires Astropy; install the optional sky dependencies.") from exc
    frame = BarycentricMeanEcliptic(equinox=Time("J2000"))
    coord = SkyCoord(lon=float(lon_deg) * u.deg, lat=float(lat_deg) * u.deg, frame=frame)
    icrs = coord.icrs
    return float(icrs.ra.deg % 360.0), float(icrs.dec.deg)

def build_ecliptic_grid(latitudes_deg: Iterable[float] = DEFAULT_ECLIPTIC_LATITUDES_DEG, longitudes_deg: Iterable[float] = DEFAULT_ECLIPTIC_LONGITUDES_DEG, *, transform: Callable[[float, float], tuple[float, float]] | None = None) -> list[SkyPosition]:
    transform = ecliptic_to_icrs if transform is None else transform
    latitudes = tuple(float(v) for v in latitudes_deg)
    longitudes = tuple(float(v) % 360.0 for v in longitudes_deg)
    if not latitudes or not longitudes:
        raise ValueError("ecliptic latitude and longitude grids must be non-empty")
    if any(not np.isfinite(v) or v < -90.0 or v > 90.0 for v in latitudes):
        raise ValueError("ecliptic latitudes must be finite and in [-90, 90] deg")
    if any(not np.isfinite(v) for v in longitudes):
        raise ValueError("ecliptic longitudes must be finite")
    positions = []
    for lat in latitudes:
        lons = (0.0,) if math.isclose(abs(lat), 90.0, abs_tol=1e-12) else longitudes
        for lon in lons:
            ra, dec = transform(lon, lat)
            positions.append(SkyPosition(_sky_id(lon, lat), float(lon), float(lat), float(ra) % 360.0, float(dec)))
    return positions

def _time_angle_r2(times_yr: np.ndarray, angles_rad: np.ndarray) -> float:
    if len(times_yr) < 3 or np.ptp(times_yr) == 0.0:
        return float("nan")
    y = times_yr - np.mean(times_yr)
    x = np.column_stack((np.sin(angles_rad), np.cos(angles_rad), np.ones(len(angles_rad))))
    coef, *_ = np.linalg.lstsq(x, y, rcond=None)
    resid = y - x @ coef
    sst = float(np.dot(y, y))
    if sst == 0.0:
        return float("nan")
    return float(np.clip(1.0 - float(np.dot(resid, resid)) / sst, 0.0, 1.0))

def scan_geometry_diagnostics(schedule: GaiaScanSchedule, *, entropy_bins: int = 12) -> ScanGeometryDiagnostics:
    schedule.validate()
    if entropy_bins < 2:
        raise ValueError("entropy_bins must be >= 2")
    angles_deg = np.mod(np.asarray(schedule.scan_angle_deg, dtype=float), 360.0)
    angles_rad = np.deg2rad(angles_deg)
    resultant = float(abs(np.mean(np.exp(1j * angles_rad))))
    ordered = np.sort(angles_deg)
    max_gap = float(np.max(np.diff(np.concatenate((ordered, [ordered[0] + 360.0])))))
    counts, _ = np.histogram(angles_deg, bins=np.linspace(0.0, 360.0, entropy_bins + 1))
    probabilities = counts[counts > 0].astype(float) / float(len(angles_deg))
    entropy = float(-np.sum(probabilities * np.log(probabilities)))
    return ScanGeometryDiagnostics(schedule.n_transits, schedule.duration_yr, schedule.mission_span_yr, resultant, 1.0 - resultant, max_gap, int(entropy_bins), entropy / math.log(entropy_bins), float(math.exp(entropy)), int(np.count_nonzero(counts)), _time_angle_r2(np.asarray(schedule.times_yr, dtype=float), angles_rad))

def position_record(position: SkyPosition, diagnostics: ScanGeometryDiagnostics) -> dict[str, float | int | str]:
    row = asdict(position)
    row.update({f"scan_{key}": value for key, value in asdict(diagnostics).items()})
    return row

def augment_bias_rows(rows: list[dict], position: SkyPosition, diagnostics: ScanGeometryDiagnostics) -> list[dict]:
    metadata = position_record(position, diagnostics)
    for row in rows:
        row.update(metadata)
    return rows
