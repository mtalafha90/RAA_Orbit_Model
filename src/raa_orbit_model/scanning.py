from __future__ import annotations

from dataclasses import dataclass
import csv
from pathlib import Path
import numpy as np


_RELEASES = ("dr1", "dr2", "dr3", "dr4", "dr5")


@dataclass(frozen=True)
class GaiaScanSchedule:
    """Mission-relative Gaia field-of-view transit times and scan angles.

    ``times_yr`` are years from the scanning-law mission start used by the
    provider. ``scan_angle_deg`` is the *directional* Gaia along-scan position
    angle: 0 deg points North and 90 deg points East.  Angles are deliberately
    not reduced modulo 180 deg because reversing the AL direction reverses the
    sign of a one-dimensional AL coordinate.
    """

    times_yr: np.ndarray
    scan_angle_deg: np.ndarray
    ra_deg: float
    dec_deg: float
    release: str
    source: str
    mission_start_decimalyear: float | None = None

    def validate(self) -> None:
        t = np.asarray(self.times_yr, dtype=float)
        psi = np.asarray(self.scan_angle_deg, dtype=float)
        if t.ndim != 1 or psi.ndim != 1 or len(t) != len(psi):
            raise ValueError("times_yr and scan_angle_deg must be 1-D arrays of equal length")
        if len(t) == 0:
            raise ValueError("Gaia scan schedule contains no transits")
        if not np.all(np.isfinite(t)) or not np.all(np.isfinite(psi)):
            raise ValueError("Gaia scan schedule contains non-finite values")
        if np.any(np.diff(t) < 0):
            raise ValueError("Gaia scan times must be sorted")
        if not (0.0 <= float(self.ra_deg) < 360.0):
            raise ValueError("ra_deg must satisfy 0 <= RA < 360")
        if not (-90.0 <= float(self.dec_deg) <= 90.0):
            raise ValueError("dec_deg must satisfy -90 <= Dec <= 90")
        if self.release not in _RELEASES and self.release != "custom":
            raise ValueError(f"release must be one of {_RELEASES} or 'custom'")

    @property
    def n_transits(self) -> int:
        return len(self.times_yr)

    @property
    def duration_yr(self) -> float:
        t = np.asarray(self.times_yr, dtype=float)
        return float(t[-1] - t[0]) if len(t) > 1 else 0.0

    @property
    def mission_span_yr(self) -> float:
        """Time from provider mission start to the final retained transit."""
        return float(np.asarray(self.times_yr, dtype=float)[-1])


def schedule_from_arrays(
    times_yr,
    scan_angle_deg,
    *,
    ra_deg: float,
    dec_deg: float,
    release: str = "custom",
    source: str = "user-supplied",
    mission_start_decimalyear: float | None = None,
) -> GaiaScanSchedule:
    t = np.asarray(times_yr, dtype=float)
    psi = np.asarray(scan_angle_deg, dtype=float)
    order = np.argsort(t)
    schedule = GaiaScanSchedule(
        times_yr=t[order],
        scan_angle_deg=psi[order],
        ra_deg=float(ra_deg),
        dec_deg=float(dec_deg),
        release=str(release).lower(),
        source=str(source),
        mission_start_decimalyear=(None if mission_start_decimalyear is None else float(mission_start_decimalyear)),
    )
    schedule.validate()
    return schedule


def schedule_from_gaiascanlaw(
    ra_deg: float,
    dec_deg: float,
    *,
    release: str = "dr4",
    obstype: str | None = "astrometry",
    _module=None,
) -> GaiaScanSchedule:
    """Return a GOST-derived nominal Gaia schedule using ``gaiascanlaw``.

    The external ``gaiascanlaw`` package is an optional research dependency.
    It returns decimal-year transit times and scan angles in radians from a
    GOST-derived level-6 HEALPix schedule.  We convert times to years from the
    provider's mission start and angles to degrees, preserving the scan
    direction.

    Parameters
    ----------
    ra_deg, dec_deg
        ICRS sky position in degrees.
    release
        One of dr1, dr2, dr3, dr4, dr5.  In this project dr4 means the
        nominal 66-month DR4 baseline and dr5 the full science-operation
        nominal schedule ending in January 2025, as encoded by gaiascanlaw.
    obstype
        Passed to gaiascanlaw. ``'astrometry'`` applies its published Gaia
        astrometric gap mask; use ``None`` for the uninterrupted nominal law.
    """
    release = str(release).lower()
    if release not in _RELEASES:
        raise ValueError(f"release must be one of {_RELEASES}")
    if obstype not in (None, "astrometry", "photometry", "spectroscopy"):
        raise ValueError("obstype must be None, 'astrometry', 'photometry', or 'spectroscopy'")

    if _module is None:
        try:
            import gaiascanlaw as gsl
        except ImportError as exc:
            raise ImportError(
                "Mission-grounded schedule generation requires the optional "
                "'gaiascanlaw' package. Install it separately or with the "
                "project's scanlaw extra, or load a previously exported CSV schedule."
            ) from exc
    else:
        gsl = _module

    end_name = f"t{release}"
    if not hasattr(gsl, "tstart") or not hasattr(gsl, end_name) or not hasattr(gsl, "scanlaw"):
        raise RuntimeError("installed gaiascanlaw module does not expose the expected scanlaw/tstart/release API")

    t_start = float(gsl.tstart)
    t_end = float(getattr(gsl, end_name))
    times_decimalyear, angles_rad = gsl.scanlaw(
        float(ra_deg),
        float(dec_deg),
        tstart=t_start,
        tend=t_end,
        obstype=obstype,
    )
    times_decimalyear = np.asarray(times_decimalyear, dtype=float)
    angles_rad = np.asarray(angles_rad, dtype=float)
    source = "gaiascanlaw/GOST nominal scanning law"
    if obstype is not None:
        source += f"; {obstype} gaps applied"
    return schedule_from_arrays(
        times_decimalyear - t_start,
        np.rad2deg(angles_rad),
        ra_deg=ra_deg,
        dec_deg=dec_deg,
        release=release,
        source=source,
        mission_start_decimalyear=t_start,
    )


def write_schedule_csv(schedule: GaiaScanSchedule, path: str | Path) -> None:
    """Archive the exact schedule used by an experiment in a portable CSV."""
    schedule.validate()
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=(
                "mission_time_yr", "scan_angle_deg", "ra_deg", "dec_deg",
                "release", "source", "mission_start_decimalyear",
            ),
        )
        writer.writeheader()
        for t, psi in zip(schedule.times_yr, schedule.scan_angle_deg):
            writer.writerow(
                {
                    "mission_time_yr": f"{float(t):.12g}",
                    "scan_angle_deg": f"{float(psi):.12g}",
                    "ra_deg": f"{schedule.ra_deg:.12g}",
                    "dec_deg": f"{schedule.dec_deg:.12g}",
                    "release": schedule.release,
                    "source": schedule.source,
                    "mission_start_decimalyear": (
                        "" if schedule.mission_start_decimalyear is None
                        else f"{schedule.mission_start_decimalyear:.12g}"
                    ),
                }
            )


def schedule_from_csv(path: str | Path) -> GaiaScanSchedule:
    """Load an exact schedule previously archived with ``write_schedule_csv``."""
    path = Path(path)
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError(f"schedule CSV is empty: {path}")

    first = rows[0]
    required = {"mission_time_yr", "scan_angle_deg", "ra_deg", "dec_deg", "release", "source"}
    missing = required - set(first)
    if missing:
        raise ValueError(f"schedule CSV is missing columns: {sorted(missing)}")
    mission_start = first.get("mission_start_decimalyear", "").strip()
    return schedule_from_arrays(
        [float(row["mission_time_yr"]) for row in rows],
        [float(row["scan_angle_deg"]) for row in rows],
        ra_deg=float(first["ra_deg"]),
        dec_deg=float(first["dec_deg"]),
        release=first["release"],
        source=first["source"],
        mission_start_decimalyear=(float(mission_start) if mission_start else None),
    )
