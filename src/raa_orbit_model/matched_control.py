from __future__ import annotations

import numpy as np

from .scanning import GaiaScanSchedule, schedule_from_arrays


def stratified_transit_indices(n_total: int, n_keep: int) -> np.ndarray:
    """Return deterministic equal-count-stratum indices across an ordered schedule."""
    n_total = int(n_total)
    n_keep = int(n_keep)
    if n_total <= 0 or n_keep <= 0:
        raise ValueError("n_total and n_keep must be > 0")
    if n_keep > n_total:
        raise ValueError("n_keep cannot exceed n_total")
    if n_keep == n_total:
        return np.arange(n_total, dtype=int)
    indices = np.floor((np.arange(n_keep, dtype=float) + 0.5) * n_total / n_keep).astype(int)
    if len(np.unique(indices)) != n_keep:
        raise RuntimeError("matched-N selection did not produce unique indices")
    return indices


def matched_n_schedule(schedule: GaiaScanSchedule, n_keep: int) -> tuple[GaiaScanSchedule, np.ndarray]:
    """Return a rank-stratified matched-N subschedule and its native indices."""
    schedule.validate()
    indices = stratified_transit_indices(schedule.n_transits, n_keep)
    selected = schedule_from_arrays(
        np.asarray(schedule.times_yr, dtype=float)[indices],
        np.asarray(schedule.scan_angle_deg, dtype=float)[indices],
        ra_deg=schedule.ra_deg,
        dec_deg=schedule.dec_deg,
        release=schedule.release,
        source=f"{schedule.source}; matched-N rank-stratified control ({n_keep}/{schedule.n_transits})",
        mission_start_decimalyear=schedule.mission_start_decimalyear,
    )
    return selected, indices
