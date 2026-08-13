import numpy as np
import pytest

from raa_orbit_model.matched_control import matched_n_schedule, stratified_transit_indices
from raa_orbit_model.scanning import schedule_from_arrays


def test_stratified_indices_are_unique_and_span_ordered_schedule():
    idx = stratified_transit_indices(100, 53)
    assert len(idx) == 53
    assert len(np.unique(idx)) == 53
    assert np.all(np.diff(idx) > 0)
    assert idx[0] <= 1
    assert idx[-1] >= 98


def test_equal_counts_keep_every_native_transit():
    assert np.array_equal(stratified_transit_indices(53, 53), np.arange(53))
    with pytest.raises(ValueError):
        stratified_transit_indices(52, 53)


def test_matched_schedule_is_an_exact_native_subset():
    schedule = schedule_from_arrays(
        np.linspace(0.0, 5.0, 80),
        np.linspace(0.0, 359.0, 80),
        ra_deg=12.0,
        dec_deg=-30.0,
        release="dr4",
        source="unit-test",
    )
    matched, idx = matched_n_schedule(schedule, 53)
    assert matched.n_transits == 53
    assert np.allclose(matched.times_yr, schedule.times_yr[idx])
    assert np.allclose(matched.scan_angle_deg, schedule.scan_angle_deg[idx])
    assert "rank-stratified control" in matched.source
