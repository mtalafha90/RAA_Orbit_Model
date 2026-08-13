import pytest

from raa_orbit_model.kepler import BinaryParams
from raa_orbit_model.scanning import schedule_from_arrays
from raa_orbit_model.sky_experiments import sky_resolution_bias_scan
from raa_orbit_model.sky_study import augment_bias_rows, build_ecliptic_grid, scan_geometry_diagnostics


def _identity_transform(lon, lat):
    return lon, lat


def _truth():
    return BinaryParams(period_yr=2.0, t_peri_yr=0.15, eccentricity=0.25, inclination_deg=72.0, omega_deg=55.0, node_deg=120.0, m1_msun=1.25, m2_msun=0.85, parallax_mas=20.0, gamma_kms=7.0, beta_g=0.20)


def test_default_ecliptic_grid_has_25_unique_positions_and_one_pole():
    grid = build_ecliptic_grid(transform=_identity_transform)
    assert len(grid) == 25
    assert len({p.sky_id for p in grid}) == 25
    pole = [p for p in grid if p.ecliptic_lat_deg == 90.0]
    assert len(pole) == 1
    assert pole[0].ecliptic_lon_deg == pytest.approx(0.0)


def test_custom_grid_supports_negative_latitudes_and_deduplicates_both_poles():
    grid = build_ecliptic_grid(latitudes_deg=(-90.0, -30.0, 90.0), longitudes_deg=(0.0, 90.0, 180.0, 270.0), transform=_identity_transform)
    assert len(grid) == 6
    assert sum(abs(p.ecliptic_lat_deg) == 90.0 for p in grid) == 2


def test_directional_scan_diagnostics_do_not_fold_180_degrees():
    schedule = schedule_from_arrays([0.0, 1.0, 2.0, 3.0], [0.0, 90.0, 180.0, 270.0], ra_deg=10.0, dec_deg=20.0)
    d = scan_geometry_diagnostics(schedule, entropy_bins=4)
    assert d.n_transits == 4
    assert d.directional_resultant_length == pytest.approx(0.0, abs=1e-12)
    assert d.directional_circular_variance == pytest.approx(1.0)
    assert d.directional_max_gap_deg == pytest.approx(90.0)
    assert d.directional_entropy_normalized == pytest.approx(1.0)
    assert d.directional_effective_bins == pytest.approx(4.0)
    assert d.directional_occupied_bins == 4
    assert 0.0 <= d.time_angle_r2 <= 1.0


def test_single_direction_has_full_circular_gap_and_one_effective_bin():
    schedule = schedule_from_arrays([0.0, 1.0, 2.0, 3.0], [15.0, 15.0, 15.0, 15.0], ra_deg=10.0, dec_deg=20.0)
    d = scan_geometry_diagnostics(schedule, entropy_bins=12)
    assert d.directional_resultant_length == pytest.approx(1.0)
    assert d.directional_circular_variance == pytest.approx(0.0)
    assert d.directional_max_gap_deg == pytest.approx(360.0)
    assert d.directional_entropy_normalized == pytest.approx(0.0)
    assert d.directional_effective_bins == pytest.approx(1.0)
    assert d.directional_occupied_bins == 1


def test_bias_rows_receive_sky_and_scan_metadata():
    grid = build_ecliptic_grid(latitudes_deg=(30.0,), longitudes_deg=(90.0,), transform=_identity_transform)
    schedule = schedule_from_arrays([0.0, 1.0, 2.0], [0.0, 120.0, 240.0], ra_deg=90.0, dec_deg=30.0)
    diagnostics = scan_geometry_diagnostics(schedule, entropy_bins=6)
    rows = augment_bias_rows([{"model": "photocentre"}], grid[0], diagnostics)
    assert rows[0]["sky_id"] == grid[0].sky_id
    assert rows[0]["ecliptic_lat_deg"] == pytest.approx(30.0)
    assert rows[0]["scan_n_transits"] == 3


def test_external_baseline_must_be_positive_before_fit():
    schedule = schedule_from_arrays([0.1, 0.5, 1.0], [0.0, 90.0, 180.0], ra_deg=10.0, dec_deg=20.0)
    with pytest.raises(ValueError, match="external_baseline_yr"):
        sky_resolution_bias_scan(_truth(), schedule, a_over_sigma_values=(0.5,), beta_values=(0.2,), seeds=(0,), external_baseline_yr=0.0)
