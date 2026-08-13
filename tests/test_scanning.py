import numpy as np

from raa_orbit_model.scanning import (
    schedule_from_arrays,
    schedule_from_csv,
    schedule_from_gaiascanlaw,
    write_schedule_csv,
)


class FakeGaiaScanLaw:
    tstart = 2014.5
    tdr1 = 2015.5
    tdr2 = 2016.3
    tdr3 = 2017.4
    tdr4 = 2020.0
    tdr5 = 2025.0

    @staticmethod
    def scanlaw(ra, dec, tstart, tend, obstype=None):
        assert ra == 20.0
        assert dec == 45.0
        assert tstart == FakeGaiaScanLaw.tstart
        assert tend == FakeGaiaScanLaw.tdr4
        assert obstype == "astrometry"
        # Preserve directional angles: 200 deg must not be reduced modulo 180.
        return np.array([2014.75, 2015.10, 2015.12, 2019.80]), np.deg2rad([5.0, 200.0, 72.0, 319.0])


def test_gaiascanlaw_adapter_converts_to_mission_relative_years_and_degrees():
    schedule = schedule_from_gaiascanlaw(
        20.0,
        45.0,
        release="dr4",
        obstype="astrometry",
        _module=FakeGaiaScanLaw,
    )
    assert np.allclose(schedule.times_yr, [0.25, 0.60, 0.62, 5.30])
    assert np.allclose(schedule.scan_angle_deg, [5.0, 200.0, 72.0, 319.0])
    assert schedule.n_transits == 4
    assert schedule.release == "dr4"
    assert "GOST" in schedule.source


def test_schedule_csv_roundtrip(tmp_path):
    original = schedule_from_arrays(
        [0.2, 0.8, 1.4],
        [10.0, 185.0, 270.0],
        ra_deg=123.4,
        dec_deg=-21.5,
        release="dr3",
        source="test schedule",
        mission_start_decimalyear=2014.56,
    )
    path = tmp_path / "schedule.csv"
    write_schedule_csv(original, path)
    loaded = schedule_from_csv(path)
    assert np.allclose(loaded.times_yr, original.times_yr)
    assert np.allclose(loaded.scan_angle_deg, original.scan_angle_deg)
    assert loaded.ra_deg == original.ra_deg
    assert loaded.dec_deg == original.dec_deg
    assert loaded.release == original.release
    assert loaded.source == original.source
    assert loaded.mission_start_decimalyear == original.mission_start_decimalyear
