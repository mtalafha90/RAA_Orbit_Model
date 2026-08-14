"""Absolute astrometric motion: parallax factors and linear proper motion.

The along-scan channel originally modelled the orbital wobble alone. That is
sufficient to expose a measurement-model mismatch in isolation, but it omits
the mechanism by which such a mismatch reaches published Gaia parameters: the
orbital signature is fitted *simultaneously* with parallax and proper motion,
and it is their covariance — worst when the period approaches one year — that
converts a wrong measurement model into a wrong parallax.

This module supplies the standard five-parameter terms so the orbit can be
fitted jointly with them.

The observer position uses a low-precision analytic Earth ephemeris. Gaia
orbits the Sun-Earth L2 point roughly 0.01 au from Earth, so this reproduces
the parallax factors to about one per cent of the parallax amplitude. That is
appropriate for a controlled synthetic study and is *not* adequate for fitting
real Gaia epoch astrometry, which must use the actual spacecraft ephemeris
distributed with the epoch data.
"""

from __future__ import annotations

import numpy as np

OBLIQUITY_J2000_DEG = 23.439291111
J2000_DECIMALYEAR = 2000.0
DAYS_PER_YEAR = 365.25


def earth_barycentric_position_au(decimalyear):
    """Low-precision barycentric position of the Earth, equatorial, in au.

    Uses the standard low-precision solar-coordinates series. The Earth's
    heliocentric position is the negative of the Sun's geocentric position.
    """
    t = np.asarray(decimalyear, dtype=float)
    days = (t - J2000_DECIMALYEAR) * DAYS_PER_YEAR

    mean_anomaly = np.deg2rad(357.528 + 0.9856003 * days)
    mean_longitude = 280.460 + 0.9856474 * days
    ecliptic_longitude = np.deg2rad(
        mean_longitude
        + 1.915 * np.sin(mean_anomaly)
        + 0.020 * np.sin(2.0 * mean_anomaly)
    )
    radius_au = (
        1.00014
        - 0.01671 * np.cos(mean_anomaly)
        - 0.00014 * np.cos(2.0 * mean_anomaly)
    )

    obliquity = np.deg2rad(OBLIQUITY_J2000_DEG)
    x_ecliptic = -radius_au * np.cos(ecliptic_longitude)
    y_ecliptic = -radius_au * np.sin(ecliptic_longitude)
    return np.column_stack((
        x_ecliptic,
        y_ecliptic * np.cos(obliquity),
        y_ecliptic * np.sin(obliquity),
    ))


def parallax_factors(decimalyear, ra_deg: float, dec_deg: float):
    """Return the (alpha*, delta) parallax factors for a sky position.

    The parallactic displacement of a star of parallax ``varpi`` is
    ``varpi * P_alpha`` in ``alpha*`` and ``varpi * P_delta`` in ``delta``,
    with the observer's barycentric position ``(X, Y, Z)`` in au:

        P_alpha = X sin(alpha) - Y cos(alpha)
        P_delta = X cos(alpha) sin(delta) + Y sin(alpha) sin(delta) - Z cos(delta)
    """
    position = earth_barycentric_position_au(decimalyear)
    X, Y, Z = position[:, 0], position[:, 1], position[:, 2]
    alpha = np.deg2rad(float(ra_deg))
    delta = np.deg2rad(float(dec_deg))
    p_alpha_star = X * np.sin(alpha) - Y * np.cos(alpha)
    p_delta = (
        X * np.cos(alpha) * np.sin(delta)
        + Y * np.sin(alpha) * np.sin(delta)
        - Z * np.cos(delta)
    )
    return p_alpha_star, p_delta


def absolute_offsets_mas(
    times_yr,
    params,
    *,
    ra_deg: float,
    dec_deg: float,
    mission_start_decimalyear: float,
    reference_time_yr: float = 0.0,
):
    """Barycentric position offset from position, proper motion and parallax.

    ``times_yr`` are mission-relative years, matching the scan schedule.
    Returns ``(delta_alpha_star_mas, delta_delta_mas)`` excluding the orbit.
    """
    t = np.asarray(times_yr, dtype=float)
    dt = t - float(reference_time_yr)
    p_alpha_star, p_delta = parallax_factors(
        float(mission_start_decimalyear) + t, ra_deg, dec_deg
    )
    alpha_star = (
        float(params.delta_alpha_star_mas)
        + float(params.pmra_mas_yr) * dt
        + float(params.parallax_mas) * p_alpha_star
    )
    delta = (
        float(params.delta_delta_mas)
        + float(params.pmdec_mas_yr) * dt
        + float(params.parallax_mas) * p_delta
    )
    return alpha_star, delta
