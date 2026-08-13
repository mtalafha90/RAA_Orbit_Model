from __future__ import annotations
from dataclasses import dataclass
import numpy as np

AU_PER_YR_TO_KMS = 4.740470463533349
TWOPI = 2.0 * np.pi

@dataclass(frozen=True)
class BinaryParams:
    period_yr: float
    t_peri_yr: float
    eccentricity: float
    inclination_deg: float
    omega_deg: float
    node_deg: float
    m1_msun: float
    m2_msun: float
    parallax_mas: float
    gamma_kms: float = 0.0
    beta_g: float = 0.0

    def validate(self) -> None:
        if self.period_yr <= 0:
            raise ValueError("period_yr must be > 0")
        if not (0.0 <= self.eccentricity < 1.0):
            raise ValueError("eccentricity must satisfy 0 <= e < 1")
        if self.m1_msun <= 0 or self.m2_msun <= 0:
            raise ValueError("component masses must be > 0")
        if self.parallax_mas <= 0:
            raise ValueError("parallax_mas must be > 0")
        if not (0.0 <= self.beta_g <= 1.0):
            raise ValueError("beta_g must satisfy 0 <= beta_g <= 1")

    @property
    def q(self) -> float:
        return self.m2_msun / self.m1_msun

    @property
    def mass_fraction_secondary(self) -> float:
        return self.m2_msun / (self.m1_msun + self.m2_msun)

    @property
    def a_rel_au(self) -> float:
        return ((self.m1_msun + self.m2_msun) * self.period_yr**2) ** (1.0 / 3.0)

    @property
    def a1_au(self) -> float:
        return self.a_rel_au * self.m2_msun / (self.m1_msun + self.m2_msun)

    @property
    def a2_au(self) -> float:
        return self.a_rel_au * self.m1_msun / (self.m1_msun + self.m2_msun)

    @property
    def omega_relative_deg(self) -> float:
        return (self.omega_deg + 180.0) % 360.0

def solve_kepler(mean_anomaly_rad, eccentricity: float, tol: float = 1e-13, max_iter: int = 100):
    if not (0.0 <= eccentricity < 1.0):
        raise ValueError("elliptic Kepler solver requires 0 <= e < 1")
    M = np.asarray(mean_anomaly_rad, dtype=float)
    Mwrapped = (M + np.pi) % TWOPI - np.pi
    E = Mwrapped.copy()
    if eccentricity >= 0.8:
        E = np.where(Mwrapped >= 0, np.pi, -np.pi)
    for _ in range(max_iter):
        f = E - eccentricity * np.sin(E) - Mwrapped
        fp = 1.0 - eccentricity * np.cos(E)
        step = f / fp
        E -= step
        if np.all(np.abs(step) < tol):
            return E + (M - Mwrapped)
    raise RuntimeError("Kepler solver did not converge")

def anomalies(times_yr, p: BinaryParams):
    p.validate()
    t = np.asarray(times_yr, dtype=float)
    M = TWOPI * (t - p.t_peri_yr) / p.period_yr
    E = solve_kepler(M, p.eccentricity)
    e = p.eccentricity
    nu = 2.0 * np.arctan2(
        np.sqrt(1 + e) * np.sin(E / 2.0),
        np.sqrt(1 - e) * np.cos(E / 2.0),
    )
    return M, E, nu

def relative_position_au(times_yr, p: BinaryParams):
    """r2-r1 as (East, North, line-of-sight away) in AU.

    omega_deg is the primary RV argument omega1. Relative astrometry uses
    omega1+180 deg. node_deg is position angle from North through East.
    """
    _, E, nu = anomalies(times_yr, p)
    r = p.a_rel_au * (1.0 - p.eccentricity * np.cos(E))
    u = nu + np.deg2rad(p.omega_deg + 180.0)
    inc = np.deg2rad(p.inclination_deg)
    Om = np.deg2rad(p.node_deg)
    cu, su = np.cos(u), np.sin(u)
    cO, sO = np.cos(Om), np.sin(Om)
    ci, si = np.cos(inc), np.sin(inc)
    north = r * (cO * cu - sO * su * ci)
    east = r * (sO * cu + cO * su * ci)
    los = r * su * si
    return np.column_stack((east, north, los))

def relative_astrometry_mas(times_yr, p: BinaryParams):
    """(Delta alpha*, Delta delta) of component 2 relative to 1 in mas."""
    return relative_position_au(times_yr, p)[:, :2] * p.parallax_mas

def rv_semiamplitudes_kms(p: BinaryParams):
    p.validate()
    sini = np.sin(np.deg2rad(p.inclination_deg))
    fac = TWOPI * sini / (p.period_yr * np.sqrt(1.0 - p.eccentricity**2))
    return p.a1_au * fac * AU_PER_YR_TO_KMS, p.a2_au * fac * AU_PER_YR_TO_KMS

def radial_velocities_kms(times_yr, p: BinaryParams):
    """Primary/secondary RVs; positive is receding."""
    _, _, nu = anomalies(times_yr, p)
    K1, K2 = rv_semiamplitudes_kms(p)
    omega1 = np.deg2rad(p.omega_deg)
    f = np.cos(nu + omega1) + p.eccentricity * np.cos(omega1)
    return np.column_stack((p.gamma_kms + K1 * f, p.gamma_kms - K2 * f))
