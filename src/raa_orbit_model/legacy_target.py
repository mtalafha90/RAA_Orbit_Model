"""Reproduce the real-binary V6a fit from a legacy visual + SB2 input file.

`docs/real_data_validation.md` records V6a as executed, and the manuscript
reports a full parameter table for GJ 765.2 (`tab:gl765`). The measurements
behind it are **not in this repository**: there is no data file, and no runner
that regenerates the table. Every synthetic result here has a frozen product
and a documented runner; the one real-stellar result did not.

This module closes the code half of that gap. Supply the legacy file and
``scripts/run_legacy_target_fit.py`` regenerates the fit and writes a frozen
CSV alongside the synthetic products.

Input format
------------
Reconstructed from the description in the manuscript, which states that the
astrometry is given as position angle and separation converted through

    delta_alpha* = rho sin(theta),   delta_delta = rho cos(theta)

with ``theta`` measured North through East, and that "the fourth astrometric
column is treated as an isotropic one-sigma tangent-plane positional
uncertainty". The reader is therefore whitespace- or comma-separated with

    astrometry:  epoch_yr   theta_deg   rho_arcsec_or_mas   sigma
    velocities:  epoch_yr   rv1_kms     rv2_kms             sigma1  [sigma2]

Blocks are introduced by a line containing ``ASTROMETRY`` or ``VELOCITIES``
(case-insensitive); ``#`` begins a comment. Header ``key = value`` pairs are
captured but never used as priors, matching the manuscript's statement that
header values initialise only. In particular the legacy header parallax of
54.27 mas and the header coordinates are known to be wrong and must not be
trusted.

**If your file differs, adjust `parse_legacy_file` rather than editing the
data.** The separation unit is explicit (``--separation-unit``) because a
legacy file may quote either arcseconds or milliarcseconds, and guessing would
silently rescale the orbit.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .kepler import BinaryParams
from .synthetic import GaiaALData, JointData, RelativeAstrometryData, SB2RVData

_NUMBER = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][-+]?\d+)?")


@dataclass(frozen=True)
class LegacyTargetData:
    """Parsed legacy visual-plus-SB2 measurements for one real binary."""

    header: dict[str, str]
    astrometry_epoch_yr: np.ndarray
    position_angle_deg: np.ndarray
    separation_mas: np.ndarray
    astrometry_sigma_mas: np.ndarray
    rv_epoch_yr: np.ndarray
    rv_primary_kms: np.ndarray
    rv_secondary_kms: np.ndarray
    rv_primary_sigma_kms: np.ndarray
    rv_secondary_sigma_kms: np.ndarray

    @property
    def n_astrometry(self) -> int:
        return len(self.astrometry_epoch_yr)

    @property
    def n_rv(self) -> int:
        return len(self.rv_epoch_yr)

    @property
    def n_constraints(self) -> int:
        """Scalar constraints: two per relative position, two per RV epoch."""
        return 2 * self.n_astrometry + 2 * self.n_rv

    def header_float(self, key: str) -> float | None:
        for name, value in self.header.items():
            if name.lower() == key.lower():
                match = _NUMBER.search(value)
                if match:
                    return float(match.group().replace("d", "e").replace("D", "e"))
        return None


def _numbers(line: str) -> list[float]:
    return [float(m.replace("d", "e").replace("D", "e")) for m in _NUMBER.findall(line)]


def parse_legacy_file(
    path: str | Path,
    *,
    separation_unit: str = "arcsec",
) -> LegacyTargetData:
    """Parse a legacy visual + SB2 orbit input file.

    ``separation_unit`` must be ``"arcsec"`` or ``"mas"`` and is applied to both
    the separation and its uncertainty.
    """
    if separation_unit not in ("arcsec", "mas"):
        raise ValueError("separation_unit must be 'arcsec' or 'mas'")
    scale = 1000.0 if separation_unit == "arcsec" else 1.0

    header: dict[str, str] = {}
    astrometry: list[list[float]] = []
    velocities: list[list[float]] = []
    block = None

    for raw in Path(path).read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        upper = line.upper()
        if "ASTROMETR" in upper and not _numbers(line):
            block = "astrometry"
            continue
        if ("VELOCIT" in upper or "RADIAL" in upper) and not _numbers(line):
            block = "velocities"
            continue
        if "=" in line and block is None:
            key, value = line.split("=", 1)
            header[key.strip()] = value.strip()
            continue
        values = _numbers(line)
        if not values:
            continue
        if block == "astrometry":
            if len(values) < 4:
                raise ValueError(f"astrometry row needs 4 columns, got {len(values)}: {raw!r}")
            astrometry.append(values)
        elif block == "velocities":
            if len(values) < 4:
                raise ValueError(f"velocity row needs at least 4 columns, got {len(values)}: {raw!r}")
            velocities.append(values)

    if not astrometry:
        raise ValueError(f"no relative astrometry found in {path}")
    if not velocities:
        raise ValueError(f"no radial velocities found in {path}")

    a = np.array([row[:4] for row in astrometry], dtype=float)
    v_primary_sigma = np.array([row[3] for row in velocities], dtype=float)
    v_secondary_sigma = np.array(
        [row[4] if len(row) > 4 else row[3] for row in velocities], dtype=float
    )
    v = np.array([row[:3] for row in velocities], dtype=float)

    return LegacyTargetData(
        header=header,
        astrometry_epoch_yr=a[:, 0],
        position_angle_deg=a[:, 1],
        separation_mas=a[:, 2] * scale,
        astrometry_sigma_mas=a[:, 3] * scale,
        rv_epoch_yr=v[:, 0],
        rv_primary_kms=v[:, 1],
        rv_secondary_kms=v[:, 2],
        rv_primary_sigma_kms=v_primary_sigma,
        rv_secondary_sigma_kms=v_secondary_sigma,
    )


def tangent_plane_offsets_mas(position_angle_deg, separation_mas):
    """Convert polar (theta, rho) to (delta_alpha*, delta_delta) in mas.

    ``theta`` is measured from North through East, matching the convention
    asserted throughout this project.
    """
    theta = np.deg2rad(np.asarray(position_angle_deg, dtype=float))
    rho = np.asarray(separation_mas, dtype=float)
    return rho * np.sin(theta), rho * np.cos(theta)


def legacy_joint_data(data: LegacyTargetData, *, reference_epoch_yr: float | None = None) -> JointData:
    """Build a :class:`JointData` with an empty Gaia channel.

    The legacy set contains no Gaia measurements, so the Gaia likelihood is
    disabled entirely and ``beta_g`` is irrelevant to the fit.
    """
    epoch0 = (
        float(np.min(data.astrometry_epoch_yr)) if reference_epoch_yr is None
        else float(reference_epoch_yr)
    )
    east, north = tangent_plane_offsets_mas(data.position_angle_deg, data.separation_mas)
    values = np.column_stack((east, north))

    # The quoted uncertainty is an isotropic one-sigma tangent-plane error.
    covariance = np.zeros((len(values), 2, 2), dtype=float)
    variance = data.astrometry_sigma_mas.astype(float) ** 2
    covariance[:, 0, 0] = variance
    covariance[:, 1, 1] = variance

    rv_values = np.column_stack((data.rv_primary_kms, data.rv_secondary_kms))
    rv_sigma = np.column_stack((data.rv_primary_sigma_kms, data.rv_secondary_sigma_kms))

    empty = np.empty(0, dtype=float)
    return JointData(
        relative_astrometry=RelativeAstrometryData(
            data.astrometry_epoch_yr - epoch0, values, covariance
        ),
        sb2_rv=SB2RVData(data.rv_epoch_yr - epoch0, rv_values, rv_sigma),
        gaia_al=GaiaALData(empty, empty, empty, empty),
    )


def initial_guess(
    data: LegacyTargetData,
    *,
    reference_epoch_yr: float | None = None,
    period_yr: float,
    eccentricity: float,
    inclination_deg: float,
    node_deg: float,
    omega_relative_deg: float,
    m1_msun: float,
    m2_msun: float,
    parallax_mas: float,
    t_peri_yr: float | None = None,
    gamma_kms: float | None = None,
) -> BinaryParams:
    """Starting point for the legacy fit. Used only to initialise, never as a prior."""
    epoch0 = (
        float(np.min(data.astrometry_epoch_yr)) if reference_epoch_yr is None
        else float(reference_epoch_yr)
    )
    if gamma_kms is None:
        gamma_kms = float(np.median(np.concatenate([data.rv_primary_kms, data.rv_secondary_kms])))
    return BinaryParams(
        period_yr=float(period_yr),
        t_peri_yr=0.0 if t_peri_yr is None else float(t_peri_yr) - epoch0,
        eccentricity=float(eccentricity),
        inclination_deg=float(inclination_deg),
        # The project stores the primary argument of periastron; relative
        # astrometry uses omega + 180 deg.
        omega_deg=(float(omega_relative_deg) - 180.0) % 360.0,
        node_deg=float(node_deg) % 360.0,
        m1_msun=float(m1_msun),
        m2_msun=float(m2_msun),
        parallax_mas=float(parallax_mas),
        gamma_kms=float(gamma_kms),
        beta_g=0.0,
    )


def node_branches(params: BinaryParams) -> tuple[BinaryParams, BinaryParams]:
    """Both visual-orbit node branches.

    Relative astrometry alone leaves the ascending node ambiguous by 180
    degrees. The manuscript states that both branches are initialised and the
    lower-chi-square joint solution retained, so the caller must fit both.
    """
    from dataclasses import replace

    flipped = replace(
        params,
        node_deg=(params.node_deg + 180.0) % 360.0,
        omega_deg=(params.omega_deg + 180.0) % 360.0,
    )
    return params, flipped


def formal_covariance(residual_jacobian: np.ndarray) -> np.ndarray:
    """Local parameter covariance ``(J^T J)^-1`` from whitened residuals.

    The residuals are already whitened, so no further variance scaling applies.
    A singular normal matrix yields an all-infinite covariance rather than an
    exception, so an unconstrained direction degrades gracefully.
    """
    J = np.asarray(residual_jacobian, dtype=float)
    try:
        return np.linalg.inv(J.T @ J)
    except np.linalg.LinAlgError:
        return np.full((J.shape[1], J.shape[1]), np.inf)


def formal_uncertainties(residual_jacobian: np.ndarray) -> np.ndarray:
    """Local one-sigma errors from the weighted least-squares Jacobian."""
    diagonal = np.diag(formal_covariance(residual_jacobian))
    return np.where(diagonal > 0, np.sqrt(np.abs(diagonal)), np.inf)


def total_mass_uncertainty(covariance: np.ndarray, names) -> float:
    """Propagate the correlated component-mass errors onto the total mass.

    The two masses are strongly correlated, so quoting
    ``sqrt(sigma_1^2 + sigma_2^2)`` would misstate the total-mass error. This
    adds the covariance term: ``var(M1+M2) = v11 + v22 + 2 v12``.
    """
    names = tuple(names)
    if "m1_msun" not in names or "m2_msun" not in names:
        return float("nan")
    i, j = names.index("m1_msun"), names.index("m2_msun")
    C = np.asarray(covariance, dtype=float)
    variance = C[i, i] + C[j, j] + 2.0 * C[i, j]
    if not np.isfinite(variance):
        return float("nan")
    if variance < 0.0:
        # Exactly anti-correlated masses give zero variance, which is a valid
        # result: the total is perfectly determined even though neither
        # component is. Only a materially negative value indicates a covariance
        # that cannot be used.
        return 0.0 if variance > -1e-12 else float("nan")
    return float(np.sqrt(variance))


def total_mass_msun(params: BinaryParams) -> float:
    """Label-invariant total mass, which is what the manuscript compares."""
    return float(params.m1_msun + params.m2_msun)


def summarise_fit(result, data: LegacyTargetData) -> dict:
    """Compact record of a legacy-target fit, for freezing alongside the paper."""
    params = result.params
    return {
        "n_astrometry": data.n_astrometry,
        "n_rv_epochs": data.n_rv,
        "n_constraints": data.n_constraints,
        "n_free": result.n_free,
        "chi2": float(result.chi2),
        "dof": int(result.dof),
        "reduced_chi2": float(result.reduced_chi2),
        "success": bool(result.success),
        "period_yr": float(params.period_yr),
        "eccentricity": float(params.eccentricity),
        "inclination_deg": float(params.inclination_deg),
        "node_deg": float(params.node_deg % 360.0),
        "omega_relative_deg": float(params.omega_relative_deg),
        "m1_msun": float(params.m1_msun),
        "m2_msun": float(params.m2_msun),
        "total_mass_msun": total_mass_msun(params),
        "parallax_mas": float(params.parallax_mas),
        "gamma_kms": float(params.gamma_kms),
    }
