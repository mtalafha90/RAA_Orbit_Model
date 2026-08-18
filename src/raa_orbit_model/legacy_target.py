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
The project's own legacy CSV layout. Lines beginning with ``#`` or ``C`` are
ignored; comment lines also introduce the measurement blocks. Fields are
comma-separated::

    #Objectinfo
    Object, GL765.2
    RA,19.404
    Dec,76.1812
    par,54.27

    #Orbital elements
    P,11.769
    T,1993.513
    ...

    # RV1 Measurements
    45533.4644,-10.69,0.51,Va,COR      # time, RV km/s, sigma, label, source

    # RV2 Measurements
    45533.4644,2.81,0.66,Vb,COR

    # Visual Measurements
    1971.6,275.6,0.21,0.04,I1,M1       # epoch, PA deg, rho arcsec, sigma arcsec

**The two measurement blocks use different time systems.** Radial velocities
are Julian dates minus 2 400 000, while visual measurements are decimal years.
Mixing them would displace the velocities by millennia relative to an orbit of
about twelve years, so the reader converts the velocities and verifies that the
two blocks end up overlapping in time.

Separations and their uncertainties are in arcseconds. Position angle is
measured North through East and converted through

    delta_alpha* = rho sin(theta),   delta_delta = rho cos(theta)

with the fourth column treated as an isotropic one-sigma tangent-plane
uncertainty, matching the polar weighting of the legacy orbit code to first
order.

Header values initialise the search only and are never priors. For GJ 765.2
the header parallax of 54.27 mas and the header coordinates are both known to
be wrong, so neither may be trusted.

Component labelling: in this file the ``Va`` curve has the *larger* velocity
amplitude, so it belongs to the *less massive* component. The labels therefore
do not map onto the later A/B convention, and only the label-invariant total
mass should be compared with published solutions.
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
        """Look up a numeric header value.

        The match is **case sensitive** first, because this format uses case to
        distinguish different elements: ``W`` is the node and ``w`` the
        argument of periastron. Folding case would silently conflate them. A
        case-insensitive fallback applies only when no exact key exists and the
        fold is unambiguous.
        """
        def as_float(text):
            match = _NUMBER.search(text)
            return float(match.group()) if match else None

        if key in self.header:
            return as_float(self.header[key])
        folded = [n for n in self.header if n.lower() == key.lower()]
        if len(folded) == 1:
            return as_float(self.header[folded[0]])
        return None


#: Julian date of the epoch J2000.0, used to convert JD - 2400000 to years.
_JD_J2000 = 2451545.0
_JD_OFFSET = 2400000.0
_DAYS_PER_YEAR = 365.25


def jd2400000_to_decimalyear(value):
    """Convert a Julian date minus 2 400 000 to a decimal year."""
    x = np.asarray(value, dtype=float)
    return 2000.0 + (x + _JD_OFFSET - _JD_J2000) / _DAYS_PER_YEAR


def _fields(line: str) -> list[str]:
    return [f.strip() for f in line.split(",")]


def _is_comment(line: str) -> bool:
    stripped = line.strip()
    return (not stripped) or stripped.startswith("#") or stripped.upper().startswith("C,")


def parse_legacy_file(
    path: str | Path,
    *,
    separation_unit: str = "arcsec",
    rv_time_system: str = "jd2400000",
) -> LegacyTargetData:
    """Parse the project's legacy visual + SB2 CSV.

    ``separation_unit`` applies to the separation and its uncertainty.
    ``rv_time_system`` is ``"jd2400000"`` (the file convention) or
    ``"decimalyear"``; it is declared rather than guessed because an
    unconverted velocity epoch would sit millennia from the visual orbit.
    """
    if separation_unit not in ("arcsec", "mas"):
        raise ValueError("separation_unit must be 'arcsec' or 'mas'")
    if rv_time_system not in ("jd2400000", "decimalyear"):
        raise ValueError("rv_time_system must be 'jd2400000' or 'decimalyear'")
    scale = 1000.0 if separation_unit == "arcsec" else 1.0

    header: dict[str, str] = {}
    rv1: list[list[float]] = []
    rv2: list[list[float]] = []
    visual: list[list[float]] = []
    block = None

    for raw in Path(path).read_text().splitlines():
        if _is_comment(raw):
            marker = raw.upper()
            if "RV1" in marker:
                block = "rv1"
            elif "RV2" in marker:
                block = "rv2"
            elif "VISUAL" in marker:
                block = "visual"
            continue
        fields = _fields(raw)
        if len(fields) < 2:
            continue
        # Header rows are "key,value" with a non-numeric key.
        if block is None:
            try:
                float(fields[0])
            except ValueError:
                header[fields[0]] = fields[1]
                continue
        try:
            numbers = [float(f) for f in fields if _NUMBER.fullmatch(f)]
        except ValueError:  # pragma: no cover - guarded by fullmatch
            continue
        if block == "visual":
            if len(numbers) < 4:
                raise ValueError(f"visual row needs 4 numeric columns: {raw!r}")
            visual.append(numbers[:4])
        elif block == "rv1":
            if len(numbers) < 3:
                raise ValueError(f"RV1 row needs 3 numeric columns: {raw!r}")
            rv1.append(numbers[:3])
        elif block == "rv2":
            if len(numbers) < 3:
                raise ValueError(f"RV2 row needs 3 numeric columns: {raw!r}")
            rv2.append(numbers[:3])

    if not visual:
        raise ValueError(f"no relative astrometry found in {path}")
    if not rv1 or not rv2:
        raise ValueError(f"no radial velocities found in {path}")
    if len(rv1) != len(rv2):
        raise ValueError(
            f"RV1 and RV2 blocks must be paired: got {len(rv1)} and {len(rv2)} rows"
        )

    a = np.array(visual, dtype=float)
    v1 = np.array(rv1, dtype=float)
    v2 = np.array(rv2, dtype=float)
    if not np.allclose(v1[:, 0], v2[:, 0], rtol=0, atol=1e-6):
        raise ValueError("RV1 and RV2 epochs differ; the two curves must be paired")

    rv_epoch = (
        jd2400000_to_decimalyear(v1[:, 0]) if rv_time_system == "jd2400000" else v1[:, 0]
    )

    # Guard against an unconverted or wrongly declared time system: an orbit of
    # a few years cannot be constrained by blocks that do not overlap at all.
    visual_epoch = a[:, 0]
    gap = max(visual_epoch.min(), rv_epoch.min()) - min(visual_epoch.max(), rv_epoch.max())
    if gap > 200.0:
        raise ValueError(
            "visual and radial-velocity epochs do not overlap "
            f"(visual {visual_epoch.min():.1f}-{visual_epoch.max():.1f}, "
            f"RV {rv_epoch.min():.1f}-{rv_epoch.max():.1f}); "
            "check rv_time_system"
        )

    return LegacyTargetData(
        header=header,
        astrometry_epoch_yr=visual_epoch,
        position_angle_deg=a[:, 1],
        separation_mas=a[:, 2] * scale,
        astrometry_sigma_mas=a[:, 3] * scale,
        rv_epoch_yr=rv_epoch,
        rv_primary_kms=v1[:, 1],
        rv_secondary_kms=v2[:, 1],
        rv_primary_sigma_kms=v1[:, 2],
        rv_secondary_sigma_kms=v2[:, 2],
    )


def initial_guess_from_header(data: LegacyTargetData) -> BinaryParams:
    """Starting orbit built from the file's own header elements.

    The header supplies ``P``, ``T``, ``e``, ``a`` (arcsec), ``W`` (node),
    ``w``, ``i``, ``K1``, ``K2``, ``V0`` and ``par``. Total mass follows from
    Kepler's third law using the header semi-major axis and parallax, and the
    mass ratio from ``M2/M1 = K1/K2``. These values initialise only.
    """
    def need(key):
        value = data.header_float(key)
        if value is None:
            raise ValueError(f"legacy header is missing '{key}'")
        return value

    period = need("P")
    parallax = need("par")
    a_arcsec = need("a")
    k1, k2 = need("K1"), need("K2")

    a_rel_au = (a_arcsec * 1000.0) / parallax
    total_mass = a_rel_au**3 / period**2
    q = k1 / k2                      # M2/M1, so Va is the less massive component
    m1 = total_mass / (1.0 + q)

    epoch0 = float(np.min(data.astrometry_epoch_yr))
    return BinaryParams(
        period_yr=period,
        t_peri_yr=need("T") - epoch0,
        eccentricity=need("e"),
        inclination_deg=need("i"),
        # The file's w is the relative-orbit argument; the model stores the
        # primary's, which differs by 180 degrees.
        omega_deg=(need("w") - 180.0) % 360.0,
        node_deg=need("W") % 360.0,
        m1_msun=m1,
        m2_msun=q * m1,
        parallax_mas=parallax,
        gamma_kms=need("V0"),
        beta_g=0.0,
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
