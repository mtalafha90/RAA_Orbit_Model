"""Command-line configuration for the physical and noise axes of an experiment.

`docs/methodology.md` section 12 names seven axes along which the photocentre
versus resolution-aware comparison should be mapped:

    rho, Delta G (or beta_G), q, P, e, scan-angle coverage, S/N.

Until now only three were reachable from the command line. Separation entered
through ``--a-over-sigma-values``, light fraction through ``--beta-values``,
and scan coverage through the choice of sky position, but the orbit itself and
every noise level were fixed in the body of each runner. Mass ratio, period,
eccentricity and signal-to-noise therefore could not be varied without editing
source, which is why the published experiment holds them constant.

This module centralises those axes so both experiment runners expose the same
options and the same defaults. The defaults reproduce the orbit and noise
prescription used for the frozen results, so adding these switches does not
change any existing experiment.
"""

from __future__ import annotations

import argparse

from .kepler import BinaryParams

# The orbit and noise prescription behind every frozen result.
DEFAULT_ORBIT = dict(
    period_yr=2.0,
    t_peri_yr=0.15,
    eccentricity=0.25,
    inclination_deg=72.0,
    omega_deg=55.0,
    node_deg=120.0,
    m1_msun=1.25,
    m2_msun=0.85,
    parallax_mas=20.0,
    gamma_kms=7.0,
    beta_g=0.20,
)
DEFAULT_NOISE = dict(
    n_ast=24,
    n_rv=48,
    ast_sigma_mas=0.20,
    rv_sigma_kms=0.10,
    gaia_sigma_mas=0.10,
)


def add_orbit_arguments(parser: argparse.ArgumentParser) -> None:
    """Expose the orbital axes that were previously hardcoded."""
    group = parser.add_argument_group(
        "orbit",
        "Physical orbit of the injected binary. Defaults reproduce the frozen results.",
    )
    group.add_argument("--period-yr", type=float, default=DEFAULT_ORBIT["period_yr"],
                       help="orbital period; near 1 yr it competes with the parallax signal")
    group.add_argument("--t-peri-yr", type=float, default=DEFAULT_ORBIT["t_peri_yr"])
    group.add_argument("--eccentricity", type=float, default=DEFAULT_ORBIT["eccentricity"])
    group.add_argument("--inclination-deg", type=float, default=DEFAULT_ORBIT["inclination_deg"])
    group.add_argument("--omega-deg", type=float, default=DEFAULT_ORBIT["omega_deg"])
    group.add_argument("--node-deg", type=float, default=DEFAULT_ORBIT["node_deg"])
    group.add_argument("--m1-msun", type=float, default=DEFAULT_ORBIT["m1_msun"])
    group.add_argument("--m2-msun", type=float, default=DEFAULT_ORBIT["m2_msun"],
                       help="ignored when --mass-ratio is given")
    group.add_argument("--mass-ratio", type=float, default=None,
                       help="q = M2/M1; sets M2 from M1, overriding --m2-msun")
    group.add_argument("--gamma-kms", type=float, default=DEFAULT_ORBIT["gamma_kms"])


def add_noise_arguments(parser: argparse.ArgumentParser) -> None:
    """Expose the signal-to-noise axis for all three data channels."""
    group = parser.add_argument_group(
        "noise",
        "Per-channel measurement precision and epoch counts. Defaults reproduce "
        "the frozen results.",
    )
    group.add_argument("--n-ast", type=int, default=DEFAULT_NOISE["n_ast"],
                       help="number of resolved relative-astrometry epochs")
    group.add_argument("--n-rv", type=int, default=DEFAULT_NOISE["n_rv"],
                       help="number of SB2 radial-velocity epochs")
    group.add_argument("--ast-sigma-mas", type=float, default=DEFAULT_NOISE["ast_sigma_mas"])
    group.add_argument("--rv-sigma-kms", type=float, default=DEFAULT_NOISE["rv_sigma_kms"])
    group.add_argument("--gaia-sigma-mas", type=float, default=DEFAULT_NOISE["gaia_sigma_mas"])


def truth_from_args(args, parser: argparse.ArgumentParser | None = None) -> BinaryParams:
    """Build the injected binary from parsed orbit arguments."""

    def fail(message: str):
        if parser is not None:
            parser.error(message)
        raise ValueError(message)

    m1 = float(args.m1_msun)
    if m1 <= 0:
        fail("--m1-msun must be > 0")

    if getattr(args, "mass_ratio", None) is not None:
        q = float(args.mass_ratio)
        if not (0.0 < q <= 1.0):
            fail("--mass-ratio must satisfy 0 < q <= 1 so component 1 is the primary")
        m2 = q * m1
    else:
        m2 = float(args.m2_msun)
        if m2 <= 0:
            fail("--m2-msun must be > 0")
    if m2 > m1:
        fail("component 1 must be the more massive: require M2 <= M1")
    if args.period_yr <= 0:
        fail("--period-yr must be > 0")
    if not (0.0 <= args.eccentricity < 1.0):
        fail("--eccentricity must satisfy 0 <= e < 1")

    params = BinaryParams(
        period_yr=float(args.period_yr),
        t_peri_yr=float(args.t_peri_yr),
        eccentricity=float(args.eccentricity),
        inclination_deg=float(args.inclination_deg),
        omega_deg=float(args.omega_deg),
        node_deg=float(args.node_deg),
        m1_msun=m1,
        m2_msun=m2,
        parallax_mas=float(DEFAULT_ORBIT["parallax_mas"]),
        gamma_kms=float(args.gamma_kms),
        beta_g=float(DEFAULT_ORBIT["beta_g"]),
    )
    params.validate()
    return params


def noise_kwargs_from_args(args, parser: argparse.ArgumentParser | None = None) -> dict:
    """Collect the per-channel noise settings for the experiment runners."""

    def fail(message: str):
        if parser is not None:
            parser.error(message)
        raise ValueError(message)

    for name, value in (
        ("--ast-sigma-mas", args.ast_sigma_mas),
        ("--rv-sigma-kms", args.rv_sigma_kms),
        ("--gaia-sigma-mas", args.gaia_sigma_mas),
    ):
        if value <= 0:
            fail(f"{name} must be > 0")
    for name, value in (("--n-ast", args.n_ast), ("--n-rv", args.n_rv)):
        if value < 0:
            fail(f"{name} must be >= 0")

    return dict(
        n_ast=int(args.n_ast),
        n_rv=int(args.n_rv),
        ast_sigma_mas=float(args.ast_sigma_mas),
        rv_sigma_kms=float(args.rv_sigma_kms),
        gaia_sigma_mas=float(args.gaia_sigma_mas),
    )
