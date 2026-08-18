#!/usr/bin/env python
"""Regenerate the real-binary V6a fit and the manuscript table tab:gl765.

The measurements themselves are not in this repository. Supply the legacy
visual + SB2 input file and this reproduces the fit, writing a frozen CSV
beside the synthetic products so the real-data result meets the same
reproducibility standard as everything else.

    python scripts/run_legacy_target_fit.py gl765.dat --preset gj765
    python scripts/run_legacy_target_fit.py --describe-format

Both visual-orbit node branches are fitted and the lower-chi-square solution is
retained, because relative astrometry alone leaves the ascending node ambiguous
by 180 degrees. Header values in the legacy file initialise only; they are
never used as priors. For GJ 765.2 the header parallax (54.27 mas) and the
header coordinates are known to be wrong.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402

from raa_orbit_model.fit import fit_joint, joint_residuals  # noqa: E402
from raa_orbit_model.legacy_target import (  # noqa: E402
    formal_covariance,
    formal_uncertainties,
    initial_guess,
    legacy_joint_data,
    node_branches,
    parse_legacy_file,
    summarise_fit,
    total_mass_uncertainty,
)
from raa_orbit_model.model import GaiaResponseConfig  # noqa: E402

FORMAT_HELP = """\
Expected legacy input file
--------------------------
Whitespace- or comma-separated, '#' begins a comment. Blocks are introduced by
a line containing ASTROMETRY or VELOCITIES.

    # optional header, used only to initialise, never as a prior
    parallax = 54.27

    ASTROMETRY
    # epoch_yr   theta_deg   rho         sigma
    1978.5460    103.4       0.2140      0.010
    ...

    VELOCITIES
    # epoch_yr   rv1_kms   rv2_kms   sigma1   [sigma2]
    1980.2100    -12.44     18.03     0.35
    ...

theta is the position angle measured North through East, converted through
delta_alpha* = rho sin(theta), delta_delta = rho cos(theta). The fourth
astrometric column is an isotropic one-sigma tangent-plane uncertainty.
Separation units are declared with --separation-unit; they are never guessed.

If your file differs, adjust parse_legacy_file in
src/raa_orbit_model/legacy_target.py rather than editing the measurements.
"""

# Balega et al. (2007), A&A 464, 635. Used only to initialise the search and as
# the external comparison printed alongside the fit.
GJ765_PRESET = dict(
    period_yr=11.919,
    eccentricity=0.240,
    inclination_deg=80.2,
    node_deg=293.0,
    omega_relative_deg=250.0,
    m1_msun=0.831,
    m2_msun=0.763,
    parallax_mas=31.0,
)

FREE_NAMES = (
    "period_yr", "t_peri_yr", "eccentricity", "inclination_deg",
    "omega_deg", "node_deg", "m1_msun", "m2_msun", "parallax_mas", "gamma_kms",
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("legacy_file", nargs="?", help="legacy visual + SB2 input file")
    parser.add_argument("--describe-format", action="store_true",
                        help="print the expected input format and exit")
    parser.add_argument("--preset", choices=("gj765",), default="gj765",
                        help="starting orbit for the search")
    parser.add_argument("--separation-unit", choices=("arcsec", "mas"), default="arcsec")
    parser.add_argument("--output", default="results/frozen/legacy_target_fit.csv")
    parser.add_argument("--fixed-parallax-mas", nargs="*", type=float, default=(),
                        help="also refit with the parallax held at each of these values")
    args = parser.parse_args()

    if args.describe_format:
        print(FORMAT_HELP)
        return
    if not args.legacy_file:
        parser.error("provide a legacy input file, or use --describe-format")

    data = parse_legacy_file(args.legacy_file, separation_unit=args.separation_unit)
    print(f"loaded {data.n_astrometry} relative positions and {data.n_rv} paired RV epochs")
    print(f"scalar constraints: {data.n_constraints}")
    if data.header:
        print(f"header (initialisation only, never a prior): {data.header}")

    joint = legacy_joint_data(data)
    response = GaiaResponseConfig("photocentre")   # Gaia channel is empty
    start = initial_guess(data, **GJ765_PRESET)

    best = None
    for label, branch in zip(("node", "node+180"), node_branches(start)):
        result = fit_joint(joint, branch, response, free_names=FREE_NAMES)
        print(f"  branch {label:9s}: chi2={result.chi2:10.3f} success={result.success}")
        if best is None or result.chi2 < best[1].chi2:
            best = (label, result)

    label, result = best
    print(f"\nretained branch: {label}")

    summary = summarise_fit(result, data)
    step = 1e-6
    jacobian = np.column_stack([
        (joint_residuals(_bump(result.params, name, step), joint, response)
         - joint_residuals(result.params, joint, response)) / step
        for name in FREE_NAMES
    ])
    covariance = formal_covariance(jacobian)
    errors = dict(zip(FREE_NAMES, formal_uncertainties(jacobian)))
    for name, value in errors.items():
        summary[f"sigma_{name}"] = float(value)
    # The component masses are strongly correlated, so the total-mass error
    # needs the covariance term rather than a quadrature sum.
    errors["total_mass_msun"] = total_mass_uncertainty(covariance, FREE_NAMES)
    # omega_relative differs from the stored primary omega by a constant.
    errors["omega_relative_deg"] = errors["omega_deg"]
    summary["sigma_total_mass_msun"] = float(errors["total_mass_msun"])

    print(f"\nchi2={summary['chi2']:.3f}  dof={summary['dof']}  "
          f"reduced={summary['reduced_chi2']:.3f}")
    print(f"{'quantity':>22} | {'fit':>12} | {'formal 1-sigma':>14}")
    print("-" * 56)
    for name in ("period_yr", "eccentricity", "inclination_deg", "node_deg",
                 "omega_relative_deg", "total_mass_msun", "parallax_mas"):
        print(f"{name:>22} | {summary[name]:12.4f} | {errors[name]:14.4f}")

    for fixed in args.fixed_parallax_mas:
        from dataclasses import replace
        held = fit_joint(
            joint, replace(result.params, parallax_mas=float(fixed)), response,
            free_names=tuple(n for n in FREE_NAMES if n != "parallax_mas"),
        )
        print(f"  parallax fixed at {fixed:7.2f} mas -> chi2={held.chi2:9.3f} "
              f"reduced={held.reduced_chi2:.3f} (dof {held.dof})")
        summary[f"chi2_parallax_fixed_{fixed:g}"] = float(held.chi2)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=sorted(summary))
        writer.writeheader()
        writer.writerow(summary)
    print(f"\nwrote frozen summary to {output}")


def _bump(params, name, step):
    from dataclasses import replace
    return replace(params, **{name: float(getattr(params, name)) + step})


if __name__ == "__main__":
    main()
