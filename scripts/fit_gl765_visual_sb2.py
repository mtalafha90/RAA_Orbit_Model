#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from raa_orbit_model.kepler import radial_velocities_kms, relative_astrometry_mas, rv_semiamplitudes_kms
from raa_orbit_model.real_data import fit_visual_sb2, parse_legacy_binary_csv


def _write_summary(path: Path, fit, label: str) -> None:
    p = fit.params
    k1, k2 = rv_semiamplitudes_kms(p)
    rows = [
        ("fit", label, ""), ("n_constraints", 110, "scalar"), ("n_free", len(fit.free_names), "parameters"),
        ("dof", fit.dof, ""), ("chi2", fit.chi2, ""), ("reduced_chi2", fit.reduced_chi2, ""),
        ("chi2_astrometry", fit.chi2_astrometry, ""), ("chi2_rv1", fit.chi2_rv1, ""), ("chi2_rv2", fit.chi2_rv2, ""),
        ("period_yr", p.period_yr, fit.uncertainties.get("period_yr", 0.0)),
        ("t_peri_yr", p.t_peri_yr, fit.uncertainties.get("t_peri_yr", 0.0)),
        ("eccentricity", p.eccentricity, fit.uncertainties.get("eccentricity", 0.0)),
        ("inclination_deg", p.inclination_deg, fit.uncertainties.get("inclination_deg", 0.0)),
        ("omega1_deg", p.omega_deg, fit.uncertainties.get("omega_deg", 0.0)),
        ("omega_relative_deg", p.omega_relative_deg, fit.uncertainties.get("omega_deg", 0.0)),
        ("node_deg", p.node_deg, fit.uncertainties.get("node_deg", 0.0)),
        ("m1_msun", p.m1_msun, fit.uncertainties.get("m1_msun", 0.0)),
        ("m2_msun", p.m2_msun, fit.uncertainties.get("m2_msun", 0.0)),
        ("m_total_msun", p.m1_msun + p.m2_msun, "derived"),
        ("parallax_mas", p.parallax_mas, fit.uncertainties.get("parallax_mas", 0.0)),
        ("gamma_kms", p.gamma_kms, fit.uncertainties.get("gamma_kms", 0.0)),
        ("k1_kms", k1, "derived"), ("k2_kms", k2, "derived"),
        ("a_rel_au", p.a_rel_au, "derived"), ("a_rel_ang_mas", p.a_rel_au * p.parallax_mas, "derived"),
    ]
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh); w.writerow(["quantity", "value", "uncertainty_or_note"]); w.writerows(rows)


def _plots(outdir: Path, data, fit) -> None:
    import matplotlib.pyplot as plt
    p = fit.params
    t = np.linspace(min(data.visual_time_yr.min(), 1970), max(data.visual_time_yr.max(), 1996), 1500)
    model = relative_astrometry_mas(t, p)
    fig, ax = plt.subplots(figsize=(5.6, 5.2))
    ax.plot(model[:, 0], model[:, 1], lw=1.6, label="best-fit relative orbit")
    ax.errorbar(data.visual_east_mas, data.visual_north_mas, xerr=data.visual_sigma_mas, yerr=data.visual_sigma_mas,
                fmt="o", ms=4, capsize=2, label="legacy visual/speckle")
    ax.scatter([0], [0], marker="+", s=80, label="primary")
    ax.set_xlabel(r"$\Delta\alpha^*$ (mas)"); ax.set_ylabel(r"$\Delta\delta$ (mas)")
    ax.set_aspect("equal", adjustable="datalim"); ax.legend(fontsize=8); fig.tight_layout()
    fig.savefig(outdir / "05_gl765_visual_orbit.svg"); plt.close(fig)

    tr = np.linspace(min(data.rv1_time_yr.min(), data.rv2_time_yr.min()), max(data.rv1_time_yr.max(), data.rv2_time_yr.max()), 1500)
    rv = radial_velocities_kms(tr, p)
    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    ax.plot(tr, rv[:, 0], lw=1.5, label="primary model"); ax.plot(tr, rv[:, 1], lw=1.5, label="secondary model")
    ax.errorbar(data.rv1_time_yr, data.rv1_kms, yerr=data.rv1_sigma_kms, fmt="o", ms=3, capsize=1.5, label="Va")
    ax.errorbar(data.rv2_time_yr, data.rv2_kms, yerr=data.rv2_sigma_kms, fmt="o", ms=3, capsize=1.5, label="Vb")
    ax.set_xlabel("Decimal year"); ax.set_ylabel(r"Radial velocity (km s$^{-1}$)")
    ax.legend(fontsize=8, ncol=2); fig.tight_layout(); fig.savefig(outdir / "06_gl765_sb2_rv.svg"); plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description="Reproduce the GJ 765.2 V6a visual+SB2 validation fit")
    ap.add_argument("--input", type=Path, default=Path("data/real/gj7652/GL765_Test1.csv"))
    ap.add_argument("--outdir", type=Path, default=Path("results/real/gj7652"))
    ap.add_argument("--figures-dir", type=Path, default=Path("figures"))
    args = ap.parse_args()
    data = parse_legacy_binary_csv(args.input)
    fit = fit_visual_sb2(data); fit_header = fit_visual_sb2(data, fixed_parallax_mas=54.27); fit_balega = fit_visual_sb2(data, fixed_parallax_mas=31.0)
    args.outdir.mkdir(parents=True, exist_ok=True); args.figures_dir.mkdir(parents=True, exist_ok=True)
    _write_summary(args.outdir / "gj765_v6a_fit_summary.csv", fit, "free_parallax")
    _write_summary(args.outdir / "gj765_fixed_parallax_54p27.csv", fit_header, "fixed_parallax_54.27")
    _write_summary(args.outdir / "gj765_fixed_parallax_31p0.csv", fit_balega, "fixed_parallax_31.0")
    _plots(args.figures_dir, data, fit)
    print(f"{data.object_name}: chi2={fit.chi2:.12f}, dof={fit.dof}, redchi2={fit.reduced_chi2:.12f}")
    print(f"Mtot={fit.params.m1_msun + fit.params.m2_msun:.9f} Msun, parallax={fit.params.parallax_mas:.9f} mas")
    print(f"fixed 54.27 mas: chi2={fit_header.chi2:.6f}, redchi2={fit_header.reduced_chi2:.6f}")
    print(f"fixed 31.0 mas: chi2={fit_balega.chi2:.6f}, redchi2={fit_balega.reduced_chi2:.6f}")


if __name__ == "__main__": main()
