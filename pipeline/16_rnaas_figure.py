"""Figure 1 of the RNAAS note: phase-folded TESS-SPOC light curves of both candidates with transit models.

Data: transit-masked detrended PDCSAP flux (vet.load_vet) from the 200-s sectors (sharpest shape),
15-minute bins with standard errors. Model: best stellar-prior fit from 07_fit (instantaneous).

Writes docs/rnaas/figure1.pdf and .png. Usage: python pipeline/16_rnaas_figure.py
"""
import json

import env  # noqa: F401
import setuptools  # noqa: F401  distutils shim for batman

import batman
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import vet as V
from env import ROOT

OUT = ROOT / "docs" / "rnaas"
PANELS = ((370612437, "TIC 370612437.01"), (252721481, "TIC 252721481.01 (HIP 114417)"))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True)
    for ax, (tic, label) in zip(axes, PANELS):
        fit = json.loads((V.OUT / f"{tic}_fit.json").read_text())
        s, best = fit["stellar"], fit["stellar"]["best"]
        P, t0 = s["P"]["median"], s["t0"]["median"]
        c = V.load_vet(tic)
        m = c["exptime"] < 300
        ph = V.phase_of(c["time"][m], P, t0) * 24
        y = (c["flux"][m] - 1) * 1e6
        ax.plot(ph, y, ".", color="0.85", ms=1, rasterized=True, zorder=0)
        edges = np.arange(-4, 4.01, 0.25)
        idx = np.digitize(ph, edges)
        xb = 0.5 * (edges[1:] + edges[:-1])
        yb = np.array([np.mean(y[idx == i]) for i in range(1, len(edges))])
        eb = np.array([np.std(y[idx == i]) / np.sqrt((idx == i).sum()) for i in range(1, len(edges))])
        ax.errorbar(xb, yb, eb, fmt="o", ms=3, color="k", lw=0.8, zorder=2)
        th = np.linspace(-4, 4, 2000)
        p = batman.TransitParams()
        a = (6.674e-8 * best["rho_cgs"] * (best["P"] * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)
        p.t0, p.per, p.rp, p.a = 0.0, best["P"], best["k"], a
        p.inc = np.degrees(np.arccos(best["b"] / a))
        p.ecc, p.w, p.limb_dark, p.u = 0.0, 90.0, "quadratic", fit["limb_darkening"]
        ax.plot(th, (batman.TransitModel(p, th / 24).light_curve(p) - 1) * 1e6, color="C3", lw=1.5, zorder=3)
        ax.set_xlim(-4, 4)
        ax.set_ylim(-550, 350)
        ax.set_xlabel("Hours from mid-transit")
        ax.set_title(f"{label}\nP = {P:.5f} d, depth {s['depth_ppm']['median']:.0f} ppm, "
                     f"$R_p$ = {s['rp_earth']['median']:.2f} R$_\\oplus$", fontsize=8)
    axes[0].set_ylabel("Relative flux (ppm)")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"figure1.{ext}", dpi=200)
    print("wrote", OUT / "figure1.pdf")


if __name__ == "__main__":
    main()
