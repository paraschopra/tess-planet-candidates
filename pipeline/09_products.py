"""Phase 3d: does the dip appear in an independently processed light curve?

The search used NASA's TESS-SPOC light curves. MIT's Quick-Look Pipeline (QLP) extracts its own light
curves from the same full-frame images with different apertures, background subtraction and
systematics correction. A processing artifact in one pipeline should not reappear in the other.

Both products get identical cleaning (transit-masked detrending) and the depth is measured at the fixed
candidate ephemeris, per sector and in total, with empirical errors.

Caveat: TESS-SPOC PDCSAP removes neighbours' light (CROWDSAP); QLP's KSPSAP flux may handle dilution
differently, so a modest depth offset between products is not by itself a problem.

Writes results/phase3/<TIC>_products.json and results/phase3/plots/<TIC>_products.png.
Usage: python pipeline/09_products.py [TIC ...]
"""
import importlib
import json
import sys

import env  # noqa: F401

import lightkurve as lk
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2 as chi2_dist

import lcsearch as L
import vet as V

dl = importlib.import_module("02_download")


def qlp_lightcurve(tic):
    sr = lk.search_lightcurve(f"TIC {tic}", author="QLP")
    parts = {"time": [], "flux": [], "flux_err": [], "sector": [], "exptime": []}
    used = []
    for i in range(len(sr)):
        row = sr[i]
        sec = int(str(row.table["mission"][0]).split()[-1])
        lc = dl.download_product(row)
        t, f, e, col = dl.flux_arrays(lc, "QLP")
        parts["time"].append(t)
        parts["flux"].append(f)
        parts["flux_err"].append(e)
        parts["sector"].append(np.full(len(t), sec))
        parts["exptime"].append(np.full(len(t), float(row.table["exptime"][0])))
        used.append(f"{sec}:{col}")
    a = {k: np.concatenate(v) for k, v in parts.items()}
    o = np.argsort(a["time"])
    return {k: v[o] for k, v in a.items()}, used


def measure(c, raw, P, t0, dur):
    """Per-sector depths (quick empirical errors) and the total with matched-processing errors."""
    rows = {}
    for s in np.unique(c["sector"]):
        m = c["sector"] == s
        if np.sum(m & (np.abs(V.phase_of(c["time"], P, t0)) < dur / 2)) < 3:
            continue
        d, e = V.depth_emp(c, P, t0, dur, sel=m)
        rows[int(s)] = (d * 1e6, e * 1e6)
    # total: every control box processed like the transit box (as in 07_fit; the Codex review found the
    # products comparison still used the quicker, slightly optimistic depth_emp)
    d, e, _ = V.depth_emp_matched(raw, P, t0, dur)
    return rows, (d * 1e6, e * 1e6)


def run(tic):
    P, t0, dur = V.ephemeris(tic)
    spoc_raw, _, _ = L.reject_bad_sectors(L.load_star(tic))
    qlp_raw, used = qlp_lightcurve(tic)
    qlp_raw, bad_q, _ = L.reject_bad_sectors(qlp_raw)
    common = np.intersect1d(np.unique(spoc_raw["sector"]), np.unique(qlp_raw["sector"]))
    res = {"tic": tic, "qlp_products": used, "qlp_rejected_sectors": bad_q, "common_sectors": common.tolist()}
    cleaned = {}
    for name, raw in (("TESS-SPOC", spoc_raw), ("QLP", qlp_raw)):
        keep = np.isin(raw["sector"], common)
        raw = {k: (v[keep] if isinstance(v, np.ndarray) and v.shape == raw["time"].shape else v)
               for k, v in raw.items()}
        c = V.masked_clean(raw, P, t0, dur)
        cleaned[name] = c
        per, tot = measure(c, raw, P, t0, dur)
        d = np.array([v[0] for v in per.values()])
        e = np.array([v[1] for v in per.values()])
        mean = np.sum(d / e ** 2) / np.sum(1 / e ** 2)
        chi2 = float(np.sum(((d - mean) / e) ** 2))
        res[name] = {"total_depth_ppm": tot[0], "total_err_ppm": tot[1], "snr": tot[0] / tot[1],
                     "per_sector": {str(k): v for k, v in per.items()},
                     "sector_chi2": chi2, "sector_consistency_p": float(chi2_dist.sf(chi2, len(d) - 1))}
    (V.OUT / f"{tic}_products.json").write_text(json.dumps(res, indent=1))
    plot(tic, P, t0, dur, cleaned, res)
    return res


def plot(tic, P, t0, dur, cleaned, res):
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.2))
    for ax, name, col in ((axes[0], "TESS-SPOC", "C0"), (axes[1], "QLP", "C2")):
        c = cleaned[name]
        ph = V.phase_of(c["time"], P, t0) * 24
        sel = np.abs(ph) < 4 * dur * 24
        edges = np.linspace(-4 * dur * 24, 4 * dur * 24, 33)
        idx = np.digitize(ph[sel], edges)
        y = (c["flux"][sel] - 1) * 1e6
        w = 1 / c["flux_err"][sel] ** 2
        xb = 0.5 * (edges[1:] + edges[:-1])
        yb = np.array([np.sum(y[idx == i] * w[idx == i]) / np.sum(w[idx == i]) for i in range(1, len(edges))])
        eb = np.array([1 / np.sqrt(np.sum(w[idx == i])) for i in range(1, len(edges))])
        ax.errorbar(xb, yb, eb, fmt="o", ms=3, color=col)
        ax.axvspan(-dur * 12, dur * 12, color="C1", alpha=0.12)
        ax.axhline(0, color="k", lw=0.5)
        r = res[name]
        ax.set_title(f"{name}: {r['total_depth_ppm']:.0f} ± {r['total_err_ppm']:.0f} ppm ({r['snr']:.1f}σ)", fontsize=10)
        ax.set_xlabel("hours from mid-transit")
        ax.set_ylabel("ppm")
    ax = axes[2]
    secs = res["common_sectors"]
    for j, (name, col) in enumerate((("TESS-SPOC", "C0"), ("QLP", "C2"))):
        ps = res[name]["per_sector"]
        x = [i + (j - 0.5) * 0.2 for i, s in enumerate(secs) if str(s) in ps]
        y = [ps[str(s)][0] for s in secs if str(s) in ps]
        e = [ps[str(s)][1] for s in secs if str(s) in ps]
        ax.errorbar(x, y, e, fmt="o", color=col, label=name)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xticks(range(len(secs)), [f"S{s}" for s in secs])
    ax.set_ylabel("depth (ppm)")
    ax.legend(fontsize=8)
    ax.set_title("per-sector depth, both pipelines", fontsize=10)
    fig.suptitle(f"TIC {tic}: same dip in two independently processed light curves?", fontsize=11)
    fig.tight_layout()
    fig.savefig(V.OUT / "plots" / f"{tic}_products.png", dpi=80)
    plt.close(fig)


def main():
    (V.OUT / "plots").mkdir(parents=True, exist_ok=True)
    for tic in [int(a) for a in sys.argv[1:]] or list(V.CANDIDATES):
        r = run(tic)
        print(f"TIC {tic}  QLP products {r['qlp_products']}  rejected {r['qlp_rejected_sectors']}")
        for name in ("TESS-SPOC", "QLP"):
            x = r[name]
            per = "  ".join(f"S{s}:{d:.0f}±{e:.0f}" for s, (d, e) in x["per_sector"].items())
            print(f"  {name:9s} total {x['total_depth_ppm']:.0f} ± {x['total_err_ppm']:.0f} ppm ({x['snr']:.1f}σ), "
                  f"sectors agree p = {x['sector_consistency_p']:.3f}   {per}")


if __name__ == "__main__":
    main()
