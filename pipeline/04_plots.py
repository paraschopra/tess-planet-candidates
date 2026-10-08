"""Phase 2c: one diagnostic page per signal that passed the search thresholds (run after 05_triage).

Panels: full light curve with predicted transits, periodogram, folded transit, odd vs even transits,
phase 0.5 (secondary eclipse), and depth measured separately in each sector.

Usage: python pipeline/04_plots.py [MAX_PLOTS]
"""
import sys
from multiprocessing import Pool

import env  # noqa: F401

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.timeseries import BoxLeastSquares

import lcsearch as L
from env import RESULTS

OUT = RESULTS / "phase2"
PLOTS = OUT / "plots"


def binned(x, y, lo, hi, n):
    edges = np.linspace(lo, hi, n + 1)
    idx = np.digitize(x, edges) - 1
    ok = (idx >= 0) & (idx < n)
    c = np.bincount(idx[ok], minlength=n)
    m = np.bincount(idx[ok], y[ok], minlength=n) / np.maximum(c, 1)
    s = np.sqrt(np.bincount(idx[ok], (y[ok] - m[idx[ok]]) ** 2, minlength=n) / np.maximum(c, 1))
    centers = 0.5 * (edges[1:] + edges[:-1])
    good = c >= 3
    return centers[good], m[good], (s / np.sqrt(np.maximum(c, 1)))[good]


def fold_panel(ax, phase_hr, flux, dur_hr, title, depth_ppm, color="C0"):
    w = 3 * dur_hr
    sel = np.abs(phase_hr) < w
    ax.plot(phase_hr[sel], (flux[sel] - 1) * 1e6, ".", color="0.75", ms=1.5, zorder=0)
    x, m, e = binned(phase_hr[sel], (flux[sel] - 1) * 1e6, -w, w, 40)
    ax.errorbar(x, m, e, fmt="o", color=color, ms=3, zorder=2)
    ax.axvspan(-dur_hr / 2, dur_hr / 2, color="C1", alpha=0.12)
    ax.axhline(-depth_ppm, color="C1", lw=0.8, ls="--")
    lim = max(3 * depth_ppm, 3 * np.nanstd(m) if len(m) else 100)
    ax.set_ylim(-lim, lim * 0.7)
    ax.set_xlabel("hours from mid-transit")
    ax.set_ylabel("ppm")
    ax.set_title(title, fontsize=9)


def L_rho(row):
    return float(row["rho"]) if "rho" in row and np.isfinite(row["rho"]) else 1.0


def plot_signal(row):
    tic, k = int(row["tic"]), int(row["signal"])
    lc, bad, _ = L.reject_bad_sectors(L.load_star(tic))
    c = L.bin_lc(L.clean(lc))
    P, t0, dur = row["period"], row["t0"], row["duration_hr"] / 24
    t, f = c["time"], c["flux"]
    phase = ((t - t0 + 0.5 * P) % P) - 0.5 * P
    epoch = np.round((t - t0) / P).astype(int)

    fig = plt.figure(figsize=(15, 11))
    gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 1])
    flags = row["flags"] if isinstance(row["flags"], str) else ""
    fig.suptitle(f"TIC {tic}  signal {k}:  P = {P:.5f} d   depth = {row['depth_ppm']:.0f} ppm   "
                 f"dur = {row['duration_hr']:.2f} h   Rp ≈ {row['rp_earth']:.2f} R⊕   SNR = {row['snr']:.1f}   "
                 f"SDE = {row['sde']:.1f}   N_tr = {row['n_transits']}\n"
                 f"Tmag {row['tmag']:.2f}  Teff {row['teff']:.0f} K  R* {row['rstar']:.2f}  d {row['dist_pc']:.1f} pc   "
                 f"sectors used {row['n_sectors_used']}, rejected {row['sectors_rejected']}   flags: {flags or 'none'}",
                 fontsize=10)

    # 1. full light curve, sectors compressed side by side
    ax = fig.add_subplot(gs[0, :])
    x0 = 0
    for s in np.unique(c["sector"]):
        m = c["sector"] == s
        ts = t[m]
        ax.plot(ts - ts.min() + x0, (f[m] - 1) * 1e6, ",", color="k" if int(s) % 2 else "0.4")
        intr = np.abs(phase[m]) < dur / 2
        ax.plot(ts[intr] - ts.min() + x0, (f[m][intr] - 1) * 1e6, ".", color="C3", ms=2)
        ax.text(x0 + 1, 0, f"S{int(s)}", fontsize=7, va="bottom", transform=ax.get_xaxis_transform())
        x0 += np.ptp(ts) + 2
    lim = 5 * L.robust_std(f) * 1e6
    ax.set_ylim(-lim, lim)
    ax.set_ylabel("ppm")
    ax.set_title("cleaned light curve, sectors placed side by side (red = predicted in-transit points)", fontsize=9)

    # 2. periodogram
    ax = fig.add_subplot(gs[1, 0])
    pg = np.load(OUT / "periodograms" / f"{tic}.npz")
    pg_p, pg_s = pg["period"], pg["sde"]
    if not np.isfinite(pg_s).any():  # saved by an early buggy run: recompute (signal-1 search only)
        _, (pg_p, pg_s) = L.bls_search(c, L_rho(row), max_signals=1)
    ax.plot(pg_p, pg_s, lw=0.6, color="k")
    ax.axvline(P, color="C3", alpha=0.4, lw=3)
    ax.set_xscale("log")
    ax.set_xlabel("period (d)")
    ax.set_ylabel("SDE")
    ax.set_title("periodogram (signal 1 search)", fontsize=9)

    # 3. folded transit
    fold_panel(fig.add_subplot(gs[1, 1]), phase * 24, f, dur * 24, "all transits folded", row["depth_ppm"])

    # 4. secondary eclipse at phase 0.5
    ph2 = ((t - t0) % P) - 0.5 * P
    fold_panel(fig.add_subplot(gs[1, 2]), ph2 * 24, f, dur * 24,
               f"phase 0.5 (secondary eclipse?) BLS depth {row['depth_half_ppm']:.0f} ppm", row["depth_ppm"],
               color="C2")

    # 5. odd vs even
    ax = fig.add_subplot(gs[2, 0])
    for parity, col in ((0, "C0"), (1, "C3")):
        m = (epoch % 2) == parity
        x, mm, e = binned(phase[m] * 24, (f[m] - 1) * 1e6, -3 * dur * 24, 3 * dur * 24, 30)
        ax.errorbar(x, mm, e, fmt="o", ms=3, color=col, label="even" if parity == 0 else "odd")
    ax.axvspan(-dur * 12, dur * 12, color="C1", alpha=0.12)
    ax.legend(fontsize=8)
    ax.set_xlabel("hours from mid-transit")
    ax.set_title(f"odd {row['depth_odd_ppm']:.0f}±{row['depth_odd_err_ppm']:.0f} vs even "
                 f"{row['depth_even_ppm']:.0f}±{row['depth_even_err_ppm']:.0f} ppm", fontsize=9)

    # 6. per-sector depth
    ax = fig.add_subplot(gs[2, 1:])
    secs, d, de = [], [], []
    for s in np.unique(c["sector"]):
        m = c["sector"] == s
        if np.sum(np.abs(phase[m]) < dur / 2) < 2:
            continue
        st = BoxLeastSquares(t[m], f[m], c["flux_err"][m]).compute_stats(P, dur, t0)
        secs.append(int(s))
        d.append(st["depth"][0] * 1e6)
        de.append(st["depth"][1] * 1e6)
    ax.errorbar(range(len(secs)), d, de, fmt="o")
    ax.axhline(row["depth_ppm"], color="C1", ls="--", lw=0.8)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xticks(range(len(secs)), [f"S{s}" for s in secs], fontsize=7, rotation=90)
    ax.set_ylabel("depth (ppm)")
    ax.set_title("depth measured in each sector separately (should agree if real)", fontsize=9)

    fig.tight_layout()
    out = PLOTS / f"{tic}_s{k}.png"
    fig.savefig(out, dpi=80)
    plt.close(fig)
    return str(out)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    PLOTS.mkdir(parents=True, exist_ok=True)
    sig = pd.read_csv(OUT / "triage.csv")  # sorted: unflagged signals first
    sig = sig[sig["passes"]].head(n)
    with Pool(12) as pool:
        out = pool.map(plot_signal, [r for _, r in sig.iterrows()])
    print(f"wrote {len(out)} plots to {PLOTS}")


if __name__ == "__main__":
    main()
