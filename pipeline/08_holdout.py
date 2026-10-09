"""Phase 3c: leave-one-season-out prediction test.

For each observing season of a candidate:
  1. hide it, and run the full blind phase 2 search (same cleaning, whole period grid) on the rest;
  2. check the strongest period found is still the candidate's period;
  3. use that blind period and epoch to predict when transits occur in the hidden season;
  4. in the hidden season, slide a transit-shaped box across one full orbit and measure the dip at every
     offset. A real planet gives the deepest dip at offset 0 (the predicted time); noise peaks anywhere.

Two predictions are reported per season: `as_designed` uses the search's own box width, and
`refined_width` first refits width and epoch on the training seasons with a fine grid of widths (added
after the first run, where one fold's search picked a box twice the dip's width and washed it out).

Everything uses the phase 2 search cleaning (no transit mask), so no knowledge of the full-data
ephemeris leaks into the hidden season.

Writes results/phase3/<TIC>_holdout.json and results/phase3/plots/<TIC>_holdout.png.
Usage: python pipeline/08_holdout.py [WORKERS] [TIC ...]
"""
import json
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
import vet as V
from env import RESULTS

SCAN_STEP_DUR = 0.25  # phase-scan step in units of the transit duration
FINE_DURATIONS_HR = np.geomspace(0.4, 4.0, 25)


def rho_of(tic):
    tg = pd.read_csv(RESULTS / "targets_ranked.csv").set_index("TIC").loc[tic]
    return float(np.clip(tg["Mass"] / tg["Rad"] ** 3, 0.2, 50))


def predict(tst, P, t0, dur):
    """Phase scan of the hidden season at a predicted ephemeris."""
    step = SCAN_STEP_DUR * dur
    n = int(0.5 * P / step)
    offsets = np.arange(-n, n + 1) * step  # symmetric, so offset 0 is on the grid
    scan = np.array([V.box_depth(tst, np.abs(V.phase_of(tst["time"], P, t0 + o)) < dur / 2) for o in offsets])
    d0, sig = V.depth_emp(tst, P, t0, dur)
    # fraction of the rest of the orbit (more than one duration from the prediction) with a dip at
    # least as deep as at the predicted time: a chance-coincidence probability for this season
    away = np.abs(offsets) > dur
    frac = float(np.mean(scan[away] >= scan[n]))
    n_tr = int(len(np.unique(np.round((tst["time"][np.abs(V.phase_of(tst["time"], P, t0)) < dur / 2] - t0) / P))))
    return {"pred_t0": t0, "pred_dur_hr": dur * 24, "depth_ppm": d0 * 1e6, "err_emp_ppm": sig * 1e6,
            "snr": d0 / sig if sig > 0 else np.nan, "n_transits": n_tr, "frac_orbit_as_deep": frac,
            "scan_peak_offset_hr": float(offsets[np.nanargmax(scan)] * 24),
            "_scan": {"offset_hr": (offsets * 24).tolist(), "depth_ppm": (scan * 1e6).tolist()}}


def fold(args):
    tic, hold = args
    P_full, t0_full, _ = V.ephemeris(tic)
    lc, _, _ = L.reject_bad_sectors(L.load_star(tic))
    c = L.bin_lc(L.clean(lc))
    test = np.isin(c["sector"], hold)
    train = {k: v[~test] for k, v in c.items()}
    tst = {k: v[test] for k, v in c.items()}
    # same iterative search as phase 2 (up to 3 signals); predict with the one matching the candidate,
    # else with the strongest (which then tests nothing about the candidate)
    signals, (periods, sde) = L.bls_search(train, rho_of(tic), max_signals=V.CANDIDATES.get(tic, 1))
    hits = [x for x in signals if abs(x["period"] / P_full - 1) < 1e-3]
    s = hits[0] if hits else signals[0]
    P, t0, dur = s["period"], s["t0"], s["duration_hr"] / 24
    match = bool(hits)
    # SDE of the candidate period in the training periodogram, whatever won
    near = np.abs(periods / P_full - 1) < 2e-4
    sde_at_true = float(np.nanmax(sde[near])) if near.any() else np.nan
    # Refinement (added after the first run, still blind): the search tries only 7 box widths,
    # ~1.45x apart, so refit width and epoch on the TRAINING data at the training period.
    bls = BoxLeastSquares(train["time"], train["flux"], train["flux_err"])
    r = bls.power([P], FINE_DURATIONS_HR / 24, objective="snr")
    t0_ref, dur_ref = float(r.transit_time[0]), float(r.duration[0])
    # how far the blind prediction drifts from the full-data ephemeris inside the hidden season
    t_mid = np.median(tst["time"])
    pred_mid = t0 + np.round((t_mid - t0) / P) * P
    full_mid = t0_full + np.round((t_mid - t0_full) / P_full) * P_full
    return {"tic": tic, "held_out": [int(x) for x in hold],
            "held_out_dates": V.btjd_to_date(tst["time"].min()) + " to " + V.btjd_to_date(tst["time"].max()),
            "train_best_period": P, "train_best_sde": s["sde"], "train_best_snr": s["snr"],
            "train_signal_number": int(s["signal"]),
            "recovers_period": bool(match), "train_sde_at_candidate_period": sde_at_true,
            "prediction_drift_vs_full_hr": float((pred_mid - full_mid) * 24),
            "as_designed": predict(tst, P, t0, dur), "refined_width": predict(tst, P, t0_ref, dur_ref),
            "_pg": {"period": periods[::20].tolist(), "sde": np.nan_to_num(sde[::20], nan=0).tolist()}}


def plot(tic, folds):
    P_full = V.ephemeris(tic)[0]
    fig, axes = plt.subplots(len(folds), 2, figsize=(14, 2.8 * len(folds)), squeeze=False)
    for (a1, a2), f in zip(axes, folds):
        a1.plot(f["_pg"]["period"], f["_pg"]["sde"], lw=0.5, color="k")
        a1.axvline(P_full, color="C3", alpha=0.4, lw=3)
        a1.set_xscale("log")
        a1.set_ylabel("SDE")
        a1.set_title(f"blind search WITHOUT {f['held_out_dates']}: best P = {f['train_best_period']:.4f} d "
                     f"(SDE {f['train_best_sde']:.1f}) {'= candidate' if f['recovers_period'] else '≠ candidate'}",
                     fontsize=9)
        for key, col, lab in (("as_designed", "0.6", "search box"), ("refined_width", "C0", "refined box")):
            p = f[key]
            a2.plot(p["_scan"]["offset_hr"], p["_scan"]["depth_ppm"], "-", color=col, lw=1,
                    label=f"{lab} {p['pred_dur_hr']:.2f} h: {p['depth_ppm']:.0f}±{p['err_emp_ppm']:.0f} ppm")
        p = f["refined_width"]
        a2.axhspan(-p["err_emp_ppm"], p["err_emp_ppm"], color="0.9")
        a2.axvline(0, color="C3", ls="--", lw=1)
        a2.axhline(0, color="k", lw=0.5)
        a2.set_ylabel("dip depth (ppm)")
        a2.legend(fontsize=7, loc="upper right")
        a2.set_title(f"hidden season, refined box: {p['snr']:.1f}σ at the predicted time (red), "
                     f"{p['n_transits']} transits; as deep elsewhere in orbit: {100 * p['frac_orbit_as_deep']:.1f}%",
                     fontsize=9)
    axes[-1, 0].set_xlabel("period (d)")
    axes[-1, 1].set_xlabel("hours from predicted transit time (one full orbit)")
    fig.suptitle(f"TIC {tic}: hide one season, predict it from the others", fontsize=11)
    fig.tight_layout()
    fig.savefig(V.OUT / "plots" / f"{tic}_holdout.png", dpi=80)
    plt.close(fig)


def strip(d):
    return {k: (strip(v) if isinstance(v, dict) else v) for k, v in d.items() if not k.startswith("_")}


def main():
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    (V.OUT / "plots").mkdir(parents=True, exist_ok=True)
    tics = [int(a) for a in sys.argv[2:]] or list(V.CANDIDATES)
    jobs = []
    for tic in tics:
        lc, _, _ = L.reject_bad_sectors(L.load_star(tic))
        jobs += [(tic, g) for g in V.seasons(lc)]
    with Pool(workers) as pool:
        res = pool.map(fold, jobs)
    for tic in tics:
        folds = [r for r in res if r["tic"] == tic]
        plot(tic, folds)
        (V.OUT / f"{tic}_holdout.json").write_text(json.dumps([strip(f) for f in folds], indent=1))
        print(f"TIC {tic}")
        for f in folds:
            print(f"  hide {f['held_out_dates']}: blind best P {f['train_best_period']:.5f} (SDE {f['train_best_sde']:.1f}, "
                  f"{'match' if f['recovers_period'] else 'NO MATCH, SDE at candidate %.1f' % f['train_sde_at_candidate_period']})"
                  f", drift {f['prediction_drift_vs_full_hr']:+.2f} h")
            for key in ("as_designed", "refined_width"):
                p = f[key]
                print(f"      {key:13s} box {p['pred_dur_hr']:.2f} h: {p['depth_ppm']:.0f}±{p['err_emp_ppm']:.0f} ppm "
                      f"({p['snr']:.1f}σ, {p['n_transits']} tr), as deep elsewhere {100 * p['frac_orbit_as_deep']:.1f}%, "
                      f"scan peak {p['scan_peak_offset_hr']:+.1f} h")


if __name__ == "__main__":
    main()
