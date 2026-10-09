"""Phase 2d: triage search signals with automatic red flags, so humans only look at plausible ones.

Reads the per-star search results (works on partial runs) and, for every signal that passed the
search thresholds, recomputes simple diagnostics on the cleaned light curve:

  variability_repeat   same period (or 2x / 0.5x) found again after masking -> smooth stellar
                       variability, not a transit (a masked transit disappears; a sinusoid doesn't)
  long_duration        dip lasts > 1.6x the longest a planet could on this period around this star
  tess_orbit_alias     period within 3% of a multiple of the 13.7 d spacecraft orbit
  edge_dominated       > 50% of in-transit points within 0.5 d of a data gap (ramps, scattered light)
  one_event_dominated  a single transit contributes > 50% of the total signal (one glitch, not periodic)
  few_events_dominate  dropping the two strongest transits more than halves the SNR
  odd_even_mismatch    alternate dips differ by > 3 sigma (two stars eclipsing)
  secondary_eclipse    a > 3 sigma dip at phase 0.5 (an eclipsing star, not a planet)
  low_snr_red_noise    depth < 6 sigma, where sigma is measured empirically by placing the same box
                       at 40 off-transit phases (white-noise errors are 2-3x too small here)
  inconsistent_sectors per-sector depths disagree (chi-square p < 1e-4, per-sector empirical noise)
  broad_depression     flux 1-3 durations either side of the dip is > 40% as low as the dip itself
                       (part of a slower wiggle, not a crisp transit)
  too_deep             implied radius > 25 R_earth (~2.2 R_jup) or depth > 3% (a star, not a planet)
  few_transits         fewer than 3 transits with data
  common_period        same period (0.1%) in >= 3 other stars (instrument systematic)
  crowded_period       > 3x the typical number of other stars' signals within +-3% in period
                       (soft: TESS-orbit pile-up at 12-25 d, variability at 0.5-1 d)

Hard flags are direct evidence of an impostor; soft flags (SOFT_FLAGS) are context.

Writes results/phase2/triage.csv.
Usage: python pipeline/05_triage.py [WORKERS] [STARS_DIR OUT_CSV]
  (phase 3 runs it on the upside-down light curves: results/phase3/inverted/stars -> inverted_triage.csv)
"""
import json
import sys
from pathlib import Path
from multiprocessing import Pool

import env  # noqa: F401

import numpy as np
import pandas as pd
from astropy.timeseries import BoxLeastSquares
from scipy.stats import chi2 as chi2_dist

import lcsearch as L
from env import RESULTS

OUT = RESULTS / "phase2"
TESS_ORBIT_D = 13.7
# Context that makes a signal less likely but is not evidence against it on its own.
SOFT_FLAGS = {"crowded_period", "tess_orbit_alias", "low_snr_red_noise"}
R_EARTH_PER_R_SUN = 0.009168


def same_period(p1, p2, tol=0.005):
    return any(abs(p1 / p2 - h) < tol * h for h in (0.5, 1.0, 2.0))


def diagnostics(args):
    """Light-curve based checks for all signals of one star."""
    rec, _ = args
    lc, _, _ = L.reject_bad_sectors(L.load_star(rec["tic"]))
    if rec.get("flipped"):  # false-alarm baseline: dips become bumps (see 06_inverted.py)
        lc["flux"] = 2 - lc["flux"]
    return diagnostics_on(L.bin_lc(L.clean(lc)), rec)


def diagnostics_on(c, rec):
    t, f, e = c["time"], c["flux"], c["flux_err"]
    gaps = np.where(np.diff(t) > 0.5)[0]
    edges = np.concatenate([[t[0]], t[gaps], t[gaps + 1], [t[-1]]])
    near_edge = np.min(np.abs(t[:, None] - edges[None, :]), axis=1) < 0.5
    bls = BoxLeastSquares(t, f, e)
    out = []
    for s in rec["signals"]:
        P, t0, dur = s["period"], s["t0"], s["duration_hr"] / 24
        phase = ((t - t0 + 0.5 * P) % P) - 0.5 * P
        intr = np.abs(phase) < dur / 2
        st = bls.compute_stats(P, dur, t0)
        # contribution of each transit to the total signal: sum over in-transit points of (1-f)/e^2
        epoch = np.round((t - t0) / P).astype(int)
        w = (1 - f) / e ** 2
        per_event = pd.Series(w[intr]).groupby(epoch[intr]).sum()
        total = per_event.clip(lower=0).sum()
        top_frac = float(per_event.max() / total) if total > 0 else 1.0
        # does the signal survive without its two strongest events? (two glitches can fake a period)
        def box_snr(sel):
            if sel.sum() <= 3:
                return 0.0
            wk = 1 / e[sel] ** 2
            return float(np.sum((1 - f[sel]) * wk) / wk.sum() * np.sqrt(wk.sum()))

        snr_all = box_snr(intr)
        snr_wo = box_snr(intr & ~np.isin(epoch, per_event.nlargest(2).index))

        # Empirical noise: the same box placed at 40 phases away from transit and secondary.
        # Correlated (red) noise makes the white-noise BLS errors 2-3x too small (see flag_calibration).
        def box_depth(sel):
            wk = 1 / e[sel] ** 2
            return np.sum((1 - f[sel]) * wk) / wk.sum() if sel.sum() > 3 else np.nan

        offsets = np.concatenate([np.linspace(0.1, 0.4, 20), np.linspace(0.6, 0.9, 20)]) * P
        null = np.array([box_depth(np.abs(((t - t0 - o + 0.5 * P) % P) - 0.5 * P) < dur / 2) for o in offsets])
        null = null[np.isfinite(null)]
        sigma_emp = 1.4826 * np.median(np.abs(null - np.median(null))) if len(null) > 10 else np.nan

        # Per-sector depth consistency, each sector with its own empirical noise.
        def robust_sigma(x):
            x = x[np.isfinite(x)]
            return 1.4826 * np.median(np.abs(x - np.median(x))) if len(x) > 10 else np.nan

        sec_d, sec_s = [], []
        for sec in np.unique(c["sector"]):
            ms = c["sector"] == sec
            if np.sum(intr & ms) < 3:
                continue
            sec_d.append(box_depth(intr & ms))
            sec_s.append(robust_sigma(np.array([
                box_depth(ms & (np.abs(((t - t0 - o + 0.5 * P) % P) - 0.5 * P) < dur / 2)) for o in offsets])))
        sec_d, sec_s = np.array(sec_d), np.array(sec_s)
        ok = np.isfinite(sec_d) & np.isfinite(sec_s) & (sec_s > 0)
        if ok.sum() >= 3:
            w_s = 1 / sec_s[ok] ** 2
            mean_d = np.sum(sec_d[ok] * w_s) / w_s.sum()
            chi2 = np.sum(((sec_d[ok] - mean_d) / sec_s[ok]) ** 2)
            sector_p = float(chi2_dist.sf(chi2, ok.sum() - 1))
        else:
            sector_p = np.nan

        # Broad depression: flux 1-3 durations either side of transit also dragged down
        wings = (np.abs(phase) > dur) & (np.abs(phase) < 3 * dur)
        wing_depth = box_depth(wings)
        rho = rec["rho_sun"]
        out.append({
            "tic": rec["tic"], "signal": s["signal"], "rho": rho,
            "dur_ratio": dur / L.duration_days(P, rho),
            "edge_frac": float(near_edge[intr].mean()) if intr.any() else 1.0,
            "top_event_frac": top_frac,
            "snr_all_events": snr_all,
            "snr_without_top2": snr_wo,
            "sigma_emp_ppm": float(sigma_emp * 1e6),
            "sector_consistency_p": sector_p,
            "wing_depth_ppm": float(wing_depth * 1e6),
            "snr_emp": float(s["depth_ppm"] / (sigma_emp * 1e6)) if sigma_emp > 0 else np.nan,
            "n_events": int(per_event.size),
            "depth_secondary_ppm": float(st["depth_phased"][0] * 1e6),
            "depth_secondary_err_ppm": float(st["depth_phased"][1] * 1e6),
        })
    return out


def main():
    stars_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else OUT / "stars"
    out_csv = Path(sys.argv[3]) if len(sys.argv) > 3 else OUT / "triage.csv"
    recs = [json.loads(p.read_text()) for p in sorted(stars_dir.glob("*.json"))]
    recs = [r for r in recs if r["status"] == "ok"]
    tg = pd.read_csv(RESULTS / "targets_ranked.csv").set_index("TIC")

    rows = []
    for r in recs:
        m = tg.loc[r["tic"]]
        for s in r["signals"]:
            rcont = m["Rcont"] if np.isfinite(m["Rcont"]) else 0
            depth = max(s["depth_ppm"], 0) * 1e-6 * (1 + rcont)
            rows.append({"tic": r["tic"], **s, "rp_earth": np.sqrt(depth) * m["Rad"] / R_EARTH_PER_R_SUN,
                         "tmag": m["Tmag"], "teff": m["Teff"], "rstar": m["Rad"], "dist_pc": m["Dist"],
                         "n_sectors_used": len(r["sectors_used"]), "sectors_rejected": len(r["rejected"]),
                         "baseline_d": r["baseline_d"]})
    sig = pd.DataFrame(rows)
    sig["passes"] = (sig["sde"] >= L.SDE_THRESHOLD) & (sig["snr"] >= L.SNR_THRESHOLD) & \
                    (sig["n_transits"] >= L.MIN_TRANSITS)

    # light-curve diagnostics only for stars with passing signals
    todo = [(r, None) for r in recs if r["tic"] in set(sig.loc[sig["passes"], "tic"])]
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 12) as pool:
        diag = pd.DataFrame([d for ds in pool.map(diagnostics, todo) for d in ds])
    sig = sig.merge(diag, on=["tic", "signal"], how="left")

    # Systematics show up as periods shared across stars. Two measures:
    #  - exact: other stars within 0.1% in period (chance expectation ~0.4 for ~750 signals)
    #  - crowding: other stars' signals within +-3% vs the median density over log period
    passing = sig[sig["passes"]]
    pp, pt = passing["period"].to_numpy(), passing["tic"].to_numpy()
    sig["n_other_stars_same_period"] = [int(np.sum((np.abs(pp / p - 1) < 0.001) & (pt != t)))
                                        for p, t in zip(sig["period"], sig["tic"])]
    local = np.array([np.sum((np.abs(np.log(pp / p)) < 0.03) & (pt != t)) for p, t in zip(sig["period"], sig["tic"])])
    grid = np.geomspace(L.PERIOD_MIN_D, L.PERIOD_MAX_D, 200)
    typical = np.median([np.sum(np.abs(np.log(pp / g)) < 0.03) for g in grid])
    sig["period_crowding"] = local / max(typical, 1)

    def flags(s):
        f = []
        others = sig[(sig["tic"] == s["tic"]) & (sig["signal"] != s["signal"])]
        if any(same_period(s["period"], p) for p in others["period"]):
            f.append("variability_repeat")
        if s.get("dur_ratio", 0) > 1.6:
            f.append("long_duration")
        if any(abs(s["period"] / (TESS_ORBIT_D * h) - 1) < 0.03 for h in (0.5, 1, 1.5, 2)):
            f.append("tess_orbit_alias")
        if s.get("edge_frac", 0) > 0.5:
            f.append("edge_dominated")
        if s.get("top_event_frac", 0) > 0.5:
            f.append("one_event_dominated")
        if s.get("snr_without_top2", 99) < 0.5 * s.get("snr_all_events", 0):
            f.append("few_events_dominate")
        # empirical noise: each half (odd/even) has sqrt(2)x the noise, so the difference has 2x
        sig_e = s.get("sigma_emp_ppm", np.nan)
        if np.isfinite(sig_e) and sig_e > 0:
            if abs(s["depth_odd_ppm"] - s["depth_even_ppm"]) / (2 * sig_e) > 3:
                f.append("odd_even_mismatch")
            if s["depth_secondary_ppm"] / sig_e > 3:
                f.append("secondary_eclipse")
            if s["snr_emp"] < L.SNR_THRESHOLD:
                f.append("low_snr_red_noise")
        if s.get("sector_consistency_p", 1) < 1e-4:  # 1e-3 flagged 2/20 injected planets, 1e-4 only 1 (on a glitch)
            f.append("inconsistent_sectors")
        if s.get("wing_depth_ppm", 0) > 0.4 * s["depth_ppm"]:
            f.append("broad_depression")
        if s["depth_ppm"] > 30000 or s["rp_earth"] > 25:
            f.append("too_deep")
        if s["n_transits"] < L.MIN_TRANSITS:
            f.append("few_transits")
        if s["n_other_stars_same_period"] >= 3:
            f.append("common_period")
        if s["period_crowding"] > 3:
            f.append("crowded_period")
        return ";".join(f)

    sig["flags"] = [flags(s) for _, s in sig.iterrows()]
    split = [fl.split(";") if fl else [] for fl in sig["flags"]]
    sig["hard_flags"] = [";".join(x for x in fl if x not in SOFT_FLAGS) for fl in split]
    sig["soft_flags"] = [";".join(x for x in fl if x in SOFT_FLAGS) for fl in split]
    sig["clean"] = sig["passes"] & (sig["flags"] == "")
    sig["survivor"] = sig["passes"] & (sig["hard_flags"] == "")
    sig = sig.sort_values(["survivor", "clean", "sde"], ascending=False)
    sig.to_csv(out_csv, index=False)

    p = sig[sig["passes"]]
    print(f"{len(recs)} stars searched; {len(p)} signals pass thresholds; "
          f"{int(p['survivor'].sum())} have no hard flags; {int(p['clean'].sum())} have no flags at all")
    print("\nred-flag counts among passing signals:")
    print(pd.Series([x for fl in p["flags"] for x in fl.split(";") if x]).value_counts().to_string())
    cols = ["tic", "signal", "period", "depth_ppm", "duration_hr", "dur_ratio", "rp_earth", "snr", "sde",
            "n_transits", "snr_emp", "snr_without_top2"]
    print("\nsignals with no hard flags (soft flags shown):")
    print(sig[sig["survivor"]][cols + ["soft_flags"]].head(40).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
