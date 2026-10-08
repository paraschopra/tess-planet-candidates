"""Phase 2b: iterative BLS transit search over every downloaded target star.

Per star: reject bad sectors -> detrend/clip -> bin -> BLS (mask & repeat, up to 3 signals).
Results are cached per star in results/phase2/stars/<TIC>.json, so the script can be rerun while
downloads are still in progress and only new stars get searched.

Writes results/phase2/signals.csv (one row per signal) with simple red flags for later vetting.

Usage: python pipeline/03_search.py [WORKERS]
"""
import json
import sys
import time
import traceback
from multiprocessing import Pool

import env  # noqa: F401

import numpy as np
import pandas as pd
from tqdm import tqdm

import lcsearch as L
from env import LIGHTCURVES, RESULTS

OUT = RESULTS / "phase2"
STARS = OUT / "stars"
PERIODOGRAMS = OUT / "periodograms"
TESS_ORBIT_D = 13.7  # spacecraft orbit; its harmonics are a classic source of fake periods


def search_star(args):
    tic, rho = args
    out = STARS / f"{tic}.json"
    if out.exists():
        return json.loads(out.read_text())
    t0 = time.time()
    try:
        lc, bad, noise = L.reject_bad_sectors(L.load_star(tic))
        if len(lc["time"]) == 0:
            rec = {"tic": tic, "status": "all_sectors_rejected", "rejected": bad}
        else:
            c = L.bin_lc(L.clean(lc))
            signals, (periods, sde) = L.bls_search(c, rho)
            # keep a downsampled periodogram (max SDE per log-period bin) for plots
            edges = np.geomspace(periods.min(), periods.max(), 4001)
            idx = np.clip(np.digitize(periods, edges) - 1, 0, 3999)
            pg = np.full(4000, -np.inf)  # not NaN: max(NaN, x) stays NaN
            np.maximum.at(pg, idx, sde)
            pg[~np.isfinite(pg)] = np.nan
            np.savez_compressed(PERIODOGRAMS / f"{tic}.npz", period=np.sqrt(edges[1:] * edges[:-1]), sde=pg)
            rec = {"tic": tic, "status": "ok", "rejected": bad,
                   "sectors_used": sorted(int(s) for s in np.unique(c["sector"])),
                   "noise_1hr_ppm": {str(k): round(v * 1e6, 1) for k, v in noise.items()},
                   "n_points": len(c["time"]), "baseline_d": float(np.ptp(c["time"])),
                   "rho_sun": rho, "signals": signals}
    except Exception as e:
        rec = {"tic": tic, "status": "error", "error": f"{type(e).__name__}: {e}",
               "trace": traceback.format_exc(limit=3)}
    rec["seconds"] = round(time.time() - t0, 1)
    out.write_text(json.dumps(rec))
    return rec


def red_flags(s):
    flags = []
    ratio = s["period"] / TESS_ORBIT_D
    for h in (0.25, 1 / 3, 0.5, 1, 2):
        if abs(ratio - h) < 0.01 * h:
            flags.append("tess_orbit_alias")
    odd, even = s["depth_odd_ppm"], s["depth_even_ppm"]
    oe_err = np.hypot(s["depth_odd_err_ppm"], s["depth_even_err_ppm"])
    if oe_err > 0 and abs(odd - even) / oe_err > 3:
        flags.append("odd_even_mismatch")
    if s["depth_ppm"] > 30000:
        flags.append("too_deep_for_planet")
    if s["period"] < 0.55:
        flags.append("at_min_period")
    if s["n_transits"] < L.MIN_TRANSITS:
        flags.append("few_transits")
    return flags


def main():
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    for d in (STARS, PERIODOGRAMS):
        d.mkdir(parents=True, exist_ok=True)
    tg = pd.read_csv(RESULTS / "targets_ranked.csv")
    tg = tg[(~tg["has_2min"]) & (tg["snr_ref"] >= 10)]
    tg = tg[[(LIGHTCURVES / f"{t}.npz").exists() for t in tg["TIC"]]]
    rho = (tg["Mass"] / tg["Rad"] ** 3).fillna(1.0).clip(0.2, 50)
    jobs = list(zip(tg["TIC"].astype(int), rho.astype(float)))
    print(f"Searching {len(jobs)} stars with {workers} workers")
    with Pool(workers) as pool:
        recs = list(tqdm(pool.imap_unordered(search_star, jobs), total=len(jobs)))

    status = pd.Series([r["status"] for r in recs]).value_counts()
    print(status.to_string())
    rows = []
    meta = tg.set_index("TIC")
    for r in recs:
        if r["status"] != "ok":
            continue
        m = meta.loc[r["tic"]]
        for s in r["signals"]:
            # planet radius implied by the depth (diluted depth corrected by Rcont)
            depth = max(s["depth_ppm"], 0) * 1e-6 * (1 + (m["Rcont"] if np.isfinite(m["Rcont"]) else 0))
            rp = np.sqrt(depth) * m["Rad"] / 0.009168
            rows.append({"tic": r["tic"], **s, "rp_earth": rp, "tmag": m["Tmag"], "teff": m["Teff"],
                         "rstar": m["Rad"], "dist_pc": m["Dist"], "n_sectors_used": len(r["sectors_used"]),
                         "sectors_rejected": len(r["rejected"]), "flags": ";".join(red_flags(s))})
    sig = pd.DataFrame(rows)
    sig["passes"] = (sig["sde"] >= L.SDE_THRESHOLD) & (sig["snr"] >= L.SNR_THRESHOLD) & \
                    (sig["n_transits"] >= L.MIN_TRANSITS)
    sig = sig.sort_values(["passes", "sde"], ascending=False)
    sig.to_csv(OUT / "signals.csv", index=False)
    p = sig[sig["passes"]]
    print(f"\n{len(p)} signals pass SDE>={L.SDE_THRESHOLD}, SNR>={L.SNR_THRESHOLD}, "
          f">={L.MIN_TRANSITS} transits ({(p['flags'] == '').sum()} with no red flags)")
    cols = ["tic", "signal", "period", "depth_ppm", "duration_hr", "rp_earth", "snr", "sde",
            "n_transits", "flags"]
    print(p[cols].head(40).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
