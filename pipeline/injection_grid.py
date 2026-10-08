"""Injection-recovery test: plant box transits in real light curves, count how many the search finds.

Used to choose search settings (e.g. detrending window) on evidence rather than single examples.
Usage: python pipeline/injection_grid.py WINDOW_D [N_STARS] [RP_EARTH]
"""
import sys
from multiprocessing import Pool

import env  # noqa: F401
import numpy as np
import pandas as pd

import lcsearch as L
from env import LIGHTCURVES, RESULTS


def trial(args):
    tic, rho, rstar, rcont, rp, period, phase, window = args
    depth = (rp * 0.009168 / rstar) ** 2 / (1 + rcont)
    dur = L.duration_days(period, rho)
    lc, _, _ = L.reject_bad_sectors(L.load_star(tic))
    t0 = lc["time"].min() + phase * period
    lc["flux"] = lc["flux"] - depth * L.transit_mask(lc["time"], period, t0, dur, factor=1.0)
    s = L.bls_search(L.bin_lc(L.clean(lc, window)), rho, max_signals=1)[0][0]
    dt = abs(((s["t0"] - t0 + period / 2) % period) - period / 2)
    hit = abs(s["period"] - period) / period < 1e-3 and dt < dur and s["sde"] >= L.SDE_THRESHOLD
    return {"tic": tic, "window": window, "rp": rp, "period": period, "depth_in": depth * 1e6,
            "depth_out": s["depth_ppm"], "sde": s["sde"], "snr": s["snr"], "recovered": hit}


def main():
    window = float(sys.argv[1])
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    rp = float(sys.argv[3]) if len(sys.argv) > 3 else 1.5
    tg = pd.read_csv(RESULTS / "targets_ranked.csv")
    tg = tg[(~tg["has_2min"]) & (tg["snr_ref"] >= 10)]
    tg = tg[[(LIGHTCURVES / f"{t}.npz").exists() for t in tg["TIC"]]]
    rng = np.random.default_rng(42)  # same stars, periods and phases for every window
    pick = tg.sample(n, random_state=1)
    jobs = []
    for _, r in pick.iterrows():
        for _ in range(2):
            jobs.append((int(r["TIC"]), float(r["Mass"] / r["Rad"] ** 3), float(r["Rad"]),
                         float(r["Rcont"]), rp, float(rng.uniform(1, 15)), float(rng.uniform()), window))
    with Pool(12) as pool:
        res = pd.DataFrame(pool.map(trial, jobs))
    res.to_csv(RESULTS / f"injection_w{window}_rp{rp}.csv", index=False)
    print(f"window {window} d, Rp {rp}: recovered {res.recovered.sum()}/{len(res)}; "
          f"median depth retained {np.median(res.depth_out / res.depth_in):.2f}; median SDE {res.sde.median():.1f}")


if __name__ == "__main__":
    main()
