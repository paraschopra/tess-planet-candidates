"""Calibrate triage red flags on injected planets: a good flag should almost never fire on them.

Injects a box transit (default 1.5 R_earth) into real light curves, runs the search, and applies the
triage diagnostics to the recovered signal.
Usage: python pipeline/flag_calibration.py [N_STARS] [WORKERS]
"""
import importlib
import sys
from multiprocessing import Pool

import env  # noqa: F401
import numpy as np
import pandas as pd

import lcsearch as L
from env import LIGHTCURVES, RESULTS

triage = importlib.import_module("05_triage")


def trial(args):
    tic, rho, rstar, rcont, period, phase = args
    depth = (1.5 * 0.009168 / rstar) ** 2 / (1 + rcont)
    dur = L.duration_days(period, rho)
    lc, _, _ = L.reject_bad_sectors(L.load_star(tic))
    t0 = lc["time"].min() + phase * period
    lc["flux"] = lc["flux"] - depth * L.transit_mask(lc["time"], period, t0, dur, factor=1.0)
    c = L.bin_lc(L.clean(lc))
    # fit at the injected period only (we are testing the flags, not detection)
    sig = L.bls_search_at(c, rho, period)
    # diagnostics() reloads the star from disk, so pass the injected light curve through instead
    d = triage.diagnostics_on(c, {"tic": tic, "rho_sun": rho, "signals": [sig]})[0]
    oe = np.hypot(sig["depth_odd_err_ppm"], sig["depth_even_err_ppm"])
    se = d["sigma_emp_ppm"]
    return {"tic": tic, "period": period, "snr": sig["snr"], "snr_emp": d["snr_emp"],
            "odd_even_sigma_emp": abs(sig["depth_odd_ppm"] - sig["depth_even_ppm"]) / (2 * se),
            "secondary_sigma_emp": d["depth_secondary_ppm"] / se,
            "odd_even_sigma": abs(sig["depth_odd_ppm"] - sig["depth_even_ppm"]) / oe,
            "secondary_sigma": d["depth_secondary_ppm"] / d["depth_secondary_err_ppm"],
            "dur_ratio": d["dur_ratio"], "top_event_frac": d["top_event_frac"],
            "top2_ratio": d["snr_without_top2"] / max(d["snr_all_events"], 1e-9),
            "edge_frac": d["edge_frac"], "sector_p": d["sector_consistency_p"],
            "wing_ratio": d["wing_depth_ppm"] / max(sig["depth_ppm"], 1e-9)}


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    tg = pd.read_csv(RESULTS / "targets_ranked.csv")
    tg = tg[(~tg["has_2min"]) & (tg["snr_ref"] >= 10)]
    tg = tg[[(LIGHTCURVES / f"{t}.npz").exists() for t in tg["TIC"]]].sample(n, random_state=7)
    rng = np.random.default_rng(7)
    jobs = [(int(r.TIC), float(r.Mass / r.Rad ** 3), float(r.Rad), float(r.Rcont),
             float(rng.uniform(1, 15)), float(rng.uniform())) for r in tg.itertuples() for _ in range(2)]
    with Pool(workers) as pool:
        res = pd.DataFrame(pool.map(trial, jobs))
    res.to_csv(RESULTS / "flag_calibration.csv", index=False)
    print(res.round(2).to_string(index=False))
    print("\nfraction of injected planets that would be flagged:")
    print(f"  odd_even (>3 sigma)       {np.mean(res.odd_even_sigma > 3):.2f}   empirical noise: "
          f"{np.mean(res.odd_even_sigma_emp > 3):.2f}")
    print(f"  secondary (>3 sigma)      {np.mean(res.secondary_sigma > 3):.2f}   empirical noise: "
          f"{np.mean(res.secondary_sigma_emp > 3):.2f}")
    print(f"  low_snr_red_noise (<6)    {np.mean(res.snr_emp < 6):.2f}")
    print(f"  long_duration (>1.6)      {np.mean(res.dur_ratio > 1.6):.2f}")
    print(f"  one_event (>0.5)          {np.mean(res.top_event_frac > 0.5):.2f}")
    print(f"  few_events (ratio <0.5)   {np.mean(res.top2_ratio < 0.5):.2f}")
    print(f"  edge (>0.5)               {np.mean(res.edge_frac > 0.5):.2f}")
    print(f"  inconsistent_sectors      {np.mean(res.sector_p < 1e-3):.2f}")
    print(f"  broad_depression (>0.4)   {np.mean(res.wing_ratio > 0.4):.2f}")


if __name__ == "__main__":
    main()
