"""Phase 3a: false-alarm baseline from upside-down light curves.

A planet only ever makes dips. Flipping each raw light curve (flux -> 2 - flux) turns its dips into
bumps, so the best "transit" BLS finds in flipped data comes purely from noise, stellar variability and
instrument artifacts, run through exactly the same cleaning and search. The distribution of the best
SDE per flipped star tells us how strong a signal chance alone produces across this sample.

Flipping before cleaning keeps the comparison fair: clean() clips upward outliers (now the original
dips) at 3 sigma and downward ones (the original flares) at 10 sigma, just as for the real search.
Up to 3 signals per star and the same record format as 03_search, so the flipped survey then goes
through the identical triage: python pipeline/05_triage.py 12 results/phase3/inverted/stars
results/phase3/inverted_triage.csv

Writes results/phase3/inverted/stars/<TIC>.json. Usage: python pipeline/06_inverted.py [WORKERS]
"""
import json
import sys
from multiprocessing import Pool

import env  # noqa: F401

import numpy as np
import pandas as pd
from tqdm import tqdm

import lcsearch as L
from env import RESULTS

OUT = RESULTS / "phase3"
CACHE = OUT / "inverted" / "stars"


def search_flipped(args):
    tic, rho = args
    out = CACHE / f"{tic}.json"
    if out.exists():
        return json.loads(out.read_text())
    try:
        lc, bad, _ = L.reject_bad_sectors(L.load_star(tic))
        lc["flux"] = 2 - lc["flux"]
        c = L.bin_lc(L.clean(lc))
        signals, _ = L.bls_search(c, rho)
        rec = {"tic": tic, "status": "ok", "flipped": True, "rejected": bad,
               "sectors_used": sorted(int(s) for s in np.unique(c["sector"])),
               "n_points": len(c["time"]), "baseline_d": float(np.ptp(c["time"])), "rho_sun": rho,
               "signals": signals}
    except Exception as e:
        rec = {"tic": tic, "status": "error", "error": f"{type(e).__name__}: {e}"}
    out.write_text(json.dumps(rec))
    return rec


def main():
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    CACHE.mkdir(parents=True, exist_ok=True)
    stars = pd.read_csv(RESULTS / "phase2" / "triage.csv").drop_duplicates("tic")
    done = {int(p.stem) for p in (RESULTS / "phase2" / "stars").glob("*.json") if p.stem.isdigit()}
    # same star list and densities as the real search
    tg = pd.read_csv(RESULTS / "targets_ranked.csv").set_index("TIC")
    tics = sorted(done)
    rho = (tg.loc[tics, "Mass"] / tg.loc[tics, "Rad"] ** 3).fillna(1.0).clip(0.2, 50)
    jobs = list(zip(tics, rho.astype(float)))
    print(f"Searching {len(jobs)} flipped light curves ({len(stars)} stars had signals)")
    with Pool(workers) as pool:
        recs = list(tqdm(pool.imap_unordered(search_flipped, jobs), total=len(jobs)))
    print(pd.Series([r["status"] for r in recs]).value_counts().to_string())


if __name__ == "__main__":
    main()
