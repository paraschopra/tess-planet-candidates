"""Phase 2a: download full-frame-image light curves for the target stars.

For every sector we keep one light curve: TESS-SPOC (NASA) if it exists, else QLP (MIT).
Each star is saved as data/lightcurves/<TIC>.npz with per-point arrays:
  time (BTJD = BJD - 2457000), flux (normalised to 1 per sector), flux_err, sector, exptime (s)
plus `authors`, the light-curve product used for each sector.

Bad sectors are kept here and rejected at search time (QLP can be junk in crowded fields).

Usage: python pipeline/02_download.py [N_STARS]   (default: all FFI-only stars with snr_ref >= 10)
"""
import re
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor

from pathlib import Path

import env  # noqa: F401

import lightkurve as lk
import numpy as np
import pandas as pd
from tqdm import tqdm

from env import CACHE, LIGHTCURVES, RESULTS

AUTHOR_PREFERENCE = ["TESS-SPOC", "QLP"]
# QLP files changed format over time; take the first systematics-corrected flux that exists.
QLP_FLUX_COLUMNS = ["kspsap_flux", "det_flux", "sap_flux"]
MIN_SNR_REF = 10
WORKERS = 8


def select_targets():
    df = pd.read_csv(RESULTS / "targets_ranked.csv")
    return df[(~df["has_2min"]) & (df["snr_ref"] >= MIN_SNR_REF)]


def pick_products(sr):
    """One product per sector: preferred author first, then shortest cadence."""
    t = sr.table.to_pandas()
    t["row"] = np.arange(len(t))
    t["sector"] = t["mission"].str.extract(r"Sector (\d+)").astype(float)
    t = t[t["author"].isin(AUTHOR_PREFERENCE)].dropna(subset=["sector"])
    t["rank"] = t["author"].map({a: i for i, a in enumerate(AUTHOR_PREFERENCE)})
    t = t.sort_values(["sector", "rank", "exptime"]).drop_duplicates("sector")
    return t


def flux_arrays(lc, author):
    if author == "QLP":
        col = next(c for c in QLP_FLUX_COLUMNS if c in lc.colnames)
    else:
        col = "pdcsap_flux"
    flux = np.asarray(lc[col].value, dtype=float)
    err_col = col + "_err" if col + "_err" in lc.colnames else "flux_err"
    err = np.asarray(lc[err_col].value, dtype=float)
    time = np.asarray(lc.time.value, dtype=float)
    ok = np.isfinite(time) & np.isfinite(flux) & np.isfinite(err) & (err > 0)
    time, flux, err = time[ok], flux[ok], err[ok]
    med = np.nanmedian(flux)
    return time, flux / med, err / med, col


def download_product(row):
    """Download one light curve; if the cached file is corrupt (interrupted download), refetch once."""
    for attempt in range(2):
        try:
            return row.download(download_dir=str(CACHE / "lightkurve"), quality_bitmask="default")
        except lk.LightkurveError as e:
            bad = re.search(r"product (\S+\.fits)", str(e))
            if attempt == 1 or not bad:
                raise
            Path(bad.group(1)).unlink(missing_ok=True)


def download_star(tic):
    out = LIGHTCURVES / f"{tic}.npz"
    if out.exists():
        return {"tic": tic, "status": "cached"}
    try:
        sr = lk.search_lightcurve(f"TIC {tic}", author=AUTHOR_PREFERENCE)
        if len(sr) == 0:
            return {"tic": tic, "status": "no_products"}
        picks = pick_products(sr)
        parts = {"time": [], "flux": [], "flux_err": [], "sector": [], "exptime": []}
        authors, failed = [], []
        for _, p in picks.iterrows():
            try:
                lc = download_product(sr[int(p["row"])])
                if lc is None:
                    continue
                t, f, e, col = flux_arrays(lc, p["author"])
            except Exception as e:  # one bad sector shouldn't lose the whole star
                failed.append(f"{int(p['sector'])}:{type(e).__name__}")
                continue
            if len(t) < 100:
                continue
            parts["time"].append(t)
            parts["flux"].append(f)
            parts["flux_err"].append(e)
            parts["sector"].append(np.full(len(t), int(p["sector"])))
            parts["exptime"].append(np.full(len(t), float(p["exptime"])))
            authors.append(f"{int(p['sector'])}:{p['author']}:{col}")
        if not parts["time"]:
            return {"tic": tic, "status": "empty", "failed_sectors": ";".join(failed)}
        arrays = {k: np.concatenate(v) for k, v in parts.items()}
        order = np.argsort(arrays["time"])
        np.savez_compressed(out, **{k: v[order] for k, v in arrays.items()},
                            authors=np.array(authors))
        return {"tic": tic, "status": "ok", "n_sectors": len(authors), "n_points": len(order),
                "failed_sectors": ";".join(failed)}
    except Exception as e:
        return {"tic": tic, "status": "error", "error": f"{type(e).__name__}: {e}",
                "trace": traceback.format_exc(limit=2)}


def main():
    LIGHTCURVES.mkdir(parents=True, exist_ok=True)
    targets = select_targets()
    if len(sys.argv) > 1:
        targets = targets.head(int(sys.argv[1]))
    print(f"Downloading {len(targets)} stars")
    # Processes, not threads: astropy's FITS reader is not thread-safe ("I/O operation on closed file").
    with ProcessPoolExecutor(WORKERS) as ex:
        log = list(tqdm(ex.map(download_star, targets["TIC"].astype(int)), total=len(targets)))
    log = pd.DataFrame(log)
    log.to_csv(RESULTS / "02_download_log.csv", index=False)
    print(log["status"].value_counts().to_string())
    errs = log[log["status"] == "error"]
    if len(errs):
        print(errs["error"].value_counts().head(10).to_string())


if __name__ == "__main__":
    main()
