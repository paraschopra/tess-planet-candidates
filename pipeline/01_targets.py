"""Phase 1: build a ranked list of nearby stars worth searching for new small planets.

Steps (each cached to data/catalogs so reruns are cheap):
  1. Gaia nearby stars (GCNS) within MAX_DIST_PC, matched to TIC v8.2, Tmag < MAX_TMAG.
  2. Drop known planet hosts, TOIs and CTOIs.
  3. Ask MAST which TESS light curves exist (SPOC 2-min, TESS-SPOC FFI, QLP FFI).
  4. Score each star by expected SNR of a REF_RADIUS_EARTH planet at REF_PERIOD_D.

Usage: python pipeline/01_targets.py [--refresh]
"""
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import env  # noqa: F401  (sets cache dirs)

import numpy as np
import pandas as pd
from astroquery.mast import Catalogs, Observations
from astroquery.vizier import Vizier
from tqdm import tqdm

from env import CATALOGS, RESULTS

MAX_DIST_PC = 60
MAX_TMAG = 13.0
MAX_STAR_RADIUS = 1.5  # R_sun; small stars make small planets detectable
MIN_SECTORS = 3
MIN_YEARS_SPANNED = 2.0
REF_RADIUS_EARTH = 1.5
REF_PERIOD_D = 10.0
LC_AUTHORS = ["SPOC", "TESS-SPOC", "QLP"]

TIC_CSV = CATALOGS / "tic_nearby.csv"
COVERAGE_CSV = CATALOGS / "tess_coverage.csv"
OUT_CSV = RESULTS / "targets_ranked.csv"

R_EARTH_PER_R_SUN = 0.009168

# Rough 1-hour TESS photometric precision (ppm) vs Tmag, used only for ranking;
# real noise gets measured from each light curve in phase 2.
NOISE_TMAG = np.array([6, 7, 8, 9, 10, 11, 12, 13, 14])
NOISE_PPM = np.array([40, 45, 60, 110, 230, 480, 1000, 2200, 5000])


def run_chunks(func, items, chunk, name, workers=8, retries=4):
    """Apply func to chunks of items in parallel; each chunk is cached as CSV so reruns resume."""
    cache = CATALOGS / f"_chunks_{name}"
    cache.mkdir(exist_ok=True)
    chunks = [items[i:i + chunk] for i in range(0, len(items), chunk)]

    def work(k):
        f = cache / f"{k:05d}.csv"
        if f.exists():
            return pd.read_csv(f)
        for attempt in range(retries):
            try:
                df = func(chunks[k])
                df.to_csv(f, index=False)
                return df
            except Exception as e:  # MAST drops connections under load; back off and retry
                if attempt == retries - 1:
                    print(f"chunk {k} failed: {e}")
                    return pd.DataFrame()
                time.sleep(10 * (attempt + 1))

    with ThreadPoolExecutor(workers) as ex:
        out = list(tqdm(ex.map(work, range(len(chunks))), total=len(chunks), desc=name))
    return pd.concat([d for d in out if len(d)], ignore_index=True)


def fetch_tic(refresh=False):
    """Nearby stars from the Gaia Catalogue of Nearby Stars, matched to TIC v8.2 by Gaia ID.

    (Filtering TIC directly on distance is a full-table scan and times out.)
    """
    if TIC_CSV.exists() and not refresh:
        return pd.read_csv(TIC_CSV)
    v = Vizier(columns=["GaiaEDR3", "Plx", "Gmag"], row_limit=-1, timeout=600)
    gcns = v.query_constraints(catalog="J/A+A/649/A6/table1c", Plx=f">{1000 / MAX_DIST_PC:.3f}",
                               Gmag=f"<{MAX_TMAG + 1}")[0].to_pandas()
    print(f"GCNS stars within {MAX_DIST_PC} pc, G<{MAX_TMAG + 1}: {len(gcns)}")
    cols = ["ID", "ra", "dec", "GAIA", "HIP", "Tmag", "Teff", "logg", "rad", "mass", "d",
            "lumclass", "disposition", "contratio", "objType"]

    def lookup(gaia_ids):
        return Catalogs.query_criteria(catalog="Tic", GAIA=list(gaia_ids)).to_pandas()[cols]

    tic = run_chunks(lookup, [str(g) for g in gcns["GaiaEDR3"]], 500, "tic")
    tic = tic.rename(columns={"ID": "TIC", "rad": "Rad", "mass": "Mass", "d": "Dist",
                              "lumclass": "LClass", "disposition": "Disp", "contratio": "Rcont"})
    tic = tic.drop_duplicates("TIC")
    tic = tic[(tic["Tmag"] < MAX_TMAG) & (tic["Dist"] < MAX_DIST_PC)]
    tic.to_csv(TIC_CSV, index=False)
    return tic


def known_hosts():
    toi = pd.read_csv(CATALOGS / "toi.csv", usecols=["TIC ID"])["TIC ID"]
    ctoi = pd.read_csv(CATALOGS / "ctoi.csv", usecols=["TIC ID"])["TIC ID"]
    ps = pd.read_csv(CATALOGS / "nasa_ps_tic.csv")["tic_id"].dropna()
    ps = ps.str.replace("TIC ", "", regex=False).astype(np.int64)
    return set(toi.astype(np.int64)) | set(ctoi.astype(np.int64)) | set(ps)


def fetch_coverage(tic_ids, refresh=False):
    if COVERAGE_CSV.exists() and not refresh:
        return pd.read_csv(COVERAGE_CSV)
    keep = ["target_name", "provenance_name", "sequence_number", "t_exptime", "t_min", "t_max",
            "obs_collection"]

    def lookup(ids):
        obs = Observations.query_criteria(target_name=list(ids), dataproduct_type="timeseries",
                                          provenance_name=LC_AUTHORS)
        return obs.to_pandas()[keep] if len(obs) else pd.DataFrame(columns=keep)

    cov = run_chunks(lookup, [str(int(t)) for t in tic_ids], 300, "coverage")
    cov["tic"] = cov["target_name"].astype(np.int64)
    cov.to_csv(COVERAGE_CSV, index=False)
    return cov


def summarize_coverage(cov):
    cov = cov.dropna(subset=["sequence_number"])
    g = cov.groupby("tic")
    s = pd.DataFrame({
        "n_sectors": g["sequence_number"].nunique(),
        "sectors": g["sequence_number"].apply(lambda x: ",".join(str(int(v)) for v in sorted(set(x)))),
        "first_mjd": g["t_min"].min(),
        "last_mjd": g["t_max"].max(),
        "has_2min": g.apply(lambda d: bool(((d.provenance_name == "SPOC") & (d.t_exptime <= 120)).any()),
                            include_groups=False),
        "authors": g["provenance_name"].apply(lambda x: ",".join(sorted(set(x)))),
    })
    s["years_spanned"] = (s["last_mjd"] - s["first_mjd"]) / 365.25
    return s.reset_index()


def expected_snr(df):
    """SNR of REF_RADIUS_EARTH planet at REF_PERIOD_D, central transit, ~27 d per sector."""
    rs = df["Rad"].fillna(1.0).clip(0.1, MAX_STAR_RADIUS)
    ms = df["Mass"].fillna(rs).clip(0.08, 2.0)
    depth_ppm = (REF_RADIUS_EARTH * R_EARTH_PER_R_SUN / rs) ** 2 * 1e6
    # Duration (hr) for circular orbit, b=0: T ~ 13 hr * (P/365)^(1/3) * rho^(-1/3)
    rho = ms / rs ** 3
    dur_hr = 13.0 * (REF_PERIOD_D / 365.25) ** (1 / 3) * rho ** (-1 / 3)
    sigma_1hr = np.interp(df["Tmag"], NOISE_TMAG, NOISE_PPM)
    n_transits = df["n_sectors"] * 27.0 / REF_PERIOD_D
    snr = depth_ppm / sigma_1hr * np.sqrt(n_transits * dur_hr)
    return depth_ppm, dur_hr, snr


def main():
    refresh = "--refresh" in sys.argv
    tic = fetch_tic(refresh)
    print(f"TIC stars <{MAX_DIST_PC} pc, Tmag<{MAX_TMAG}: {len(tic)}")

    tic = tic[~tic["Disp"].isin(["SPLIT", "DUPLICATE", "ARTIFACT"])]
    tic = tic[(tic["LClass"] == "DWARF") & (tic["Rad"] < MAX_STAR_RADIUS)]
    print(f"  dwarfs with R < {MAX_STAR_RADIUS} R_sun: {len(tic)}")

    known = known_hosts()
    tic = tic[~tic["TIC"].isin(known)]
    print(f"  not a known planet host / TOI / CTOI: {len(tic)}")

    cov = summarize_coverage(fetch_coverage(tic["TIC"], refresh))
    df = tic.merge(cov, left_on="TIC", right_on="tic", how="inner").drop(columns="tic")
    print(f"  with any TESS light curve: {len(df)}")

    df = df[(df["n_sectors"] >= MIN_SECTORS) & (df["years_spanned"] >= MIN_YEARS_SPANNED)]
    print(f"  >= {MIN_SECTORS} sectors spanning >= {MIN_YEARS_SPANNED} yr: {len(df)}")

    df["depth_ppm"], df["dur_hr"], df["snr_ref"] = expected_snr(df)
    df = df.sort_values("snr_ref", ascending=False)
    RESULTS.mkdir(exist_ok=True)
    df.to_csv(OUT_CSV, index=False)
    print(f"\nWrote {OUT_CSV} ({len(df)} stars). Expected SNR of a {REF_RADIUS_EARTH} R_earth, "
          f"{REF_PERIOD_D} d planet (7.1 is the usual detection threshold):")
    print(f"  SNR>=7: {(df.snr_ref >= 7).sum()}   SNR>=10: {(df.snr_ref >= 10).sum()}")
    show = ["TIC", "Tmag", "Teff", "Rad", "Dist", "n_sectors", "years_spanned", "has_2min",
            "authors", "depth_ppm", "snr_ref"]
    print(df[show].head(20).round(2).to_string(index=False))


if __name__ == "__main__":
    main()
