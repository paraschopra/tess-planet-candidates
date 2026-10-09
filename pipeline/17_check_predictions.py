"""Phase 4: test the preregistered predictions when TESS re-observes the candidates (2027-28).

Run any time; it reports "no data yet" until the new sectors are public on MAST.

  1. Verifies the local prediction CSVs are byte-identical to the preregistered Zenodo copies (MD5),
     doi:10.5281/zenodo.23257784, so the test uses exactly what was published in advance.
  2. Takes the ephemeris (P, T0, duration, expected depth) straight from those CSVs.
  3. Downloads only the NEW sectors (121, 122, 129): SPOC 2-min if granted, else TESS-SPOC, else QLP.
  4. Measures the dip at the predicted times (transit-masked detrending, matched-processing empirical
     error), and scans the whole orbit to see where the deepest dip actually is.

DECISION RULE, fixed on 2026-10-09 before any new data exist (per candidate, all new sectors combined):
  RECURS     depth at the predicted times >= 3 sigma, AND the deepest point of the orbit scan lies within
             the predicted window (|offset| <= T14/2 + 3 sigma_t), AND depth within 2 sigma of the
             preregistered depth.
  NOT SEEN   depth < 1 sigma while the preregistered depth would have been >= 3 sigma in these data
             (the data were sensitive enough, and the dip is absent).
  INCONCLUSIVE  anything else (e.g. too few transits, or a weak dip).
RECURS supports a real periodic signal, but says nothing about WHICH star it is on (blends remain).

Usage: python pipeline/17_check_predictions.py            (new sectors 121, 122, 129)
       python pipeline/17_check_predictions.py --test 83 84 (dry run on old sectors, to check the method)
Writes results/prediction_check/<TIC>.json and <TIC>.png.
"""
import hashlib
import importlib
import json
import sys
import urllib.request

import env  # noqa: F401

import lightkurve as lk
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import lcsearch as L
import vet as V
from env import RESULTS

dl = importlib.import_module("02_download")

ZENODO_RECORD = "23257784"
NEW_SECTORS = [121, 122, 129]
CANDIDATES = {370612437: "Aakanksha", 252721481: "Orbit"}
OUT = RESULTS / "prediction_check"
AUTHORS = ["SPOC", "TESS-SPOC", "QLP"]  # preference order: 2-min SPOC exists only if a DDT/GI request is granted


def verify_preregistration():
    """MD5 of local prediction CSVs vs the published Zenodo files."""
    try:
        with urllib.request.urlopen(f"https://zenodo.org/api/records/{ZENODO_RECORD}", timeout=60) as r:
            files = {f["key"]: f["checksum"].split(":")[1] for f in json.load(r)["files"]}
    except Exception as e:
        print(f"  (could not reach Zenodo to verify checksums: {e})")
        return None
    ok = True
    for tic in CANDIDATES:
        name = f"{tic}_predicted_transits.csv"
        local = hashlib.md5((V.OUT / name).read_bytes()).hexdigest()
        match = files.get(name) == local
        ok &= match
        print(f"  {name}: {'matches Zenodo' if match else 'DOES NOT MATCH Zenodo'}")
    return ok


def ephemeris_from_csv(tic):
    d = pd.read_csv(V.OUT / f"{tic}_predicted_transits.csv")
    P, t0_bjd = np.polyfit(d.epoch.values, d.mid_BJD_TDB.values, 1)  # exact linear ephemeris
    return {"P": P, "t0": t0_bjd - 2457000, "dur": d.duration_hr.iloc[0] / 24, "depth_ppm": float(d.depth_ppm.iloc[0]),
            "sigma_t_d": d.sigma_min.max() / 1440}


def fetch(tic, sectors):
    sr = lk.search_lightcurve(f"TIC {tic}", sector=sectors, author=AUTHORS)
    if len(sr) == 0:
        return None, []
    t = sr.table.to_pandas()
    t["row"] = np.arange(len(t))
    t["sector"] = t["mission"].str.extract(r"Sector (\d+)").astype(float)
    t["rank"] = t["author"].map({a: i for i, a in enumerate(AUTHORS)})
    t = t.sort_values(["sector", "rank", "exptime"]).drop_duplicates("sector")
    parts = {k: [] for k in ("time", "flux", "flux_err", "sector", "exptime")}
    used = []
    for _, p in t.iterrows():
        lc = dl.download_product(sr[int(p["row"])])
        author = p["author"] if p["author"] != "SPOC" else "TESS-SPOC"  # SPOC 2-min uses pdcsap_flux too
        tt, f, e, col = dl.flux_arrays(lc, author)
        parts["time"].append(tt)
        parts["flux"].append(f)
        parts["flux_err"].append(e)
        parts["sector"].append(np.full(len(tt), int(p["sector"])))
        parts["exptime"].append(np.full(len(tt), float(p["exptime"])))
        used.append(f"S{int(p['sector'])}:{p['author']}:{int(p['exptime'])}s")
    a = {k: np.concatenate(v) for k, v in parts.items()}
    o = np.argsort(a["time"])
    return {k: v[o] for k, v in a.items()}, used


def check(tic, sectors, tag=""):
    e = ephemeris_from_csv(tic)
    raw, used = fetch(tic, sectors)
    if raw is None:
        return {"tic": tic, "status": "no data yet", "sectors_requested": sectors}
    raw, rejected, _ = L.reject_bad_sectors(raw)
    P, t0, dur = e["P"], e["t0"], e["dur"]
    c = V.masked_clean(raw, P, t0, dur)
    n_tr = int(len(np.unique(np.round((c["time"][np.abs(V.phase_of(c["time"], P, t0)) < dur / 2] - t0) / P))))
    depth, sig, _ = V.depth_emp_matched(raw, P, t0, dur)
    # whole-orbit scan on unmasked (search-style) cleaning, so the predicted phase gets no special treatment
    cs = L.bin_lc(L.clean(raw))
    step = dur / 4
    n = int(0.5 * P / step)
    offs = np.arange(-n, n + 1) * step
    scan = np.array([V.box_depth(cs, np.abs(V.phase_of(cs["time"], P, t0 + o)) < dur / 2) for o in offs])
    peak = float(offs[np.nanargmax(scan)])
    window = dur / 2 + 3 * e["sigma_t_d"]
    snr = depth / sig
    expected_snr = e["depth_ppm"] * 1e-6 / sig
    consistent = abs(depth * 1e6 - e["depth_ppm"]) <= 2 * sig * 1e6
    if snr >= 3 and abs(peak) <= window and consistent:
        verdict = "RECURS"
    elif snr < 1 and expected_snr >= 3:
        verdict = "NOT SEEN"
    else:
        verdict = "INCONCLUSIVE"
    res = {"tic": tic, "status": "measured", "verdict": verdict, "products": used, "rejected_sectors": rejected,
           "n_transits": n_tr, "depth_ppm": depth * 1e6, "err_ppm": sig * 1e6, "snr": snr,
           "expected_depth_ppm": e["depth_ppm"], "expected_snr": expected_snr, "depth_consistent_2sigma": bool(consistent),
           "scan_peak_offset_hr": peak * 24, "allowed_window_hr": window * 24,
           "frac_orbit_as_deep": float(np.mean(scan[np.abs(offs) > dur] >= scan[n]))}
    plot(tic, c, P, t0, dur, e, offs, scan, window, res, tag)
    return res


def plot(tic, c, P, t0, dur, e, offs, scan, window, res, tag=""):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 4.2))
    ph = V.phase_of(c["time"], P, t0) * 24
    sel = np.abs(ph) < 4 * dur * 24
    edges = np.linspace(-4 * dur * 24, 4 * dur * 24, 33)
    idx = np.digitize(ph[sel], edges)
    y = (c["flux"][sel] - 1) * 1e6
    xb = 0.5 * (edges[1:] + edges[:-1])
    yb = [np.mean(y[idx == i]) if (idx == i).sum() else np.nan for i in range(1, len(edges))]
    eb = [np.std(y[idx == i]) / np.sqrt(max((idx == i).sum(), 1)) for i in range(1, len(edges))]
    a1.plot(ph[sel], y, ".", color="0.85", ms=1)
    a1.errorbar(xb, yb, eb, fmt="o", color="k", ms=3)
    a1.axvspan(-dur * 12, dur * 12, color="C1", alpha=0.15, label="predicted transit")
    a1.axhline(-e["depth_ppm"], color="C3", ls="--", lw=1, label=f"preregistered depth {e['depth_ppm']:.0f} ppm")
    a1.set_ylim(-4 * e["depth_ppm"], 3 * e["depth_ppm"])
    a1.set_xlabel("hours from PREDICTED mid-transit")
    a1.set_ylabel("ppm")
    a1.legend(fontsize=8)
    a1.set_title(f"new data folded on the preregistered ephemeris: {res['depth_ppm']:.0f} ± {res['err_ppm']:.0f} ppm "
                 f"({res['snr']:.1f}σ, {res['n_transits']} transits)", fontsize=9)
    a2.plot(offs * 24, scan * 1e6, color="C0")
    a2.axvspan(-window * 24, window * 24, color="C1", alpha=0.15, label="allowed window")
    a2.axhline(0, color="k", lw=0.5)
    a2.set_xlabel("offset from predicted time (hours, whole orbit)")
    a2.set_ylabel("dip depth (ppm)")
    a2.legend(fontsize=8)
    a2.set_title(f"where is the deepest dip? peak at {res['scan_peak_offset_hr']:+.1f} h", fontsize=9)
    fig.suptitle(f"TIC {tic} \"{CANDIDATES[tic]}\": prediction test → {res['verdict']}   ({', '.join(res['products'])})",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / f"{tic}{tag}.png", dpi=90)
    plt.close(fig)


def main():
    test = "--test" in sys.argv
    sectors = [int(a) for a in sys.argv[sys.argv.index("--test") + 1:]] if test else NEW_SECTORS
    OUT.mkdir(parents=True, exist_ok=True)
    print("Preregistration check (Zenodo doi:10.5281/zenodo.23257784):")
    verify_preregistration()
    if test:
        print(f"\nDRY RUN on old sectors {sectors}: these were used to MAKE the predictions, so this only checks the method.")
    for tic in CANDIDATES:
        r = check(tic, sectors, tag="_dryrun" if test else "")
        if not test:
            (OUT / f"{tic}.json").write_text(json.dumps(r, indent=1))
        if r["status"] != "measured":
            print(f"\nTIC {tic}: {r['status']} for sectors {sectors}. Data appear on MAST weeks to months after each sector.")
            continue
        print(f"\nTIC {tic} \"{CANDIDATES[tic]}\": {r['verdict']}")
        print(f"  data: {', '.join(r['products'])}; {r['n_transits']} predicted transits covered")
        print(f"  depth at predicted times: {r['depth_ppm']:.0f} ± {r['err_ppm']:.0f} ppm ({r['snr']:.1f}σ); "
              f"preregistered {r['expected_depth_ppm']:.0f} ppm (would be {r['expected_snr']:.1f}σ in these data)")
        print(f"  deepest dip in the orbit at {r['scan_peak_offset_hr']:+.1f} h (allowed ±{r['allowed_window_hr']:.1f} h); "
              f"{100 * r['frac_orbit_as_deep']:.1f}% of the orbit as deep")
        print(f"  plot: {OUT / (str(tic) + ('_dryrun' if test else '') + '.png')}")


if __name__ == "__main__":
    main()
