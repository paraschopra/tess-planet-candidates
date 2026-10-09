"""Phase 3f: statistical validation with TRICERATOPS (Giacalone et al. 2021).

TRICERATOPS weighs "a planet transits the target" against every false-positive scenario, using the real
TIC/Gaia neighbours, which of them fall in TESS-SPOC's photometric apertures, and the transit shape:
  TP / EB / EBx2P   target: planet / eclipsing binary / EB at twice the period
  PTP / PEB / ...   an unresolved bound companion of the target hosts the transit or eclipse
  STP / SEB / ...   an unresolved chance-aligned background/foreground star does
  DTP / DEB / ...   same, diluted by an unresolved companion
  NTP / NEB / ...   a known nearby star resolved by Gaia hosts it
Outputs FPP (probability it is any false positive) and NFPP (probability the signal is on a resolved
neighbour). Validation customarily needs FPP < 0.015 and NFPP < 0.001 (Giacalone et al. 2021), and
in practice high-resolution imaging; we have none, so unresolved-companion scenarios stay open.

TRICERATOPS pins old numpy/scipy, so it runs in its own environment (.venv-tri, see
requirements-triceratops.txt). Two steps:
  python pipeline/12_triceratops.py export          (main .venv: folded light curves + apertures)
  REQUESTS_CA_BUNDLE=cache/certs/bundle.pem .venv-tri/bin/python pipeline/12_triceratops.py run [N_RUNS]
The TRILEGAL server (stev.oapd.inaf.it, background-star model) omits its intermediate certificate;
cache/certs/bundle.pem is certifi plus that intermediate (ZeroSSL RSA DV SSL CA 2, fetched from the
leaf's AIA URL), so the connection is still fully verified.

Writes results/phase3/<TIC>_triceratops.json and results/phase3/plots/<TIC>_triceratops.png.
"""
import glob
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "phase3"
CANDIDATES = [370612437, 252721481]
WINDOW_DUR = 3  # light curve passed to TRICERATOPS: within 3 durations of transit
N_SIM = 1_000_000  # samples per scenario (the TRICERATOPS default)


def export():
    import env  # noqa: F401
    import numpy as np
    from astropy.io import fits
    import vet as V
    for tic in CANDIDATES:
        fit = json.loads((OUT / f"{tic}_fit.json").read_text())["stellar"]
        P, t0 = fit["P"]["median"], fit["t0"]["median"]
        dur = fit["t14_hr"]["median"] / 24
        c = V.load_vet(tic)
        # 200-s sectors only: one exposure time for the model, and the sharpest transit shape
        m = (c["exptime"] < 300) & (np.abs(V.phase_of(c["time"], P, t0)) < WINDOW_DUR * dur)
        ph = V.phase_of(c["time"][m], P, t0)
        o = np.argsort(ph)
        apertures, sectors = [], []
        for f in sorted(glob.glob(str(ROOT / "cache" / "**" / f"*{tic:016d}-s00*_lc.fits"), recursive=True)):
            with fits.open(f) as h:
                ap = h["APERTURE"]
                rows, cols = np.where((ap.data & 2) > 0)
                apertures.append([[int(ap.header["CRVAL1P"] + cc), int(ap.header["CRVAL2P"] + rr)]
                                  for rr, cc in zip(rows, cols)])
                sectors.append(int(h[0].header["SECTOR"]))
        np.savez(OUT / f"{tic}_triceratops_input.npz", time=ph[o], flux=c["flux"][m][o],
                 flux_err=float(np.median(c["flux_err"][m])), P=P, depth=fit["depth_ppm"]["median"] * 1e-6,
                 sectors=np.array(sectors), apertures=np.array(apertures, dtype=object), allow_pickle=True)
        print(f"TIC {tic}: {m.sum()} points, sectors {sectors}, aperture sizes {[len(a) for a in apertures]}")


def run(n_runs):
    import os
    import time
    os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "cache" / "mpl"))
    import warnings
    warnings.filterwarnings("ignore")
    import numpy as np
    if not hasattr(np, "asscalar"):  # removed in numpy 1.23; triceratops 1.0.8 still calls it
        np.asscalar = lambda a: np.asarray(a).item()
    import triceratops.triceratops as tr
    # triceratops 1.0.8 passes the coordinate positionally; astroquery >= 0.4.8 made it keyword-only
    _get_cutouts = tr.Tesscut.get_cutouts

    from astropy.io import fits
    cache = ROOT / "cache" / "triceratops_cutouts"
    cache.mkdir(parents=True, exist_ok=True)

    def get_cutouts_compat(*args, **kwargs):
        if args:
            kwargs["coordinates"] = args[0]
        c = kwargs["coordinates"]
        f = cache / f"{c.ra.deg:.5f}_{c.dec.deg:.5f}_s{int(kwargs.get('sector', 0))}_{kwargs.get('size')}.fits"
        if f.exists():  # cached from an earlier run (these 200-s cutouts take minutes each)
            return [fits.open(f)]
        for attempt in range(6):  # MAST drops long TESScut connections now and then
            try:
                hdus = _get_cutouts(**kwargs)
                hdus[0].writeto(f, overwrite=True)
                return hdus
            except Exception as e:
                if attempt == 5:
                    raise
                print(f"TESScut retry {attempt + 1}: {type(e).__name__}", flush=True)
                time.sleep(30 * (attempt + 1))

    tr.Tesscut.get_cutouts = get_cutouts_compat
    for tic in CANDIDATES:
        d = np.load(OUT / f"{tic}_triceratops_input.npz", allow_pickle=True)
        sectors = d["sectors"]
        target = tr.target(ID=tic, sectors=sectors, search_radius=10)
        target.calc_depths(tdepth=float(d["depth"]), all_ap_pixels=[np.array(a) for a in d["apertures"]])
        stars = target.stars
        runs = []
        for i in range(n_runs):
            target.calc_probs(time=d["time"], flux_0=d["flux"], flux_err_0=float(d["flux_err"]),
                              P_orb=float(d["P"]), N=N_SIM, exptime=200 / 86400, parallel=True, verbose=0)
            runs.append({"FPP": float(target.FPP), "NFPP": float(target.NFPP),
                         "scenarios": target.probs.groupby("scenario")["prob"].sum().to_dict()})
            print(f"TIC {tic} run {i + 1}: FPP {target.FPP:.4f}  NFPP {target.NFPP:.2e}", flush=True)
        fpp = np.array([r["FPP"] for r in runs])
        nfpp = np.array([r["NFPP"] for r in runs])
        scen = {k: float(np.mean([r["scenarios"].get(k, 0) for r in runs])) for k in runs[0]["scenarios"]}
        # neighbours that could produce the signal at all (needed eclipse depth below 100%)
        cols = [c for c in ("ID", "Tmag", "sep (arcsec)", "tdepth") if c in stars]
        possible = stars[(stars.get("tdepth", 0) > 0) & (stars.get("tdepth", 0) < 1)][cols]
        out = {"tic": tic, "n_runs": n_runs, "N_per_scenario": N_SIM,
               "FPP_mean": float(fpp.mean()), "FPP_std": float(fpp.std()),
               "NFPP_mean": float(nfpp.mean()), "NFPP_std": float(nfpp.std()),
               "scenario_probs_mean": dict(sorted(scen.items(), key=lambda kv: -kv[1])),
               "n_stars_within_10px": int(len(stars) - 1),
               "n_neighbours_able_to_produce_signal": int(len(possible) - 1 if len(possible) else 0),
               "possible_sources": possible.head(20).to_dict("records"),
               "validated": bool(fpp.mean() < 0.015 and nfpp.mean() < 0.001),
               "runs": runs}
        (OUT / f"{tic}_triceratops.json").write_text(json.dumps(out, indent=1, default=str))
        try:
            target.plot_fits(time=d["time"], flux_0=d["flux"], flux_err_0=float(d["flux_err"]),
                             save=True, fname=str(OUT / "plots" / f"{tic}_triceratops"))
        except Exception as e:
            print("plot failed:", e)
        print(f"TIC {tic}: FPP = {fpp.mean():.4f} ± {fpp.std():.4f}, NFPP = {nfpp.mean():.2e} ± {nfpp.std():.1e}; "
              f"top scenarios: " + ", ".join(f"{k} {v:.3f}" for k, v in list(out['scenario_probs_mean'].items())[:5]))


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).parent))
    if sys.argv[1] == "export":
        export()
    else:
        run(int(sys.argv[2]) if len(sys.argv) > 2 else 5)
