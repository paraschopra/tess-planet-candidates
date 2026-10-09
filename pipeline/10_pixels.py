"""Phase 3e: which star is actually dimming? Difference imaging on full-frame-image cutouts.

TESS pixels are 21" wide, so light from neighbours lands in the target's aperture. A faint eclipsing
binary next to the target can mimic a tiny planet. Method:

  1. download a 15x15-pixel cutout of the calibrated full-frame images (TESScut) per sector, subtract sky;
  2. for every transit: mean of frames in transit, and of frames 1-3 durations BEFORE and AFTER (both
     sides required, so a linear drift cancels); difference = (before + after)/2 - in transit;
  3. average per sector. Where light went missing shows as a positive spot. Per-pixel noise is measured
     empirically from the same difference images at 40 off-transit phases;
  4. joint fit over ALL sectors of one PSF-shaped spot (Gaussian, width from the target's own image)
     at a common sky position, with free per-sector amplitude and background. Its significance is the
     chi-square improvement over no spot; its uncertainty comes from bootstrapping transits per sector;
  5. source comparison: the same fit with the spot pinned on the target and, in turn, on each Gaia
     neighbour (positions moved to each sector's date with proper motion). The true source fits best;
     delta-chi2 says how much worse each neighbour does.

Also: depth in growing apertures (flat if the target is the source) and, per neighbour, the minimum
eclipse depth it would need to produce the signal (TESS-band magnitudes from Gaia G and BP-RP, Stassun
et al. 2019; a lower bound, since less of a neighbour's light falls in the target aperture).

Caveats: Gaussian PSF rather than the TESS PRF; TESScut timestamps are barycentred for the FFI centre,
offsets up to ~40 s (frames are 200-1800 s, the dips ~1.4 h).

Writes results/phase3/<TIC>_pixels.json and results/phase3/plots/<TIC>_pixels.png.
Usage: python pipeline/10_pixels.py [TIC ...]
"""
import json
import sys
import time

import env  # noqa: F401

import astropy.units as u
import lightkurve as lk
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astroquery.gaia import Gaia
from scipy.optimize import least_squares
from scipy.stats import chi2 as chi2_dist

import vet as V
from env import CACHE, RESULTS

CUTOUT = 15
N_BOOT = 200
FIT_RADIUS_PX = 4.5
NEIGHBOUR_RADIUS_ARCSEC = 150
NEIGHBOUR_MAX_DMAG = 9.5  # 9.5 mag fainter = 1/6300 of the light: cannot make 200 ppm even if fully eclipsed
COMPARE_RADIUS_ARCSEC = 90  # neighbours tested as the possible source
PLOT_MAX_DT = 7.0  # only circle neighbours that could make the dip with a <~30% eclipse


def tess_mag(g, bp_rp):
    """TESS magnitude from Gaia G and BP-RP (Stassun et al. 2019, TIC v8)."""
    if not np.isfinite(bp_rp):
        return g - 0.43
    x = np.clip(bp_rp, -0.2, 3.0)
    return g - 0.00522555 * x ** 3 + 0.0891337 * x ** 2 - 0.633923 * x + 0.0324473


def gaia_info(tic):
    tg = pd.read_csv(RESULTS / "targets_ranked.csv").set_index("TIC").loc[tic]
    gid = int(tg["GAIA"])
    q = f"""SELECT source_id, ra, dec, pmra, pmdec, phot_g_mean_mag, bp_rp FROM gaiadr3.gaia_source
            WHERE source_id = {gid}"""
    me = Gaia.launch_job(q).get_results().to_pandas().iloc[0]
    me["tic_ra"], me["tic_dec"] = tg["ra"], tg["dec"]
    me["tmag"] = tess_mag(me.phot_g_mean_mag, me.bp_rp)
    q = f"""SELECT source_id, ra, dec, pmra, pmdec, phot_g_mean_mag, bp_rp, ruwe FROM gaiadr3.gaia_source
            WHERE 1 = CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {me.ra}, {me.dec}, {NEIGHBOUR_RADIUS_ARCSEC / 3600}))
            AND phot_g_mean_mag < {me.phot_g_mean_mag + NEIGHBOUR_MAX_DMAG}"""
    nb = Gaia.launch_job(q).get_results().to_pandas()
    nb = nb[nb.source_id != gid].copy()
    nb["tmag"] = [tess_mag(g, c) for g, c in zip(nb.phot_g_mean_mag, nb.bp_rp)]
    c0 = SkyCoord(me.ra, me.dec, unit="deg")
    nb["sep_arcsec"] = c0.separation(SkyCoord(nb.ra.values, nb.dec.values, unit="deg")).arcsec
    return me, nb.sort_values("sep_arcsec").reset_index(drop=True)


def at_epoch(row, jd):
    """Gaia DR3 (J2016.0) position moved to a given date with proper motion."""
    dt = (jd - Time(2016.0, format="jyear").jd) / 365.25
    pmra = row["pmra"] if np.isfinite(row["pmra"]) else 0.0
    pmdec = row["pmdec"] if np.isfinite(row["pmdec"]) else 0.0
    dec = row["dec"] + pmdec * dt / 3.6e6
    ra = row["ra"] + pmra * dt / 3.6e6 / np.cos(np.radians(dec))
    return ra, dec


def gauss(yy, xx, x0, y0, sigma):
    return np.exp(-((xx - x0) ** 2 + (yy - y0) ** 2) / (2 * sigma ** 2))


def cutout(tic, sector, me, retries=5):
    """TESScut cutout for one sector, reusing a previous download (MAST drops connections)."""
    hits = sorted((CACHE / "tesscut").rglob(f"tess-s{int(sector):04d}-*_{me.tic_ra:.6f}_{me.tic_dec:.6f}_*.fits"))
    if hits:
        return lk.read(str(hits[0]), quality_bitmask="hard")
    for attempt in range(retries):
        try:
            sr = lk.search_tesscut(f"TIC {tic}", sector=int(sector))
            return sr.download(cutout_size=CUTOUT, download_dir=str(CACHE / "tesscut"), quality_bitmask="hard")
        except Exception:
            if attempt == retries - 1:
                raise
            time.sleep(30 * (attempt + 1))


def sector_cube(tic, sector, me):
    tpf = cutout(tic, sector, me)
    t = np.asarray(tpf.time.value, float)
    cube = np.asarray(tpf.flux.value, float)
    ok = np.isfinite(t) & np.all(np.isfinite(cube.reshape(len(t), -1)), axis=1)
    t, cube = t[ok], cube[ok]
    # sky: median of the faintest 30% of pixels in each frame
    med = np.nanmedian(cube, axis=0)
    faint = med < np.nanpercentile(med, 30)
    cube = cube - np.median(cube[:, faint], axis=1)[:, None, None]
    return t, cube, tpf


def diff_images(t, cube, P, t0, dur):
    """Per-transit (before+after)/2 - in-transit images; both sides of each transit required."""
    ph = V.phase_of(t, P, t0)
    ep = np.round((t - t0) / P).astype(int)
    intr = np.abs(ph) < dur / 2
    before = (ph < -dur) & (ph > -3 * dur)
    after = (ph > dur) & (ph < 3 * dur)
    imgs = []
    for e in np.unique(ep[intr]):
        a, b1, b2 = intr & (ep == e), before & (ep == e), after & (ep == e)
        if a.sum() >= 1 and b1.sum() >= a.sum() and b2.sum() >= a.sum():
            imgs.append(0.5 * (cube[b1].mean(0) + cube[b2].mean(0)) - cube[a].mean(0))
    return np.array(imgs)


def load_sector(tic, sector, P, t0, dur, me, nb):
    t, cube, tpf = sector_cube(tic, sector, me)
    jd = np.median(t) + 2457000
    wcs = tpf.wcs
    tc = SkyCoord(*at_epoch(me, jd), unit="deg")
    tx, ty = map(float, wcs.world_to_pixel(tc))
    # pixel shift per arcsec east / north at the target (local linear WCS)
    jac = np.zeros((2, 2))
    for k, (pa, name) in enumerate(((90, "E"), (0, "N"))):
        x1, y1 = map(float, wcs.world_to_pixel(tc.directional_offset_by(pa * u.deg, 10 * u.arcsec)))
        jac[:, k] = [(x1 - tx) / 10, (y1 - ty) / 10]
    imgs = diff_images(t, cube, P, t0, dur)
    nulls = [diff_images(t, cube, P, t0 + o * P, dur) for o in V.NULL_OFFSETS]
    nulls = np.array([n.mean(0) for n in nulls if len(n) >= max(2, len(imgs) // 2)])
    noise = 1.4826 * np.median(np.abs(nulls - np.median(nulls, 0)), axis=0)
    direct = np.median(cube, axis=0)
    yy, xx = np.mgrid[:CUTOUT, :CUTOUT]
    near = np.hypot(xx - tx, yy - ty) < FIT_RADIUS_PX
    # PSF width from the target's own image
    r = least_squares(lambda p: (p[0] * gauss(yy, xx, p[1], p[2], p[3]) + p[4] - direct)[near],
                      [direct[near].max(), tx, ty, 0.8, 0], bounds=([0, tx - 2, ty - 2, 0.3, -np.inf],
                                                                    [np.inf, tx + 2, ty + 2, 3, np.inf]))
    sigma = float(r.x[3])
    # every neighbour's position in this sector's pixels, at this sector's date
    nbpx = np.array([list(map(float, wcs.world_to_pixel(SkyCoord(*at_epoch(n, jd), unit="deg")))) for _, n in nb.iterrows()]) \
        if len(nb) else np.zeros((0, 2))
    # depth in growing circular apertures around the target
    rr = np.hypot(xx - tx, yy - ty)
    diff = imgs.mean(0)
    ap = {rad: float(diff[rr <= rad].sum() / direct[rr <= rad].sum() * 1e6) for rad in (1.0, 1.5, 2.0, 2.5, 3.0)}
    return {"sector": int(sector), "imgs": imgs, "diff": diff, "noise": noise, "direct": direct, "near": near,
            "target_px": (tx, ty), "jac": jac, "sigma": sigma, "nbpx": nbpx, "aperture_depth_ppm": ap,
            "yy": yy, "xx": xx, "nb_dT": (nb["tmag"] - me["tmag"]).to_numpy()}


def design(s, x0, y0):
    """Per-sector columns: PSF at (x0, y0) and a constant background, over the fit pixels."""
    return np.column_stack([gauss(s["yy"], s["xx"], x0, y0, s["sigma"])[s["near"]], np.ones(s["near"].sum())])


def chi2_at(secs, diffs, positions):
    """Chi-square with the spot pinned at given per-sector pixel positions (amplitudes, bg linear)."""
    total, amps = 0.0, []
    for s, d, (x0, y0) in zip(secs, diffs, positions):
        w = 1 / s["noise"][s["near"]]
        A = design(s, x0, y0) * w[:, None]
        b = d[s["near"]] * w
        coef, *_ = np.linalg.lstsq(A, b, rcond=None)
        total += float(np.sum((A @ coef - b) ** 2))
        amps.append(float(coef[0]))
    return total, amps


def chi2_none(secs, diffs):
    total = 0.0
    for s, d in zip(secs, diffs):
        w = 1 / s["noise"][s["near"]]
        b = d[s["near"]] * w
        total += float(np.sum((b - w * np.sum(b * w) / np.sum(w ** 2)) ** 2))
    return total


def offset_px(s, dEN):
    return s["target_px"][0] + s["jac"][0] @ dEN, s["target_px"][1] + s["jac"][1] @ dEN


def joint_fit(secs, diffs):
    """Common sky offset (arcsec E, N from the target) of the missing-light spot across all sectors."""
    f = lambda dEN: chi2_at(secs, diffs, [offset_px(s, dEN) for s in secs])[0]
    grid = np.arange(-60, 61, 6.0)
    best = min(((f(np.array([e, n])), e, n) for e in grid for n in grid))
    r = least_squares(lambda p: [np.sqrt(f(p))], [best[1], best[2]], diff_step=1e-3)
    return r.x, f(r.x)


def run(tic):
    P, t0, dur = V.ephemeris(tic)
    me, nb = gaia_info(tic)
    lc = V.load_vet(tic, P, t0, dur)
    depth, _ = V.depth_emp(lc, P, t0, dur)
    secs = [load_sector(tic, s, P, t0, dur, me, nb) for s in np.unique(lc["sector"])]
    secs = [s for s in secs if len(s["imgs"]) >= 2]
    diffs = [s["diff"] for s in secs]
    dEN, chi_best = joint_fit(secs, diffs)
    chi_target, amps_target = chi2_at(secs, diffs, [s["target_px"] for s in secs])
    chi0 = chi2_none(secs, diffs)
    # bootstrap transits within each sector for the position uncertainty
    rng = np.random.default_rng(tic % 997)
    boots = []
    for _ in range(N_BOOT):
        bd = [s["imgs"][rng.integers(0, len(s["imgs"]), len(s["imgs"]))].mean(0) for s in secs]
        boots.append(joint_fit(secs, bd)[0])
    boots = np.array(boots)
    cov = np.cov(boots.T)
    off_sigma = float(np.sqrt(dEN @ np.linalg.inv(cov) @ dEN))
    # each nearby star in turn as the source
    near_nb = nb[nb.sep_arcsec < COMPARE_RADIUS_ARCSEC].copy()
    c0 = SkyCoord(me.ra, me.dec, unit="deg")  # offsets at J2016 for the plot
    de, dn = c0.spherical_offsets_to(SkyCoord(near_nb.ra.values, near_nb.dec.values, unit="deg"))
    near_nb["dE"], near_nb["dN"] = de.to(u.arcsec).value, dn.to(u.arcsec).value
    rows = []
    for i, n in near_nb.iterrows():
        chi_n, amps_n = chi2_at(secs, diffs, [tuple(s["nbpx"][i]) for s in secs])
        dT = n.tmag - me.tmag
        rows.append({"source_id": int(n.source_id), "sep_arcsec": float(n.sep_arcsec), "G": float(n.phot_g_mean_mag),
                     "Tmag": float(n.tmag), "dT": float(dT), "delta_chi2_vs_target": chi_n - chi_target,
                     "mean_amp_sign": float(np.sign(np.mean(amps_n))),
                     "min_eclipse_depth_needed": float(min(1.0, depth * 10 ** (0.4 * dT))),
                     "impossible": bool(depth * 10 ** (0.4 * dT) > 1)})
    nbt = pd.DataFrame(rows)
    out = {"tic": tic, "depth_ppm": depth * 1e6, "target_tmag_from_gaia": float(me.tmag),
           "n_sectors": len(secs), "n_transits": int(sum(len(s["imgs"]) for s in secs)),
           "spot_detection_delta_chi2": chi0 - chi_target,
           "spot_detection_p": float(chi2_dist.sf(chi0 - chi_target, len(secs))),
           "target_amplitudes_positive": int(sum(a > 0 for a in amps_target)),
           "joint_offset_arcsec_EN": dEN.tolist(), "joint_offset_cov": cov.tolist(),
           "joint_offset_sep_arcsec": float(np.hypot(*dEN)), "joint_offset_sigma": off_sigma,
           "chi2_free_position_minus_target": chi_best - chi_target,
           "sectors": [{"sector": s["sector"], "n_transits": len(s["imgs"]), "psf_sigma_px": s["sigma"],
                        "aperture_depth_ppm": s["aperture_depth_ppm"],
                        "amp_at_target": a} for s, a in zip(secs, amps_target)],
           "neighbours_tested": nbt.to_dict("records"),
           "neighbours_within_150arcsec": int(len(nb))}
    (V.OUT / f"{tic}_pixels.json").write_text(json.dumps(out, indent=1))
    plot(tic, secs, out, dEN, boots, near_nb)
    return out


def plot(tic, secs, out, dEN, boots, near_nb):
    n = len(secs)
    fig = plt.figure(figsize=(2.4 * n + 5, 7.5))
    gs = fig.add_gridspec(3, n + 2)
    for j, s in enumerate(secs):
        for i, key in enumerate(("direct", "snr")):
            ax = fig.add_subplot(gs[i, j])
            if key == "direct":
                im = s["direct"]
                ax.imshow(np.log10(np.clip(im, np.percentile(im, 5), None)), origin="lower", cmap="gray")
            else:
                ax.imshow(s["diff"] / s["noise"], origin="lower", cmap="RdBu_r", vmin=-4, vmax=4)
            for (x, y), dT in zip(s["nbpx"], s["nb_dT"]):
                if -0.5 < x < CUTOUT - 0.5 and -0.5 < y < CUTOUT - 0.5 and dT < PLOT_MAX_DT:
                    ax.plot(x, y, "o", mfc="none", mec="orange", ms=4)
            ax.plot(*s["target_px"], "+", color="lime", ms=10, mew=1.5)
            if key == "snr":
                ax.plot(*offset_px(s, dEN), "x", color="k", ms=8, mew=1.5)
            ax.set_title(f"S{s['sector']} {'image' if key == 'direct' else 'missing light / noise'}", fontsize=8)
            ax.set_xticks([])
            ax.set_yticks([])
        ax = fig.add_subplot(gs[2, j])
        ap = s["aperture_depth_ppm"]
        ax.plot(list(ap.keys()), list(ap.values()), "o-")
        ax.axhline(out["depth_ppm"], color="C1", ls="--", lw=0.8)
        ax.set_xlabel("aperture radius (px)", fontsize=7)
        ax.set_title(f"S{s['sector']} depth vs aperture", fontsize=8)
    ax = fig.add_subplot(gs[:2, n:])
    ax.plot(boots[:, 0], boots[:, 1], ".", color="0.6", ms=2, label="bootstrap")
    ax.plot(*dEN, "x", color="k", ms=10, mew=2, label="missing-light position")
    ax.plot(0, 0, "+", color="lime", ms=14, mew=2, label="target")
    for (_, g), nrow in zip(near_nb.iterrows(), out["neighbours_tested"]):
        if nrow["dT"] >= PLOT_MAX_DT and nrow["delta_chi2_vs_target"] > 0:
            continue
        ax.plot(g["dE"], g["dN"], "o", mfc="none", mec="orange", ms=max(3, 12 - nrow["dT"]))
        ax.annotate(f"Δχ² {nrow['delta_chi2_vs_target']:+.0f}", (g["dE"], g["dN"]), fontsize=6, color="C1")
    ax.set_xlim(70, -70)
    ax.set_ylim(-70, 70)
    ax.set_aspect("equal")
    ax.set_xlabel("arcsec east of target")
    ax.set_ylabel("arcsec north")
    ax.legend(fontsize=7, loc="lower left")
    ax.set_title("where the light went missing (all sectors jointly);\norange = Gaia neighbours (ΔT<7 or fitting better), "
                 "Δχ² = how much worse if it were the source", fontsize=8)
    fig.suptitle(f"TIC {tic}: missing-light spot detected at Δχ² = {out['spot_detection_delta_chi2']:.0f} "
                 f"({out['n_transits']} transits, {out['n_sectors']} sectors); offset from target "
                 f"{out['joint_offset_sep_arcsec']:.1f}\" ({out['joint_offset_sigma']:.1f}σ)", fontsize=10)
    fig.tight_layout()
    fig.savefig(V.OUT / "plots" / f"{tic}_pixels.png", dpi=80)
    plt.close(fig)


def main():
    tics = [int(a) for a in sys.argv[1:]] or V.CANDIDATES
    (V.OUT / "plots").mkdir(parents=True, exist_ok=True)
    for tic in tics:
        o = run(tic)
        print(f"TIC {tic}: depth {o['depth_ppm']:.0f} ppm; spot at target detected with Δχ² = "
              f"{o['spot_detection_delta_chi2']:.1f} (p = {o['spot_detection_p']:.1e}), "
              f"{o['target_amplitudes_positive']}/{o['n_sectors']} sectors positive")
        print(f"  free-position spot offset {o['joint_offset_sep_arcsec']:.1f}\" ({o['joint_offset_sigma']:.1f}σ), "
              f"E/N {o['joint_offset_arcsec_EN'][0]:+.1f}, {o['joint_offset_arcsec_EN'][1]:+.1f}\"")
        for s in o["sectors"]:
            print(f"  S{s['sector']}: {s['n_transits']} tr, amp {s['amp_at_target']:+.1f}, aperture depths "
                  f"{[round(v) for v in s['aperture_depth_ppm'].values()]}")
        nbt = pd.DataFrame(o["neighbours_tested"])
        print(f"  {o['neighbours_within_150arcsec']} Gaia neighbours within 150\"; tested as source (within 90\"):")
        if len(nbt):
            print(nbt.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
