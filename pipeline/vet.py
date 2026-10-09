"""Shared helpers for phase 3 vetting of individual candidates.

Differences from the phase 2 search cleaning:
  - the trend is fitted with the candidate's transits masked out, so the dip no longer pulls its own
    baseline down (phase 2 lost ~duration/window of the depth this way);
  - data are kept at native cadence (no binning) for model fitting.

Depths at a fixed ephemeris are quoted with an empirical error: the same box placed at 40 phases away
from the transit and secondary eclipse, so correlated noise is included (as in phase 2 triage).
"""
import numpy as np
import pandas as pd
from wotan import flatten

import lcsearch as L
from env import RESULTS

# candidate TIC -> which search signal it is (350883383's first signal was a glitch at 23.8 d)
CANDIDATES = {370612437: 1, 252721481: 1, 350883383: 2}
OUT = RESULTS / "phase3"
NULL_OFFSETS = np.concatenate([np.linspace(0.1, 0.4, 20), np.linspace(0.6, 0.9, 20)])


def ephemeris(tic):
    """Phase 2 search ephemeris (P, t0, duration in days) for a candidate's signal."""
    t = pd.read_csv(RESULTS / "phase2" / "triage.csv")
    r = t[(t.tic == tic) & (t.signal == CANDIDATES.get(tic, 1))].iloc[0]
    return float(r.period), float(r.t0), float(r.duration_hr) / 24


def phase_of(t, P, t0):
    return ((t - t0 + 0.5 * P) % P) - 0.5 * P


def masked_clean(lc, P, t0, dur, window=L.DETREND_WINDOW_D, mask_factor=2.0):
    """Detrend each sector with in-transit points excluded from the trend, then clip as in L.clean."""
    t, f, sec = lc["time"], lc["flux"], lc["sector"]
    intr = np.abs(phase_of(t, P, t0)) < mask_factor * dur / 2
    trend = np.full_like(f, np.nan)
    for s in np.unique(sec):
        m = sec == s
        o = m & ~intr
        _, tr = flatten(t[o], f[o], method="biweight", window_length=window,
                        break_tolerance=0.5, edge_cutoff=0.0, return_trend=True)
        ok = np.isfinite(tr)
        # trend under the masked transits comes from its neighbours on either side
        trend[m] = np.interp(t[m], t[o][ok], tr[ok], left=np.nan, right=np.nan)
    flat = f / trend
    keep = np.isfinite(flat)
    err = np.zeros_like(f)
    for s in np.unique(sec):
        m = (sec == s) & keep
        oot = m & ~intr
        sd = L.robust_std(flat[oot])
        resid = flat[m] - 1
        # flares clipped everywhere; dropouts only outside transit, so the transit itself is untouched
        bad = (resid > L.FLARE_CLIP_SIGMA * sd) | ((resid < -L.DROPOUT_CLIP_SIGMA * sd) & ~intr[m])
        idx = np.where(m)[0]
        keep[idx[bad]] = False
        err[m] = sd
    return {"time": t[keep], "flux": flat[keep], "flux_err": err[keep], "sector": sec[keep],
            "exptime": lc["exptime"][keep]}


def load_vet(tic, P=None, t0=None, dur=None):
    """Raw light curve -> bad sectors rejected -> transit-masked detrending, native cadence."""
    if P is None:
        P, t0, dur = ephemeris(tic)
    lc, _, _ = L.reject_bad_sectors(L.load_star(tic))
    return masked_clean(lc, P, t0, dur)


def box_depth(c, sel):
    e = c["flux_err"][sel]
    if sel.sum() <= 3:
        return np.nan
    w = 1 / e ** 2
    return float(np.sum((1 - c["flux"][sel]) * w) / w.sum())


def depth_emp(c, P, t0, dur, sel=None):
    """Box depth at a fixed ephemeris and its empirical error (box at 40 off-transit phases)."""
    sel = np.ones(len(c["time"]), bool) if sel is None else sel
    t = c["time"]
    d = box_depth(c, sel & (np.abs(phase_of(t, P, t0)) < dur / 2))
    null = np.array([box_depth(c, sel & (np.abs(phase_of(t, P, t0 + o * P)) < dur / 2))
                     for o in NULL_OFFSETS])
    null = null[np.isfinite(null)]
    sig = 1.4826 * np.median(np.abs(null - np.median(null))) if len(null) > 10 else np.nan
    return d, float(sig)


def seasons(c, gap_d=100):
    """Group sectors into observing seasons (TESS revisits a field roughly yearly)."""
    secs = sorted(np.unique(c["sector"]).astype(int))
    start = {s: c["time"][c["sector"] == s].min() for s in secs}
    groups, cur = [], [secs[0]]
    for a, b in zip(secs, secs[1:]):
        if start[b] - start[a] > gap_d:
            groups.append(cur)
            cur = []
        cur.append(b)
    groups.append(cur)
    return groups


def btjd_to_date(t):
    from astropy.time import Time
    return Time(t + 2457000, format="jd").to_value("iso", subfmt="date")


def depth_emp_matched(lc_raw, P, t0, dur):
    """Depth and empirical error with every control box processed exactly like the transit box.

    depth_emp() measures control boxes on data detrended with only the real transits masked, so the
    controls helped fit their own trend and come out slightly attenuated (noise ~10% low; caught by
    the Codex audit). Here each control phase gets its own masked detrend.
    """
    c = masked_clean(lc_raw, P, t0, dur)
    d = box_depth(c, np.abs(phase_of(c["time"], P, t0)) < dur / 2)
    null = []
    for o in NULL_OFFSETS:
        tc = t0 + o * P
        co = masked_clean(lc_raw, P, tc, dur)
        null.append(box_depth(co, np.abs(phase_of(co["time"], P, tc)) < dur / 2))
    null = np.array(null)
    null = null[np.isfinite(null)]
    sig = 1.4826 * np.median(np.abs(null - np.median(null)))
    return d, float(sig), null
