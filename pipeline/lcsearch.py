"""Light-curve cleaning and transit search, shared by the phase 2 search and later vetting.

Pipeline for one star:
  load_star -> reject_bad_sectors -> clean (detrend, clip flares) -> bin -> bls_search (iterative)
"""
import numpy as np
from astropy.timeseries import BoxLeastSquares
from wotan import flatten

from env import LIGHTCURVES

DETREND_WINDOW_D = 0.75  # biweight window. In-transit points pull the trend down by ~duration/window
# (~18% for a 2.9 h dip), but wider windows leave stellar variability: injection tests
# (8 stars x 2, Rp=1.5) recovered 3/16, 10/16, 8/16, 5/16 for windows 0.5, 0.75, 1.0, 1.5 d.
MAX_SECTOR_NOISE_1HR = 0.01  # sectors noisier than 1% per hour can't show a ~500 ppm dip
SECTOR_NOISE_RATIO = 3.0  # ...or 3x noisier than the star's typical sector
FLARE_CLIP_SIGMA = 3.0  # upward outliers only; downward points could be transits
DROPOUT_CLIP_SIGMA = 10.0  # single deep points are artifacts; a small planet is a few sigma
BIN_S = 600  # bin fast cadences to 10 min for search speed
PERIOD_MIN_D = 0.5
PERIOD_MAX_D = 25.0
MIN_TRANSITS = 3
OVERSAMPLE = 3  # transit-time drift across the baseline stays below duration/OVERSAMPLE
SDE_THRESHOLD = 7.0
SNR_THRESHOLD = 6.0
MAX_SIGNALS = 3


def robust_std(x):
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def ptp_noise(flux):
    """Point-to-point scatter: insensitive to slow stellar variability."""
    return 1.4826 * np.median(np.abs(np.diff(flux))) / np.sqrt(2)


def load_star(tic):
    d = np.load(LIGHTCURVES / f"{tic}.npz")
    return {k: d[k] for k in d.files}


def reject_bad_sectors(lc):
    """Drop sectors whose noise (scaled to 1 hr) is absurd or far above the star's typical sector."""
    sectors = np.unique(lc["sector"])
    noise = {}
    for s in sectors:
        m = lc["sector"] == s
        noise[s] = ptp_noise(lc["flux"][m]) * np.sqrt(np.median(lc["exptime"][m]) / 3600)
    typical = np.median(list(noise.values()))
    bad = [s for s in sectors
           if noise[s] > MAX_SECTOR_NOISE_1HR or (len(sectors) >= 3 and noise[s] > SECTOR_NOISE_RATIO * typical)]
    keep = ~np.isin(lc["sector"], bad)
    out = {k: (v[keep] if isinstance(v, np.ndarray) and v.shape == lc["time"].shape else v)
           for k, v in lc.items()}
    return out, [int(s) for s in bad], {int(s): float(n) for s, n in noise.items()}


def clean(lc, window=DETREND_WINDOW_D):
    """Detrend each sector, clip flares and dropouts, set per-sector errors from measured scatter."""
    t, f, sec = lc["time"], lc["flux"], lc["sector"]
    flat = np.full_like(f, np.nan)
    for s in np.unique(sec):
        m = sec == s
        flat[m] = flatten(t[m], f[m], method="biweight", window_length=window,
                          break_tolerance=0.5, edge_cutoff=0.0)
    keep = np.isfinite(flat)
    err = np.zeros_like(f)
    for s in np.unique(sec):
        m = (sec == s) & keep
        sd = robust_std(flat[m])
        resid = flat[m] - 1
        good = (resid < FLARE_CLIP_SIGMA * sd) & (resid > -DROPOUT_CLIP_SIGMA * sd)
        idx = np.where(m)[0]
        keep[idx[~good]] = False
        err[m] = sd
    return {"time": t[keep], "flux": flat[keep], "flux_err": err[keep], "sector": sec[keep],
            "exptime": lc["exptime"][keep]}


def bin_lc(lc, bin_s=BIN_S):
    """Bin cadences shorter than bin_s; longer cadences pass through unchanged."""
    fast = lc["exptime"] < bin_s
    out = {k: [v[~fast]] for k, v in lc.items()}
    if fast.any():
        t, f, e, s = (lc[k][fast] for k in ("time", "flux", "flux_err", "sector"))
        width = bin_s / 86400
        key = np.floor((t - t.min()) / width).astype(np.int64)
        # keep bins within a sector so per-sector errors stay meaningful
        key = key * 1000 + s
        uniq, inv, counts = np.unique(key, return_inverse=True, return_counts=True)
        w = 1 / e ** 2
        out["time"].append(np.bincount(inv, t * w) / np.bincount(inv, w))
        out["flux"].append(np.bincount(inv, f * w) / np.bincount(inv, w))
        out["flux_err"].append(1 / np.sqrt(np.bincount(inv, w)))
        out["sector"].append(np.bincount(inv, s) / counts)
        out["exptime"].append(np.full(len(uniq), float(bin_s)))
    out = {k: np.concatenate(v) for k, v in out.items()}
    order = np.argsort(out["time"])
    return {k: v[order] for k, v in out.items()}


def duration_days(period, rho_sun):
    """Central-transit duration for a circular orbit around a star of density rho (solar units)."""
    return 13.0 / 24 * (period / 365.25) ** (1 / 3) * rho_sun ** (-1 / 3)


def period_grid(baseline, rho_sun, pmin=PERIOD_MIN_D, pmax=PERIOD_MAX_D, oversample=OVERSAMPLE):
    """Frequency grid fine enough that transits stay aligned across the whole baseline.

    A frequency error df shifts the last transit by ~baseline*df*P; keep that below D/oversample,
    using a conservative (grazing, 0.5x) duration.
    """
    pmax = min(pmax, baseline / MIN_TRANSITS)
    f, fmin = 1 / pmin, 1 / pmax
    freqs = []
    while f > fmin:
        freqs.append(f)
        p = 1 / f
        f -= 0.5 * duration_days(p, rho_sun) / (oversample * baseline * p)
    return np.sort(1 / np.array(freqs))


def duration_grid(rho_sun, pmin=PERIOD_MIN_D, pmax=PERIOD_MAX_D, n=7):
    lo = 0.35 * duration_days(pmin, rho_sun)
    hi = min(1.3 * duration_days(pmax, rho_sun), 0.4 * pmin)
    return np.geomspace(max(lo, 0.02), hi, n)


def sde(power, periods):
    """Signal detection efficiency after removing the slow trend of BLS power with period."""
    logp = np.log(periods)
    edges = np.linspace(logp.min(), logp.max(), 60)
    idx = np.clip(np.digitize(logp, edges) - 1, 0, len(edges) - 2)
    med = np.array([np.median(power[idx == i]) if np.any(idx == i) else np.nan
                    for i in range(len(edges) - 1)])
    centers = 0.5 * (edges[1:] + edges[:-1])
    ok = np.isfinite(med)
    trend = np.interp(logp, centers[ok], med[ok])
    p = power - trend
    return (p - np.median(p)) / np.std(p), p


def transit_mask(t, period, t0, duration, factor=1.5):
    phase = ((t - t0 + 0.5 * period) % period) - 0.5 * period
    return np.abs(phase) < factor * duration / 2


def signal_stats(bls, res, k, n_periods):
    """Signal dict for entry k of a BLS result."""
    P, t0, dur = res.period[k], res.transit_time[k], res.duration[k]
    stats = bls.compute_stats(P, dur, t0)
    return {
        "period": float(P), "t0": float(t0), "duration_hr": float(dur * 24),
        "depth_ppm": float(res.depth[k] * 1e6), "depth_err_ppm": float(res.depth_err[k] * 1e6),
        "snr": float(res.depth[k] / res.depth_err[k]),
        "n_transits": int(np.sum(stats["per_transit_count"] > 0)),
        "depth_odd_ppm": float(stats["depth_odd"][0] * 1e6),
        "depth_odd_err_ppm": float(stats["depth_odd"][1] * 1e6),
        "depth_even_ppm": float(stats["depth_even"][0] * 1e6),
        "depth_even_err_ppm": float(stats["depth_even"][1] * 1e6),
        "depth_half_ppm": float(stats["depth_half"][0] * 1e6),
        "harmonic_dll": float(stats["harmonic_delta_log_likelihood"]),
        "n_periods": n_periods,
    }


def bls_search_at(lc, rho_sun, period):
    """BLS at one known period (best duration and epoch), e.g. to test vetting on injected signals."""
    bls = BoxLeastSquares(lc["time"], lc["flux"], lc["flux_err"])
    res = bls.power([period], duration_grid(rho_sun), objective="snr")
    return {"signal": 1, **signal_stats(bls, res, 0, 1), "sde": np.nan}


def bls_search(lc, rho_sun, max_signals=MAX_SIGNALS):
    """Iterative BLS: find the best signal, record it, mask its transits, repeat."""
    t, y, dy = lc["time"], lc["flux"], lc["flux_err"]
    baseline = t.max() - t.min()
    periods = period_grid(baseline, rho_sun)
    durations = duration_grid(rho_sun, pmax=periods.max())
    signals = []
    mask = np.zeros(len(t), bool)
    for i in range(max_signals):
        if (~mask).sum() < 500:  # a very short-period signal's mask can swallow the light curve
            break
        tt, yy, ee = t[~mask], y[~mask], dy[~mask]
        bls = BoxLeastSquares(tt, yy, ee)
        res = bls.power(periods, durations, objective="snr")
        s, _ = sde(np.asarray(res.power), periods)
        k = int(np.nanargmax(s))
        P, t0, dur = res.period[k], res.transit_time[k], res.duration[k]
        stats = bls.compute_stats(P, dur, t0)
        n_with_data = int(np.sum(stats["per_transit_count"] > 0))
        sig = {
            "signal": i + 1, "period": float(P), "t0": float(t0), "duration_hr": float(dur * 24),
            "depth_ppm": float(res.depth[k] * 1e6), "depth_err_ppm": float(res.depth_err[k] * 1e6),
            "snr": float(res.depth[k] / res.depth_err[k]), "sde": float(s[k]),
            "n_transits": n_with_data,
            "depth_odd_ppm": float(stats["depth_odd"][0] * 1e6),
            "depth_odd_err_ppm": float(stats["depth_odd"][1] * 1e6),
            "depth_even_ppm": float(stats["depth_even"][0] * 1e6),
            "depth_even_err_ppm": float(stats["depth_even"][1] * 1e6),
            "depth_half_ppm": float(stats["depth_half"][0] * 1e6),
            "harmonic_dll": float(stats["harmonic_delta_log_likelihood"]),
            "n_periods": len(periods),
        }
        signals.append(sig)
        if sig["sde"] < SDE_THRESHOLD or sig["snr"] < SNR_THRESHOLD or n_with_data < MIN_TRANSITS:
            break
        mask |= transit_mask(t, P, t0, dur, factor=2.0)
    return signals, (periods, s)
