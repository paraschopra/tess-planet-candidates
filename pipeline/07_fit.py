"""Phase 3b: physical transit model fit (MCMC) for each candidate.

Fits a limb-darkened transit (batman) to the transit-masked-detrended light curve at native cadence,
integrating the model over each exposure (30-min frames in 2019 smear a 1-hour dip badly).

Parameters: t0, P, k = Rp/R*, impact parameter b, stellar density rho (g/cc).
Per-point errors are inflated by a red-noise factor beta so that a transit-length box has the
empirical (matched-processing) error, and the TIC density prior is Gaussian in rho (Jacobian included).
Two fits:
  free     - rho has a wide log-uniform prior: the transit shape alone says how dense the star must be.
             If that disagrees with the catalog star, the dip may come from a different star.
  stellar  - rho has the TIC density as a Gaussian prior: the physically sensible fit, gives b and Rp.

PDCSAP flux from TESS-SPOC is already corrected for light from neighbours (CROWDSAP), so no extra
dilution correction is applied to the depth.

Writes results/phase3/<TIC>_fit.json and results/phase3/plots/<TIC>_fit.png.
Usage: python pipeline/07_fit.py [TIC ...]
"""
import json
import sys

import env  # noqa: F401
import setuptools  # noqa: F401  provides the distutils shim batman needs on Python 3.12

import batman
import emcee
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import vet as V
from env import RESULTS

G_CGS = 6.674e-8
R_EARTH_PER_R_SUN = 0.009168
# Quadratic limb darkening for a ~4350 K, logg 4.65 dwarf in the TESS band (Claret 2017 tables).
LD = [0.50, 0.17]
WINDOW_DUR = 4  # fit points within 4 durations of transit centre; the rest are flat
SUPERSAMPLE = {1800.0: 15, 600.0: 5, 200.0: 3, 120.0: 1}
N_WALKERS, N_BURN, N_STEPS = 40, 4000, 12000  # autocorrelation ~200 steps: >50 tau


def a_over_rs(rho_cgs, P_d):
    return (G_CGS * rho_cgs * (P_d * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)


class Model:
    def __init__(self, c):
        self.groups = []
        for exp in np.unique(c["exptime"]):
            m = c["exptime"] == exp
            params = batman.TransitParams()
            params.t0, params.per, params.rp, params.a, params.inc = 0, 1, 0.02, 15, 89
            params.ecc, params.w, params.limb_dark, params.u = 0, 90, "quadratic", LD
            ss = SUPERSAMPLE.get(float(exp), max(1, int(exp // 120)))
            model = batman.TransitModel(params, c["time"][m], supersample_factor=ss, exp_time=exp / 86400)
            self.groups.append((m, params, model))
        self.n = len(c["time"])

    def __call__(self, t0, P, k, b, rho):
        out = np.ones(self.n)
        a = a_over_rs(rho, P)
        if b >= a:
            return None
        for m, params, model in self.groups:
            params.t0, params.per, params.rp, params.a = t0, P, k, a
            params.inc = np.degrees(np.arccos(b / a))
            out[m] = model.light_curve(params)
        return out


def log_prob(theta, model, y, e, rho_prior):
    t0, P, k, b, log_rho = theta
    if not (0.001 < k < 0.2 and 0 <= b < 1 + k and np.log(0.05) < log_rho < np.log(80)):
        return -np.inf
    rho = np.exp(log_rho)
    lp = 0.0
    if rho_prior is not None:
        mu, sd = rho_prior
        # Gaussian in rho while sampling log(rho): include the Jacobian d rho / d log rho = rho
        lp = -0.5 * ((rho - mu) / sd) ** 2 + log_rho
    mod = model(t0, P, k, b, rho)
    if mod is None:
        return -np.inf
    return lp - 0.5 * np.sum(((y - mod) / e) ** 2)


def run_fit(c, P, t0, dur, depth, rho_prior, seed):
    rng = np.random.default_rng(seed)
    p0c = np.array([t0, P, np.sqrt(depth), 0.5, np.log(rho_prior[0] if rho_prior else 3.0)])
    scale = np.array([dur / 20, P * 1e-5, 0.1 * np.sqrt(depth), 0.2, 0.2])
    model = Model(c)
    y, e = c["flux"], c["flux_err"]
    p0 = []
    while len(p0) < N_WALKERS:
        p = p0c + scale * rng.standard_normal(5)
        p[3] = abs(p[3])
        if np.isfinite(log_prob(p, model, y, e, rho_prior)):
            p0.append(p)
    sampler = emcee.EnsembleSampler(N_WALKERS, 5, log_prob, args=(model, y, e, rho_prior))
    state = sampler.run_mcmc(np.array(p0), N_BURN, progress=False)
    sampler.reset()
    sampler.run_mcmc(state, N_STEPS, progress=False)
    chain = sampler.get_chain(flat=True, thin=5)
    lp = sampler.get_log_prob(flat=True, thin=5)
    best = chain[np.argmax(lp)]
    chi2_flat = np.sum(((y - 1) / e) ** 2)
    chi2_best = np.sum(((y - model(*best[:4], np.exp(best[4]))) / e) ** 2)
    try:
        tau = float(np.max(sampler.get_autocorr_time(quiet=True)))
    except Exception:
        tau = np.nan
    return chain, best, float(chi2_flat - chi2_best), tau, model


def summarize(x):
    lo, med, hi = np.percentile(x, [16, 50, 84])
    return {"median": float(med), "minus": float(med - lo), "plus": float(hi - med)}


def fit_star(tic):
    P, t0, dur = V.ephemeris(tic)
    tg = pd.read_csv(RESULTS / "targets_ranked.csv").set_index("TIC").loc[tic]
    import lcsearch as L
    raw, _, _ = L.reject_bad_sectors(L.load_star(tic))
    full = V.masked_clean(raw, P, t0, dur)
    near = np.abs(V.phase_of(full["time"], P, t0)) < WINDOW_DUR * dur
    c = {k: v[near] for k, v in full.items()}
    depth, sig, _ = V.depth_emp_matched(raw, P, t0, dur)
    # Correlated noise: inflate per-point errors so a transit-length box has the empirical error
    # (beta factor). Without this the likelihood treats hour-long wiggles as independent points.
    intr = np.abs(V.phase_of(full["time"], P, t0)) < dur / 2
    sig_white = 1 / np.sqrt(np.sum(1 / full["flux_err"][intr] ** 2))
    beta = float(max(1.0, sig / sig_white))
    c["flux_err"] = c["flux_err"] * beta
    # TIC density in g/cc; mass and radius errors from ExoFOP (TIC v8.2), ~12% and ~9%
    rho_sun_cgs = 1.408
    rho_tic = tg["Mass"] / tg["Rad"] ** 3 * rho_sun_cgs
    rho_tic_err = rho_tic * np.hypot(0.123, 3 * 0.087)
    out = {"tic": tic, "masked_detrend_depth_ppm": depth * 1e6, "masked_detrend_depth_err_emp_ppm": sig * 1e6,
           "white_noise_box_err_ppm": sig_white * 1e6, "red_noise_beta": beta,
           "phase2_ephemeris": {"P": P, "t0": t0, "dur_hr": dur * 24},
           "rho_tic_cgs": float(rho_tic), "rho_tic_err_cgs": float(rho_tic_err), "limb_darkening": LD}
    chains = {}
    for name, prior in (("free", None), ("stellar", (rho_tic, rho_tic_err))):
        chain, best, dchi2, tau, model = run_fit(c, P, t0, dur, depth, prior, seed=tic % 1000)
        t0s, Ps, ks, bs, rhos = chain.T
        rhos = np.exp(rhos)
        rstar = tg["Rad"] * (1 + 0.087 * np.random.default_rng(1).standard_normal(len(ks)))
        a = a_over_rs(rhos, Ps)
        inc = np.arccos(bs / a)
        # total duration T14 for a circular orbit
        arg = np.sqrt(np.clip((1 + ks) ** 2 - bs ** 2, 0, None)) / (a * np.sin(inc))
        t14 = Ps / np.pi * np.arcsin(np.clip(arg, 0, 1)) * 24
        out[name] = {
            "t0": summarize(t0s), "P": summarize(Ps), "k": summarize(ks), "b": summarize(bs),
            "rho_cgs": summarize(rhos), "rp_earth": summarize(ks * rstar / R_EARTH_PER_R_SUN),
            "depth_ppm": summarize(ks ** 2 * 1e6), "t14_hr": summarize(t14),
            "a_over_rs": summarize(a), "delta_chi2_vs_flat": dchi2, "autocorr_steps": tau,
            "best": dict(zip(["t0", "P", "k", "b", "rho_cgs"], [*map(float, best[:4]), float(np.exp(best[4]))])),
        }
        chains[name] = (chain, best, model)
    # how probable is the catalog density under the transit-shape-only fit?
    free_rho = np.exp(chains["free"][0][:, 4])
    out["free_fit_P_rho_below_tic"] = float(np.mean(free_rho < rho_tic))
    (V.OUT / f"{tic}_fit.json").write_text(json.dumps(out, indent=1))
    plot_fit(tic, full, c, chains, out)
    return out


def plot_fit(tic, full, c, chains, out):
    chain, best, model = chains["stellar"]
    t0, P, k, b, log_rho = best
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    ph = V.phase_of(c["time"], P, t0) * 24
    ax = axes[0, 0]
    o = np.argsort(ph)
    for exp, col in ((1800.0, "C1"), (200.0, "C0")):
        m = c["exptime"] == exp
        if m.sum() == 0:
            continue
        edges = np.arange(-6, 6.01, 0.25)
        idx = np.digitize(ph[m], edges)
        yb = [np.mean(c["flux"][m][idx == i]) for i in range(1, len(edges))]
        eb = [np.std(c["flux"][m][idx == i]) / np.sqrt(max((idx == i).sum(), 1)) for i in range(1, len(edges))]
        ax.errorbar(0.5 * (edges[1:] + edges[:-1]), (np.array(yb) - 1) * 1e6, np.array(eb) * 1e6, fmt="o", ms=3,
                    color=col, label=f"{int(exp)} s frames (15-min bins)")
        # model evaluated on these cadences, then binned identically
        mod = model(t0, P, k, b, np.exp(log_rho))
        mb = [np.mean(mod[m][idx == i]) for i in range(1, len(edges))]
        ax.plot(0.5 * (edges[1:] + edges[:-1]), (np.array(mb) - 1) * 1e6, "-", color=col, lw=1.5)
    ax.set_xlim(-6, 6)
    ax.set_xlabel("hours from mid-transit")
    ax.set_ylabel("ppm")
    ax.legend(fontsize=7)
    ax.set_title("folded data and best model (stellar-density fit)", fontsize=9)
    labels = [("k", "Rp/R*"), ("b", "impact parameter b"), ("rho_cgs", "stellar density (g/cc)"),
              ("rp_earth", "planet radius (R⊕)"), ("t14_hr", "duration T14 (h)")]
    for ax, (key, lab) in zip(axes.flat[1:], labels):
        for name, col in (("free", "0.5"), ("stellar", "C3")):
            ch = chains[name][0]
            vals = {"k": ch[:, 2], "b": ch[:, 3], "rho_cgs": np.exp(ch[:, 4]),
                    "rp_earth": ch[:, 2] * pd.read_csv(RESULTS / "targets_ranked.csv").set_index("TIC").loc[tic, "Rad"]
                    / R_EARTH_PER_R_SUN, "t14_hr": None}[key]
            if vals is None:
                continue
            rng = (0, 15) if key == "rho_cgs" else None
            ax.hist(vals, 60, range=rng, histtype="step", color=col, density=True, label=f"{name} density prior")
        if key == "rho_cgs":
            ax.axvline(out["rho_tic_cgs"], color="k", ls="--", lw=1, label="catalog star")
        ax.set_xlabel(lab)
        ax.legend(fontsize=7)
    s = out["stellar"]
    fig.suptitle(f"TIC {tic} transit fit:  Rp = {s['rp_earth']['median']:.2f} ± {s['rp_earth']['plus']:.2f} R⊕,  "
                 f"b = {s['b']['median']:.2f},  depth = {s['depth_ppm']['median']:.0f} ppm,  "
                 f"Δχ² vs no transit = {s['delta_chi2_vs_flat']:.0f}", fontsize=10)
    fig.tight_layout()
    (V.OUT / "plots").mkdir(parents=True, exist_ok=True)
    fig.savefig(V.OUT / "plots" / f"{tic}_fit.png", dpi=80)
    plt.close(fig)


def main():
    tics = [int(a) for a in sys.argv[1:]] or V.CANDIDATES
    V.OUT.mkdir(parents=True, exist_ok=True)
    for tic in tics:
        o = fit_star(tic)
        f, s = o["free"], o["stellar"]
        print(f"TIC {tic}: masked-detrend box depth {o['masked_detrend_depth_ppm']:.0f} ± "
              f"{o['masked_detrend_depth_err_emp_ppm']:.0f} ppm (empirical, matched processing; white "
              f"{o['white_noise_box_err_ppm']:.0f}, beta {o['red_noise_beta']:.2f})")
        for name, r in (("free", f), ("stellar", s)):
            print(f"  {name:8s} Rp {r['rp_earth']['median']:.2f}+{r['rp_earth']['plus']:.2f}-{r['rp_earth']['minus']:.2f} R⊕"
                  f"  b {r['b']['median']:.2f}+{r['b']['plus']:.2f}-{r['b']['minus']:.2f}"
                  f"  rho {r['rho_cgs']['median']:.1f}+{r['rho_cgs']['plus']:.1f}-{r['rho_cgs']['minus']:.1f} g/cc"
                  f"  T14 {r['t14_hr']['median']:.2f} h  depth {r['depth_ppm']['median']:.0f} ppm"
                  f"  P {r['P']['median']:.6f}±{r['P']['plus']:.6f}  dchi2 {r['delta_chi2_vs_flat']:.1f}"
                  f"  tau {r['autocorr_steps']:.0f}")
        print(f"  catalog rho {o['rho_tic_cgs']:.2f} ± {o['rho_tic_err_cgs']:.2f};"
              f" P(free-fit rho < catalog) = {o['free_fit_P_rho_below_tic']:.3f}")


if __name__ == "__main__":
    main()
