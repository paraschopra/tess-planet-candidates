"""Combined explainer video (~50 s): both candidates up front, how we found them, a fast real-3D flight to
TIC 370612437, its transit, and what happens next. Original score from 14_score.py.

Scenes (seconds):
  0-7    hook: "Two new exoplanet candidates found", sky map, a card per candidate with a live transit
  7-17   how we found them: search funnel, "is it new?" checks, signal strength
  17-21  where candidate #1 is in the sky (Cygnus), from Earth
  21-28  fast 3D flight from the Sun through real stars (32,434 stars within 60 pc from the TESS Input
         Catalog, Gaia distances; bright stars beyond as a fixed backdrop; schematic Milky Way glow), arrival
  28-41  the transit: planet crossing on the fitted tilted path, real folded TESS data and model
  41-50  what happens next
Rendered in parallel chunks, concatenated, and muxed with the score.

Writes results/media/candidates.mp4. Usage: python pipeline/15_combined_video.py [WORKERS]
"""
import importlib
import json
import subprocess
import sys
from multiprocessing import Pool

import env  # noqa: F401
import setuptools  # noqa: F401  distutils shim for batman

import batman
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord, get_constellation
from matplotlib.animation import FFMpegWriter

import vet as V
from env import CATALOGS, RESULTS

anim = importlib.import_module("11_animation")
score_mod = importlib.import_module("14_score")

W, H, DPI, FPS, DURATION = 1920, 1080, 100, 30, 50.0
MAIN, SECOND = 370612437, 252721481
NICKS = {370612437: "Aakanksha", 252721481: "Orbit"}  # informal nicknames
NICK = NICKS[370612437]
LY = 3.2616
BG, FG, DIM, ACCENT, RED = "#05070d", "#ece9e4", "#8d93a1", "#ffb066", "#ff5a4e"
R_END_PC = 0.002  # stop ~400 AU from the star
PULSES = (35.0,)  # main transit mid-time in the video (score heartbeat)
FLIGHT = (21.0, 26.5)
FAMOUS = {"Deneb": 102098, "Vega": 91262, "Altair": 97649, "Polaris": 11767}


# ---------------------------------------------------------------- helpers
def smooth(x):
    x = np.clip(x, 0, 1)
    return x * x * (3 - 2 * x)


def ramp(t, a, b):
    return smooth((t - a) / (b - a))


def unit(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], axis=-1)


def teff_rgb(T):
    Ts = np.array([2500, 3500, 4500, 5500, 6500, 8000, 12000])
    cols = np.array([[1.00, .55, .25], [1.00, .70, .45], [1.00, .80, .62], [1.00, .92, .82],
                     [1.00, .98, .98], [.85, .90, 1.00], [.70, .80, 1.00]])
    T = np.clip(np.nan_to_num(T, nan=5500), 2500, 12000)
    return np.stack([np.interp(T, Ts, cols[:, i]) for i in range(3)], axis=-1)


def bv_teff(bv):
    bv = np.nan_to_num(bv, nan=0.6)
    return 4600 * (1 / (0.92 * bv + 1.7) + 1 / (0.92 * bv + 0.62))


def slerp(a, b, x):
    om = np.arccos(np.clip(a @ b, -1, 1))
    return a if om < 1e-9 else (np.sin((1 - x) * om) * a + np.sin(x * om) * b) / np.sin(om)


class Camera:
    def __init__(self, pos, look, fov_deg):
        self.pos, self.f = pos, look / np.linalg.norm(look)
        up = np.array([0.0, 0.0, 1.0])
        self.r = np.cross(self.f, up)
        self.r /= np.linalg.norm(self.r)
        self.u = np.cross(self.r, self.f)
        self.fpx = (W / 2) / np.tan(np.radians(fov_deg) / 2)

    def project(self, xyz, at_infinity=False):
        v = xyz if at_infinity else xyz - self.pos
        z = v @ self.f
        x = W / 2 + self.fpx * (v @ self.r) / np.where(z > 0, z, np.nan)  # r points west: east on the left
        y = H / 2 + self.fpx * (v @ self.u) / np.where(z > 0, z, np.nan)
        ok = (z > 0) & (x > -50) & (x < W + 50) & (y > -50) & (y < H + 50)
        return x, y, ok, np.linalg.norm(v, axis=-1)


def mag_to_px(m):
    return np.clip(1.4 * 10 ** (-0.16 * (m - 6)), 0.35, 26)


def mag_to_alpha(m):
    return np.clip(10 ** (-0.4 * (m - 7.0)), 0.04, 1)


# ---------------------------------------------------------------- data (loaded once per worker)
D = {}


def load():
    tg = pd.read_csv(RESULTS / "targets_ranked.csv").set_index("TIC")
    near = pd.read_csv(CATALOGS / "tic_nearby.csv").dropna(subset=["ra", "dec", "Dist", "Tmag"])
    D["near_xyz"] = unit(near.ra.values, near.dec.values) * near.Dist.values[:, None]
    D["near_M"] = near.Tmag.values - 5 * np.log10(near.Dist.values / 10)  # absolute (TESS band)
    D["near_rgb"] = teff_rgb(near.Teff.values)
    D["near_tic"] = near.TIC.values
    br = pd.read_csv(CATALOGS / "bright_stars.csv").dropna(subset=["RAICRS", "DEICRS", "Vmag"])
    D["bright_dir"] = unit(br.RAICRS.values, br.DEICRS.values)
    D["bright_m"] = br.Vmag.values
    D["bright_rgb"] = teff_rgb(bv_teff(br["B-V"].values))
    D["bright_hip"] = br.HIP.values
    rng = np.random.default_rng(7)
    n = 60000  # schematic Milky Way: points concentrated to the galactic plane, brighter towards the centre
    l = rng.uniform(0, 360, n)
    b = rng.normal(0, 7, n) * (1 + 0.6 * np.cos(np.radians(l)))
    keep = rng.uniform(0, 1, n) < 0.45 + 0.55 * (1 + np.cos(np.radians(l))) / 2
    g = SkyCoord(l=l[keep], b=b[keep], unit="deg", frame="galactic").icrs
    D["mw_dir"] = unit(g.ra.deg, g.dec.deg)
    for tic in (MAIN, SECOND):
        m = tg.loc[tic]
        fit = json.loads((V.OUT / f"{tic}_fit.json").read_text())
        s = fit["stellar"]
        info = {"ra": m.ra, "dec": m.dec, "dist": m.Dist, "teff": m.Teff, "rad": m.Rad, "tmag": m.Tmag,
                "P": s["P"]["median"], "t0": s["t0"]["median"], "k": s["k"]["median"], "b": s["b"]["median"],
                "a": s["a_over_rs"]["median"], "t14": s["t14_hr"]["median"], "rp": s["rp_earth"]["median"],
                "rp_err": s["rp_earth"]["plus"], "depth": s["depth_ppm"]["median"], "ld": fit["limb_darkening"],
                "const": get_constellation(SkyCoord(m.ra, m.dec, unit="deg"))}
        info["teq"] = m.Teff * np.sqrt(1 / (2 * info["a"]))
        # folded 200-s data and model, as in 11_animation
        c = V.load_vet(tic)
        sel = c["exptime"] < 300
        ph = V.phase_of(c["time"][sel], info["P"], info["t0"]) * 24
        y = (c["flux"][sel] - 1) * 1e6
        edges = np.arange(-4, 4.01, 0.25)
        idx = np.digitize(ph, edges)
        info["xb"] = 0.5 * (edges[1:] + edges[:-1])
        info["yb"] = np.array([np.mean(y[idx == i]) for i in range(1, len(edges))])
        info["eb"] = np.array([np.std(y[idx == i]) / np.sqrt(max((idx == i).sum(), 1)) for i in range(1, len(edges))])
        th = np.linspace(-4, 4, 1600)
        p = batman.TransitParams()
        p.t0, p.per, p.rp, p.a = 0.0, info["P"], info["k"], info["a"]
        p.inc = np.degrees(np.arccos(info["b"] / info["a"]))
        p.ecc, p.w, p.limb_dark, p.u = 0.0, 90.0, "quadratic", info["ld"]
        info["th"], info["model"] = th, (batman.TransitModel(p, th / 24).light_curve(p) - 1) * 1e6
        info["img"] = anim.star_image(500, *info["ld"], anim.STAR_COLOUR[min(anim.STAR_COLOUR, key=lambda T: abs(T - m.Teff))])
        info["rgb"] = anim.STAR_COLOUR[min(anim.STAR_COLOUR, key=lambda T: abs(T - m.Teff))]
        D[tic] = info
    D["u_target"] = unit(D[MAIN]["ra"], D[MAIN]["dec"])
    D["u_start"] = unit(312.0, 44.0)
    D["main_idx"] = np.where(D["near_tic"] == MAIN)[0]


# ---------------------------------------------------------------- scene pieces
def remaining_pc(t):
    Dm = D[MAIN]["dist"]
    f0, f1 = FLIGHT
    mid = f0 + 0.65 * (f1 - f0)
    if t < f0:
        return Dm
    if t < mid:
        return Dm - (Dm - 2.0) * smooth((t - f0) / (mid - f0))
    if t < f1:
        return 2.0 * (R_END_PC / 2.0) ** smooth((t - mid) / (f1 - mid))
    return R_END_PC


def camera_at(t):
    Dm = D[MAIN]["dist"]
    look = slerp(D["u_start"], D["u_target"], ramp(t, 17, 20))
    fov = 70 - 12 * ramp(t, 17, 22)
    return Camera((Dm - remaining_pc(t)) * D["u_target"], look, fov)


def draw_sky(ax, t, dim=1.0, labels=0.0, target_alpha=1.0):
    cam = camera_at(min(t, 26))
    # Milky Way glow and backdrop at infinity
    x, y, ok, _ = cam.project(D["mw_dir"], at_infinity=True)
    ax.scatter(x[ok], y[ok], s=4.5, color="#9fb3ff", alpha=0.075 * dim, lw=0)
    x, y, ok, _ = cam.project(D["bright_dir"], at_infinity=True)
    m = D["bright_m"][ok]
    rgba = np.concatenate([D["bright_rgb"][ok], np.full((ok.sum(), 1), dim)], axis=1)
    ax.scatter(x[ok], y[ok], s=(2 * mag_to_px(m) * 0.72) ** 2, c=rgba, lw=0, zorder=2)
    # nearby stars in 3D at their real distances
    x, y, ok, dist = cam.project(D["near_xyz"])
    ok &= dist > 1e-4
    ok[D["main_idx"]] = False  # the target is drawn separately
    m = D["near_M"][ok] + 5 * np.log10(dist[ok] / 10)
    vis = m < 12.5
    rgba = np.concatenate([D["near_rgb"][ok][vis], (mag_to_alpha(m[vis]) * dim)[:, None]], axis=1)
    ax.scatter(x[ok][vis], y[ok][vis], s=(2 * mag_to_px(m[vis]) * 0.72) ** 2, c=rgba, lw=0, zorder=3)
    bright = m[vis] < 2.5
    if bright.any():  # soft halo for stars we pass close by
        ax.scatter(x[ok][vis][bright], y[ok][vis][bright], s=(6 * mag_to_px(m[vis][bright]) * 0.72) ** 2,
                   c=D["near_rgb"][ok][vis][bright], alpha=0.12 * dim, lw=0, zorder=3)
    # the target star
    tx, ty, tok, tdist = cam.project(D["u_target"][None, :] * D[MAIN]["dist"])
    tm = D[MAIN]["tmag"] - 5 * np.log10(D[MAIN]["dist"] / tdist[0])
    if tok[0] and target_alpha > 0.01:
        r = mag_to_px(tm)
        dim_t = dim * target_alpha
        ax.scatter(tx, ty, s=(2 * r * 0.72) ** 2, color=D[MAIN]["rgb"], alpha=max(0.6, mag_to_alpha(tm)) * dim_t, zorder=5)
        for k, a in ((3, 0.18), (7, 0.07), (14, 0.03)):
            ax.scatter(tx, ty, s=(2 * r * k * 0.72) ** 2, color=D[MAIN]["rgb"], alpha=a * dim_t * min(1, 10 ** (-0.15 * (tm - 4))), lw=0, zorder=4)
    if labels > 0:
        x, y, ok, _ = cam.project(D["bright_dir"], at_infinity=True)
        for name, hip in FAMOUS.items():
            i = np.where(D["bright_hip"] == hip)[0]
            if len(i) and ok[i[0]]:
                ax.text(x[i[0]] + 12, y[i[0]] + 8, name, color=DIM, fontsize=13, alpha=labels, zorder=6)
        if tok[0]:
            pulse = 1 + 0.25 * np.sin(t * 2 * np.pi * 0.8)
            ax.scatter(tx, ty, s=900 * pulse, facecolor="none", edgecolor=RED, lw=2, alpha=labels, zorder=6)
            ax.text(tx[0] + 30, ty[0] - 40, f'Candidate #1 "{NICK}"\nTIC {MAIN} · {D[MAIN]["dist"] * LY:.0f} light-years',
                    color=FG, fontsize=16, alpha=labels, zorder=6)
    return cam


def draw_disk(ax, cx, cy, R, alpha, info, planet_t=None, zoom=4, path=True, z=6):
    for rr, a in ((1.16, 0.05), (1.08, 0.08), (1.03, 0.12)):
        ax.add_patch(plt.Circle((cx, cy), R * rr, color=info["rgb"], alpha=a * alpha, lw=0, zorder=z))
    ax.imshow(info["img"], extent=(cx - R, cx + R, cy - R, cy + R), origin="lower", alpha=alpha, zorder=z + 1)
    if planet_t is not None:
        def xy(th):
            ph = 2 * np.pi * th / 24 / info["P"]
            return info["a"] * np.sin(ph), info["b"] * np.cos(ph)
        if path:
            tt = np.linspace(-3.2, 3.2, 600)
            px, py = xy(tt)
            near = np.abs(px) < 1.15
            ax.plot(cx + R * px[near], cy + R * py[near], ls=(0, (2, 4)), color=FG, alpha=0.25 * alpha, lw=1, zorder=z + 2)
        px, py = xy(planet_t)
        fade = np.clip((1.45 - abs(px)) / 0.3, 0, 1)  # planet fades in/out beyond the limb
        if fade > 0:
            ax.add_patch(plt.Circle((cx + R * px, cy + R * py), R * info["k"] * zoom, color="#030303",
                                    alpha=alpha * fade, zorder=z + 3))


def draw_curve(fig, rect, info, t_hr, alpha, small=False):
    ax = fig.add_axes(rect, facecolor="none")
    ax.errorbar(info["xb"], info["yb"], info["eb"], fmt="o", ms=3 if small else 4, color="#6f7a91",
                ecolor="#3d4456", alpha=alpha)
    ax.plot(info["th"], info["model"], color=ACCENT, lw=2, alpha=alpha)
    ax.plot([t_hr], [np.interp(t_hr, info["th"], info["model"])], "o", color=RED, ms=9 if not small else 7, alpha=alpha)
    ax.axvline(t_hr, color=RED, lw=1, alpha=0.5 * alpha)
    lo = min(info["model"].min(), np.nanmin(info["yb"] - info["eb"]))
    ax.set_xlim(-3.2, 3.2)
    ax.set_ylim(lo * 1.6, -lo * 0.9)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(DIM)
        ax.spines[sp].set_alpha(alpha)
    ax.tick_params(colors=DIM, labelsize=9 if small else 11)
    for lab in ax.get_xticklabels() + ax.get_yticklabels():
        lab.set_alpha(alpha)
    ax.set_xlabel("hours from mid-transit", color=FG, alpha=alpha, fontsize=10 if small else 12)
    ax.set_ylabel("brightness change (ppm)", color=FG, alpha=alpha, fontsize=10 if small else 12)
    return ax


def text(ax, x, y, s, size, alpha, color=FG, weight="normal", ha="left"):
    if alpha > 0.01:
        ax.text(x, y, s, fontsize=size, color=color, alpha=alpha, weight=weight, ha=ha, va="center", zorder=20)


# ---------------------------------------------------------------- frame
def box(ax, x, y, w, h, alpha, big, small, colour="#141b2d", edge=ACCENT):
    from matplotlib.patches import FancyBboxPatch
    if alpha <= 0.01:
        return
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=18", fc=colour, ec=edge,
                                lw=1.5, alpha=0.95 * alpha, zorder=15))
    text(ax, x + w / 2, y + h * 0.64, big, 30, alpha, ACCENT, weight="bold", ha="center")
    for j, line in enumerate(small.split("\n")):
        text(ax, x + w / 2, y + h * 0.33 - j * 22, line, 13.5, alpha, FG, ha="center")


def arrow(ax, x0, x1, y, alpha):
    if alpha > 0.01:
        ax.annotate("", xy=(x1, y), xytext=(x0, y), zorder=15,
                    arrowprops=dict(arrowstyle="-|>", color=DIM, lw=2, alpha=alpha, mutation_scale=18))


def sky_map(fig, rect, t, a):
    sky = fig.add_axes(rect, projection="mollweide", facecolor="#0e1324")
    br = D["bright_dir"]
    ra = np.degrees(np.arctan2(br[:, 1], br[:, 0])) % 360
    dec = np.degrees(np.arcsin(br[:, 2]))
    sky.scatter(-np.radians(((ra + 180) % 360) - 180), np.radians(dec), s=np.clip(5.5 - D["bright_m"], 0.2, 5) ** 1.7,
                color="white", alpha=0.7 * a, lw=0)
    mw = D["mw_dir"][::6]
    mra = np.degrees(np.arctan2(mw[:, 1], mw[:, 0])) % 360
    sky.scatter(-np.radians(((mra + 180) % 360) - 180), np.arcsin(mw[:, 2]), s=3, color="#7b8cff", alpha=0.14 * a, lw=0)
    for k, tic in enumerate((MAIN, SECOND)):
        i = D[tic]
        x0 = -np.radians(((i["ra"] + 180) % 360) - 180)
        sky.scatter([x0], [np.radians(i["dec"])], s=150 + 60 * np.sin(t * 5 + k), facecolor="none", edgecolor=RED, lw=2, alpha=a)
        # candidate #2 (higher RA) sits left of #1 on an east-left chart
        sky.text(x0 + (-0.12 if tic == MAIN else 0.12), np.radians(i["dec"]) + (-0.17 if tic == MAIN else 0.1),
                 f"#{k + 1}", color=FG, fontsize=12, weight="bold", ha="left" if tic == MAIN else "right", alpha=a)
    sky.set_xticklabels([])
    sky.set_yticklabels([])
    sky.grid(alpha=0.12 * a)
    for sp in sky.spines.values():
        sp.set_alpha(0.4 * a)


def render(fig, t):
    fig.clf()
    ax = fig.add_axes([0, 0, 1, 1], facecolor=BG)
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis("off")
    main = D[MAIN]

    if t < 7:  # hook: both candidates, alive from the first frame
        draw_sky(ax, 17, dim=0.18, target_alpha=0)
        a = 1 - ramp(t, 6.4, 7)
        text(ax, W / 2, H - 70, "Two new exoplanet candidates found", 50, a, weight="bold", ha="center")
        sky_map(fig, [0.345, 0.50, 0.31, 0.28], t, a)
        cards = ((MAIN, f'#1  "{NICKS[MAIN]}"', f"TIC {MAIN}"), (SECOND, f'#2  "{NICKS[SECOND]}"', "HIP 114417  (TIC 252721481)"))
        for k, (tic, title, name) in enumerate(cards):
            i = D[tic]
            x0 = 70 + k * W / 2
            ax.add_patch(plt.Rectangle((x0 - 30, 50), W / 2 - 100, 440, color="#121827", alpha=0.88 * a, zorder=10))
            text(ax, x0, 448, title, 30, a, ACCENT, weight="bold")
            text(ax, x0, 398, name, 18, a, FG)
            # "constellation" spelled out: a direction on the sky, not a place (Andromeda is not the galaxy)
            lines = [f"{i['const']} constellation, {i['dist'] * LY:.0f} light-years away",
                     f"about {i['rp']:.1f}× Earth's size",
                     f"orbits its star every {i['P']:.2f} days",
                     f"dims it by {i['depth'] / 1e4:.3f}% for {i['t14']:.1f} hours",
                     f"~{i['teq']:.0f} K: far too hot for life"]
            for j, ln in enumerate(lines):
                text(ax, x0, 345 - j * 50, ln, 18, a, FG if j < 2 else DIM)
            loop = (t / 4.0 + 0.5 - 0.12 * k) % 1.0  # a crossing every 4 s; both planets on their stars at frame one
            draw_disk(ax, x0 + 715, 250, 115, a, i, planet_t=-1.7 + 3.4 * loop, zoom=4, path=True, z=12)
        return

    if t < 17:  # how we found them
        draw_sky(ax, 17, dim=0.15, target_alpha=0)
        out = 1 - ramp(t, 16.4, 17)
        text(ax, W / 2, H - 70, "How we found them", 42, ramp(t, 7, 7.6) * out, weight="bold", ha="center")
        steps = [("526", "nearby stars, never in\nTESS's sharpest mode"), ("5 years", "of NASA TESS images,\nup to 44 sectors per star"),
                 ("~1 million", "possible orbits\ntested per star"), ("743", "repeating dips\nfound"),
                 ("31", "pass all major\nred-flag tests"), ("2", "stronger than any\nclean fake signal")]
        bw, gap, x = 262, 48, 70
        for j, (big, small) in enumerate(steps):
            al = ramp(t, 7.6 + 0.55 * j, 8.1 + 0.55 * j) * out
            box(ax, x + j * (bw + gap), 690, bw, 170, al, big, small, edge=RED if j == 5 else ACCENT)
            if j:
                arrow(ax, x + j * (bw + gap) - gap + 6, x + j * (bw + gap) - 6, 775, al)
        groups = [(70, "Are they new?", ["not among 341,351 NASA automatic planet alerts", "no TOI or community candidate on ExoFOP",
                                          "no planet papers on either star (SIMBAD)"]),
                  (W / 2 + 10, "Are they real?", ["local signal-to-noise about 7–8σ",
                                                   "same dip in NASA's and MIT's separate processing",
                                                   "hid a season, predicted it blind: 3 of 4 hits",
                                                   "in 526 upside-down light curves, chance never got this strong"])]
        for g, (x0, head, items) in enumerate(groups):
            al = ramp(t, 11.3 + 1.6 * g, 11.9 + 1.6 * g) * out
            ax.add_patch(plt.Rectangle((x0, 170), W / 2 - 80, 440, color="#121827", alpha=0.85 * al, zorder=10))
            text(ax, x0 + 30, 565, head, 26, al, ACCENT, weight="bold")
            for j, it in enumerate(items):
                aj = ramp(t, 11.6 + 1.6 * g + 0.35 * j, 12.1 + 1.6 * g + 0.35 * j) * out
                text(ax, x0 + 30, 495 - j * 78, "✓", 24, aj, "#7ee08a", weight="bold")
                text(ax, x0 + 75, 495 - j * 78, it, 17, aj, FG)
        text(ax, W / 2, 105, "Not yet confirmed: a faint neighbouring pair of eclipsing stars could still mimic the dip",
             19, ramp(t, 14.6, 15.2) * out, ACCENT, ha="center")
        return

    if t < 28.5:  # where is it, fast flight, arrival
        dim = ramp(t, 17, 17.6)
        labels = ramp(t, 17.6, 18.4) * (1 - ramp(t, 20.6, 21.3))
        draw_sky(ax, t, dim=dim, labels=labels, target_alpha=1 - ramp(t, 25.6, 27.5))
        text(ax, W / 2, H - 70, f'Where is candidate #1, "{NICK}"?', 36, ramp(t, 17.2, 17.9) * (1 - ramp(t, 20.6, 21.3)),
             weight="bold", ha="center")
        if t >= FLIGHT[0]:
            traveled = (main["dist"] - remaining_pc(t)) * LY
            a = ramp(t, 21, 21.5) * (1 - ramp(t, 26.3, 27))
            text(ax, 60, 80, f"{traveled:5.1f} light-years from Earth", 26, a, ACCENT, weight="bold")
            text(ax, 60, 40, "every dot is a real star, at its real distance (ESA Gaia)", 15, a, DIM)
            text(ax, W / 2, H - 70, f"Its light takes {main['dist'] * LY:.0f} years to reach us", 24,
                 ramp(t, 21.5, 22) * (1 - ramp(t, 25, 25.6)), FG, ha="center")
        if t >= 25.5:  # the star resolves into its disk
            g = ramp(t, 25.5, 28.3)
            draw_disk(ax, W / 2, H / 2, 30 + 330 * g, g, main)
            text(ax, W / 2, 90, f"An orange dwarf, {main['rad']:.2f}× the Sun's width, {main['teff']:.0f} K", 22,
                 ramp(t, 26.8, 27.8), DIM, ha="center")
        return

    if t < 41:  # the transit
        draw_sky(ax, 26.5, dim=0.25, target_alpha=0)
        s = ramp(t, 28.5, 29.8)
        cx, R = W / 2 + (0.27 * W - W / 2) * s, 360
        planet_hr = -3.2 + 6.4 * np.clip((t - 29.5) / 11.0, 0, 1)
        draw_disk(ax, cx, H / 2 - 30, R, 1.0, main, planet_t=planet_hr if t >= 29.2 else None)
        a = ramp(t, 28.8, 30)
        text(ax, 60, H - 60, f'Candidate #1 "{NICK}"', 40, a, weight="bold")
        text(ax, 60, H - 110, f"around TIC {MAIN}  ·  about {main['rp']:.2f} ± {main['rp_err']:.2f} Earth radii  ·  "
             f"one orbit every {main['P']:.2f} days  ·  ~{main['teq']:.0f} K", 19, a, DIM)
        draw_curve(fig, [0.56, 0.28, 0.40, 0.46], main, planet_hr, a)
        text(ax, 0.56 * W, 0.78 * H, "The star's brightness, all transits stacked", 18, a, FG)
        text(ax, 0.56 * W, 0.22 * H, "grey: real TESS measurements (2022–2024)  ·  orange: transit model", 14, a, DIM)
        text(ax, 60, 40, "Planet drawn 4× larger than true size. Nicknames are informal, not official IAU names.", 13, a, DIM)
        if t > 40.2:
            ax.add_patch(plt.Rectangle((0, 0), W, H, color=BG, alpha=ramp(t, 40.2, 41), zorder=50))
        return

    # what next
    draw_sky(ax, 26.5, dim=0.3, target_alpha=0)
    lines = [("What happens next", 38, FG, "bold"),
             ("Report both to NASA's ExoFOP, where astronomers coordinate follow-up", 22, FG, "normal"),
             ("Ground telescopes check that no faint neighbouring star is the real source", 22, FG, "normal"),
             ("TESS watches both stars again in 2027–28: the dips should return on schedule", 22, FG, "normal"),
             ("Until then: tentative planet candidates, not confirmed planets", 22, ACCENT, "normal")]
    for j, (txt, size, col, wt) in enumerate(lines):
        a = ramp(t, 41.3 + 0.9 * j, 42.1 + 0.9 * j) * (1 - ramp(t, 48.8, 50))
        text(ax, W / 2, H * 0.70 - j * 80, txt, size, a, col, weight=wt, ha="center")
    a = ramp(t, 45, 46) * (1 - ramp(t, 48.8, 50))
    text(ax, W / 2, 70, "Made with Claude Code  ·  data: NASA TESS, ESA Gaia, TESS Input Catalog  ·  original score generated in code",
         14, a, DIM, ha="center")


# ---------------------------------------------------------------- driver
def render_chunk(args):
    k, f0, f1, path = args
    load()
    plt.rcParams.update({"font.family": "DejaVu Sans"})
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI, facecolor=BG)
    writer = FFMpegWriter(fps=FPS, bitrate=9000, extra_args=["-pix_fmt", "yuv420p"])
    with writer.saving(fig, path, dpi=DPI):
        for f in range(f0, f1):
            render(fig, f / FPS)
            writer.grab_frame(facecolor=BG)
    plt.close(fig)
    return path


def main():
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    out_dir = RESULTS / "media"
    seg_dir = out_dir / "segments"
    seg_dir.mkdir(parents=True, exist_ok=True)
    n = int(DURATION * FPS)
    bounds = np.linspace(0, n, workers * 2 + 1).astype(int)
    jobs = [(k, bounds[k], bounds[k + 1], str(seg_dir / f"seg_{k:03d}.mp4")) for k in range(len(bounds) - 1)]
    with Pool(workers) as pool:
        paths = pool.map(render_chunk, jobs)
    (seg_dir / "list.txt").write_text("".join(f"file '{p}'\n" for p in paths))
    # music: same scene timings as the video
    audio = score_mod.score(DURATION, flight=FLIGHT, transits=PULSES, sparkle_spans=((0.5, 17.0), (41.0, 47.0)))
    score_mod.write_wav(out_dir / "score.wav", audio)
    out = out_dir / "candidates.mp4"
    subprocess.run(["ffmpeg", "-loglevel", "error", "-y", "-f", "concat", "-safe", "0", "-i", str(seg_dir / "list.txt"),
                    "-i", str(out_dir / "score.wav"), "-c:v", "copy", "-c:a", "aac", "-b:a", "192k", "-shortest",
                    "-movflags", "+faststart", str(out)], check=True)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
