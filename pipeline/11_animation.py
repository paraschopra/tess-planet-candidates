"""Explainer animation for a candidate: the star, the planet crossing it, the real dimming data, where the
star is in the sky, and how far away it is.

Everything shown comes from the analysis: the star's colour, size and limb darkening; the planet's size
and tilted path (impact parameter) from the 07_fit transit model; the dimming curve is that model over
real TESS 200-s data folded on the period. The planet is drawn enlarged (labelled) because at true scale
it is a speck 1.4% of the star's width.

Writes results/phase3/<TIC>_animation.mp4 (needs ffmpeg).
Usage: python pipeline/11_animation.py TIC ["NICKNAME"]
"""
import json
import sys

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
from matplotlib.patches import Circle

import vet as V
from env import CATALOGS, RESULTS

FPS = 30
HOLD_S, CROSS_S = 1.5, 11.0
PLANET_ZOOM = 4  # planet drawn this many times its true size
LIGHT_YEARS_PER_PC = 3.2616
STAR_COLOUR = {4300: (1.0, 0.69, 0.40), 4400: (1.0, 0.71, 0.43)}  # blackbody RGB, ~K5 dwarf
FAMOUS = {"Polaris": 11767, "Vega": 91262, "Deneb": 102098, "Altair": 97649, "Sirius": 32349,
          "Betelgeuse": 27989, "Capella": 24608}
BG, FG, DIM = "#0b0e17", "#e8e6e3", "#8a8f9c"


def star_image(n, u1, u2, rgb):
    y, x = np.mgrid[-1:1:n * 1j, -1:1:n * 1j]
    r2 = x ** 2 + y ** 2
    mu = np.sqrt(np.clip(1 - r2, 0, 1))
    inten = 1 - u1 * (1 - mu) - u2 * (1 - mu) ** 2
    img = np.zeros((n, n, 4))
    img[..., :3] = np.array(rgb)[None, None, :] * inten[..., None]
    img[..., 3] = (r2 <= 1).astype(float)
    return img


def main():
    tic = int(sys.argv[1])
    nickname = sys.argv[2] if len(sys.argv) > 2 else None
    fit = json.loads((V.OUT / f"{tic}_fit.json").read_text())
    s = fit["stellar"]
    P, t0, k, b = s["P"]["median"], s["t0"]["median"], s["k"]["median"], s["b"]["median"]
    a, t14 = s["a_over_rs"]["median"], s["t14_hr"]["median"]
    rp, rp_err = s["rp_earth"]["median"], s["rp_earth"]["plus"]
    tg = pd.read_csv(RESULTS / "targets_ranked.csv").set_index("TIC").loc[tic]
    ly = tg["Dist"] * LIGHT_YEARS_PER_PC
    coord = SkyCoord(tg["ra"], tg["dec"], unit="deg")
    constellation = get_constellation(coord)
    rgb = STAR_COLOUR[min(STAR_COLOUR, key=lambda T: abs(T - tg["Teff"]))]

    # real data: transit-masked detrended TESS-SPOC light curve, 200-s sectors only (sharp shape), folded
    c = V.load_vet(tic)
    m = c["exptime"] < 300
    ph = V.phase_of(c["time"][m], P, t0) * 24
    y = (c["flux"][m] - 1) * 1e6
    edges = np.arange(-4, 4.01, 0.25)
    idx = np.digitize(ph, edges)
    xb = 0.5 * (edges[1:] + edges[:-1])
    yb = np.array([np.mean(y[idx == i]) for i in range(1, len(edges))])
    eb = np.array([np.std(y[idx == i]) / np.sqrt(max((idx == i).sum(), 1)) for i in range(1, len(edges))])
    # model: instantaneous transit with the fitted parameters
    th = np.linspace(-4, 4, 2000)
    params = batman.TransitParams()
    params.t0, params.per, params.rp, params.a = 0.0, P, k, a
    params.inc = np.degrees(np.arccos(b / a))
    params.ecc, params.w, params.limb_dark, params.u = 0.0, 90.0, "quadratic", fit["limb_darkening"]
    model = (batman.TransitModel(params, th / 24).light_curve(params) - 1) * 1e6

    # planet position on the sky-projected disk (stellar radii) as a function of hours from mid-transit
    def planet_xy(t_hr):
        phase = 2 * np.pi * t_hr / 24 / P
        return a * np.sin(phase), b * np.cos(phase)

    span = 3.2  # hours either side shown in the animation
    n_hold, n_cross = int(HOLD_S * FPS), int(CROSS_S * FPS)
    times = np.concatenate([np.full(n_hold, -span), np.linspace(-span, span, n_cross), np.full(n_hold, span)])

    plt.rcParams.update({"text.color": FG, "axes.labelcolor": FG, "xtick.color": DIM, "ytick.color": DIM,
                         "axes.edgecolor": DIM, "font.size": 11})
    fig = plt.figure(figsize=(16, 9), dpi=120, facecolor=BG)
    title = f"Planet candidate around TIC {tic}"
    if nickname:
        title += f'   ·   nicknamed "{nickname}"'
    fig.text(0.03, 0.95, title, fontsize=22, weight="bold")
    fig.text(0.03, 0.915, f"Found in NASA TESS data · about {rp:.2f} ± {rp_err:.2f} Earth radii · one orbit every "
             f"{P * 24:.1f} hours ({P:.2f} days) · dims its star by {s['depth_ppm']['median'] / 1e4:.2f}% for "
             f"{t14:.1f} hours", fontsize=13, color=DIM)
    fig.text(0.03, 0.02, "Unconfirmed candidate: the dip could still come from a faint neighbouring star. "
             + ("Nickname is informal, not an official IAU name. " if nickname else "")
             + "Planet drawn 4× larger than true size.", fontsize=10, color=DIM)

    # 1. the star and the planet
    ax = fig.add_axes([0.02, 0.08, 0.46, 0.80], facecolor=BG)
    ax.set_xlim(-1.45, 1.45)
    ax.set_ylim(-1.45, 1.45)
    ax.set_aspect("equal")
    ax.axis("off")
    for r, al in ((1.18, 0.05), (1.10, 0.08), (1.04, 0.12)):  # soft glow
        ax.add_patch(Circle((0, 0), r, color=rgb, alpha=al, lw=0))
    ax.imshow(star_image(600, *fit["limb_darkening"], rgb), extent=(-1, 1, -1, 1), origin="lower", zorder=2)
    tt = np.linspace(-span, span, 400)
    px, py = planet_xy(tt)
    ax.plot(px, py, ls=(0, (2, 4)), color=FG, alpha=0.25, lw=1, zorder=3)
    planet = Circle(planet_xy(-span), k * PLANET_ZOOM, color="#050505", zorder=4)
    ax.add_patch(planet)
    ax.text(0, -1.33, f"{'K' if tg['Teff'] < 5200 else 'G'}-type dwarf star · {tg['Teff']:.0f} K · "
            f"{tg['Rad']:.2f}× the Sun's width", ha="center", color=DIM, fontsize=11)

    # 2. the dimming graph
    ax2 = fig.add_axes([0.54, 0.50, 0.43, 0.36], facecolor=BG)
    ax2.errorbar(xb, yb, eb, fmt="o", ms=4, color="#6f7a91", ecolor="#3d4456", label="TESS data, 2022–2024 (15-min bins)")
    ax2.plot(th, model, color="#ffb066", lw=2, label="transit model")
    dot, = ax2.plot([], [], "o", color="#ff5a4e", ms=10, zorder=5)
    vline = ax2.axvline(-span, color="#ff5a4e", lw=1, alpha=0.6)
    ax2.set_xlim(-span, span)
    lo = min(model.min(), np.nanmin(yb - eb))
    ax2.set_ylim(lo * 1.6, -lo * 0.9)
    ax2.set_xlabel("hours from the middle of the transit")
    ax2.set_ylabel("brightness change (ppm)")
    ax2.set_title("The star's brightness, all transits stacked", color=FG, fontsize=13, loc="left")
    ax2.legend(loc="lower left", fontsize=9, frameon=False)
    for sp in ("top", "right"):
        ax2.spines[sp].set_visible(False)

    # 3. where in the sky
    ax3 = fig.add_axes([0.53, 0.07, 0.26, 0.33], projection="mollweide", facecolor="#111626")
    bright = pd.read_csv(CATALOGS / "bright_stars.csv").dropna(subset=["RAICRS", "DEICRS", "Vmag"])

    def sky(ra, dec):  # sky-chart orientation: RA increases to the left
        return -np.radians(((np.asarray(ra) + 180) % 360) - 180), np.radians(dec)

    bx, by = sky(bright.RAICRS, bright.DEICRS)
    ax3.scatter(bx, by, s=np.clip(6 - bright.Vmag, 0.2, 6) ** 1.8, color="white", alpha=0.8, lw=0)
    gal = SkyCoord(l=np.linspace(0, 360, 1441), b=np.zeros(1441), unit="deg", frame="galactic").icrs
    gx, gy = sky(gal.ra.deg, gal.dec.deg)
    ax3.scatter(gx, gy, s=40, color="#7b8cff", alpha=0.08, lw=0)
    for name, hip in FAMOUS.items():
        row = bright[bright.HIP == hip]
        if len(row):
            x, y_ = sky(row.RAICRS.iloc[0], row.DEICRS.iloc[0])
            ax3.text(x, y_ + 0.06, name, fontsize=7, color=DIM, ha="center")
    tx, ty = sky(tg["ra"], tg["dec"])
    ring = ax3.scatter([tx], [ty], s=200, facecolor="none", edgecolor="#ff5a4e", lw=2)
    ax3.scatter([tx], [ty], s=20, color="#ff5a4e")
    ax3.set_xticklabels([])
    ax3.set_yticklabels([])
    ax3.grid(alpha=0.15)
    ax3.set_title(f"Where: in the constellation {constellation}\n(blue band = Milky Way)", color=FG, fontsize=11)

    # 4. how far
    ax4 = fig.add_axes([0.80, 0.07, 0.18, 0.33], facecolor=BG)
    ax4.axis("off")
    ax4.text(0.0, 0.92, "How far", fontsize=13, weight="bold")
    ax4.text(0.0, 0.72, f"{ly:.0f} light-years", fontsize=20, color="#ffb066")
    ax4.text(0.0, 0.60, f"({tg['Dist']:.1f} parsecs, measured by Gaia)", fontsize=9, color=DIM)
    ax4.text(0.0, 0.40, f"The light TESS caught left\nthis star in the {int((2022 - ly) // 10 * 10)}s.", fontsize=11)
    # V - T is about 1.2 mag for a K dwarf (HIP 114417: V 11.0, T 9.73)
    ax4.text(0.0, 0.15, f"Visible brightness: magnitude ~{tg['Tmag'] + 1.2:.0f}\n"
             "(too faint for the naked eye;\nbinoculars or a small telescope)", fontsize=9, color=DIM)

    out = V.OUT / f"{tic}_animation.mp4"
    writer = FFMpegWriter(fps=FPS, bitrate=4000, extra_args=["-pix_fmt", "yuv420p"])
    with writer.saving(fig, str(out), dpi=120):
        for i, t in enumerate(times):
            planet.center = planet_xy(t)
            mod_t = np.interp(t, th, model)
            dot.set_data([t], [mod_t])
            vline.set_xdata([t, t])
            ring.set_sizes([200 + 120 * np.sin(i / FPS * 2 * np.pi * 0.8) ** 2])
            writer.grab_frame(facecolor=BG)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
