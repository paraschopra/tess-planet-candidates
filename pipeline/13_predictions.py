"""Phase 3g: predicted transit times for the next TESS observations (for preregistration).

tess-point says TESS returns to TIC 370612437 in sectors 121, 122, 129 and to TIC 252721481 (HIP 114417)
in sectors 122, 129. Sector 121 is 2027-08-23 to 2027-09-19 (MIT TESS schedule); 122 and 129 dates are
estimates (27.3-day sectors after 121) and must be checked against the official schedule.

Uncertainty: sigma_t(n) = sqrt(sigma_t0^2 + (n sigma_P)^2) from the 07_fit posterior, ignoring the
(negative) P-t0 covariance at the 2019 reference epoch, so it is conservative.

Times: predictions are BJD_TDB (barycentric). The UTC columns are the corresponding GEOCENTRIC UTC
times: the barycentric light-travel correction (up to +-8.3 min, depends on the star's direction and
the date) is inverted iteratively, then TDB -> UTC. An earlier version formatted BJD_TDB directly as UTC,
off by up to ~6.5 min (caught by the independent Codex review).

Writes results/phase3/<TIC>_predicted_transits.csv. Usage: python pipeline/13_predictions.py
"""
import json

import env  # noqa: F401

import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time

import vet as V

SECTORS = {121: ("2027-08-23", "2027-09-19", "official"),
           122: ("2027-09-19", "2027-10-16", "estimated"),
           129: ("2028-03-28", "2028-04-24", "estimated")}
OBSERVED_IN = {370612437: [121, 122, 129], 252721481: [122, 129]}


GEOCENTRE = EarthLocation.from_geocentric(0, 0, 0, unit=u.m)


def bjd_to_geocentric_utc(bjd_tdb, coord):
    """Invert the barycentric correction: find the geocentric time t with t + ltt(t) = BJD_TDB."""
    bary = Time(bjd_tdb, format="jd", scale="tdb", location=GEOCENTRE)
    t = bary
    for _ in range(4):  # converges to << 1 s
        t = bary - t.light_travel_time(coord, kind="barycentric")
    return t.utc


def main():
    tg = pd.read_csv(V.RESULTS / "targets_ranked.csv").set_index("TIC")
    for tic, secs in OBSERVED_IN.items():
        coord = SkyCoord(tg.loc[tic, "ra"], tg.loc[tic, "dec"], unit="deg")
        s = json.loads((V.OUT / f"{tic}_fit.json").read_text())["stellar"]
        P, t0 = s["P"]["median"], s["t0"]["median"] + 2457000
        sP, st0 = max(s["P"]["plus"], s["P"]["minus"]), max(s["t0"]["plus"], s["t0"]["minus"])
        t14 = s["t14_hr"]["median"] / 24
        rows = []
        for sec in secs:
            start, end, status = SECTORS[sec]
            j0, j1 = Time(start).jd, Time(end).jd
            for n in range(int(np.ceil((j0 - t0) / P)), int(np.floor((j1 - t0) / P)) + 1):
                tc = t0 + n * P
                sig = np.sqrt(st0 ** 2 + (n * sP) ** 2)
                utc = bjd_to_geocentric_utc([tc, tc - t14 / 2 - 3 * sig, tc + t14 / 2 + 3 * sig], coord)
                rows.append({"sector": sec, "sector_dates": status, "epoch": n, "mid_BJD_TDB": round(tc, 5),
                             "mid_UTC_geocentric": utc[0].iso[:16], "sigma_min": round(sig * 1440, 1),
                             "window_start_UTC": utc[1].iso[:16], "window_end_UTC": utc[2].iso[:16],
                             "depth_ppm": round(s["depth_ppm"]["median"]), "duration_hr": round(t14 * 24, 2)})
        df = pd.DataFrame(rows)
        df.to_csv(V.OUT / f"{tic}_predicted_transits.csv", index=False)
        print(f"TIC {tic}: {len(df)} predicted transits in sectors {secs}; timing uncertainty "
              f"{df.sigma_min.min():.0f}-{df.sigma_min.max():.0f} min (1σ)")
        print(df.head(4).to_string(index=False))


if __name__ == "__main__":
    main()
