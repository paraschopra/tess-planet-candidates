"""Sanity check: inject a box-shaped transit into a real light curve and see if bls_search recovers it.

Usage: python pipeline/inject_test.py TIC RP_EARTH PERIOD_D
"""
import sys
import env  # noqa: F401
import numpy as np
import pandas as pd
import lcsearch as L
from env import RESULTS

tic, rp, period = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
m = pd.read_csv(RESULTS / "targets_ranked.csv").set_index("TIC").loc[tic]
rho = float(m["Mass"] / m["Rad"] ** 3)
depth = (rp * 0.009168 / m["Rad"]) ** 2 / (1 + m["Rcont"])
dur = L.duration_days(period, rho)
lc, bad, _ = L.reject_bad_sectors(L.load_star(tic))
t0 = lc["time"].min() + 0.37 * period
lc["flux"] = lc["flux"] - depth * L.transit_mask(lc["time"], period, t0, dur, factor=1.0)
c = L.bin_lc(L.clean(lc))
signals, _ = L.bls_search(c, rho, max_signals=1)
s = signals[0]
hit = abs(s["period"] - period) / period < 1e-3 and abs(((s["t0"] - t0 + period / 2) % period) - period / 2) < dur
print(f"injected Rp={rp} P={period} depth={depth*1e6:.0f}ppm dur={dur*24:.2f}h -> found P={s['period']:.5f} "
      f"depth={s['depth_ppm']:.0f} SNR={s['snr']:.1f} SDE={s['sde']:.1f}  {'RECOVERED' if hit else 'MISSED'}")
