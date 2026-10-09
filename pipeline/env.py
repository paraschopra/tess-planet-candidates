"""Import first in every script: points all caches into the project (home dirs aren't writable here)."""
import os
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "cache"
DATA = ROOT / "data"
CATALOGS = DATA / "catalogs"
LIGHTCURVES = DATA / "lightcurves"
RESULTS = ROOT / "results"

for sub in ("mpl", "astropy", "lightkurve", "xdg"):
    (CACHE / sub).mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(CACHE / "mpl"))
os.environ.setdefault("XDG_CACHE_HOME", str(CACHE / "xdg"))
os.environ.setdefault("XDG_CONFIG_HOME", str(CACHE / "xdg"))
os.environ.setdefault("ASTROPY_CACHE_DIR", str(CACHE / "astropy"))
os.environ.setdefault("LIGHTKURVE_CACHE_DIR", str(CACHE / "lightkurve"))
warnings.filterwarnings("ignore", message=".*oktopus.*")
# Parallelism comes from process pools; keep numerical libraries single-threaded inside each worker
# (otherwise 12 workers x 14 threads oversubscribe the CPU badly).
for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
            "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(var, "1")
