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
