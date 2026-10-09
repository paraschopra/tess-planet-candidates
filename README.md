# Two Earth-sized planet candidates in NASA TESS data

A search of five years of NASA TESS full-frame-image data on **526 nearby stars that were never observed at TESS's 2-minute cadence**. It found two Earth-sized transiting-planet candidates that are not on any NASA candidate or planet list (checked 2026-10-09).

> **Status: tentative, unconfirmed candidates.** The dips look real, but TESS's 21-arcsecond pixels can't rule out a faint neighbouring pair of eclipsing stars as the source. Only follow-up with other telescopes can turn these into planets.

| | **#1 "Aakanksha"** | **#2 "Orbit"** |
|---|---|---|
| Star | TIC 370612437 (Cygnus constellation, 158 ly) | HIP 114417 = TIC 252721481 (Andromeda constellation, 134 ly) |
| Orbital period | 3.564921 ± 0.000022 d | 4.933533 ± 0.000019 d |
| Transit depth / duration | 210 ppm / 1.44 h | 279 ppm / 1.59 h |
| Size | 1.04 ± 0.15 Earth radii | 1.10 ± 0.14 Earth radii |
| Equilibrium temperature | ~850 K | ~716 K |
| Local depth-to-noise | ~7.2 | ~7.7 |

Nicknames are informal, not IAU names.

- **Video (50 s):** [`docs/media/candidates.mp4`](docs/media/candidates.mp4)
- **Plain-language guide** to the concepts, methods and every result: [`docs/GUIDE.md`](docs/GUIDE.md)
- **Research note draft** (RNAAS): [`docs/rnaas/note.pdf`](docs/rnaas/note.pdf)
- **Preregistered transit predictions** for TESS sectors 121, 122 and 129 (2027–28): [Zenodo, doi:10.5281/zenodo.23257784](https://doi.org/10.5281/zenodo.23257784)
- **Draft community-candidate (CTOI) parameters:** [`docs/CTOI_DRAFT.md`](docs/CTOI_DRAFT.md)
- **Next steps and the 2027 prediction test:** [`docs/next_steps.md`](docs/next_steps.md)

## How it was done
| Step | Script | What it does |
|---|---|---|
| 1 | `pipeline/01_targets.py` | Nearby stars (Gaia GCNS, ≤ 60 pc) matched to the TESS Input Catalog; cuts on size, brightness, contamination and coverage; excludes known planet and candidate hosts |
| 2 | `pipeline/02_download.py` | One TESS-SPOC (else QLP) light curve per sector |
| 3 | `pipeline/lcsearch.py`, `03_search.py` | Biweight detrending, then an iterative box-least-squares search (0.5–25 d, up to 3 signals per star) |
| 4 | `pipeline/05_triage.py` | Automatic red flags (odd/even, secondary eclipse, variability, single-event dominance, sector consistency, …), calibrated on injected planets (`injection_grid.py`, `flag_calibration.py`) |
| 5 | `pipeline/06_inverted.py` | False-alarm baseline: the identical search and triage on upside-down light curves |
| 6 | `pipeline/07_fit.py` … `10_pixels.py` | Vetting: MCMC transit fits, retrospective hold-out seasons, NASA vs MIT reductions, pixel-level difference imaging |
| 7 | `pipeline/13_predictions.py` | Predicted transit times for future TESS sectors |
| 8 | `pipeline/11_animation.py`, `14_score.py`, `15_combined_video.py`, `16_rnaas_figure.py` | Videos, original score and the note's figure |

`pipeline/12_triceratops.py` (statistical validation) was set up but not completed; see guide Part 4.

## Independent reviews
- **Two AI audits** (OpenAI Codex and a fresh Claude agent) tried to break the analysis: `results/phase3/audit_*.md`.
- **A second Codex review** (`results/review/REVIEW_CODEX.md`) found overstated claims and two code errors. Every finding was addressed: `results/review/RESPONSE_TO_CODEX.md`.

## Reproducing
```bash
python -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python pipeline/01_targets.py      # then 02_download.py, 03_search.py, 05_triage.py, ...
```
Light curves (~440 MB) and download caches are not in the repo. `02_download.py` re-fetches them from MAST. Steps 3 and 5 take about 1.5–3 hours on 12 cores. TRICERATOPS needs its own environment (`requirements-triceratops.txt`).

## How this was made
The code was written with the AI coding agent Claude Code (Anthropic). Independent reviews used OpenAI Codex and a separate Claude session. All data are public: NASA TESS via MAST, ESA Gaia, the TESS Input Catalog, ExoFOP and the NASA Exoplanet Archive. The project was inspired by a similar r/ClaudeAI post.

## License
Code: MIT (see `LICENSE`). Text, figures and videos: CC BY 4.0.
