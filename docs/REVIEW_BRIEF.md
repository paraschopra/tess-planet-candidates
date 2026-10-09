# Independent review brief: two TESS planet candidates

**For:** an independent coding/science agent with no prior context on this project.
**Goal:** check whether our claims hold up **before** they go public (tweet) and are submitted to NASA's ExoFOP as Community TESS Objects of Interest (CTOIs). Please try hard to **break** the analysis; don't just confirm it. A finding that kills a candidate is a success.

---

## 1. Ground rules
- **Read-only.** Do not modify, move or delete anything in the repo except your own report (see §8). Do not commit, push, submit to ExoFOP, or post anything.
- Put scratch scripts and outputs in a new folder **outside** the repo, or under `results/review/scratch/`.
- Use the project environment:
  - Python is `.venv/bin/python` (3.12).
  - Every pipeline script starts with `import env`, which sets cache paths and single-threaded BLAS. Run scripts from the repo root as `.venv/bin/python pipeline/<script>.py`. For ad-hoc code, run from `pipeline/` so that `import env, lcsearch, vet` works.
  - TRICERATOPS lives in a separate env, `.venv-tri`, because it pins old numpy/scipy.
- **Don't rerun the full 526-star searches** (`03_search.py`, `06_inverted.py`); they take 1.5–3 hours. Their outputs are saved (see §4). Targeted re-checks on the candidates take seconds to minutes.
- Network: MAST, Gaia, VizieR and SIMBAD are used, and are slow or flaky at times. `vizier.cfa.harvard.edu` worked when the main VizieR mirror refused connections.

## 2. What the project did, in brief
1. **Targets** (`pipeline/01_targets.py`):
   - Started from Gaia's catalogue of nearby stars (≤ 60 pc) and matched them to the TESS Input Catalog (TIC).
   - Kept stars that are small (R ≤ 1.5 R☉), bright (Tmag ≤ 13), lightly contaminated (TIC contamination ratio < 0.2), observed in ≥ 3 sectors spanning ≥ 2 years, and **never observed at 2-minute cadence**.
   - Excluded known TOI, CTOI and confirmed-planet hosts (snapshot of 2026-10-08).
   - Result: **526 stars** with full-frame-image (FFI) light curves only. These were never searched by NASA's SPOC planet-search pipeline, which runs only on 2-minute targets.
2. **Download** (`02_download.py`): one light curve per sector, TESS-SPOC PDCSAP preferred, else QLP.
3. **Search** (`lcsearch.py`, `03_search.py`):
   - Reject bad sectors; detrend with a biweight filter (wotan, 0.75 d window); clip flares; bin to 10 min.
   - Run box least squares (BLS) on a baseline-aware grid of ~0.5–1 million periods (0.5–25 d), using SNR as the objective.
   - SDE is computed after removing the trend of power against period.
   - Up to 3 signals per star, masking each one before searching again.
   - Thresholds: SDE ≥ 7, SNR ≥ 6, ≥ 3 transits.
4. **Triage** (`05_triage.py`): automatic red flags (odd/even, secondary eclipse, variability, duration, edge, few-event dominance, sector consistency, broad depression, shared period, crowded period…). Thresholds were calibrated on injected planets (`flag_calibration.py`, `injection_grid.py`).
   - **743 → 31 (no hard flags) → 10 (no flags at all) → 2 that stand out from chance.**
5. **Vetting** (Phase 3, `vet.py`, `06`–`13`): described in §4.

Plain-language background and full results are in `docs/GUIDE.md`, Parts 1–5.

## 3. The claims to verify

| # | Claim | Where the evidence is |
|---|---|---|
| C1 | **Both candidates are new.** No TOI, CTOI, confirmed planet, SPOC TCE or planet paper exists for either star (as of 2026-10-09). | `data/catalogs/toi.csv`, `ctoi.csv`, `nasa_ps_tic.csv`; `data/catalogs/tce/*.csv` (132 files, 341,351 TCEs); GUIDE Part 4 table |
| C2 | **TIC 370612437**: P = 3.564921 ± 0.000022 d, epoch BJD 2458742.1409 ± 0.0069. Depth 208 ± 29 ppm with an empirical error, i.e. **7.2σ**. SDE 11.26. Rp = 1.04 ± 0.15 R⊕, b ≈ 0.73, T14 ≈ 1.44 h. 8 sectors (16, 17, 56, 57, 76, 77, 83, 84), 51 transits. | `results/phase3/370612437_fit.json`, `results/phase2/triage.csv` |
| C3 | **TIC 252721481 = HIP 114417**: P = 4.933533 ± 0.000019 d, epoch BJD 2458740.6616 ± 0.0060. Depth 263 ± 34 ppm, **7.7σ**. SDE 10.26. Rp = 1.10 ± 0.14 R⊕, b ≈ 0.68. 5 sectors (16, 17, 57, 77, 84), 19 transits. | `results/phase3/252721481_fit.json` |
| C4 | **The signals are not chance.** We flipped all 526 light curves upside down and ran the same search and triage: 17 flag-free fake signals, the strongest at SDE 8.38. None reached 9. The real data has 10 flag-free signals; the top two are at 11.26 and 10.26. | `results/phase3/inverted/stars/*.json`, `results/phase3/inverted_triage.csv` vs `results/phase2/triage.csv` |
| C5 | **The signals are not a processing artifact.** The same dip at the same ephemeris appears in MIT's QLP light curves: 225 ± 29 and 238 ± 52 ppm. | `results/phase3/*_products.json`, `pipeline/09_products.py` |
| C6 | **Hold-out prediction.** Blind search of the other seasons, then prediction of the hidden season: TIC 370612437 hits 3 of 4 (2022 fails); TIC 252721481 hits 3 of 4 (hiding 2019 leaves too little data to recover the period). | `results/phase3/*_holdout.json`, `pipeline/08_holdout.py` |
| C7 | **Physically consistent.** The stellar density implied by the transit shape alone agrees with the TIC density. No secondary eclipse. Odd/even depths agree within 1.4σ and 2.0σ. | `*_fit.json`, `triage.csv`, `results/phase3/audit_claude.md` |
| C8 | **Not localized: a blend is NOT excluded.** We say this openly. Difference imaging detects the dimming near the target but cannot pin the source (TIC 370612437 is at galactic latitude −1.3°, with 471 Gaia stars within 150"). For TIC 252721481, SPOC centroids show a ~3.5σ in-transit shift. | `results/phase3/*_pixels.json`, `pipeline/10_pixels.py`, `audit_claude.md` |
| C9 | **A third signal was set aside.** HIP 93326 / TIC 350883383, P = 14.93 d: no pixel detection, hold-out mostly fails, period in the TESS-orbit artifact pile-up. | `results/phase3/350883383_*.json` |
| C10 | **Predicted transits** for future TESS sectors 121, 122 and 129 (2027–28). | `results/phase3/*_predicted_transits.csv`, `pipeline/13_predictions.py` |
| C11 | **Public wording is accurate:** "tentative", "candidates", "Earth-sized", "not yet confirmed", "7.2σ / 7.7σ", "3 of 4 hits", "chance never got this strong". | `docs/TWEET_DRAFT.md`, `docs/CTOI_DRAFT.md`, `docs/media/candidates.mp4` |

## 4. Data and outputs (what exists on disk)
- **Light curves:** `data/lightcurves/<TIC>.npz` (gitignored, 704 files, 441 MB). Each holds `time` (BTJD), `flux` (normalised per sector), `flux_err`, `sector`, `exptime`, and `authors` (product per sector). Load with `lcsearch.load_star(tic)`.
  - If you work from a fresh clone without this folder, regenerate the candidates' files with `pipeline/02_download.py`. It downloads only the stars it's given; adapt it to pass the TICs.
- **Raw FITS cache:** `cache/` (gitignored). It contains TESS-SPOC/QLP light-curve FITS, including CROWDSAP and aperture masks, and TESScut cutouts.
- **Phase 2:** `results/phase2/triage.csv` (all signals, with flags and diagnostics) and `results/phase2/stars/*.json` (per-star search records).
- **Phase 3:** in `results/phase3/`:
  - `*_fit.json`: batman + emcee fit with red-noise inflation (β) and the TIC density prior.
  - `*_holdout.json`, `*_products.json`, `*_pixels.json`.
  - `inverted_triage.csv`.
  - `plots/*.png`: fits, hold-out, products, pixels.
- **Previous audits:** `results/phase3/audit_codex.md` (OpenAI Codex) and `results/phase3/audit_claude.md` (a fresh Claude agent). Their criticisms led to fixes: matched-processing noise estimate, Jacobian in the prior, red-noise β, joint multi-sector difference imaging, and the full-triage flipped survey. **Please don't just repeat them.** Check whether the fixes are correct and look for what they missed.

## 5. Known weak spots (probe these, don't take our word)
1. **Blend.** This is the main risk. Is there any additional test we could run with existing data to discriminate further? For example: per-pixel depth maps; PRF-based fitting instead of a Gaussian; comparing aperture sizes using the SPOC aperture masks; TESS-band flux of the Gaia neighbours inside the actual SPOC apertures.
2. **TIC 370612437's 2022 season.** The season depths are 186 ± 50, 103 ± 48, 235 ± 42 and 243 ± 41 ppm (χ² p ≈ 0.12), and the blind hide-2022 prediction found ~40 ± 42 ppm. Is 2022 consistent with noise, or a sign of something else?
3. **TIC 252721481's centroid shift** (~3.5σ in MOM_CENTR, measured by the earlier auditor). Real, or systematics? What does it imply?
4. **Look-elsewhere and selection effects.**
   - The triage thresholds were calibrated on injected planets. They were partly designed *after* looking at real signals (e.g. the "broad depression" and "few events dominate" flags were added after visual review).
   - Does the flipped survey fairly calibrate the *whole* selection procedure?
   - Is "SDE > 10 and flag-free" a legitimate post-hoc criterion?
   - Is the upside-down null biased? The cleaning clips upward outliers at 3σ and downward ones at 10σ, so flipping swaps which side gets clipped hard.
5. **Empirical error bars** (`vet.depth_emp_matched`: the same box at 40 phase offsets, each with its own masked detrend, and a MAD-based σ). Are the 7.2σ / 7.7σ figures defensible? Is the 40-offset sample adequate?
6. **The hold-out test was modified after a fold failed:** a finer box-width refit was added. Both versions are reported. Is the remaining test still meaningful?
7. **Dilution.** TESS-SPOC PDCSAP is already crowding-corrected (CROWDSAP 0.86–0.95). Phase 2's `rp_earth` additionally multiplied by (1 + Rcont), a known double count; Phase 3's fit doesn't. Check that nothing user-facing uses the double-counted radius.
8. **Stellar parameters.** Radius and density come from TIC v8.2. Is the planet radius error budget reasonable?
9. **Time systems.** The light curves are BTJD (BJD_TDB − 2457000). The predicted transits add 2457000. Verify the epochs and the predicted times end to end.
10. **Code correctness.** Spot-check these:
    - `lcsearch.py`: `period_grid`, `sde`, `bls_search` masking;
    - `vet.py`: `masked_clean`, `depth_emp_matched`;
    - `07_fit.py`: a/R*, inclination, exposure supersampling, β, prior Jacobian;
    - `08_holdout.py`: leakage between training and test seasons;
    - `10_pixels.py`: WCS and proper-motion handling, the joint fit, neighbour Δχ².

## 6. Suggested independent checks (pick what you think matters most)
- **Recompute both candidates from scratch with your own simple method:** load the npz, detrend your own way, fold at the claimed period, and measure the depth and an honest error. Do you get ~200–260 ppm at 7–8σ?
- **Run TLS** (`transitleastsquares`, installed in `.venv`) or your own BLS on each candidate's light curve. Does the same period come out on top?
- **Per-transit depths** with momentum-dump and quality-flag timestamps overlaid. The flags are available in the cached FITS files under `cache/lightkurve/`.
- **Re-check novelty live:** ExoFOP target pages `https://exofop.ipac.caltech.edu/tess/target.php?id=<TIC>&json`, the NASA Exoplanet Archive TOI table, SIMBAD references, and NASA ADS full text. Report anything filed since 2026-10-08.
- **Estimate a false-alarm probability your own way** for each candidate, for example with per-star bootstrap or scramble searches, and compare it with ours.
- **Check the wording** of `docs/TWEET_DRAFT.md`, `docs/CTOI_DRAFT.md` and the video text (`pipeline/15_combined_video.py`, the `render()` text) against the evidence. Flag anything overstated.

## 7. Facts that are easy to get wrong
- The data span is **Sep 2019 – Oct 2024**, not 2025.
- "Sharpest mode" means TESS's 2-minute cadence. These stars only have FFI data: 30-min in 2019, 200-s from 2022.
- Both would be hot: Teq ≈ 850 K and 716 K. Neither is habitable.
- The nicknames "Aakanksha" (TIC 370612437) and "Orbit" (HIP 114417) are informal.
- TESS revisits: TIC 370612437 in S121 (official dates 2027-08-23 to 09-19), S122 and S129. HIP 114417 in S122 and S129. The S122 and S129 dates are estimates.

## 8. What to send back
Write your report to **`results/review/REVIEW_<your-name>.md`**, structured as follows:

1. **Verdict per candidate:** one of (a) strong candidate, fine to submit as a CTOI; (b) weak candidate, submit with caveats; (c) likely false positive, don't submit. Give a one-paragraph justification for each.
2. **Findings**, most severe first. For each:
   - severity (critical / moderate / minor);
   - the claim affected (C1–C11);
   - file:line where relevant;
   - what you did to check it, with commands or code;
   - what you found, with numbers;
   - a concrete fix or extra test.
3. **Claims you verified independently,** and how. This matters as much as the problems.
4. **Wording changes** needed in the tweet, the CTOI draft or the video, with exact before → after text.
5. **Anything you could not check,** and why.

Keep it under ~1,500 words, plus any tables. Number every finding (F1, F2, …) so they can be answered point by point.
