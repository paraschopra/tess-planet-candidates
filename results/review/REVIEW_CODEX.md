# Independent review — Codex, 2026-10-09

## 1. Verdicts

**TIC 370612437: (b) weak candidate, submit with caveats.** Independent detrending and local transit measurements recover the signal. The weaker 2022 season is not independently decisive against it. However, the source remains unlocalized in a crowded field, and the global false-alarm claims are unsupported. Submission as an explicitly unvalidated transit-like candidate is reasonable after correcting the wording below; “Earth-sized” remains conditional on a planet transiting the catalogued star.

**TIC 252721481 / HIP 114417: (b) weak candidate, submit with caveats.** Independent photometry and period searching support repeatable dimming. QLP corroboration is weaker than advertised, and pixel detection is poorly calibrated. Conversely, the previously reported centroid shift does **not** establish a consistent sky displacement. Neither candidate is demonstrably a likely false positive, but neither deserves a strong-candidate verdict.

## 2. Findings, most severe first

### F1 — Critical: discovery significance remains uncalibrated (C4, C11)

**Locations:** `docs/CTOI_DRAFT.md:31`, `:60`; `pipeline/15_combined_video.py:362`; `pipeline/06_inverted.py:8`.

I counted both sets of 526 successful search records and recomputed CSV selections: real **743 → 31 → 10**, inverted **752 → 39 → 17**; maximum flag-free inverted SDE **8.3765**. This verifies the revised full-triage comparison, not the quoted **0.01/0.2** survey false alarms. Those numbers extrapolate an earlier auditor’s 38 scrambles into unobserved tails; their realizations are unavailable. Even a valid, fixed, independent zero-in-526 null gives a 95% upper expected count of approximately **2.99 per survey**, using `526*(1-.05**(1/526))`. Adaptive selection and asymmetric clipping further prevent a calibrated interpretation here.

**Fix:** remove numerical false-alarm claims; freeze selection and calibrate repeated complete null surveys, including a symmetric-clipping sensitivity test. Local depth/error ratios are not discovery significances.

### F2 — Moderate: centroid coherence disappears in sky coordinates (C8, C11)

**Location:** `docs/CTOI_DRAFT.md:61`.

`centroid.py` fits an exposure-integrated box plus local intercept/slope to cached MOM_CENTR1/2, then transforms shifts using each aperture WCS. HIP 114417’s northward shifts in S17/57/77/84 are **+6.61, +3.33, −7.11, +2.44 mas**: detector-row agreement crosses a spacecraft rotation. Combining sectors gives east/north **−0.17 ± 1.39 / −0.28 ± 1.90 mas**, with errors from 40 control phases. Including POS_CORR1/2 as nuisance regressors also gives no significant common shift.

**Fix:** describe the detector-coordinate anomaly without inferring a fixed one-pixel source; fit a common sky-source model with sector-specific apertures/PRFs. This does not clear blends.

### F3 — Moderate: revised pixel probabilities fail a direct control check (C8, C11)

**Locations:** `pipeline/10_pixels.py:184`, `:257`; `docs/GUIDE.md:302`.

`pixel_null.py` reproduces the target-position Δχ² values **23.240/20.176**, then applies precisely the same statistic to the 40 existing control phases. For HIP 114417, **2/40** exceed the real score; the largest is **636.47**. The add-one rank estimate is **3/41 = 0.073**, inconsistent with treating the analytic **p=0.00116** as calibrated. These controls are finite and correlated, so 0.073 is diagnostic, not a precise replacement probability. For TIC 370612437, **0/40** exceed its score, insufficient to verify p=0.0031. Positive amplitudes in 8/8 or 5/5 sectors are not individual detections.

**Fix:** calibrate the entire pixel statistic on more control phases and account for correlated pixels; use PRF/aperture-aware localization.

### F4 — Moderate: noise correction was not propagated to QLP (C2, C3, C5, C11)

**Locations:** `pipeline/09_products.py:62`, `:64`; `pipeline/vet.py:112`.

`check.py` reruns `depth_emp_matched` on cached QLP products: **225.37 ± 34.89 ppm (6.46)** and **237.83 ± 71.25 ppm (3.34)**, versus published errors 29/52 ppm. The same 40 controls reproduce the fit’s SPOC errors exactly. Adding 80 interleaved controls changes the combined MAD errors to **27.94/37.55 ppm**, giving local ratios **7.46/7.00**. Overlapping controls do not constitute 120 independent samples.

**Fix:** update product errors and describe the 7–8 figures as approximate local SNRs. Different reductions of the same exposures cannot exclude shared instrumental artifacts.

### F5 — Moderate: “never searched by SPOC” is false (C1, C11)

**Locations:** `docs/REVIEW_BRIEF.md:23`; `docs/GUIDE.md:297`.

[MAST’s TESS-SPOC documentation](https://archive.stsci.edu/hlsp/tess-spoc) explicitly dates FFI planet searching to May 2021/S36 and multi-sector searching to S56 data. Two-minute exclusion does not imply search exclusion. I additionally searched the public DV download manifests for S56/57/76/77/83/84 and S56–69: neither TIC appears.

**Fix:** distinguish “no matching alert found in the checked products” from “never searched”; replace “sharpest mode” with “two-minute cadence” because TESS also has 20-second targets.

### F6 — Moderate: predicted UTC columns are incorrect (C10)

**Location:** `pipeline/13_predictions.py:42`.

`timing.py` verifies BJD predictions against `t0+n*P` to the **4.32-second rounding tolerance**, but UTC columns merely format BJD as UTC. They omit TDB→UTC and barycentric light-travel corrections. Inverting these with Astropy at Earth’s center gives listed-minus-correct UTC ranges **−3.79 to +5.13 minutes** and **−5.18 to +6.51 minutes**. Example: HIP 114417’s first prediction should be approximately **2027-09-22 16:13:12 UTC**, not 16:19.

**Fix:** use an explicit observer location and invert the barycentric correction before `.utc`; fix window endpoints too. [Astropy’s timing documentation](https://docs.astropy.org/en/stable/time/index.html#barycentric-and-heliocentric-light-travel-time-corrections) distinguishes these operations.

### F7 — Moderate: retrospective hold-outs are overstated (C6, C11)

**Locations:** `pipeline/08_holdout.py:65`, `:83`; `docs/GUIDE.md:301`.

Code inspection finds training-only width refinement and no direct hidden-flux fit for these two first-signal candidates. Nevertheless, discovery used all seasons and the test was subsequently changed. Moreover, HIP 114417’s refined S77 scan peaks **21.66 hours away**, with **0.662%** of comparison phases deeper; it is not “the deepest point in the orbit.” “Three supportive seasons” is descriptive, not independent blind validation.

**Fix:** preserve both versions and reserve confirmation language for preregistered future observations.

### F8 — Minor: density agreement is weak and convergence incomplete (C7)

**Locations:** `pipeline/07_fit.py:44`, `:114`; both `*_fit.json`.

Reading saved posteriors gives free-density 16–84% intervals **1.39–12.51** and **0.29–12.12 g/cm³**. Compatibility is unsurprising across such broad ranges. `12000/autocorr_steps` is only **27.7/14.5**, below the code’s stated >50τ criterion. The stellar-prior fits do exceed 50τ; this is not evidence that their radii are wrong.

**Fix:** extend/check independent chains, save samples, and describe density as “not inconsistent, weakly constrained.”

## 3. Independently verified evidence

`search.py` uses a one-day running median, symmetric 5σ clipping, excludes 0.5-day gap edges, and searches 0.5–25 days with a different duration grid. It retains the project’s period grid/SDE definition, so this is a processing robustness check, not a wholly independent algorithm. Both periods maximize **raw BLS SNR as well as SDE**.

| Check | TIC 370612437 | TIC 252721481 |
|---|---:|---:|
| Recovered period, days | 3.56492074 | 4.93355523 |
| Recomputed SDE | 10.318 | 10.046 |
| Independent local-box depth, ppm | 194.6 ± 25.5 | 274.6 ± 34.0 |
| Local season consistency p | 0.087 | 0.477 |
| Secondary depth, ppm | −14 ± 25 | 15 ± 37 |
| Odd/even discrepancy | 1.63σ | 1.17σ |

Local fits use raw NPZ flux, ±0.3-day windows, exposure-integrated boxes at fitted ephemerides, independent linear baselines, and 3,000 event bootstrap resamples; these are conditional uncertainties. The 2022 depth is **108 ± 41 ppm**, supporting a deficit without establishing inconsistency.

Live [ExoFOP target queries](https://exofop.ipac.caltech.edu/tess/target.php?id=370612437&json), [HIP 114417’s query](https://exofop.ipac.caltech.edu/tess/target.php?id=252721481&json), and NASA Archive TOI/PS queries returned no planets/candidates. The 132 cached TCE files contain 341,351 rows and no matching TIC. I verified the density Jacobian, a/R*, inclination, exposure integration setup, BJD offset, and Phase-3/video radius provenance; the known Phase-2 dilution double count is not used by those radii. Live TIC radius errors are 8.3%/8.8%, close to the assumed 8.7%. The third candidate’s saved pixel/hold-out failures justify setting it aside.

## 4. Exact wording changes

| Location | Before → after |
|---|---|
| Tweet 1 | “TESS's sharpest observing mode” → “two-minute cadence” |
| Tweet 1 | “nobody has reported” → “absent from the TOI/CTOI and confirmed-planet lists checked on 9 Oct 2026” |
| Tweet 3 | “hid a year of data, predicted the transit times blind, found the dip there in 3 of 4 seasons” → “retrospective season tests support each signal in three seasons; one test fails for each” |
| Tweet 3 | “in 526 upside-down light curves, chance never made a signal this strong” → “no flag-free signal in one inverted 526-star survey reached either candidate’s SDE” |
| Video | “pass every\nred-flag test” → “pass every\nhard-flag test” |
| Video | “stronger than\nchance allows” → “above this\ninverted-survey sample” |
| Video | “hid a season, predicted it blind: 3 of 4 hits” → “retrospective season tests: three supportive, one fails” |
| Video | “signal strength 7.2σ (#1) and 7.7σ (#2), red-noise aware” → “local depth SNR roughly 7–8; global false-alarm probability uncalibrated” |
| Video | “in 526 upside-down light curves, chance never got this strong” → “one inverted survey produced no flag-free signal reaching either SDE” |
| CTOI, both | “**Recovered independently:**” → “**Seen in a second reduction of the same exposures:**” |
| CTOI | “Blind leave-one-season-out predictions recover it in 3 of 4 seasons” → “Retrospective leave-one-season-out tests support it in three of four seasons”; “**Blind hold-out predictions:**” → “**Retrospective season checks:**” |
| CTOI | “225±29 ppm” / “238±52 ppm” → “225±35 ppm” / “238±71 ppm” |
| CTOI, both | “**Survey-wide false alarms:** about 0.01 expected from scrambled and inverted light curves.” / “**Survey-wide false alarms:** about 0.2 expected.” → “**Survey-wide false alarms:** not reliably calibrated; neither SDE was reached by a flag-free signal in one inverted survey.” |
| CTOI, 370612437 | “Difference imaging detects the dimming near the target (8/8 sectors) but cannot localise it.” → “Target-position difference-image amplitudes are positive in all eight sectors, but the source is not localized.” |
| CTOI, 252721481 | “SPOC flux-weighted centroids show a small consistent in-transit shift (−0.28 ± 0.08 millipixel, about 3.5σ), and difference-imaging localisation is inconclusive. A blended source about 1 pixel away is possible.” → “Detector-row centroid shifts were reported, but a check accounting for spacecraft rotation finds no significant common sky displacement. Pixel localization is inconclusive; a blend remains possible.” |
| CTOI, 370612437 | “The transit-shape stellar density agrees with the TIC.” → “The weakly constrained transit-shape density is not inconsistent with the TIC.” |

## 5. Scope and reproducibility

No existing analysis, draft, or media was changed; nothing was committed or submitted. Reproduction scripts/results are in `results/review/scratch/codex/`. From the repo root, run `.venv/bin/python results/review/scratch/codex/<script>.py`; scripts write only there. `check.py`, `pixel_null.py`, `centroid.py`, `secondary.py`, `search.py`, and `timing.py` correspond to the checks above; `novelty.json` preserves live responses and DV-manifest hashes.

I did not rerun the full surveys, MCMC, TRICERATOPS, or a new global bootstrap. I did not establish exhaustive ADS full-text/RAVEN-table novelty, independently reconstruct HGCA, verify every momentum-dump association, or decode the MP4; wording review used its renderer. The official S121 dates match the [TESS schedule](https://heasarc.gsfc.nasa.gov/docs/tess/sector.html); I did not independently recompute footprint coverage or establish S122/S129 dates. No test here excludes an unresolved eclipsing binary.
