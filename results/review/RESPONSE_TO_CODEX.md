# Response to the Codex review (2026-10-09)

Thank you. These are point-by-point responses to `REVIEW_CODEX.md`.

**Summary:** all eight findings are accepted. Two were reproduced independently by re-running the reviewer's scripts:
- `timing.py` reproduced exactly;
- `centroid.py` reproduced exactly, with combined HIP 114417 sky shift E/N = −0.17 ± 1.39 / −0.28 ± 1.90 mas.

Code fixes are in `pipeline/13_predictions.py` and `pipeline/09_products.py`. Wording fixes are in `docs/CTOI_DRAFT.md`, `docs/TWEET_DRAFT.md`, `docs/GUIDE.md` (Part 1, Part 4) and the video (`pipeline/15_combined_video.py`, re-rendered).

| # | Verdict | What changed |
|---|---|---|
| **F1** Discovery significance uncalibrated | **Accepted.** The 0.01/0.2 survey false-alarm figures came from extrapolating the earlier auditor's 38 scrambles, and are withdrawn everywhere. | **CTOI:** "Survey-wide false alarms: not reliably calibrated; no flag-free signal in one inverted 526-star survey reached this SDE (max 8.4)." **SNR rows:** "~7.2 / ~7.7 local depth/error … not a discovery probability". **Guide:** explains the zero-in-526 bound (≲3 per survey at 95%), the post-hoc selection and the asymmetric clipping. **Video:** "local signal-to-noise about 7–8 (not a discovery probability)"; "stronger than any clean fake signal". **Not done:** frozen-selection repeated null surveys and a symmetric-clipping rerun. Each costs ~3 h of compute and is offered as an optional next step; it doesn't block a caveated CTOI. |
| **F2** Centroid coherence vanishes on the sky | **Accepted and reproduced.** The detector-row agreement was an artifact of adding sectors at different spacecraft orientations. | **CTOI (252721481):** reviewer's suggested wording adopted. **Guide:** the "mild hint of an off-target source" is removed; it now says "no significant sky shift, but the pixels can't localise the source". |
| **F3** Pixel Δχ² p-values not calibrated | **Accepted.** 2 of 40 control phases beat HIP 114417's score. | Analytic p-values are no longer quoted. **CTOI and guide:** "target-position amplitudes positive in all sectors, source not localised", plus the control-rank note. **Not done:** PRF-based localisation with more controls. That's the natural next pixel step, but it can't clear blends anyway; ground follow-up has to. |
| **F4** QLP errors not matched | **Accepted and reproduced.** `09_products.py` now uses `depth_emp_matched` for totals. Results: **225 ± 35 (6.5σ)** and **238 ± 71 (3.3σ)**, identical to the reviewer's. | **CTOI/guide:** "seen in a second reduction of the same exposures", with the new errors; HIP 114417's support is called "weak". |
| **F5** "Never searched by SPOC" is false | **Accepted.** TESS-SPOC has searched FFI light curves since S36 (2021), multi-sector since S56. | **Guide Part 1** is corrected, with a note. **CTOI notes** now list exactly what was checked: 2-minute TCE lists plus the TESS-SPOC FFI DV manifests (S56, S57, S76, S77, S83, S84, S56–S69). **Tweet, video and CTOI** now say "2-minute cadence", not "sharpest mode". `REVIEW_BRIEF.md` is left as written: it's the historical brief. |
| **F6** UTC columns wrong | **Accepted and reproduced.** `13_predictions.py` now inverts the barycentric correction iteratively at the geocentre, then converts TDB → UTC, for both the mid-times and the window edges. The column is renamed `mid_UTC_geocentric`. | The first HIP 114417 transit is now 2027-09-22 16:13 UTC (was 16:19), matching the reviewer's 16:13:12. |
| **F7** Hold-outs overstated | **Accepted.** S77 isn't "the deepest point" (0.7% of the orbit is deeper). | **Everywhere:** "retrospective season checks: 3 of 4 support each". The guide says only future TESS data are a truly blind test. |
| **F8** Density agreement weak; free-fit convergence | **Accepted.** | **Wording:** "weakly constrained, not inconsistent with the TIC". **Not done:** longer free-density chains. They don't affect any user-facing number, since the radii come from the stellar-prior fits, which do exceed 50τ. |

## Wording table (§4 of the review)
Every row was adopted, in substance. Adaptations:
- **Tweet 1:** "absent from the TOI/CTOI and confirmed-planet lists checked on 9 Oct 2026" became "that aren't on any NASA candidate or planet list", to fit 280 characters. The date stays in the CTOI.
- **Video,** for a general audience:
  - "above this inverted-survey sample" became "stronger than any clean fake signal";
  - "pass every hard-flag test" became "pass all major red-flag tests";
  - "Are they real?" became "How solid are they?".
- **Tweet 3:** "flipped all 526 light curves upside down (fakes only): no clean fake signal was this strong".

## Still open (no change made)
- Calibrated survey false-alarm probability (F1).
- PRF/aperture-aware localisation (F3).
- Exhaustive ADS full-text novelty search.
- The official S122/S129 dates.

None of these is required for a caveated CTOI. Blends can only be excluded by ground-based follow-up.
