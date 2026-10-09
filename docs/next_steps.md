# Next steps after announcing the candidates

Written 2026-10-09. The candidates are TIC 370612437.01 ("Aakanksha") and TIC 252721481.01 (HIP 114417, "Orbit"). Predictions are preregistered at doi:10.5281/zenodo.23257784.

## Calendar (set these reminders)
| When (local time) | What |
|---|---|
| **Mon 11 Jan 2027, 10:00** | Look for NASA's **TESS Cycle 10 Guest Investigator** call (ROSES). Cycle 10 runs 19 Sep 2027 – 13 Sep 2028 and includes sectors 122 and 129. If it opens, propose 2-minute cadence for both stars. Phase-1 deadlines have been around March in recent years. |
| **Mon 14 Jun 2027, 10:00** | Submit a **TESS Director's Discretionary Time (DDT)** request for 2-minute cadence in **sector 121** (starts Mon 23 Aug 2027). The hard deadline is 6 weeks before: **Mon 12 Jul 2027**. If the Cycle 10 proposal wasn't accepted, also request sector 122 (starts about 19 Sep 2027, so the deadline is about 8 Aug 2027). |
| **Mon 6 Dec 2027, 10:00** | Run `python pipeline/17_check_predictions.py`. Sectors 121 and 122 data should be public by then. If it says "no data yet", try again monthly. |
| **Mon 10 Jul 2028, 10:00** | Run it again, now including sector 129 (about 28 Mar – 24 Apr 2028; dates still estimated). |

## The order of things
**This week**
1. **Tweet,** then watch the replies. Keep saying "candidates". Check criticism the way we handled the reviews, and fix the repo when it's right.
2. **Submit the RNAAS note** (`docs/rnaas/`) at journals.aas.org, choosing *Research Notes of the AAS*. Before submitting:
   - check the affiliation and ORCID;
   - check every reference against ADS;
   - check the UAT keywords.
3. **Request ExoFOP upload access** with the "Published Candidate Upload Request" form (exofop.ipac.caltech.edu → Candidate Guidelines).
4. **Optional:** connect GitHub to Zenodo and make a release (v1.0) to get a DOI for the code.

**Next 1–2 months**
5. **When the note is published,** upload both CTOIs to ExoFOP using `docs/CTOI_DRAFT.md` plus the paper's URL. MIT's TOI team reviews them; accepted ones become TOIs, which puts them in front of the follow-up network.
6. **Reach out to a professional collaborator** (local-only list and sample email: `docs/outreach_private.md`).
7. **Optional extra analyses,** left open after the Codex review:
   - a calibrated false-alarm probability from repeated null surveys;
   - PRF-based pixel localisation;
   - another TRICERATOPS attempt.

**2027–28**
8. **Prediction test** (see the calendar).
9. **If the dips recur and ground follow-up clears the neighbours:** a validation paper, then validated planets, then a possible name proposal in a future IAU NameExoWorlds campaign.

## Questions answered

### Can TESS be pointed at these stars before August 2027?
No.
- **The pointing is fixed in advance.** TESS follows a pre-planned schedule, one strip of sky per 27-day sector, and isn't re-pointed for individual stars.
- **Neither star is in view before then.** Its next view of them is sector 121 (23 Aug – 19 Sep 2027) for TIC 370612437, then sector 122 for both, and sector 129 in spring 2028. TESS will record them anyway in its full-frame images (one frame every 200 s).
- **What a request changes.** A DDT or Guest Investigator request asks for the stars to be added to the **2-minute target list** for those sectors. NASA would then process them with its full pipeline: sharper light curves, its own transit search, and its own pixel-level centroid tests.

### How could a neighbouring star make a *periodic* dip?
Because the most common impostor is itself periodic: an **eclipsing binary**, two stars orbiting each other.
- **They eclipse on a schedule.** Each orbit, one passes in front of the other and blocks, say, 10% of their light, as regularly as a planet.
- **TESS can mix its light into ours.** If such a pair sits behind or beside our star, within one TESS pixel (21 arcseconds), TESS measures the combined light of both.
- **Dilution shrinks the eclipse.** A pair 500 times fainter than our star, eclipsing by 10%, dims the *combined* light by 10% ÷ 500 = **0.02%**: exactly our dip, with perfect regularity.
- **Our cases.** For TIC 370612437, Gaia lists faint stars within about 12" that would need only 3–5% eclipses to do this.
- **A bound companion works too.** A pair gravitationally bound to our star (a triple system) can play the same trick.

This "blended eclipsing binary" is the classic false positive for small TESS candidates, which is why neither can be called a planet yet.

### What does "ground follow-up clears the neighbours" mean?
- **Sharper telescopes.** Ground telescopes have much finer vision than TESS: about 1–2" per pixel versus 21". They can measure each nearby star's brightness separately.
- **Watch at the predicted times.** Observers point at the field at the predicted transit times (our Zenodo list) and record every neighbour's light curve.
- **A guilty neighbour is easy to spot.** The neighbour would have to dim by 3–30%, which is easy to see from the ground, even though our star's own 0.02% dip is too small for most ground telescopes.
- **Clean neighbours clear the field.** If every neighbour stays steady, they are "cleared": the dip must come from our star, or from something extremely close to it.
- **Imaging finishes the job.** High-resolution imaging (adaptive optics or speckle cameras on large telescopes) then checks for companions closer than 1".

This is routine for the **TESS Follow-up Observing Program (TFOP)**, sub-group 1 (seeing-limited photometry), mostly for official TOIs. Hence the ExoFOP and TOI steps above.

### The prediction check (`pipeline/17_check_predictions.py`)
The script:
1. **Proves it uses the preregistered predictions,** by comparing the local CSVs with the Zenodo files (MD5 checksums).
2. **Downloads only the new sectors** (121, 122, 129): 2-minute data if a request was granted, otherwise TESS-SPOC or QLP.
3. **Measures the dip** at exactly the predicted times, with honest noise estimates.
4. **Scans the whole orbit** to see where the deepest dip actually is.

The **decision rule was fixed on 2026-10-09,** before any new data exist:
- **RECURS:** dip ≥ 3σ at the predicted times, the deepest point of the orbit inside the predicted window, and depth within 2σ of the preregistered depth.
- **NOT SEEN:** dip < 1σ, even though the predicted depth would have been ≥ 3σ in those data.
- **INCONCLUSIVE:** anything else.

A dry run on the 2024 sectors (`--test 83 84`) returns RECURS for both, with the dip exactly at the predicted time. That only checks the mechanics, since those sectors were used to make the predictions.

**RECURS would show the signal is real and periodic. It cannot show which star it's on;** that's what ground follow-up is for.

Expected sensitivity:
- **TIC 370612437:** two sectors in 2027 should detect a real 210 ppm dip at about 4–5σ.
- **HIP 114417:** sector 122 alone, about 4σ; sector 129 adds more.

### Sources
- TESS DDT program: https://heasarc.gsfc.nasa.gov/docs/tess/ddt.html and https://tess.mit.edu/science/ddt
- GI program: https://heasarc.gsfc.nasa.gov/docs/tess/proposing-investigations.html
- TFOP: https://heasarc.gsfc.nasa.gov/docs/tess/tfop.html
- ExoFOP candidate guidelines: https://exofop.ipac.caltech.edu/tess/candidate_help.php
