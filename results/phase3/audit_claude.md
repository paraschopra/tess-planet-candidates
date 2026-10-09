# Independent Claude referee audit of Phase 3, 2026-10-09

This was a fresh Claude agent with no conversation context, run read-only. It was asked to try to show the candidates are not planets.

It ran before the revised pixel analysis (joint multi-sector fit) existed. Its statement that 252721481 had no pixel analysis refers to that earlier state. The agent's scratch scripts were in the session scratchpad and are not kept in the repo.

Its report follows verbatim.

---

**Independent referee audit: TIC 370612437 and TIC 252721481**

This was a read-only review. I changed no repo files. The earlier Codex audit (results/phase3/audit_codex.md) already covers SDE calibration, empirical-error bias, hold-out post-selection and pixel-fit issues. I did not repeat those.

**Findings**

1. **CRITICAL (370612437): a blend is not excluded.**
   - The pixel analysis detects the missing-light spot only weakly (Δχ²=23, p=0.003).
   - Its position is 13.7" from the target with about 20" uncertainty, so it does not localize the source.
   - Gaia stars at 9.5", 11.1" and 12.0" (Tmag 15.3–15.8) would need only 3–5% eclipses to produce the dip. There are 471 sources within 150" at b=-1.3°.
   - SPOC centroids show no coherent in-transit shift once three epochs next to momentum dumps or gap edges are removed. For example, S76 drops from +0.82 to +0.21 millipixels without epoch 456. At this depth that cannot rule out a source within about 1 pixel (21").

2. **MODERATE (252721481): a consistent centroid shift, and no pixel analysis has been run.**
   - MOM_CENTR2 shifts in transit by -0.35, -0.26, -0.47 and -0.40 millipixels in S17/57/77/84. Combined it is -0.28±0.08 millipixels, about 3.5σ.
   - If the dip were on the target star (CROWDSAP 0.93–0.95, i.e. 5–7% contaminating light), the expected shift is ≲0.03 millipixels. The observed shift fits a source about 1 pixel away, or pointing-correlated systematics.
   - This needs difference imaging before any report.

3. **MODERATE (370612437): the 2022 deficit is real in the data but not statistically significant.**
   - Season depths at the fitted ephemeris, with empirical errors: 186±50, 103±48, 235±42, 243±41 ppm. Consistency χ²=5.8 for 3 degrees of freedom, p=0.12.
   - The 2022 transit time is unconstrained (bootstrap peak offset +0.1 h, 16–84% range -2.8 to +2.2 h). There is no evidence for the 24-minute shift, so a TTV case is unsupported.
   - Aperture contamination cannot explain it: CROWDSAP is 0.888/0.885 in S56/57 versus 0.86–0.91 in other sectors, and the uncorrected SAP flux shows the same deficit (127, 62 ppm).
   - The pipelines disagree sector by sector: QLP gives S56 348±41 vs SPOC 171±61, and QLP S76 26±90. That means per-sector errors are underestimated. The blind leave-2022-out result fails to confirm the signal but does not refute it.

4. **MODERATE (both): the quoted SNR ≈8 is not a detection significance.**
   - In 38 scrambled versions per star, 3 reached an empirical SNR at least as high as the candidate's.
   - Flipped (upside-down) light curves with the full triage applied produce flag-free signals with empirical SNR up to 8.8.
   - The flipped light curve of 252721481 itself gives a flag-free 315 ppm "signal" at P=4.66 d with empirical SNR 7.9 (SDE 7.0).
   - Only the SDE separates the candidates from chance.

5. **Favourable (per-star false-alarm rate).** I masked each candidate's transits, gave each sector a random cyclic shift, and reran the full BLS search 38 times per star.
   - 370612437: maximum null SDE 8.53 against the candidate's 11.26. A Gumbel extrapolation gives p ≈ 1.5e-5.
   - 252721481: maximum null SDE 9.33 against 10.26, p ≈ 3e-4.
   - Multiplied by 526 stars, the expected number of false alarms this strong is about 0.01 for 370612437 and about 0.2 for 252721481.
   - These figures are uncertain by roughly 10× (extrapolated from 38 samples). The scrambles also cannot test strictly periodic systematics.

6. **Favourable (survey-wide look-elsewhere).** I applied the triage flags to the 416 flipped searches finished so far.
   - 15 signals had no hard flags. Among those with empirical SNR ≥6, the highest SDE was 8.45, and none reached 10.
   - Before triage, 27% of flipped stars had SDE ≥11.3, so selection does most of the work here.
   - Zero events in 416 caps the false-alarm rate at about 3.8 per 526 stars (95%).

7. **MINOR: momentum dumps and gaps are not driving the signals.**
   - For 370612437, only 3 of 51 transits fall within about 1 h of a dump or gap edge, and none is among the top contributors. 90% of individual transits are positive (about 88% expected); dropping the top 3 gives 163 ppm.
   - For 252721481, 95% of transits are positive; dropping the top 3 gives 216 ppm.

8. **MINOR: no shared instrumental systematic.** I stacked 88 other survey stars that share ≥2 of the sectors at each ephemeris: 4±2 ppm and -2±4 ppm.

9. **MINOR: aliases.**
   - Phase 0.5: -13±26 and 13±29 ppm. Phases 1/3 and 2/3 are also empty.
   - BLS strongly prefers P over P/2 (SNR 8.5 vs 5.7).
   - Alternating transits (testing twice the period) differ by 1.4σ (370612437: 227 vs 151±37 ppm) and 2.0σ (252721481: 313 vs 199±40 ppm).
   - An eclipsing binary at twice the period with equal eclipses, in a blended star, cannot be excluded.
   - Both candidates need impact parameter b≈0.7–0.76: transit duration is 0.67–0.74 of the central-transit value.

10. **MINOR: the 370612437 depth is overstated.** My independent depths are 162–180 ppm across four detrending methods, 189±26 ppm masked, and 182±22 ppm from SPOC's corrected flux (152±26 ppm from raw SAP). The headline 208–225 ppm is at the high end; about 185 ppm (roughly 1.0 Earth radius) is more defensible. For 252721481 the depth is robust at 241–258 ppm.

**Verdicts**

- **TIC 370612437: (b), a weak candidate needing more data.** The candidate itself holds up statistically: false-alarm odds are low, there is no shared instrumental systematic, and the signal survives other detrending methods and both pipelines. But a blend in a very crowded galactic-plane field is not excluded, pixel localization is uninformative, and the 2022 season is about 2σ low with the leave-2022-out test failing. It could be submitted as a CTOI only with explicit caveats; it does not count as well vetted. Next steps: a proper TESS-band neighbour calculation with TRICERATOPS, and 2-minute data or ground photometry at a resolution of 10" or better.
- **TIC 252721481: (b), a weak candidate needing more data.** The detection is statistically decent (about 0.2 expected false alarms across the survey) and the depth is consistent across seasons. Against it: a 3.5σ consistent centroid shift, no pixel analysis yet, the leave-2019-out search cannot recover the period, and only 5 sectors. Difference imaging is the immediate next step.
- **Neither is likely a false positive (c) on current evidence.** For both, the main risk is that the signal is real but comes from a neighbouring star, not that it is noise.
