# Draft CTOI submissions (not yet submitted)

> **Cannot be uploaded yet.** Since 2026-08-19, ExoFOP accepts community candidates only after they are published: in a peer-reviewed journal, or in an RNAAS that cites previously published detection and vetting methods. Upload access must first be approved via the Published Candidate Upload Request (https://exofop.ipac.caltech.edu/tess/candidate_help.php). Each upload needs the paper's URL. See `docs/GUIDE.md` Part 5. Everything below is ready for when that is done, and can also serve as the core of the RNAAS note.

These are for submission by the project owner through an ExoFOP-TESS account. Map the values below onto the form's fields. Times are BJD_TDB, depths are in ppm and durations are full transit durations (T14).

Both signals come from TESS-SPOC full-frame-image light curves (PDCSAP, which is already corrected for crowding) and were cross-checked with QLP. The search, vetting and audit code is in this repository. Results and figures are in `results/phase3/`, and the summary is in `docs/GUIDE.md` Parts 3–4.

---

## TIC 370612437.01

| Field | Value |
|---|---|
| TIC | 370612437 (TYC 3968-348-1, Gaia DR3 1981360227180249088) |
| Period (d) | 3.564921 ± 0.000022 |
| Epoch (BJD_TDB) | 2458742.1409 ± 0.0069 |
| Depth (ppm) | 210 (+45 / −49) |
| Duration (h) | 1.44 (+0.28 / −0.14) |
| Rp/R* | 0.0144 |
| Impact parameter | 0.73 (+0.08 / −0.18) |
| Planet radius (R⊕) | 1.04 ± 0.15 |
| Equilibrium temperature (K) | ~850 (zero albedo) |
| Insolation (S⊕) | ~87 |
| SNR | ~7.2 local depth/error (empirical red-noise error; not a discovery probability); BLS SDE 11.3 |
| Sectors | 16, 17, 56, 57, 76, 77, 83, 84 (51 transits) |

**Notes (draft):**
> Small transit-like signal found in a multi-sector BLS search of TESS-SPOC FFI light curves of 526 nearby (<60 pc) dwarf stars never observed at 2-minute cadence.
>
> - **Seen in a second reduction of the same exposures:** QLP shows 225±35 ppm.
> - **Retrospective season checks:** leave-one-season-out tests support it in 3 of 4 seasons. The 2022 season (S56–57) is weaker (≈105±45 ppm) and its hold-out test fails.
> - **Clean diagnostics:** odd/even depths agree (1.4–1.6σ), no secondary eclipse, no significant centroid shift. The weakly constrained transit-shape density is not inconsistent with the TIC.
> - **Survey-wide false alarms:** not reliably calibrated. In one inverted (upside-down) 526-star survey through the same pipeline, no flag-free signal reached this SDE (max 8.4).
> - **Main caveat:** this is a crowded galactic-plane field (b = −1.3°). Target-position difference-image amplitudes are positive in all 8 sectors, but the source is not localised. Several Gaia neighbours within 12–22" (T ≈ 15–18) could produce the signal with 3–30% eclipses.
> - **TRICERATOPS:** not completed. Needed once ground-based follow-up exists.
> - **Follow-up needed:** seeing-limited photometry to clear neighbours, and high-resolution imaging.
> - **Not in checked alert products:** no SPOC TCE in the 2-minute TCE lists, and no TESS-SPOC FFI DV report in the S56, S57, S76, S77, S83, S84 and S56–S69 manifests.
> - **Next TESS coverage:** S121, S122, S129. Predicted transit times are in the repository.

## TIC 252721481.01

| Field | Value |
|---|---|
| TIC | 252721481 (HIP 114417, G 216-30, Gaia DR3 1991778271640247808) |
| Period (d) | 4.933533 ± 0.000019 |
| Epoch (BJD_TDB) | 2458740.6616 ± 0.0060 |
| Depth (ppm) | 279 (+52 / −48) |
| Duration (h) | 1.59 ± 0.15 |
| Rp/R* | 0.0167 |
| Impact parameter | 0.68 (+0.10 / −0.14) |
| Planet radius (R⊕) | 1.10 ± 0.14 |
| Equilibrium temperature (K) | ~716 (zero albedo) |
| Insolation (S⊕) | ~44 |
| SNR | ~7.7 local depth/error (empirical red-noise error; not a discovery probability); BLS SDE 10.3 |
| Sectors | 16, 17, 57, 77, 84 (19 transits) |

**Notes (draft):**
> Small transit-like signal in TESS-SPOC FFI light curves (same survey as TIC 370612437.01).
>
> - **Seen in a second reduction of the same exposures:** QLP shows 238±71 ppm (weak).
> - **Retrospective season checks:** leave-one-season-out tests support it in 3 of 4 seasons. Withholding 2019 leaves too little data to recover the period.
> - **Clean diagnostics:** no secondary eclipse; odd/even differ by 2.0σ. Gaia RUWE is 1.03, and the Hipparcos–Gaia proper-motion anomaly shows no companion (HGCA χ² = 0.35).
> - **Survey-wide false alarms:** not reliably calibrated. No flag-free signal in one inverted 526-star survey reached this SDE.
> - **Main caveat:** detector-row centroid shifts were noted, but a check that accounts for spacecraft rotation finds no significant common sky displacement. Pixel localisation is inconclusive, so a blend remains possible.
> - **TRICERATOPS:** not completed.
> - **Follow-up needed:** neighbour clearing and high-resolution imaging.
> - **Not in checked alert products:** no SPOC TCE in the 2-minute TCE lists, and no TESS-SPOC FFI DV report in the S56, S57, S76, S77, S83, S84 and S56–S69 manifests.
> - **Next TESS coverage:** S122, S129. Predicted transit times are in the repository.

---

### Before submitting
- **Publish first** (RNAAS or journal) and get upload approval (see the note at the top). Include the paper URL in the upload.
- Check that the ExoFOP CTOI form's fields and units match the values above, and that neither star has received a TOI or CTOI since 2026-10-09.
- Decide on the **"user" and "group" credit**. The form records who identified the signal.
- Optionally **preregister** first: publish `results/phase3/*_predicted_transits.csv` (for example on Zenodo, with a DOI) *before* sectors 121/122 are observed (from 2027-08-23). That makes those sectors an untouched test.
