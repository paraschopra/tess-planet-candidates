# Planet hunting in TESS data: a plain-language guide

This project tries to find a planet nobody has reported yet, using free public data from NASA's TESS space telescope. It is inspired by a Reddit post (r/ClaudeAI, Oct 2026) in which someone used Claude Code to find a ~1.4 Earth-radius planet candidate orbiting a star 116 light-years away, with a 0.05% dip in brightness every 3.18 days.

This guide explains the ideas with no astronomy background assumed, then describes each phase of the project.

---

## Part 1: the concepts

### Transits: how you find a planet you can't see
When a planet passes in front of its star (as seen from Earth), it blocks a tiny fraction of the star's light. The star looks slightly dimmer for a few hours. This is called a **transit**. Because the planet orbits, the dimming repeats at a fixed interval, like a clock ticking.

We never see the planet itself, only its shadow on the star's brightness. From the dip we can read off:
- **How often it repeats**: the planet's year (its **orbital period**).
- **How deep the dip is**: how big the planet is compared to the star. A planet covering 1% of the star's disk dims it by 1%.
- **How long the dip lasts**: tied to the orbit size and the star's size.

Small planets make very small dips. An Earth-sized planet around a Sun-like star dims it by about 0.008%. Around a small, cool star (an "M dwarf", about a third of the Sun's width) the same planet dims it ten times more. That is why small stars are good places to look for small planets.

### TESS
**TESS** (Transiting Exoplanet Survey Satellite) is a NASA space telescope launched in 2018. It is essentially a set of cameras that repeatedly measures how bright stars are, looking for these tiny repeating dips.

### Light curve
A **light curve** is a graph of one star's brightness over time: time on the horizontal axis, brightness on the vertical. A planet shows up as a small notch that repeats at regular intervals. Everything in this project works on light curves.

### Sectors
TESS does not watch the whole sky at once. It stares at one strip of sky for about **27 days**, then turns to the next strip. Each 27-day stare is a **sector**. Over the years TESS returns to the same parts of the sky, so one star might be watched in 2018, 2020 and 2025: three separate sectors. Stars near the poles of the sky get watched for many sectors; others only once or twice.

### Full-frame images vs. 2-minute targets
TESS records data in two ways:

- **Full-frame images (FFI).** TESS photographs its whole patch of sky at regular intervals: every 30 minutes in 2018–2020, every 10 minutes in 2020–2022, and every 200 seconds since. *Every* star in view is recorded this way. Light curves are extracted from these photos afterwards.
- **2-minute targets.** TESS cannot send all of its data home at high speed, so astronomers choose a shortlist of about 20,000 "interesting" stars per sector. For those, brightness is saved every 2 minutes. These light curves are sharper and these stars get the most scrutiny.

The Reddit poster's star was never on the 2-minute shortlist. It only existed in the coarser full-frame photos. After finding the signal, they asked NASA to add it to the 2-minute list, and NASA approved it.

### Who has already searched the data
NASA and MIT run automated programs that scan light curves for repeating dips:
- **SPOC** (Science Processing Operations Center) is NASA's main pipeline. It searches the 2-minute stars carefully. It also produces light curves from full-frame images for many brighter stars (called **TESS-SPOC** light curves), and since 2021 (sector 36) it has searched those too, including multi-sector searches from sector 56. (An earlier version of this guide wrongly said FFI stars were never searched by SPOC; the independent Codex review corrected this.)
- **QLP** (Quick-Look Pipeline, MIT) searches full-frame-image stars, but with stricter cut-offs, so faint signals can slip through.

When these programs flag a signal it becomes a **TCE** (Threshold Crossing Event, an automatic alert). If humans think it is promising it becomes a **TOI** (TESS Object of Interest, an official candidate). Amateurs and other teams can submit **CTOIs** (Community TOIs). A candidate becomes a confirmed planet only after other telescopes follow it up.

### Multi-sector searches
A **multi-sector search** stitches together all of a star's sectors, from different years, before looking for dips. This matters enormously for small planets. A dip too faint to notice in a single 27-day window can become obvious when a dozen transits from several years are lined up and added together. This is exactly how the Reddit poster's signal became convincing: the 2025 data looked suspicious, and the 2020 and 2018 data showed the same dip.

These combined searches are run only periodically, so the newest data (2025–26) has often not yet been combined with older years.

### Our strategy
A planet nobody has found is most likely hiding around a star that:
1. TESS has watched several times over several years, so enough transits exist to add up;
2. was **never a 2-minute target**, so it only got the quicker, less sensitive full-frame search;
3. has **recent data** not yet combined with older data.

### Signal strength: SNR and SDE
- **SNR** (signal-to-noise ratio): how big the dip is compared to the random jitter in the data. Around 7 is the usual minimum to take a signal seriously; 10+ is comfortable.
- **SDE** (signal detection efficiency): how much the best orbit stands out compared to all the other orbits the search tried. A high SDE means one period is clearly special rather than the search just picking the biggest bump in noise.

### Reading the diagnostic pictures
- **Periodogram, and "one sharp peak".** The search tries about a million possible orbit lengths. For each one it lines up the data at that rhythm and scores how much of a dip appears. The periodogram is a plot of score against orbit length. A real planet gives one tall, narrow spike at its true period, standing well above a flat "grass" of noise. Trouble looks different: several similar spikes, a broad hump, or spikes at many related periods (half, double). Those suggest noise, star spots or instrument rhythms.
- **Folded transit, and "box-shaped dip".** Once a period is chosen, every transit is stacked on top of the others ("folding"), which beats down the noise. A planet crossing a star gives a dip with a quick drop, a flat bottom and a quick rise, like a shallow box (slightly rounded, because a star's edge is dimmer than its centre). A star's slow brightness wave gives a smooth, rounded sag instead. An eclipse by a star that only grazes the edge gives a pointed "V".
- **Odd/even test, and "matching alternate dips".** Number the transits 1, 2, 3, 4… and stack the odd ones separately from the even ones. A planet makes the same dip every time, so both stacks match. The most common impostor is two stars orbiting each other, seen at twice the true period. Each star passes in front of the other once per orbit, and the two stars usually differ in size or brightness, so the eclipses alternate deep, shallow, deep, shallow. Unequal odd and even dips give it away.
- **Secondary eclipse.** Halfway between transits, the planet passes *behind* its star. A small planet gives off almost no light of its own, so nothing measurable happens. A companion *star* does glow, so when it is hidden behind its partner the total light drops a little. A small dip exactly halfway between the main dips (phase 0.5) means a star, not a planet.

---

## Part 2: the phases

### Phase 0: setup (done)
Python environment with standard astronomy tools:
- `lightkurve` downloads and handles TESS light curves.
- `astroquery` queries star catalogs.
- `wotan` removes slow trends from light curves.
- `transitleastsquares` (TLS) is a sensitive transit search.
- `batman` and `emcee` model and fit transit shapes.

All code is in `pipeline/` and versioned with git.

### Phase 1: choosing which stars to search (done)
Script: `pipeline/01_targets.py`. Output: `results/targets_ranked.csv`.

| Step | Stars |
|---|---|
| Nearby stars from Gaia (within 60 parsecs ≈ 200 light-years) bright enough for TESS | 32,434 |
| Ordinary (dwarf) stars no bigger than 1.5× the Sun | 28,089 |
| Not already a known planet host, TOI or CTOI | 27,310 |
| Has any TESS light curve | 26,706 |
| Watched in at least 3 sectors spanning at least 2 years | **20,123** |

Each star is scored by how detectable a 1.5 Earth-radius planet on a 10-day orbit would be (expected SNR). About 10,900 stars would show such a planet at SNR ≥ 10.

Being detectable is not the same as being likely undiscovered. The top-ranked stars have 2-minute data and have already been searched deeply. The better hunting ground:
- **1,012 stars** were never 2-minute targets, yet a 1.5 Earth-radius planet would still be detectable around them.
- **583** of these have 2025–26 data.
- Next tier: 2,981 stars with only 1–2 sectors of 2-minute data.

Notes:
- Filtering the TESS Input Catalog directly by distance is a full scan of 1.7 billion rows and times out. We start from the Gaia Catalogue of Nearby Stars instead and look up each star's TESS ID.
- Snapshots of the known-planet, TOI and CTOI lists (as of 2026-10-08) are in `data/catalogs/`.

### Phase 2: searching for repeating dips
Target: the 1,012 full-frame-only stars. For each star:

1. **Download** its light curves from every sector. These are brightness readings NASA already extracted from the full-frame photos. TESS-SPOC is used where available, QLP otherwise.
2. **Clean them:**
   - Stars naturally flicker and drift (starspots rotating into view, for example), and the telescope has quirks. A smoothing filter removes slow wobbles while leaving short, sharp dips intact.
   - Bad readings are dropped, such as moments when the telescope was adjusting itself.
   - Sudden upward spikes (stellar flares, common on small stars) are clipped. Downward points are never clipped, because those could be transits.
3. **Stitch all years together** into one long record. This is our own multi-sector search, including the new data.
4. **Search for repeating dips:**
   - Try every plausible orbital period, from half a day to a few tens of days, and ask whether there is a dip every X days that lines up across all years.
   - A fast method (BLS, Box Least Squares) does the first pass.
   - A slower, more precise method (TLS, Transit Least Squares) re-checks only the strongest signals.
5. **Mask and repeat:** hide any signal found and search again, to catch additional planets around the same star.
6. **Output:**
   - A table of signals (period, depth, duration, SNR, SDE, number of transits).
   - For each strong signal, a picture showing the full light curve and all the dips folded on top of each other.

Expectation: most strong signals will turn out to be impostors. Phase 2 only produces suspects.

Cost: ~3–6 GB of disk, a few hours of computer time, no money (all data is free; no paid AI API calls).

### Phase 3: trying to prove each signal is *not* a planet
The guiding principle, from the Reddit post: don't ask "is this a planet?", keep finding ways to prove it isn't.

1. **Two stars eclipsing each other**, the most common impostor. Two stars in orbit each pass in front of the other, giving two dips per orbit, often of different depths.
   - *Odd/even test:* compare alternate dips. Unequal depths mean two stars.
   - *Secondary eclipse test:* look for a faint dip exactly halfway between transits. A planet passing behind its star is nearly invisible; a star is not.
2. **Shape and timing.**
   - Planets make flat-bottomed "U" dips; grazing star eclipses make "V" dips.
   - The dip's duration must match what physics predicts for this star's size and mass.
3. **Size sanity check.** The object's size follows from dip depth and the star's size. Anything bigger than about 2× Jupiter is a star.
4. **Is it really this star?** TESS's pixels are huge (21 arcseconds), so many faint background stars can share a pixel, and a background eclipsing pair can leak its dips into our star's light.
   - Measure where on the image the dimming comes from: does the centre of light shift during the dip?
   - Check every neighbouring star in Gaia's map to see whether any could produce the dip.
5. **Telescope glitches.**
   - Do dips line up with thruster firings ("momentum dumps"), data gaps, or scattered sunlight or moonlight?
   - Does the same period show up in many unrelated stars? That would mean it is the telescope, not the star.
6. **Consistency.** The dip should appear separately in every year with the same depth, and in both independent light-curve products (TESS-SPOC and QLP).
7. **The prediction test.** Hide one year of data, fit the orbit from the remaining years, predict exactly when the hidden year's dips should occur, and check only at those times. Repeat for every year.
8. **Has someone already found it?** Check TOIs, CTOIs, confirmed planets, the automated alert (TCE) lists from SPOC and QLP, and the published literature.
9. **Independent audit.** Codex and fresh Claude sessions with read-only access try to break the analysis.

Signals that survive everything are candidates worth reporting.

### Later phases (if something survives)
- **Phase 4: careful modelling.** Fit a physical transit model to get the planet's size, orbit and uncertainties. Run statistical validation (e.g. TRICERATOPS), which estimates the probability of planet versus each impostor scenario.
- **Phase 5: reporting.** Write it up and submit a CTOI to ExoFOP, where professional astronomers coordinate follow-up. Optionally publish predicted transit times in advance (preregistration on Zenodo) and propose 2-minute observations to TESS.

Only follow-up by other telescopes, measuring the planet's mass or confirming the dips from the ground, can turn a candidate into a confirmed planet.

---

## Part 3: Phase 2 results (2026-10-08)

### What changed from the plan, and why

**1. Most "never 2-minute" stars at the top of the first ranking were fakes.** Each star has a *contamination ratio*: how much light from neighbouring stars falls into the same TESS pixels, relative to the star's own light. The top-ranked stars had ratios of 68–222. They were faint companions sitting right next to stars ~100× brighter, so their light curves were mostly the neighbour's light. That is probably *why* they were never chosen for 2-minute observation. We now:
- reject stars with contamination ≥ 0.2;
- divide each star's expected dip depth by (1 + contamination), because neighbours' light dilutes any dip.

This left **526 target stars**, mostly K dwarfs (orange stars somewhat smaller than the Sun).

**2. Some sectors are junk.** MIT's QLP light curves were 30–100× noisier than expected for some stars, for example a star next to the Large Magellanic Cloud, a crowded neighbouring galaxy. Each sector's noise is now measured, and sectors far noisier than the star's typical sector are dropped before searching.

### How we tested the search: planting fake planets

The search could look like it is working when it isn't. To check, we used a technique called **injection-recovery**: plant artificial planets into real data, then see whether the search finds them.

1. **Take a real star's real light curve,** with all its genuine noise, flicker, gaps and instrument quirks.
2. **Choose a fake planet:** a size (e.g. 1.5 Earth radii), a random orbital period (1–15 days) and a random starting time.
3. **Work out what its transits would look like:**
   - The depth follows from planet size versus the star's size, diluted by any neighbouring stars' light.
   - The duration follows from the orbital period and the star's density, using Kepler's laws.
4. **Subtract that dip** from the real brightness measurements at every moment the fake planet would be in front of the star.
5. **Run the exact same pipeline** (cleaning, smoothing, searching) without telling it anything about the plant.
6. **Check the answer.** A plant counts as recovered only if the search finds the right period (within 0.1%), the right timing, and a strong enough signal.

Because the planted planet's true properties are known, this measures the pipeline's real sensitivity on real, messy data, not idealised noise. We used it three ways:

- **Choosing the smoothing window.** The filter that removes slow stellar wobbles can also eat part of a transit. Using the same 16 plants (8 stars × 2), we tried four window lengths:

  | Window | Recovered |
  |---|---|
  | 0.5 d | 3/16 |
  | 0.75 d | **10/16** |
  | 1.0 d | 8/16 |
  | 1.5 d | 5/16 |

  We kept 0.75 days, chosen by evidence rather than guesswork.
- **Finding hidden signals.** Some plants were missed because the star already contained a *different*, stronger signal. The same unknown signal came up in both plants, which flagged those stars as interesting.
- **Calibrating the red-flag tests.** A test that rejects impostors must not reject real planets. We ran every red-flag test on planted planets, which by construction have equal odd/even dips and no secondary eclipse. Two tests fired on 15–25% of genuine plants. The cause: their error bars assumed independent noise, but real stars and instruments produce correlated wiggles lasting hours. We replaced those error bars with an **empirical** noise estimate (place the same box at 40 random off-transit times and measure the spread). Misfires dropped to 5%. Every flag threshold was checked this way.

### The search
- 526 stars, 1–44 sectors each, spanning up to 7.5 years.
- About a million trial periods per star, from 0.5 to 25 days.
- Up to 3 signals per star: find a signal, mask it, search again.
- About 1.5 hours on 12 cores.

### Triage: automatic red flags
743 signals passed the basic thresholds (SDE ≥ 7, SNR ≥ 6, ≥ 3 transits). Each was then tested for impostor signatures.

**Hard flags** are direct evidence against a planet:
- odd/even mismatch or a secondary eclipse (an eclipsing star);
- the same period re-found after masking (stellar variability);
- the dip is too long for the orbit;
- one or two events dominate (glitches);
- the dip sits at data-gap edges;
- per-sector depths disagree;
- the dip is part of a broad depression;
- the same period appears in ≥ 3 other stars (instrument systematic);
- the dip is too deep to be a planet.

**Soft flags** are context that makes a signal suspicious but isn't proof:
- the period falls in a crowded range (huge pile-ups at 12–25 d near TESS's 13.7-day orbit, and at 0.5–1 d from stellar variability);
- the period is near the TESS orbit;
- the signal is weak under realistic noise.

Results: **743 → 31 with no hard flags → 10 with no flags at all.** With ~1,500 searches, signals with SDE ≈ 7–8 are within reach of chance, so only the two with SDE > 10 stand out.

### Lessons from looking at the pictures
- `figures/306635397_s1.png`: an apparent Earth-sized planet with 42 transits. The per-sector panel shows that almost all the signal comes from two glitches, in sectors 87 and 97. This led to the "few events dominate" flag.
- `figures/258748095_s2.png`: the dip sits inside a broad, wavy depression, and its depth swings between sectors. That is stellar variability, and it led to the "broad depression" and "inconsistent sectors" flags.

### The two strongest candidates (Phase 2 view; see Part 4 for the Phase 3 vetting)

| | TIC 370612437 | TIC 252721481 |
|---|---|---|
| Period | 3.5649 d | 4.9336 d |
| Depth | ~204 ppm (0.02%) | ~252 ppm (0.025%) |
| Duration | ~1.1 h | ~1.5 h |
| Implied size | ~1.1 Earth radii | ~1.1 Earth radii |
| Star | K dwarf, 4380 K, 0.66 R☉, Tmag 9.8 | K dwarf, 4301 K, 0.61 R☉, Tmag 9.7 |
| Distance | 48 pc (~158 ly) | 41 pc (~134 ly) |
| Data | 8 sectors, Sep 2019 – Oct 2024, 55 transits | 5 sectors, Sep 2019 – Oct 2024, 21 transits |
| SDE / SNR | 11.3 / 8.4 | 10.3 / 8.1 |

Both show:
- a single clean periodogram peak;
- a box-shaped folded dip;
- consistent depth in *every* sector across five years;
- matching odd/even dips;
- no secondary eclipse.

The two stars are 12° apart on the sky, so they cannot contaminate each other. They share sectors only because TESS observed that region at the same times.

Caveats to resolve in Phase 3:
- 370612437's dip is about half as long as a central transit would be. That suggests a grazing orbit, or wrong stellar parameters.
- 370612437 has 12% light contamination from neighbours.
- Neither has had pixel-level, catalog or literature checks yet.

Figures: `figures/370612437_s1.png`, `figures/252721481_s1.png`.

---

## Part 4: Phase 3 results (2026-10-09)

**Bottom line:** both signals survived every attempt to explain them away as noise, instrument glitches or processing artifacts. Neither could be pinned to its own star. TESS's pixels are too coarse to rule out a faint neighbouring pair of eclipsing stars, so both stay **weak candidates that need more data**. They are not yet planets, and not yet ready to report as well-vetted candidates. Two independent auditors (OpenAI's Codex and a fresh Claude session) reached the same verdict.

### Why only 2 of the 10 flag-free signals?
The search tried about a million orbits per star on 526 stars, about 1,500 searches counting repeat passes. Across that many tries, chance alone produces some impressive-looking signals. We measured how impressive with an **upside-down test**:
- Flip every light curve so that dips become bumps and bumps become dips.
- Run the identical search and the identical red-flag triage.

A real planet only makes dips, so anything found in flipped data is noise, star spots or instrument behaviour, by construction.

| | Real data | Upside-down data (all fake by construction) |
|---|---|---|
| Signals passing the search | 743 | 752 |
| No hard red flags | 31 | 39 |
| No flags at all | 10 | 17 |
| Strongest flag-free signal | **SDE 11.3 and 10.3** (our two) | SDE 8.4 |

Chance produces as *many* clean-looking signals as the real data, but never one as *strong*: across all 526 flipped stars, no flag-free fake reached SDE 9.

The other 8 flag-free signals have SDE 7.0–8.1, inside the range the flipped data produce, so we can't tell them apart from chance. The two leaders, at SDE 11.3 and 10.3, stand above anything the flipped survey produced once triage is applied.

**One more signal checked, then set aside (watch list): HIP 93326 (TIC 350883383), 14.93 days.** The flipped test showed three real signals with SDE ≥ 10 that carried only soft flags. Two have dips only 2.8–3.2× the realistic noise, which is stellar wobble. The third looked good: a box-shaped 245 ppm dip, matching alternate dips, no secondary eclipse, consistent sectors, and SDE 10.7. So it got the full Phase 3 treatment, and it did not hold up:
- the dimming is **not detected in the pixels** (p = 0.35);
- with only 8 transits, the hide-a-season predictions mostly failed (one 2.7σ hit, one season too weak to search, one prediction 3 hours off);
- its period sits in the pile-up of artifacts near TESS's 13.7-day orbit.

It stays on a watch list for the next TESS data. Rejecting it is the tests doing their job.

A second check scrambled each candidate's own light curve:
- mask out the transits;
- slide each sector by a random amount;
- search again, 38 times.

The best fake signals reached SDE 8.5 and 9.3, against the real 11.3 and 10.3.

**No calibrated false-alarm probability.** An earlier draft extrapolated these 38 scrambles to about 0.01 and 0.2 expected false alarms per survey. The second independent review rejected that: 38 samples can't support a tail extrapolation. Even zero fakes in 526 flipped stars only caps the rate at about 3 per survey (95%). On top of that, the selection rules were partly designed after looking at real signals, and the clipping is asymmetric, so flipping isn't a perfect null. The honest statement: **in one upside-down survey through the identical pipeline, no flag-free fake reached either candidate's SDE.** The 7–8 figures are local depth-to-noise ratios, not discovery probabilities.

### The tests and what they found

| Test (plain language) | TIC 370612437 (P = 3.565 d) | TIC 252721481 (P = 4.934 d) |
|---|---|---|
| **Already known?** TOI, CTOI, ExoFOP, SIMBAD literature, and all 341,351 NASA automatic alerts (TCEs) | Nothing. Neither star is in the 341,351 two-minute TCEs, nor in the TESS-SPOC FFI search (DV) manifests for the relevant sectors (checked by the Codex review). | Nothing. The star (HIP 114417) is on NASA's Habitable Worlds Observatory target lists, but no planet has been reported. |
| **Hidden companion star?** Gaia wobble (RUWE) and binary flags; Hipparcos-to-Gaia motion change | None seen (RUWE 0.89). | None seen (RUWE 1.03; no change in motion between 1991 and 2016). |
| **Physical model fit:** a realistic transit shape, with 30-min exposures smoothed and correlated noise included | 1.04 ± 0.15 Earth radii, tilted orbit (b ≈ 0.73), 1.4 h duration. The shape alone gives a weakly constrained star density, not inconsistent with the catalog K dwarf. | 1.10 ± 0.14 Earth radii, b ≈ 0.68, 1.6 h duration. Density also not inconsistent (weakly constrained). |
| **Second reduction:** the same dip in MIT's QLP light curves (same images, different processing)? | Yes: 225 ± 35 ppm vs NASA's 208 ± 29. QLP's individual sectors scatter more than their error bars (p = 0.004). | Weakly: 238 ± 71 ppm vs 263 ± 34. |
| **Hide-a-season check (retrospective):** search the other seasons, predict the hidden season's transit times, look only there | 3 of 4 seasons support it (3.0–5.2σ at the predicted time). **2022 fails** (40 ± 42 ppm). | 3 of 4 support it (2.8–4.6σ; in one season, 0.7% of the orbit shows a deeper dip). Hiding 2019 leaves too little data to find the period at all. These are retrospective: the candidates were selected using all seasons, so only future TESS data are a truly blind test. |
| **Which star is dimming?** Difference images of the raw pixels | Dimming detected near the target in 8/8 sectors, but too weakly to locate. 471 Gaia stars lie within 150"; three faint ones within 12" would need only 3–5% eclipses. | Positive target-position amplitude in 5/5 sectors, location unconstrained; the pixel statistic is not well calibrated (2 of 40 off-transit control phases score higher). An apparent centroid shift seen in detector coordinates disappears once spacecraft rotation is accounted for: no significant sky shift. |
| **Instrument glitches?** Momentum dumps, data gaps, shared periods in other stars | Only 3 of 51 transits are near a dump or gap. Stacking 88 other stars at this period gives nothing. | 95% of individual transits dip. Nothing in other stars. |
| **Period aliases?** Twice or half the period, odd/even, phase 0.5 | P strongly preferred over P/2. No secondary eclipse. Odd/even differ by 1.4σ. | No secondary eclipse. Odd/even differ by 2.0σ. |

Figures: `figures/<TIC>_fit.png` (model fit), `figures/<TIC>_holdout.png` (hide-a-season), `figures/<TIC>_products.png` (NASA vs MIT), `figures/<TIC>_pixels.png` (difference images).

### The 2022 problem (TIC 370612437)
The four seasons show 186 ± 50, 103 ± 48, 235 ± 42 and 243 ± 41 ppm. The 2022 season is low:
- Statistically this is not alarming: a consistency test gives p = 0.12.
- But the blind hide-2022 prediction found essentially nothing there.

Explanations considered:
- **A blend whose light share changes between sectors.** Ruled out as a simple contamination change: the crowding numbers in 2022 match the other years.
- **Timing variations from another planet.** Unsupported: the 2022 timing is simply unconstrained.
- **Noise.** Plausible.

It stays an open question.

### What the audits changed
We asked Codex (OpenAI) and a fresh Claude agent to break the analysis, read-only. Their full reports are in `results/phase3/audit_codex.md` and `results/phase3/audit_claude.md`. What they caught, and what we did:
- **Noise was slightly underestimated.** The control boxes used to measure noise had been cleaned slightly differently from the transit box. Fixed. Significance dropped from about 8.2σ to **7.2σ and 7.7σ**.
- **The model fit had a prior-math slip** (a missing Jacobian) and ignored correlated noise. Both fixed. Sizes barely changed.
- **The first pixel test could "locate" a brightening as easily as a dimming,** and treated each sector separately. It was rewritten as a joint fit across all sectors, with balanced before/after frames, every neighbour tested as the possible source, and TESS-band brightnesses.
- **The hide-a-season test had been adjusted after one fold failed** (a finer box width). We report both versions, and the failures stand.
- **The first upside-down test** searched for only one signal and skipped triage. It was redone with the full pipeline.
- **Confirmed correct:** the orbital physics, density units, image orientation, binning and the exposure smoothing.

**Second review round (Codex, `results/review/REVIEW_CODEX.md`; response in `results/review/RESPONSE_TO_CODEX.md`):**
- **False-alarm numbers withdrawn.** The 0.01/0.2 figures were unsupported.
- **The centroid "shift" for HIP 114417 is not real in sky coordinates.**
- **The pixel Δχ² p-values aren't calibrated.**
- **QLP errors updated** to matched processing: 225 ± 35 and 238 ± 71 ppm.
- **"Never searched by SPOC" was false.** SPOC has searched FFI light curves since 2021; neither star was flagged.
- **The predicted UTC times were off by up to 6.5 minutes.** Fixed to proper geocentric UTC.
- **Hold-out checks are retrospective,** and are now described that way.
- **Density agreement is weak,** and is described that way.

We independently re-ran the reviewer's centroid and timing checks and reproduced them exactly.

### Verdict
| | TIC 370612437 | TIC 252721481 |
|---|---|---|
| Is the dip real (not noise)? | Probably: local SNR ~7, seen with independent detrending and in QLP, no clean flipped fake as strong. No calibrated false-alarm probability. | Probably: local SNR ~7.7, but QLP support is weak (3.3σ). No calibrated false-alarm probability. |
| Is it on this star? | **Unknown.** Very crowded field, so a faint background eclipsing pair can't be excluded | **Unknown.** No significant centroid shift, but the pixels can't localise the source |
| Status | Weak candidate, needs more data | Weak candidate, needs more data |

### Statistical validation (TRICERATOPS): attempted, not completed
TRICERATOPS (Giacalone et al. 2021) is the standard tool for weighing "planet on the target" against every blend scenario. Our run never completed, for reasons unrelated to the data:
- it pins old numpy/scipy versions (now isolated in `.venv-tri`, `requirements-triceratops.txt`);
- the TRILEGAL server it calls serves an incomplete certificate chain (fixed with `cache/certs/bundle.pem`);
- MAST dropped and then stalled its large image downloads during a network change.

It isn't needed for a CTOI, and without high-resolution imaging it couldn't validate these candidates anyway. The follow-up team will run it once ground-based data exist. `pipeline/12_triceratops.py` is ready to rerun.

### Predicted transits for the next TESS visits
TESS observes TIC 370612437 again in sectors 121, 122 and 129 (from 2027-08-23), and HIP 114417 in sectors 122 and 129. `pipeline/13_predictions.py` lists every predicted transit with its uncertainty (±28 and ±20 minutes, 1σ) in `results/phase3/*_predicted_transits.csv`. The dates of sectors 122 and 129 are estimates until the official schedule is published. Publishing these before the data arrive turns the next sectors into a clean, pre-announced test.

### What would settle it
1. **Ground-based telescopes with sharper vision** (pixels about 1" rather than 21"). They can't see a 0.02% dip, but they can easily check whether any nearby faint star has a deep eclipse at the predicted times. This is the standard "neighbour clearing" step that TESS follow-up teams (TFOP) do for every candidate.
2. **Statistical validation (TRICERATOPS)**, which weighs the planet scenario against every blend scenario using the real neighbours.
3. **Future TESS sectors.** More transits at the predicted times would test the period and the 2022 deficit. Requesting 2-minute data for both stars would sharpen the pixel test.
4. **Submit as CTOIs**, with every caveat above stated plainly. That puts them in front of the professionals who can do steps 1 and 2.

---

## Part 5: from candidate to planet: reporting, confirmation and naming

### What "confirmed" means
No one can confirm a planet from TESS data alone, including NASA. A TESS dip proves only that *something* periodically blocks light somewhere in those 21-arcsecond pixels. Turning a candidate into a planet takes other telescopes:

1. **Clear the neighbours (ground-based photometry).** Telescopes with ~1" pixels watch the field at the predicted transit times. They can't see a 0.02% dip, but any nearby faint star with a deep eclipse at those times shows up easily. If every neighbour is clean, the dip must come from the target or something extremely close to it.
2. **High-resolution imaging** (adaptive optics or speckle cameras on 4–10 m telescopes) rules out companion stars hiding within ~0.1–1".
3. **Reconnaissance spectra** check the star is single, quiet and as the catalog says (not a disguised binary).
4. Then either:
   - **Statistical validation**, as most small TESS planets get. A program such as TRICERATOPS combines all of the above with the light curve to show a blend is improbable, typically a false-positive probability below 1.5%.
   - **A mass measurement** with precise radial velocities. A ~1 Earth-radius planet around an 11th-magnitude K dwarf would wobble its star by under 1 m/s, which is at the limit of the best instruments today. That is why validation is the usual route.

### Who to report to
**Rule change (checked on ExoFOP on 2026-10-09):** since **19 August 2026**, ExoFOP only accepts community candidates that are **already published**. The options are:
- a peer-reviewed journal with online access; or
- a **Research Note of the AAS (RNAAS)**, provided the methods used to detect and vet the candidate were themselves published in a peer-reviewed journal and are cited.

Uploading also needs prior approval, through the **Published Candidate Upload Request** form (https://exofop.ipac.caltech.edu/tess/candidate_help.php). An earlier version of this guide said to submit CTOIs directly; that is no longer possible.

The path is therefore:
1. **Preregister now (Zenodo, free).** Upload the predicted transit times and a short description. This gives a public, dated record and a DOI immediately, before any journal step. Social posts can link to it.
2. **Publish a short note.**
   - The quickest route is an **RNAAS** (about 1,000 words, one figure). It describes the two signals and cites the published methods we used: BLS, the wotan detrending, batman transit models, standard TESS vetting tests.
   - Whether ExoFOP accepts our pipeline as "previously published methodology" is their call. If not, the route is a full paper, or a methods paper first.
   - Working with a professional astronomer, for example through TFOP, would strengthen either route.
3. **Request ExoFOP upload access,** then submit each signal as a **CTOI**, citing the paper: **TIC 370612437.01** and **TIC 252721481.01**. The values and caveats are ready in `docs/CTOI_DRAFT.md`.
4. **The TESS TOI team at MIT** reviews CTOIs. If one meets their standards it becomes an official **TOI**, which puts it on the target lists of the **TESS Follow-up Observing Program (TFOP)**: the professional and amateur observatories that do the follow-up steps above.

Two more things you can do:
- **TESS Director's Discretionary Time.** The Reddit poster got 2-minute observations scheduled for a future sector. That sharpens both the dip and the pixel test, and doesn't depend on ExoFOP.
- **Credit.** The people who identify a CTOI are credited, and are normally invited onto the first discovery paper if it is followed up. A published note plus a Zenodo record establishes that priority.

### Naming
- **The scientific name is automatic:** the star's name plus a letter, starting at "b". If confirmed these would be **TOI-XXXX b** (once a TOI number is assigned), or **TIC 370612437 b** and **HIP 114417 b**.
- **Official proper names come only from the IAU** (International Astronomical Union), through occasional public **NameExoWorlds** campaigns (2015, 2019, 2022). In these, teams in each country propose a name for a planet and its host star, chosen from a list of already-confirmed planets. Discoverers cannot name planets themselves. Under the 2022 rules:
  - names must have long-standing cultural, historical or geographical significance;
  - the star and planet names share a theme;
  - names of people (living or dead), pets and organisations are excluded.
- **"Aakanksha"** (Sanskrit/Hindi for "aspiration" or "longing") for TIC 370612437 and **"Orbit"** for HIP 114417 can be used now as informal nicknames, as the videos do. Many famous objects have informal nicknames. To become official, the planet would first need to be confirmed and picked for a future campaign. A team would then propose it as a word with cultural meaning, not as a person's name, paired with a host-star name on the same theme (for example, another Sanskrit word about striving or the sky).

Sources: IAU NameExoWorlds rules (https://nameexoworlds.iau.org/2022naming-rules), TESS Follow-up Program (https://heasarc.gsfc.nasa.gov/docs/tess/tfop.html), ExoFOP-TESS submission guide (https://tess.mit.edu/followup/exofop-tess).

### Media
`media/370612437_animation.mp4` (and `media/252721481_animation.mp4`) is a 14-second explainer showing:
- the star, and the planet crossing it on the fitted tilted path;
- the real stacked TESS dimming data with the transit model;
- where the star is in the sky (Cygnus);
- its distance (158 light-years).

Make one for any candidate with `python pipeline/11_animation.py <TIC> "<nickname>"`.

`media/candidates.mp4` is a 50-second combined video (`pipeline/15_combined_video.py`), built to hook in the first seconds:
1. both candidates up front, with live transits and info cards;
2. how we found them: the search funnel, the "is it new?" checks and signal strength;
3. where candidate #1 is in the sky;
4. a fast 3D flight through 32,434 real stars at their Gaia distances;
5. the transit with the real stacked data;
6. what happens next.

The music is an original score synthesised in code (`pipeline/14_score.py`): ambient pads, a swell during the flight, and a low pulse as each planet crosses.
