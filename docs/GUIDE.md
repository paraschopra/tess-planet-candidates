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
- **SPOC** (Science Processing Operations Center) is NASA's main pipeline. It searches the 2-minute stars carefully and also produces light curves from full-frame images for many brighter stars (called **TESS-SPOC** light curves).
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

### The two strongest candidates (not yet vetted; Phase 3 tests have not been run)

| | TIC 370612437 | TIC 252721481 |
|---|---|---|
| Period | 3.5649 d | 4.9336 d |
| Depth | ~204 ppm (0.02%) | ~252 ppm (0.025%) |
| Duration | ~1.1 h | ~1.5 h |
| Implied size | ~1.1 Earth radii | ~1.1 Earth radii |
| Star | K dwarf, 4380 K, 0.66 R☉, Tmag 9.8 | K dwarf, 4301 K, 0.61 R☉, Tmag 9.7 |
| Distance | 48 pc (~158 ly) | 41 pc (~134 ly) |
| Data | 8 sectors, 2019–2025, 55 transits | 5 sectors, 2019–2025, 21 transits |
| SDE / SNR | 11.3 / 8.4 | 10.3 / 8.1 |

Both show:
- a single clean periodogram peak;
- a box-shaped folded dip;
- consistent depth in *every* sector across six years;
- matching odd/even dips;
- no secondary eclipse.

The two stars are 12° apart on the sky, so they cannot contaminate each other. They share sectors only because TESS observed that region at the same times.

Caveats to resolve in Phase 3:
- 370612437's dip is about half as long as a central transit would be. That suggests a grazing orbit, or wrong stellar parameters.
- 370612437 has 12% light contamination from neighbours.
- Neither has had pixel-level, catalog or literature checks yet.

Figures: `figures/370612437_s1.png`, `figures/252721481_s1.png`.
