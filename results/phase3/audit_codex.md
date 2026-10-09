# Codex (OpenAI) adversarial audit of Phase 3, 2026-10-09

Run read-only with `codex exec --sandbox read-only` while some Phase 3 jobs were still running (prompt in the commit history / GUIDE Part 4). Verbatim final report:

**I would not yet report these as well-vetted planet candidates.** Both contain persistent transit-like signals, but the current analysis overstates how strongly false alarms and contaminating sources have been excluded. None of the findings below proves either signal is false.

1. **Critical — SDE is not calibrated for the search and selection procedure.**  
   [lcsearch.py:130](/Users/paraschopra/Documents/Code/random/planet/pipeline/lcsearch.py:130) standardizes the periodogram, then [lcsearch.py:188](/Users/paraschopra/Documents/Code/random/planet/pipeline/lcsearch.py:188) selects its maximum after optimizing duration and epoch. This is not a Gaussian significance or a survey-wide false-alarm probability.

   The inverted results already demonstrate the problem: in a snapshot of **405 completed inverted searches**, **95** had SDE ≥10, SNR ≥6 and ≥3 transits; **65** exceeded TIC 370612437’s SDE of 11.26 with those other cuts. These are *before triage*, so they are not candidate false-discovery rates, but they contradict the inference that SDE >10 itself makes these exceptional.

   [06_inverted.py:39](/Users/paraschopra/Documents/Code/random/planet/pipeline/06_inverted.py:39) searches only one signal and never applies the subsequent vetting/selection procedure. Also, inversion exchanges flares and dropouts under asymmetric clipping; it is not an exchangeable null for sign-asymmetric artifacts. The full search **and selection** need calibrated controls. Simply multiplying periods × durations × iterations is also inappropriate because those trials are correlated.

2. **Moderate — The empirical error compares differently processed transit and control phases; I measured an optimistic bias.**  
   [vet.py:37](/Users/paraschopra/Documents/Code/random/planet/pipeline/vet.py:37) excludes candidate transits from detrending, whereas [vet.py:85](/Users/paraschopra/Documents/Code/random/planet/pipeline/vet.py:85) measures control boxes where the data participated in their own trend estimation. Candidate-phase negative excursions are protected; control-phase excursions can be attenuated. Additionally, [vet.py:56](/Users/paraschopra/Documents/Code/random/planet/pipeline/vet.py:56) exempts candidate transits from dropout clipping.

   I repeated the 40 control measurements, independently detrending with **each control phase masked**:

   | TIC | Reported empirical error | Matched-processing error | SNR change |
   |---|---:|---:|---:|
   | 370612437 | 25.29 ppm | 28.82 ppm | 8.24 → **7.23** |
   | 252721481 | 31.81 ppm | 33.92 ppm | 8.26 → **7.75** |

   Those revised SNRs remain local measurements, not discovery significances. Forty fixed phase offsets—not random independent realizations—cannot establish an 8σ tail probability. MAD also deliberately suppresses rare excursions, precisely the tail relevant to an artifact search. Finally, `box_depth` measures relative to exactly 1, while the MAD removes the control distribution’s center without correcting the measured depth for that center.

3. **Critical for claims of independent validation — The hold-out test is post-selection and adaptively revised.**  
   The actual prediction uses training-only BLS results at [08_holdout.py:73](/Users/paraschopra/Documents/Code/random/planet/pipeline/08_holdout.py:73); I did **not** find the full-data ephemeris secretly substituted into the prediction. However, these two targets and their significance were selected using every season. Reusing those seasons afterward does not create untouched validation data.

   The width-refinement procedure was explicitly added after inspecting a failed fold: [08_holdout.py:10](/Users/paraschopra/Documents/Code/random/planet/pipeline/08_holdout.py:10). That adaptation must be included in any false-alarm calibration. Preprocessing also precedes splitting, including sector rejection based on the whole star’s sector-noise distribution, although that is a smaller issue.

   More importantly, the saved results **do not pass every fold**:

   - TIC 370612437, held-out 2022: **26.5 ±42.8 ppm**, or **39.7 ±41.8 ppm** after refinement—both below 1σ. [Result:49](/Users/paraschopra/Documents/Code/random/planet/results/phase3/370612437_holdout.json:49)
   - TIC 252721481, withholding 2019: the blind search selects **0.73090 d**, with SDE 6.12, rather than 4.93355 d. [Result:9](/Users/paraschopra/Documents/Code/random/planet/results/phase3/252721481_holdout.json:9)

   These failures are not decisive falsifications—sensitivity and ephemeris uncertainty matter—but they cannot be presented as successful independent predictions. Likewise, `frac_orbit_as_deep = 0` at [08_holdout.py:57](/Users/paraschopra/Documents/Code/random/planet/pipeline/08_holdout.py:57) means no sampled comparison exceeded the prediction, not zero chance probability.

4. **Critical — Difference imaging can return a convincing-looking position without detecting missing light.**  
   [10_pixels.py:90](/Users/paraschopra/Documents/Code/random/planet/pipeline/10_pixels.py:90) allows **negative Gaussian amplitude**. I reproduced this with a synthetic negative spot: the fitter returned amplitude −1 and its correct centroid. A brightening can therefore receive the same localization machinery as a dimming.

   The centroid is constrained to ±3 pixels around the target, and only pixels within 4.5 pixels enter the fit. All fitted positions subsequently enter the combined offset, without requiring positive amplitude, a significant source detection, acceptable fit residuals, or an interior solution: [10_pixels.py:139](/Users/paraschopra/Documents/Code/random/planet/pipeline/10_pixels.py:139), [10_pixels.py:189](/Users/paraschopra/Documents/Code/random/planet/pipeline/10_pixels.py:189).

   Consequently, “offset consistent with target” can mean **an unconstrained or misidentified feature**, not a localized transit. The bootstrap repeats the same constrained model and cannot establish that its centroid distribution represents a detected source. No candidate `*_pixels.json` results were present when inspected.

5. **Moderate — Pixel subtraction and centroid uncertainties omit important failure modes.**  
   [10_pixels.py:112](/Users/paraschopra/Documents/Code/random/planet/pipeline/10_pixels.py:112) requires enough out-of-transit frames but does not require frames on **both sides** or balanced observation times. Near gaps, a one-sided baseline does not cancel a linear drift.

   The frames receive scalar sky subtraction but no registration or explicit treatment of pointing-induced difference patterns. Localization uses a circular Gaussian, diagonal pixel weights and a bootstrap of only a handful of transits per sector. PRF mismatch, correlated pixels, WCS errors and persistent artifacts are not propagated into the combined covariance. A small bootstrap ellipse can therefore exaggerate the exclusion of neighboring stars.

6. **Moderate — Neighbor exclusions use inconsistent astrometry and inadequate flux ratios.**  
   Individual images propagate proper motion at [10_pixels.py:125](/Users/paraschopra/Documents/Code/random/planet/pipeline/10_pixels.py:125), but neighbor distances from the combined spot use **unpropagated J2016 positions** at [10_pixels.py:198](/Users/paraschopra/Documents/Code/random/planet/pipeline/10_pixels.py:198). Differential proper motion is omitted from the final exclusion calculation.

   The neighbor cutoff and required eclipse depths use **Gaia G**, rather than TESS-band flux and aperture throughput: [10_pixels.py:47](/Users/paraschopra/Documents/Code/random/planet/pipeline/10_pixels.py:47), [10_pixels.py:206](/Users/paraschopra/Documents/Code/random/planet/pipeline/10_pixels.py:206). A red neighbor can be substantially more important in TESS than its G contrast suggests. The calculation also combines a corrected PDCSAP depth with a simplified two-star dilution formula. These numbers cannot securely eliminate blended eclipsing binaries.

7. **Moderate — MCMC uncertainties and Δχ² ignore the noise model used to justify empirical vetting.**  
   [07_fit.py:86](/Users/paraschopra/Documents/Code/random/planet/pipeline/07_fit.py:86) uses independent Gaussian errors from per-sector point scatter. It includes no correlated-noise likelihood, baseline nuisance parameters or detrending uncertainty. The empirical depth error calculated at line 128 never enters the likelihood.

   Thus Δχ² ≈89/78 and the posterior intervals are conditional on substantially stronger assumptions than the empirical-noise analysis. The fit also uses only points selected around the discovery ephemeris, initializes walkers locally, and gives period/epoch no explicit bounds: [07_fit.py:74](/Users/paraschopra/Documents/Code/random/planet/pipeline/07_fit.py:74), [07_fit.py:126](/Users/paraschopra/Documents/Code/random/planet/pipeline/07_fit.py:126). It does not establish global period uniqueness or compare the planet interpretation against realistic artifact models.

8. **Moderate — The “Gaussian density prior” is not the prior actually sampled.**  
   The sampler parameter is `log_rho`, but [07_fit.py:82](/Users/paraschopra/Documents/Code/random/planet/pipeline/07_fit.py:82) adds only the Gaussian density penalty. In physical density, this gives
   \[
   p(\rho)\propto \exp[-(\rho-\mu)^2/(2\sigma^2)]/\rho,
   \]
   within the bounds. A Gaussian prior in physical density requires the corresponding `+log_rho` Jacobian.

   Mass/radius fractional uncertainties are also hard-coded for both stars, and radius is later sampled independently of density: [07_fit.py:132](/Users/paraschopra/Documents/Code/random/planet/pipeline/07_fit.py:132), [07_fit.py:141](/Users/paraschopra/Documents/Code/random/planet/pipeline/07_fit.py:141). Stellar-parameter covariance is lost. Agreement with the imposed stellar-density prior is not independent evidence that the target is the eclipsed star.

9. **Moderate — Passing the sector, odd/even and secondary tests provides weak exclusion power.**  
   [05_triage.py:202](/Users/paraschopra/Documents/Code/random/planet/pipeline/05_triage.py:202) flags odd/even differences only above `6 × sigma_emp`: approximately **163 ppm** and **186 ppm** for these candidates. Those are large fractions of their primary depths. The assumed factor of two is not calibrated using the actual parity sampling and covariance.

   The secondary test examines only **phase 0.5 at the primary duration**; eccentric binaries or different-width secondaries can escape. Its nominal thresholds are about **82/93 ppm**, so substantial weaker secondaries remain allowed.

   Sector inconsistency requires p<10⁻⁴, and missing diagnostics silently produce no flag: [05_triage.py:117](/Users/paraschopra/Documents/Code/random/planet/pipeline/05_triage.py:117), [05_triage.py:208](/Users/paraschopra/Documents/Code/random/planet/pipeline/05_triage.py:208). “Not rejected” does not establish equal depths. Calibration on injected planets checks rejection of planets, not acceptance of binaries and artifacts.

10. **Moderate — QLP agreement is useful, but is not an independent noise realization.**  
    Both products use the same photons and undergo the same candidate-conditioned detrending: [09_products.py:78](/Users/paraschopra/Documents/Code/random/planet/pipeline/09_products.py:78). Shared photon fluctuations, pointing artifacts, astrophysical blends and artifacts from the common detrending can appear in both. Their significances must not be combined as independent detections.

    TIC 370612437’s QLP sector depths also deserve attention: S56 is **348±41 ppm**, S57 **73±72 ppm**, and S76 **26±90 ppm**. [Results:76](/Users/paraschopra/Documents/Code/random/planet/results/phase3/370612437_products.json:76) Using the reported errors, I obtain χ²=20.82 for seven degrees of freedom, **p≈0.004**. This does not prove a false signal, but it weakens “consistent in every sector.”

11. **Moderate pipeline risk; not demonstrated as the cause here — Edge artifacts and downward outliers can survive cleaning.**  
    [lcsearch.py:65](/Users/paraschopra/Documents/Code/random/planet/pipeline/lcsearch.py:65) retains detrending edges and clips only after estimating the trend. One-sided windows around ramps/gaps can create depressions; negative excursions survive to 10σ while positive ones are removed at 3σ.

    The edge diagnostic measures the **fraction of points**, not their contribution to the detection. The event-dominance denominator sums only positive event contributions, potentially understating dominance when other events contribute negatively: [05_triage.py:79](/Users/paraschopra/Documents/Code/random/planet/pipeline/05_triage.py:79).

    However, my checks argue against a simple edge/window explanation for these two: windows of 0.5–1.5 days retain both dips, and excluding points within 0.5 days of major gaps leaves depths around **199/244 ppm**, with local empirical SNRs **7.45/7.20**.

12. **Minor — TESScut timing and result provenance need correction.**  
    [10_pixels.py:97](/Users/paraschopra/Documents/Code/random/planet/pipeline/10_pixels.py:97) directly uses TESScut times. Their barycentric correction can refer to the FFI reference position rather than the target, as documented by [Lightkurve](https://lightkurve.github.io/lkspacecraft/tutorials/calculating-barycentric-time-correction.html). Matching cached spacecraft exposures, I measured TESScut–SPOC offsets of roughly **−38 to +34 seconds**. This is not an explanation for an hour-long periodic dip, but it changes which frames enter the difference image.

    Separately, [03_search.py:35](/Users/paraschopra/Documents/Code/random/planet/pipeline/03_search.py:35) reuses caches without checking data or configuration changes. The returned periodogram is from the **last search iteration**, although plotting labels it signal 1. Both candidate periodogram caches currently contain only NaNs; the plotting fallback recomputes them. This is an auditability problem rather than evidence against the dips.

Several suspected basic bugs **did not materialize**: the \(a/R_\star\) formula, solar-to-cgs density conversion, circular-orbit inclination relation, pixel x/y conventions and out-minus-in difference sign are correct. Binning does not combine sectors or bridge long gaps. Exposure integration uses the [batman API correctly](https://lkreidberg.github.io/batman/docs/html/tutorial.html); increasing supersampling changed the best-fit model by only about 5–12 ppm at worst. Fixed limb darkening remains an uncertainty assumption, not an obvious signal-manufacturing error.

I changed no files. **Overall verdict: retain both as promising preliminary transit-like signals, but the present analysis is not sufficient to report them as well-vetted Earth-sized planet candidates.** The principal blockers are a false-alarm calibration covering the actual selection procedure, honest treatment of the failed/adaptive hold-out tests, and a significant, properly calibrated localization of the dimming source.
