# Pre-registration: did eRASSt J0456-20 keep flaring after May 2023?

Frozen 2026-11-21, code version 0.1.0, before any photon from an observation after
MJD 60072 was counted by us. What we had looked at: the HEASARC observation log only
(observation IDs, start times, exposures), which is needed to choose the windows.
Made with AI agents working for MV2 of Munim, Inc. Licence: Apache-2.0.

## Why this question

eRASSt J045650.3-203750 (J0456-20, z = 0.077) is one of the best candidates for a
*repeating partial tidal disruption*: a star on a tight, eccentric orbit around a
supermassive black hole that loses part of its mass at every close passage. Liu et al.
(2023, A&A 669, A75; 2024, A&A 683, L13) followed five X-ray/UV cycles up to MJD 60072
(2023 May 2) and found that:

- the recurrence time fell fast, from ~299 d to 227 d to 188 d, then 193 d, which they read
  as "an indication of constant values of ~190 days in the latest cycles";
- each flare was weaker than the one before: plateau luminosities 4.6, 1.1, 0.70, 0.27 and
  0.13 x 10^44 erg/s for cycles 1-5;
- hydrodynamic models need a star that has already lost most of its mass.

The paper made no dated forecast, and the source has not been reported on since. NASA's
public Swift archive holds ~59 observations taken after MJD 60072: a weekly campaign from
MJD 60118 to 60384 (2023 June to 2024 March) and sparser visits up to MJD 61330 (2026 October).
Those data can say whether the star went on flaring on the 190-day clock, faded on the same
trend, or stopped (a star destroyed, or an orbit changed). Nobody has published that.

## Forecast (frozen from the published numbers; jz/ephem.py)

- Rise start of cycle n: T_n = 59972 + 190.5 (n - 5) MJD (190.5 = mean of 188 and 193 d).
  Cycle 6: T_6 = 60162.5.
- High state ("ON"): T_n + 35 to T_n + 95 d. Faint state ("OFF"): T_n + 120 to T_n + 185 d.
  Cycle 6 ON = MJD 60197.5-60257.5; cycle 5 OFF = 60092-60157; cycle 6 OFF = 60282.5-60347.5.
- Cycle-6 drop: mean of T_6 + 107.0 (mean rise-to-drop delay of cycles 3 and 4) and
  59880.5 + 2 x 190.5 (cycle 4 drop mid-point plus two periods) = **MJD 60265.5**.
- Fading: per-cycle factor f = sqrt(R5 / R3) from our own cycle-3 and cycle-5 plateau rates;
  the forecast cycle-6 high-state rate is R5 x f.

## Measurement (jz/photometry.py, job.py)

Swift-XRT photon-counting cleaned event files from the HEASARC archive; 0.3-2 keV
(PI 30-200); source circle 35" around the host nucleus (RA 74.20750, Dec -20.63000);
background annulus 70"-200"; rates corrected for the PSF fraction (0.864); observations
more than 8' off-axis are dropped. Detection: Poisson tail of the source counts given the
background >= 3 sigma and >= 4 counts. Upper limits: 99.73% Bayesian (Kraft, Burrows &
Nousek 1991). Stacks sum counts and exposures.

## Validation gates (on the cycles Liu et al. 2024 already published)

All must pass or the result is reported as UNGRADED_GATE_FAILURE (with what it would
have been) and no rule is changed without a version bump and a Deviations entry here.

- G1: the cycle-3 plateau stack (MJD 59642-59699) is detected.
- G2: the cycle-4 (59817-59877) and cycle-5 (60008-60075) plateau stacks are detected.
- G3: the faint phase after cycle 4 (59886-59968) is below half the cycle-4 plateau rate
  (its rate if detected, else its 3-sigma upper limit).
- G4: the change-point method (below) run on MJD 59815-59968 finds a drop
  (rate before > rate after, TS >= 9) whose bracket overlaps MJD 59870-59892
  (published drop 59878-59883).
- G5: the cycle-3 plateau rate is above the cycle-5 plateau rate (the published fading).

## Primary outcome (MJD 60100-60400, the weekly campaign)

Let ON = stack of cycle-6 ON, OFF = stack of cycle-5 OFF plus cycle-6 OFF,
z = (rate_ON - rate_OFF) / combined error. In order:

1. **FLARE_ON_SCHEDULE**: ON detected and z >= 3.
2. **FLARE_OFF_SCHEDULE**: otherwise, some stack of 4 consecutive observations in
   60100-60400 is brighter than all the other observations in that range at z >= 4
   (stricter because of the scan over ~35 windows).
3. **STEADY_SOURCE**: otherwise, the stack of all of 60100-60400 is detected.
4. **NO_FLARE**: nothing detected.
(**NO_DATA** if no observation in 60100-60400 could be measured.)

For STEADY_SOURCE and NO_FLARE, the depth reading is **BELOW_FORECAST_LEVEL** if the
3-sigma upper limit of the ON stack is below the forecast cycle-6 rate R5 x f, else
**NOT_DEEP_ENOUGH**. Only BELOW_FORECAST_LEVEL supports "the flares stopped or collapsed
faster than the published trend".

## Secondary tests

- **T2 amplitude** (only if FLARE_ON_SCHEDULE): d = log10(rate_ON / R5) - log10(f);
  CONSISTENT if |d| <= 2 sigma (errors of both ratios in quadrature; the shared R5 is
  ignored), else FADING_FASTER (d < 0) or FADING_SLOWER (d > 0).
- **T3 timing**: one-step Poisson change-point (jz/changepoint.py: constant rate before
  and after, change between consecutive observations, at least 2 on each side) on
  MJD 60180-60350. If TS >= 9 and the rate falls, ON_TIME if the bracket overlaps
  60265.5 +- 20 d, else EARLY or LATE; we report the implied drop-to-drop recurrence
  (mid - 59880.5) / 2. Otherwise UNTIMED.
- **T4 later epochs** (descriptive): groups of observations after MJD 60400 (gaps > 20 d
  split groups), each with its stacked rate or limit and where the constant-190.5-d
  forecast puts it (RISE / ON / DROP_ZONE / OFF / PRE_RISE). Swift's 2024-2026 visits fall
  mostly in forecast faint phases, so these are not used for the primary outcome.
- **UV** (descriptive): UVOT UVM2/UVW2/UVW1 5" aperture rates (host not subtracted),
  medians in the published cycle-4 plateau and faint windows and the cycle-6 ON and OFF
  windows.

## What each outcome would mean

- FLARE_ON_SCHEDULE + ON_TIME: the 190-day clock held for a sixth cycle, the first forward
  check of this source's ephemeris; T2 then says whether the star's mass loss keeps falling.
- FLARE_OFF_SCHEDULE: still flaring, but the period changed again (Liu et al. saw it shrink
  from 299 d to 188 d; diffusive-tide models, Lau, McKeever & Yu 2026, predict irregular
  recurrence).
- NO_FLARE with BELOW_FORECAST_LEVEL: the flares stopped or collapsed after cycle 5,
  consistent with the star having been destroyed, or its orbit having changed.

## Limits stated in advance

Swift's ~1 ks weekly visits are shallow; the cycle-5 plateau was already ~40x below
cycle 1. No spectral fitting, no flux conversion: everything is in count-rate ratios. No
exposure maps. A host-galaxy X-ray floor or a nearby source in the annulus would bias all
epochs alike. The forecast windows are our reading of the published numbers, not a
forecast the authors made.

## Deviations

None yet.
