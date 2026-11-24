# Grading: eRASSt J0456-20 after May 2023 (code 0.1.0)

Graded 2026-11-24 by AI agents working for MV2 of Munim, Inc. No human expert has checked it yet.
Licence: Apache-2.0 (code), CC BY 4.0 (text and result files).

- Rules: [PREREGISTRATION.md](PREREGISTRATION.md), frozen and committed 2026-11-21 14:49 UTC
  (commit 84a94074b08ac125a12e388388e49017f4d1d874), before any post-MJD-60072 photon was counted.
- Run: compete job `j0456-clock-v1`, 2026-11-23 14:00-14:22 UTC, image python:3.12, code exactly
  as committed (unit tests passed first). Peak memory 0.41 GB.
- Raw outputs, unchanged: [results/v0.1.0/](results/v0.1.0/) (`results.json`, `summary.md`,
  `xrt_obs.csv`, `uvot_obs.csv`, `xrt_skipped.csv`, `swiftmastr_query.txt`).
- Data: 133 Swift observations in the HEASARC cone; 131 measured with XRT, 2 skipped (no
  usable PC event file; see `xrt_skipped.csv`).

## Validation gates (on the cycles Liu et al. already published)

| Gate | Rule | Measured | Result |
|---|---|---|---|
| G1 | cycle-3 plateau stack detected | 452 net counts, 40 sigma, 17.08e-3 ct/s | PASS |
| G2 | cycle-4 and cycle-5 plateau stacks detected | 7.97e-3 ct/s (27 sigma); 5.34e-3 ct/s (15 sigma) | PASS |
| G3 | faint phase after cycle 4 below half the cycle-4 plateau | 3-sigma limit 0.59e-3 vs 3.98e-3 | PASS |
| G4 | change-point finds the cycle-4 drop overlapping MJD 59870-59892 | drop between MJD 59874.8 and 59881.7, TS 176 | PASS |
| G5 | cycle-3 plateau brighter than cycle-5 plateau | 17.08e-3 > 5.34e-3 ct/s | PASS |

All five gates passed, so the outcome below is graded.

## Primary outcome: FLARE_ON_SCHEDULE

| Stack | Visits | Exposure | Source / background counts | Rate (1e-3 ct/s, 0.3-2 keV, PSF-corrected) | Detection |
|---|---|---|---|---|---|
| Cycle-6 high state (forecast window MJD 60197.5-60257.5) | 8 | 7.57 ks | 19 / 1.71 expected | 2.64 +- 0.67 | 7.5 sigma |
| Faint states (cycle-5 and cycle-6 windows) | 15 | 17.90 ks | 9 / 6.32 expected | 0.17 +- 0.20 | not detected (3-sigma limit 0.96) |
| All of MJD 60100-60400 | 38 | 39.81 ks | 39 / 12.04 expected | 0.78 +- 0.18 | 6.1 sigma |

z (ON vs OFF) = 3.55 >= 3 and the ON stack is detected, so the first pre-registered outcome
applies: **FLARE_ON_SCHEDULE**. The sixth flare happened inside the window forecast from the
published cycles, and the source was faint outside it.

Free-phase scan (reported for completeness; not needed once outcome 1 applies): the brightest
4-visit stack in 60100-60400 is MJD 60216.8-60237.5, inside the forecast ON window (z = 2.77).

## Secondary tests

- **T2 amplitude: CONSISTENT.** Per-cycle fading factor from our own cycle-3 and cycle-5
  plateaus f = sqrt(5.342/17.080) = 0.559; forecast cycle-6 rate R5 x f = 2.99e-3 ct/s.
  Measured/R5 = 0.495; d = log10(0.495) - log10(0.559) = -0.053 +- 0.127 (0.4 sigma). The star's
  flares kept fading at the same rate per cycle.
- **T3 timing: ON_TIME.** Change-point on MJD 60180-60350 (24 visits): the rate falls from
  2.02e-3 to 0.33e-3 ct/s (TS = 10.8 >= 9) between MJD 60251.1 and 60258.5 (2023 Nov 3 and 10);
  95% range 60244.3-60265.5. The forecast drop was MJD 60265.5 +- 20; mid-bracket 60254.8 is
  10.7 days early, inside the window. Implied mean drop-to-drop recurrence over cycles 5 and 6:
  (60254.8 - 59880.5) / 2 = **187.1 +- 3.1 days** (half-width from the bracket), matching the
  ~188 and ~193 days Liu et al. (2024) measured for cycles 3-5.
- **T4 later epochs (descriptive, not graded):**

| Visits (MJD) | Dates | Visits | Exposure | Rate (1e-3 ct/s) | 3-sigma limit | Forecast phase |
|---|---|---|---|---|---|---|
| 60502.0-60563.1 | 2024 Jul-Sep | 9 | 19.2 ks | 0.03 +- 0.18 | 0.80 | cycle 7 faint, cycle 8 rise |
| 60850.9-60876.3 | 2025 Jun-Jul | 8 | 35.5 ks | 0.38 +- 0.15 (3.1 sigma, 22 counts) | 0.93 | cycle 9 drop zone, then faint |
| 60962.9 | 2025 Oct | 1 | 0.5 ks | not detected | 13.6 | cycle 10 high state |
| 61049.2 | 2026 Jan | 1 | 1.4 ks | not detected | 4.9 | cycle 10 faint |
| 61282.2-61284.3 | 2026 Aug-Sep | 2 | 1.3 ks | 0.75 +- 0.92 | 7.3 | cycle 11 faint |
| 61329.7 | 2026 Oct | 1 | 0.4 ks | not detected | 15.9 | cycle 12 rise |

  The 2025-2026 single visits are far too shallow to test whether flaring continues: at the
  measured fading (f = 0.559 per cycle) a cycle-10 high state would be near
  2.64e-3 x 0.559^4 = 0.26e-3 ct/s, 20-60 times below those limits.

- **UV (descriptive):** UVOT 5" aperture rates include the host galaxy and were not
  host-subtracted. UVW1 medians are the same in the cycle-6 high state (0.0379 ct/s, 8 visits)
  and the faint windows (0.0381 ct/s, 15 visits): no UV flare is visible over the host light in
  this crude measurement. (The cycle-6 campaign has no UVM2/UVW2 visits.)

## Post-hoc notes (looked at after grading; NOT part of the pre-registered test)

1. **Hint of the seventh rise.** The four visits the frozen forecast places in cycle 7's rise
   (MJD 60356.6-60377.2; obsids 00097198036-039) hold 7 source-region counts against
   35 x 0.0349 = 1.22 expected from background (area ratio 35"^2 / (200"^2 - 70"^2) = 0.0349);
   Poisson P(>= 7 | 1.22) = 3e-4. The window was chosen after looking, so this is a hint, not a
   detection; the campaign ended one week later (MJD 60384.4, 0 counts in 0.85 ks).
2. **Mid-2025.** Extrapolating the measured fading three more cycles gives a cycle-9 high-state
   rate near 2.64e-3 x 0.559^3 = 0.46e-3 ct/s. The 2025 Jun-Jul stack (0.38 +- 0.15 e-3 ct/s,
   3.1 sigma), which the constant-period forecast places at the end of cycle 9's high state, is
   consistent with that, but the forecast phase after four extrapolated cycles is uncertain by
   more than a week and the window was not pre-registered. It is not evidence either way.
3. **Date label.** PREREGISTRATION.md and README.md call MJD 60072 "2023 May 2"; MJD 60072 is
   2023 May 8. Only the calendar label was wrong; every rule uses MJD.

## Updated ephemeris for observers (derived after grading; an extrapolation, not a test)

Switch-off of cycle n = MJD 60254.8 + 187.1 (n - 6), with an uncertainty of roughly
3.1 (n - 6) days plus the 3.7-day half-width of the cycle-6 bracket, and more if the period
drifts. The high state ends at the switch-off. At the measured fading (f = 0.559 per cycle),
the high-state rate of cycle n is about 2.64e-3 x 0.559^(n - 6) Swift-XRT ct/s.

| Cycle | Switch-off (MJD) | Date | Rough uncertainty | High-state rate (1e-3 ct/s) |
|---|---|---|---|---|
| 12 | 61377.4 | 2026 Dec 3 | +- 22 d | 0.08 |
| 13 | 61564.5 | 2027 Jun 8 | +- 25 d | 0.05 |
| 14 | 61751.6 | 2027 Dec 12 | +- 29 d | 0.03 |

These rates are far below what short Swift visits reach; only much deeper X-ray exposures could
test them. The Swift visit of 2026 Oct 17 (MJD 61329.7, 0.4 ks) falls in the cycle-12 high
state by this ephemeris; its 3-sigma limit (15.9e-3) is ~200 times above the expected rate.

## Deviations

None. The job ran the committed code unchanged.

## Credit

The setup (the five published cycles, the shrinking recurrence and the fading) is Liu, Ryu,
Goodwin, Rau, Homan, Krumpe, Merloni, Grotova, Anderson, Malyali & Miller-Jones (2024, A&A 683,
L13) and Liu et al. (2023, A&A 669, A75). The post-2023 Swift observations were requested by
others (targets 00097198 and 00098103) and are public through NASA's HEASARC.
