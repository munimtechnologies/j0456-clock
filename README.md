# j0456-clock

**Did the star in eRASSt J0456-20 keep flaring after May 2023?** A pre-registered test on
public Swift data. Licence: Apache-2.0 (code); CC BY 4.0 (text and result files). By MV2 of
Munim, Inc. (mv2@munimtech.com). The code and text were made with AI agents working in his name,
and no human expert has checked them yet.

## Result (graded 2026-11-24): yes, a sixth flare, on schedule

All five validation gates on the published cycles passed. In the window forecast from the
published cycles before any later photon was counted (MJD 60197.5-60257.5, 2023 September to
early November), Swift-XRT detects the source at **2.64 +- 0.67 x 10^-3 counts/s (7.5 sigma)**,
against 0.17 +- 0.20 x 10^-3 in the faint windows on either side (z = 3.55): outcome
**FLARE_ON_SCHEDULE**.

- **Fading (T2): CONSISTENT.** The sixth flare is 0.495 times the fifth; the trend fixed in
  advance from cycles 3 and 5 forecast 0.559.
- **Timing (T3): ON_TIME.** The flux dropped between 2023 Nov 3 and Nov 10 (MJD 60251.1-60258.5),
  10.7 days before the forecast drop and inside its +-20-day window. The mean drop-to-drop
  recurrence over cycles 5 and 6 is **187.1 +- 3.1 days**.
- **Later visits (T4, descriptive):** 2024-2026 visits are too shallow to say whether the flares
  continue at the faded level.

Details, tables and post-hoc notes: [GRADING.md](GRADING.md). Raw outputs:
[results/v0.1.0/](results/v0.1.0/). Scoop checks: [notes/](notes/).

## Background

eRASSt J045650.3-203750 is a galaxy nucleus that flared in X-rays five times between 2020 and
2023, each flare weaker than the last and the gaps shrinking from ~300 to ~190 days: the
signature of a star on a tight orbit around a supermassive black hole, losing part of
itself at each pass (Liu et al. 2023, A&A 669, A75; Liu et al. 2024, A&A 683, L13,
[arXiv:2401.14091](https://arxiv.org/abs/2401.14091)). That paper's data end on MJD 60072
(2023 May 8; the pre-registration's "May 2" label was a calendar slip). NASA's Swift kept
watching: a weekly campaign from June 2023 to March 2024 and visits up to October 2026, all
public and, as far as we can find, never analysed in print.

The repo froze, before counting any of those photons, a forecast drawn from the published
numbers (a constant ~190.5-day cycle; the sixth high state at MJD 60197.5-60257.5 and its drop
near MJD 60265.5; the flares fading by the same factor per cycle) and the rules that grade it.
See [PREREGISTRATION.md](PREREGISTRATION.md).

## Layout

- `jz/ephem.py`: published cycle times and the frozen forecast.
- `jz/photometry.py`, `jz/stats.py`, `jz/geometry.py`: Swift-XRT 0.3-2 keV aperture photometry
  and Poisson statistics (pure Python; shared with
  [at2018fyk-return](https://github.com/munimtechnologies/at2018fyk-return)).
- `jz/changepoint.py`: one-step Poisson change-point for sparse light curves.
- `jz/grade.py`: validation gates on the published cycles and the pre-registered outcome.
- `jz/uvot.py`: secondary UVOT photometry.
- `job.py`: downloads the data from HEASARC and writes `out/results.json`, `out/xrt_obs.csv`,
  `out/summary.md`.
- `test_jz.py`: unit tests (standard library only): `python -m unittest -v test_jz`.

Run: `pip install -r requirements.txt && python job.py --data-dir data --out-dir out`.

## Credit

The setup is entirely Liu, Ryu, Goodwin, Rau, Homan, Krumpe, Merloni, Grotova, Anderson,
Malyali & Miller-Jones (2024) and Liu et al. (2023); the Swift data are from NASA's HEASARC.
Theory context: Lau, McKeever & Yu 2026 (diffusive tides, arXiv:2603.27132).
