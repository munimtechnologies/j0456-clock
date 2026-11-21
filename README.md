# j0456-clock

**Did the star in eRASSt J0456-20 keep flaring after May 2023?** A pre-registered test on
public Swift data. Licence: Apache-2.0. By MV2 of Munim, Inc. (mv2@munimtech.com); the code
and text were made with AI agents working in his name, and no human expert has checked them yet.

eRASSt J045650.3-203750 is a galaxy nucleus that flared in X-rays five times between 2020 and
2023, each flare weaker than the last and the gaps shrinking from ~300 to ~190 days: the
signature of a star on a tight orbit around a supermassive black hole, losing part of
itself at each pass (Liu et al. 2023, A&A 669, A75; Liu et al. 2024, A&A 683, L13,
[arXiv:2401.14091](https://arxiv.org/abs/2401.14091)). That paper's data end on 2023 May 2.
NASA's Swift kept watching: a weekly campaign from June 2023 to March 2024 and visits up to
October 2026, all public and, as far as we can find, never analysed in print.

This repo freezes, before counting any of those photons, a forecast drawn from the published
numbers (a constant ~190.5-day cycle; the sixth high state at MJD 60197.5-60257.5 and its drop
near MJD 60265.5; the flares fading by the same factor per cycle) and the rules that grade it.
See [PREREGISTRATION.md](PREREGISTRATION.md).

The wow, either way: a sixth flare on time would be the first forward check of this star's
orbital clock; no flare where the fading trend still predicts one would mean the flares
stopped, which is what a star being destroyed would look like.

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
