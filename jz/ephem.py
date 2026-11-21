# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""Published cycle times of eRASSt J0456-20 and the forecast frozen from them.

Source: Liu et al. 2024, A&A 683, L13 (arXiv:2401.14091), Table 1 and Appendix B.
Rise-phase start times t_i (fit to all five cycles, MJD): 59065, 59363, 59591, 59778, 59972.
Drop-phase start windows (MJD): C1 59285-59291, C2 59512-59518, C3 59699-59706,
C4 59878-59883; the C5 drop was not observed (Sun constraint after MJD 60072).
Rise-phase recurrence times: 299, 227, 188 +- 3, 193 +- 3 d. The authors read the last
three cycles as "an indication of constant values of ~190 days".

Forecast used here (frozen before any post-2023 photon was counted):
  rise start of cycle n:  T_n = 59972 + 190.5 (n - 5)     (190.5 = mean of 188 and 193)
  drop of cycle n: mean of  T_n + 107.0  (mean rise-to-drop delay of cycles 3 and 4)
                   and      59880.5 + 190.5 (n - 4)  (cycle 4 drop mid-point, same period)
"""
import math

RISE_START = {1: 59065.0, 2: 59363.0, 3: 59591.0, 4: 59778.0, 5: 59972.0}
DROP_WINDOW = {1: (59285.0, 59291.0), 2: (59512.0, 59518.0), 3: (59699.0, 59706.0),
               4: (59878.0, 59883.0)}
T5 = RISE_START[5]
RECUR = (188.0 + 193.0) / 2.0
DROP_OFFSET = sum(0.5 * (DROP_WINDOW[n][0] + DROP_WINDOW[n][1]) - RISE_START[n]
                  for n in (3, 4)) / 2.0
C4_DROP_MID = 0.5 * (DROP_WINDOW[4][0] + DROP_WINDOW[4][1])

# Phase windows, days after the predicted rise start of a cycle.
ON_LO, ON_HI = 35.0, 95.0      # late rise + plateau (the "high state")
OFF_LO, OFF_HI = 120.0, 185.0  # faint phase, before the next rise


def rise_start(n):
    return T5 + RECUR * (n - 5)


def predicted_drop(n):
    rise_based = rise_start(n) + DROP_OFFSET
    drop_based = C4_DROP_MID + RECUR * (n - 4)
    return 0.5 * (rise_based + drop_based)


def on_window(n):
    t = rise_start(n)
    return (t + ON_LO, t + ON_HI)


def off_window(n):
    t = rise_start(n)
    return (t + OFF_LO, t + OFF_HI)


def cycle_of(mjd):
    return int(math.floor((mjd - T5) / RECUR)) + 5


def phase(mjd):
    n = cycle_of(mjd)
    return n, mjd - rise_start(n)


def phase_class(mjd):
    """Where the constant-recurrence forecast puts a date: RISE, ON, DROP_ZONE, OFF or PRE_RISE."""
    _, p = phase(mjd)
    if p < ON_LO:
        return "RISE"
    if p <= ON_HI:
        return "ON"
    if p < OFF_LO:
        return "DROP_ZONE"
    if p <= OFF_HI:
        return "OFF"
    return "PRE_RISE"
