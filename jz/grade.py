# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""Validation gates (published 2022-2023 cycles) and the pre-registered tests (PREREGISTRATION.md).

Rows are per-observation XRT measurements: dicts with obsid, mjd, n_src, n_bkg, exposure_s.
"""
import math

from . import changepoint, ephem, photometry

LN10 = math.log(10.0)

# Validation windows on data Liu et al. 2024 already published (MJD).
C3_PLATEAU = (59642.0, 59699.0)   # rise 59591 + 51 d  ... drop 59699
C4_PLATEAU = (59817.0, 59877.0)   # rise 59778 + 39 d  ... drop 59878
C5_PLATEAU = (60008.0, 60075.0)   # rise 59972 + 38 d  ... last published point 60072
C4_FAINT = (59886.0, 59968.0)     # after the C4 drop, before the C5 rise
G4_WINDOW = (59815.0, 59968.0)
G4_TARGET = (59870.0, 59892.0)    # published C4 drop 59878-59883, +-8 d

# Pre-registered windows on the unpublished data.
POST = (60100.0, 60400.0)         # weekly Swift campaign, 2023 June - 2024 March
ON6 = ephem.on_window(6)          # (60197.5, 60257.5)
OFF5 = ephem.off_window(5)        # (60092.0, 60157.0)
OFF6 = ephem.off_window(6)        # (60282.5, 60347.5)
CP6_WINDOW = (60180.0, 60350.0)
DROP6 = ephem.predicted_drop(6)   # 60265.5
DROP_TOL = 20.0
TS_MIN = 9.0
Z_ONOFF = 3.0
Z_SCAN = 4.0
SCAN_WIDTH = 4
LATE_FROM = 60400.0
GROUP_GAP = 20.0


def in_window(rows, w):
    return [r for r in rows if w[0] <= r["mjd"] <= w[1]]


def overlaps(a, b):
    return a[0] <= b[1] and b[0] <= a[1]


def contrast(a, b):
    """(rate_a - rate_b) / combined error, or None if either stack is missing."""
    if not a or not b:
        return None
    den = math.sqrt(a["rate_err"] ** 2 + b["rate_err"] ** 2)
    if den <= 0.0:
        return 0.0
    return (a["rate"] - b["rate"]) / den


def scan(rows, width=SCAN_WIDTH):
    """Best stack of `width` consecutive observations against all the others (free phase)."""
    rows = sorted(rows, key=lambda r: r["mjd"])
    n_win = len(rows) - width + 1
    best = None
    for i in range(max(0, n_win)):
        inside = rows[i:i + width]
        rest = rows[:i] + rows[i + width:]
        if not rest:
            continue
        a = photometry.stack(inside)
        b = photometry.stack(rest)
        z = contrast(a, b)
        if best is None or z > best["z"]:
            best = {"z": z, "mjd_first": a["mjd_first"], "mjd_last": a["mjd_last"],
                    "rate": a["rate"], "rate_err": a["rate_err"], "n_src": a["n_src"],
                    "rest_rate": b["rate"], "rest_rate_err": b["rate_err"]}
    if best is not None:
        best["n_windows"] = n_win
    return best


def _stack_w(rows, w):
    return photometry.stack(in_window(rows, w))


def _det(s):
    return bool(s and s["detected"])


def _pf(ok):
    return "PASS" if ok else "FAIL"


def _cp(rows, w):
    return changepoint.changepoint(in_window(rows, w), photometry.AREA_RATIO, photometry.EEF)


def _short(cp):
    if cp is None:
        return None
    return {k: v for k, v in cp.items() if k != "profile"}


def run_gates(rows):
    s3 = _stack_w(rows, C3_PLATEAU)
    s4 = _stack_w(rows, C4_PLATEAU)
    s5 = _stack_w(rows, C5_PLATEAU)
    sf = _stack_w(rows, C4_FAINT)
    cp = _cp(rows, G4_WINDOW)
    st = {
        "G1_cycle3_plateau_detected": _pf(_det(s3)),
        "G2_cycle4_and_cycle5_plateaus_detected": _pf(_det(s4) and _det(s5)),
        "G3_faint_phase_below_half_cycle4": _pf(bool(sf and _det(s4)
                                                     and sf["level"] < 0.5 * s4["rate"])),
        "G4_cycle4_drop_recovered": _pf(bool(cp and cp["ts"] >= TS_MIN
                                             and cp["r_before"] > cp["r_after"]
                                             and overlaps(cp["bracket"], G4_TARGET))),
        "G5_fading_cycle3_to_cycle5": _pf(bool(_det(s3) and _det(s5)
                                               and s3["rate"] > s5["rate"])),
    }
    return {"status": st, "all_passed": all(v == "PASS" for v in st.values()),
            "stacks": {"cycle3_plateau": s3, "cycle4_plateau": s4, "cycle5_plateau": s5,
                       "cycle4_faint": sf},
            "cycle4_changepoint": _short(cp)}


def fade_forecast(rows):
    """Per-cycle fading factor from our own cycle 3 and cycle 5 plateau rates."""
    s3 = _stack_w(rows, C3_PLATEAU)
    s5 = _stack_w(rows, C5_PLATEAU)
    if not (_det(s3) and _det(s5)) or s3["rate"] <= 0.0 or s5["rate"] <= 0.0:
        return None
    f = math.sqrt(s5["rate"] / s3["rate"])
    slf = 0.5 / LN10 * math.sqrt((s5["rate_err"] / s5["rate"]) ** 2
                                 + (s3["rate_err"] / s3["rate"]) ** 2)
    return {"per_cycle_factor": f, "sigma_log10": slf, "cycle5_rate": s5["rate"],
            "cycle5_rate_err": s5["rate_err"], "forecast_cycle6_rate": s5["rate"] * f}


def late_groups(rows):
    late = sorted(in_window(rows, (LATE_FROM, 1e9)), key=lambda r: r["mjd"])
    groups = []
    cur = []
    for r in late:
        if cur and r["mjd"] - cur[-1]["mjd"] > GROUP_GAP:
            groups.append(cur)
            cur = []
        cur.append(r)
    if cur:
        groups.append(cur)
    out = []
    for g in groups:
        s = photometry.stack(g)
        item = {k: s[k] for k in ("n_obs", "mjd_first", "mjd_last", "exposure_s", "n_src",
                                  "n_bkg", "rate", "rate_err", "sigma", "ul3_rate", "detected")}
        item["forecast_phases"] = [ephem.phase_class(r["mjd"]) for r in g]
        item["forecast_cycles"] = sorted(set(ephem.cycle_of(r["mjd"]) for r in g))
        out.append(item)
    return out


def grade(rows):
    fade = fade_forecast(rows)
    post = in_window(rows, POST)
    on = _stack_w(rows, ON6)
    off = photometry.stack(in_window(rows, OFF5) + in_window(rows, OFF6))
    z = contrast(on, off)
    sc = scan(post)
    allp = photometry.stack(post)

    if not post:
        outcome = "NO_DATA"
    elif _det(on) and z is not None and z >= Z_ONOFF:
        outcome = "FLARE_ON_SCHEDULE"
    elif sc is not None and sc["z"] >= Z_SCAN:
        outcome = "FLARE_OFF_SCHEDULE"
    elif _det(allp):
        outcome = "STEADY_SOURCE"
    else:
        outcome = "NO_FLARE"

    res = {"outcome": outcome, "z_on_vs_off": z, "on_stack": on, "off_stack": off,
           "post_stack": allp, "free_phase_scan": sc, "fade_forecast": fade,
           "windows": {"post": POST, "on6": ON6, "off5": OFF5, "off6": OFF6,
                       "predicted_rise6": ephem.rise_start(6), "predicted_drop6": DROP6}}

    if outcome in ("NO_FLARE", "STEADY_SOURCE") and on and fade:
        res["depth_reading"] = ("BELOW_FORECAST_LEVEL"
                                if on["ul3_rate"] < fade["forecast_cycle6_rate"]
                                else "NOT_DEEP_ENOUGH")
        res["on_ul3_over_cycle5_rate"] = on["ul3_rate"] / fade["cycle5_rate"]

    # T2: amplitude of cycle 6 against the fading forecast.
    t2 = {"status": "NOT_APPLICABLE"}
    if outcome == "FLARE_ON_SCHEDULE" and fade and on["rate"] > 0.0:
        q = on["rate"] / fade["cycle5_rate"]
        slq = math.sqrt((on["rate_err"] / on["rate"]) ** 2
                        + (fade["cycle5_rate_err"] / fade["cycle5_rate"]) ** 2) / LN10
        d = math.log10(q) - math.log10(fade["per_cycle_factor"])
        sig = math.sqrt(slq ** 2 + fade["sigma_log10"] ** 2)
        if abs(d) <= 2.0 * sig:
            cls = "CONSISTENT"
        elif d < 0.0:
            cls = "FADING_FASTER"
        else:
            cls = "FADING_SLOWER"
        t2 = {"status": cls, "ratio_cycle6_to_cycle5": q, "sigma_log10_ratio": slq,
              "forecast_ratio": fade["per_cycle_factor"], "dlog10": d, "sigma_dlog10": sig}
    res["T2_amplitude"] = t2

    # T3: timing of the cycle-6 drop.
    cp6 = _cp(rows, CP6_WINDOW)
    t3 = {"status": "UNTIMED", "changepoint": _short(cp6)}
    if cp6 and cp6["ts"] >= TS_MIN and cp6["r_before"] > cp6["r_after"]:
        lo, hi = cp6["bracket"]
        if overlaps((lo, hi), (DROP6 - DROP_TOL, DROP6 + DROP_TOL)):
            t3["status"] = "ON_TIME"
        elif hi < DROP6 - DROP_TOL:
            t3["status"] = "EARLY"
        else:
            t3["status"] = "LATE"
        mid = 0.5 * (lo + hi)
        t3["drop_mid_mjd"] = mid
        t3["minus_forecast_days"] = mid - DROP6
        t3["drop_to_drop_recurrence_days"] = (mid - ephem.C4_DROP_MID) / 2.0
        t3["recurrence_halfwidth_days"] = (0.5 * (hi - lo) + 2.5) / 2.0
    res["T3_timing"] = t3

    res["T4_late_epochs"] = late_groups(rows)
    return res
