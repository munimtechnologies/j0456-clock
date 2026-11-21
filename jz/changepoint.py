# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""One-step Poisson change-point for a sparse X-ray light curve (pure Python).

Observation i has n_i source-aperture counts with expected background b_i = alpha * n_bkg_i
(floored at 1e-3 counts) and effective exposure e_i = exposure_i * eef. The model is a
constant source rate r1 before the change and r2 after it (each fitted >= 0 by maximum
likelihood), with the change between two consecutive observations and at least two
observations on each side. TS = 2 (lnL_best_split - lnL_one_rate).
"""
import math

B_FLOOR = 1e-3


def _points(rows, alpha, eef):
    pts = []
    for r in sorted(rows, key=lambda r: r["mjd"]):
        pts.append((float(r["mjd"]), int(r["n_src"]), max(r["n_bkg"] * alpha, B_FLOOR),
                    float(r["exposure_s"]) * eef))
    return pts


def _loglik(seg, rate):
    total = 0.0
    for _, n, b, e in seg:
        mu = b + rate * e
        if n > 0:
            total += n * math.log(mu)
        total -= mu
    return total


def fit_rate(seg):
    """Maximum-likelihood source rate (>= 0) for one segment; returns (rate, lnL)."""
    e_sum = sum(p[3] for p in seg)
    n_sum = sum(p[1] for p in seg)
    if e_sum <= 0.0:
        return 0.0, _loglik(seg, 0.0)

    def deriv(r):
        return sum(n * e / (b + r * e) for _, n, b, e in seg) - e_sum

    if deriv(0.0) <= 0.0:
        rate = 0.0
    else:
        lo = 0.0
        hi = n_sum / e_sum + 1e-12
        while deriv(hi) > 0.0:
            hi *= 2.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if deriv(mid) > 0.0:
                lo = mid
            else:
                hi = mid
        rate = 0.5 * (lo + hi)
    return rate, _loglik(seg, rate)


def changepoint(rows, alpha, eef, min_side=2):
    pts = _points(rows, alpha, eef)
    n = len(pts)
    if n < 2 * min_side:
        return None
    _, l0 = fit_rate(pts)
    profile = []
    for k in range(min_side, n - min_side + 1):
        r1, l1 = fit_rate(pts[:k])
        r2, l2 = fit_rate(pts[k:])
        profile.append({"k": k, "after_mjd": pts[k - 1][0], "before_mjd": pts[k][0],
                        "lnL": l1 + l2, "r_before": r1, "r_after": r2})
    best = max(profile, key=lambda p: p["lnL"])
    conf = [p for p in profile if p["lnL"] >= best["lnL"] - 1.92]
    return {
        "n_obs": n,
        "bracket": [best["after_mjd"], best["before_mjd"]],
        "r_before": best["r_before"],
        "r_after": best["r_after"],
        "ts": 2.0 * (best["lnL"] - l0),
        "conf95_range": [min(p["after_mjd"] for p in conf), max(p["before_mjd"] for p in conf)],
        "profile": profile,
    }
