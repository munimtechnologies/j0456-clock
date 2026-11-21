# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""Swift-XRT aperture photometry on cleaned photon-counting event lists (pure Python).

J0456-20 is a super-soft source (disk temperatures ~50-100 eV, Liu et al. 2023, 2024), so
the band is 0.3-2 keV (PI 30-200), which keeps most of its photons and drops most of the
hard background. Source circle 35" (King-profile PSF fraction ~0.86), background annulus
70"-200". No exposure-map (bad column) or vignetting correction: every comparison in this
test is a ratio of rates measured the same way, so those cancel to first order.
"""
import math

from . import geometry, stats

SRC_R = 35.0
BKG_IN = 70.0
BKG_OUT = 200.0
PI_MIN = 30
PI_MAX = 200
AREA_RATIO = SRC_R ** 2 / (BKG_OUT ** 2 - BKG_IN ** 2)
EEF = geometry.xrt_eef_king(SRC_R)
DETECT_SIGMA = 3.0
DETECT_MIN_COUNTS = 4
UL_CL = 0.9973


def count_aperture(xs, ys, pis, x0, y0, pixscale_arcsec):
    """Counts in the source circle and background annulus around pixel (x0, y0)."""
    n_src = 0
    n_bkg = 0
    r_src2 = SRC_R ** 2
    r_in2 = BKG_IN ** 2
    r_out2 = BKG_OUT ** 2
    s2 = pixscale_arcsec ** 2
    for x, y, pi in zip(xs, ys, pis):
        if pi < PI_MIN or pi > PI_MAX:
            continue
        r2 = ((x - x0) ** 2 + (y - y0) ** 2) * s2
        if r2 <= r_src2:
            n_src += 1
        elif r_in2 <= r2 <= r_out2:
            n_bkg += 1
    return n_src, n_bkg


def measure(n_src, n_bkg, exposure, alpha=AREA_RATIO, eef=EEF):
    """Net rate, error, detection significance and 3-sigma upper limit for one aperture."""
    if exposure <= 0.0:
        raise ValueError("exposure must be positive")
    b = n_bkg * alpha
    mu_b = max(n_bkg, 1) * alpha  # conservative background for the significance
    net = n_src - b
    scale = 1.0 / (exposure * eef)
    p = stats.poisson_tail(n_src, mu_b)
    sig = stats.sigma_from_p(p)
    ul_counts = stats.source_upper_limit(n_src, b, UL_CL)
    detected = (sig >= DETECT_SIGMA) and (n_src >= DETECT_MIN_COUNTS)
    rate = net * scale
    return {
        "n_src": int(n_src),
        "n_bkg": int(n_bkg),
        "exposure_s": float(exposure),
        "bkg_in_src": b,
        "net_counts": net,
        "rate": rate,
        "rate_err": math.sqrt(n_src + alpha ** 2 * n_bkg) * scale,
        "p_bkg": p,
        "sigma": sig,
        "ul3_rate": ul_counts * scale,
        "detected": bool(detected),
        "level": rate if detected else ul_counts * scale,
    }


def stack(rows):
    """Measure the summed counts and exposure of several observations (None if empty)."""
    rows = list(rows)
    if not rows:
        return None
    n_src = sum(r["n_src"] for r in rows)
    n_bkg = sum(r["n_bkg"] for r in rows)
    exp = sum(r["exposure_s"] for r in rows)
    m = measure(n_src, n_bkg, exp)
    m["n_obs"] = len(rows)
    m["mjd_first"] = min(r["mjd"] for r in rows)
    m["mjd_last"] = max(r["mjd"] for r in rows)
    m["obsids"] = [r["obsid"] for r in rows]
    return m
