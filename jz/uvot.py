# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""Secondary: Swift-UVOT 5" aperture photometry on archive sky images (host + transient).

Same code as munimtechnologies/at2018fyk-return (fyk/uvot.py). No host-galaxy subtraction:
the host is constant, so epoch-to-epoch changes in one filter are what we use.
Needs numpy and astropy; the unit tests never import it.
"""
import math

# AB zero points (Breeveld et al. 2011).
ZP_AB = {"w2": 19.11, "m2": 18.54, "w1": 18.95, "uu": 19.36, "bb": 18.98, "vv": 17.88}
FRAME_TIME = 0.0110329
SRC_R = 5.0
BKG_IN = 25.0
BKG_OUT = 45.0


def phot_sky_image(path, ra, dec):
    import numpy as np
    from astropy.io import fits
    from astropy.stats import sigma_clipped_stats
    from astropy.wcs import WCS
    from astropy.wcs.utils import proj_plane_pixel_scales

    src = 0.0
    bkg_counts = 0.0
    bkg_var = 0.0
    exposure = 0.0
    n_ext = 0
    with fits.open(path, memmap=False) as h:
        for hdu in h[1:]:
            if hdu.data is None or hdu.data.ndim != 2:
                continue
            hdr = hdu.header
            exp = float(hdr.get("EXPOSURE", 0.0) or 0.0)
            if exp <= 0.0:
                continue
            w = WCS(hdr).celestial
            x0, y0 = w.all_world2pix([[ra, dec]], 0)[0]
            scale = float(abs(proj_plane_pixel_scales(w)[0])) * 3600.0
            data = np.asarray(hdu.data, dtype=float)
            ny, nx = data.shape
            m = int(math.ceil(BKG_OUT / scale)) + 2
            if not (m <= x0 < nx - m and m <= y0 < ny - m):
                continue
            xi0, yi0 = int(round(x0)), int(round(y0))
            cut = data[yi0 - m:yi0 + m + 1, xi0 - m:xi0 + m + 1]
            yy, xx = np.mgrid[yi0 - m:yi0 + m + 1, xi0 - m:xi0 + m + 1]
            r = np.hypot(xx - x0, yy - y0) * scale
            s_mask = r <= SRC_R
            b_mask = (r >= BKG_IN) & (r <= BKG_OUT)
            if s_mask.sum() == 0 or b_mask.sum() < 50:
                continue
            mean, _, std = sigma_clipped_stats(cut[b_mask], sigma=3.0)
            n_s = int(s_mask.sum())
            src += float(cut[s_mask].sum())
            bkg_counts += float(mean) * n_s
            bkg_var += (float(std) ** 2 / max(1, int(b_mask.sum()))) * n_s ** 2
            exposure += exp
            n_ext += 1
    if n_ext == 0 or exposure <= 0.0:
        return None
    raw = src / exposure
    coi = 1.0
    if raw > 0.0 and raw * FRAME_TIME < 0.99:
        coi = (-math.log(1.0 - raw * FRAME_TIME) / FRAME_TIME) / raw
    net = (src - bkg_counts) / exposure * coi
    err = math.sqrt(max(src, 1.0) + bkg_var) / exposure * coi
    return {"net_rate": net, "net_rate_err": err, "exposure_s": exposure, "n_ext": n_ext,
            "coi_factor": coi, "src_counts": src, "bkg_counts_in_src": bkg_counts}


def ab_mag(rate, filt):
    if rate is None or rate <= 0.0:
        return None
    return ZP_AB[filt] - 2.5 * math.log10(rate)
