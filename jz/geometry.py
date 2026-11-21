# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""Sky geometry: gnomonic (TAN) projection, angular separation, XRT PSF fraction, dates.

Same code as munimtechnologies/at2018fyk-return (fyk/geometry.py).
"""
import datetime
import math

MJD_ZERO = datetime.datetime(1858, 11, 17)


def tan_world_to_pixel(ra, dec, crval1, crval2, crpix1, crpix2, cdelt1, cdelt2):
    """Sky (deg) -> FITS pixel for a RA---TAN / DEC--TAN axis pair with no rotation.

    The intermediate coordinate xi is positive toward east (increasing RA); with the
    usual negative cdelt1, east is toward smaller x.
    """
    a = math.radians(ra)
    d = math.radians(dec)
    a0 = math.radians(crval1)
    d0 = math.radians(crval2)
    cosc = math.sin(d0) * math.sin(d) + math.cos(d0) * math.cos(d) * math.cos(a - a0)
    if cosc <= 0.0:
        raise ValueError("position is more than 90 degrees from the projection centre")
    xi = math.cos(d) * math.sin(a - a0) / cosc
    eta = (math.cos(d0) * math.sin(d) - math.sin(d0) * math.cos(d) * math.cos(a - a0)) / cosc
    x = crpix1 + math.degrees(xi) / cdelt1
    y = crpix2 + math.degrees(eta) / cdelt2
    return x, y


def angsep_arcsec(ra1, dec1, ra2, dec2):
    """Great-circle separation in arcseconds (haversine)."""
    p1 = math.radians(dec1)
    p2 = math.radians(dec2)
    dl = math.radians(ra2 - ra1)
    h = math.sin((p2 - p1) / 2.0) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2.0) ** 2
    return math.degrees(2.0 * math.asin(min(1.0, math.sqrt(h)))) * 3600.0


def xrt_eef_king(r_arcsec, rc=5.8, beta=1.55):
    """Fraction of the Swift-XRT PSF inside radius r (King profile, Moretti et al. 2005, ~1.5 keV)."""
    if r_arcsec <= 0.0:
        return 0.0
    return 1.0 - (1.0 + (r_arcsec / rc) ** 2) ** (1.0 - beta)


def mjd_to_datetime(mjd):
    return MJD_ZERO + datetime.timedelta(days=float(mjd))


def mjd_to_archive_month(mjd, shift=0):
    """Swift archive folder 'YYYY_MM' for the month of the observation start, optionally shifted."""
    t = mjd_to_datetime(mjd)
    m = t.year * 12 + (t.month - 1) + shift
    return "%04d_%02d" % (m // 12, m % 12 + 1)
