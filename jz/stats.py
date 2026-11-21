# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""Dependency-free Poisson statistics for low-count X-ray photometry.

Same code as munimtechnologies/at2018fyk-return (fyk/stats.py):
- regularized incomplete gamma functions (Numerical Recipes style),
- Poisson tail probability of the source-aperture counts given background,
- one-sided Gaussian-equivalent significance,
- Bayesian upper limit on net source counts with a flat prior on s >= 0
  (the Kraft, Burrows & Nousek 1991 posterior), solved by bisection.
"""
import math

_EPS = 1e-15
_TINY = 1e-300


def _gser(a, x):
    """Series for the regularized lower incomplete gamma P(a, x); use for x < a + 1."""
    if x <= 0.0:
        return 0.0
    ap = a
    total = 1.0 / a
    term = total
    for _ in range(100000):
        ap += 1.0
        term *= x / ap
        total += term
        if abs(term) < abs(total) * _EPS:
            break
    return total * math.exp(-x + a * math.log(x) - math.lgamma(a))


def _gcf(a, x):
    """Continued fraction for the regularized upper incomplete gamma Q(a, x); use for x >= a + 1."""
    b = x + 1.0 - a
    c = 1.0 / _TINY
    d = 1.0 / b
    h = d
    for i in range(1, 100000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < _TINY:
            d = _TINY
        c = b + an / c
        if abs(c) < _TINY:
            c = _TINY
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < _EPS:
            break
    return math.exp(-x + a * math.log(x) - math.lgamma(a)) * h


def gammainc_lower(a, x):
    """Regularized lower incomplete gamma P(a, x)."""
    if a <= 0.0:
        raise ValueError("a must be positive")
    if x <= 0.0:
        return 0.0
    if x < a + 1.0:
        return _gser(a, x)
    return 1.0 - _gcf(a, x)


def gammainc_upper(a, x):
    """Regularized upper incomplete gamma Q(a, x) = 1 - P(a, x)."""
    if a <= 0.0:
        raise ValueError("a must be positive")
    if x <= 0.0:
        return 1.0
    if x < a + 1.0:
        return 1.0 - _gser(a, x)
    return _gcf(a, x)


def poisson_tail(n, mu):
    """P(X >= n) for X ~ Poisson(mu)."""
    if n <= 0:
        return 1.0
    if mu <= 0.0:
        return 0.0
    return gammainc_lower(float(n), float(mu))


def sigma_from_p(p):
    """One-sided Gaussian-equivalent significance z with P(Z > z) = p (capped at 40)."""
    if p >= 0.5:
        return 0.0
    if p <= 0.0:
        return 40.0
    lo, hi = 0.0, 40.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if 0.5 * math.erfc(mid / math.sqrt(2.0)) > p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def source_upper_limit(n, b, cl=0.9973):
    """Upper limit S on net source counts: posterior p(s) ~ (s+b)^n exp(-(s+b)), s >= 0.

    The posterior CDF is [Q(n+1, b) - Q(n+1, S+b)] / Q(n+1, b).
    """
    if n < 0 or b < 0.0:
        raise ValueError("n and b must be non-negative")
    a = float(n) + 1.0
    qb = gammainc_upper(a, b)
    if qb <= 0.0:
        qb = _TINY

    def cdf(s):
        return (qb - gammainc_upper(a, s + b)) / qb

    lo = 0.0
    hi = max(10.0, n + 20.0 * math.sqrt(n + 1.0) + 50.0)
    for _ in range(60):
        if cdf(hi) >= cl:
            break
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if cdf(mid) < cl:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)
