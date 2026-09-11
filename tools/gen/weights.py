"""Edge-weight models.

Weights are deliberately not ``round(1e6 * euclidean_length)``: that would make
the graph metric a scaled copy of the Euclidean metric -- see families/road2d.py.
"""

import math


def uniform(rng, lo, hi):
    return rng.randint(lo, hi)


def loguniform(rng, lo, hi):
    """Integer drawn uniformly in log space, so every decade is equally likely.

    Spanning several decades of weight is a property the instances vary on
    purpose; see the per-instance ranges in tools/manifest/public.json.
    """
    return max(lo, min(hi, int(round(math.exp(rng.uniform(math.log(lo), math.log(hi)))))))


def bounded_lognormal(rng, sigma, clamp_sigmas=2.0):
    """Multiplicative noise in [exp(-c*sigma), exp(+c*sigma)].

    Clamped rather than open-ended so that ``min_e w(e)/||u-v||`` is
    predictable from the parameters instead of being set by the single
    unluckiest edge in a two-million-edge draw.
    """
    hi = math.exp(clamp_sigmas * sigma)
    return min(hi, max(1.0 / hi, math.exp(sigma * rng.gauss(0.0, 1.0))))


def make(spec, rng_factory):
    """Return f(rng) -> weight for a spec like ``uniform:1:1000``."""
    parts = spec.split(":")
    kind = parts[0]
    if kind == "uniform":
        lo, hi = int(parts[1]), int(parts[2])
        return lambda rng: uniform(rng, lo, hi)
    if kind == "loguniform":
        lo, hi = int(parts[1]), int(parts[2])
        return lambda rng: loguniform(rng, lo, hi)
    raise ValueError(f"unknown weight spec: {spec}")
