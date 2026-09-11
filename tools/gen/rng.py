"""Reproducible random streams.

Every random decision in the generator draws from a *named* stream derived from
(seed, salt, name).  Streams are independent, so adding a new query type does
not perturb the graph, and regenerating queries for an existing graph reproduces
it exactly.

The stream name is hashed with BLAKE2b rather than Python's ``hash()``, which is
salted per process and therefore not reproducible across runs.  ``random.Random``
seeded from an int uses Mersenne Twister, whose sequence is stable across CPython
versions -- which is why the generator deliberately avoids numpy for anything
that consumes randomness.
"""

import hashlib
import random


def derive_seed(seed: int, salt: str, name: str) -> int:
    payload = b"\x00".join((str(int(seed)).encode(), salt.encode(), name.encode()))
    return int.from_bytes(hashlib.blake2b(payload, digest_size=32).digest(), "big")


def stream(seed: int, salt: str, name: str) -> random.Random:
    """An independent, reproducible RNG for the named stream."""
    return random.Random(derive_seed(seed, salt, name))
