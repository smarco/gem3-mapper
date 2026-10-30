# GEM-Mapper v3 (GEM3) test suite.
# Copyright (c) 2026 Jacob Schreiber <jacob.schreiber@umassmed.edu>
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of GEM-Mapper v3 (GEM3), distributed under the GNU General
# Public License version 3 or (at your option) any later version; see LICENSE.
"""Seeded pseudo-random numbers that are identical on every Python version.

Python's ``random`` module only promises that ``random()`` keeps its sequence
across versions; ``randrange``, ``choice``, ``shuffle`` and ``gauss`` may
change. Fixtures whose bytes are pinned must not depend on that, so every
helper here is a SplitMix64 generator written out in full.
"""

import math

_MASK64 = (1 << 64) - 1


class Rng(object):
	"""SplitMix64. ``Rng(seed)`` gives the same stream on every platform."""

	def __init__(self, seed):
		self.state = int(seed) & _MASK64

	def next_u64(self):
		self.state = (self.state + 0x9E3779B97F4A7C15) & _MASK64
		z = self.state
		z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & _MASK64
		z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & _MASK64
		return z ^ (z >> 31)

	def random(self):
		"""Uniform float in [0, 1) with 53 random bits."""
		return (self.next_u64() >> 11) * (1.0 / (1 << 53))

	def randint(self, lo, hi):
		"""Uniform integer in [lo, hi], both ends included."""
		if hi < lo:
			raise ValueError("empty range [%d, %d]" % (lo, hi))
		span = hi - lo + 1
		# Rejection sampling keeps the distribution exactly uniform.
		limit = (_MASK64 + 1) - ((_MASK64 + 1) % span)
		while True:
			x = self.next_u64()
			if x < limit:
				return lo + (x % span)

	def choice(self, seq):
		return seq[self.randint(0, len(seq) - 1)]

	def bernoulli(self, p):
		return self.random() < p

	def normal(self, mean, sd):
		"""Box-Muller; uses two uniforms per call."""
		u1 = self.random()
		u2 = self.random()
		if u1 <= 0.0:
			u1 = 1.0 / (1 << 53)
		return mean + sd * math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * u2)

	def dna(self, length, alphabet="ACGT"):
		k = len(alphabet)
		return "".join(alphabet[self.next_u64() % k] for _ in range(length))

	def fork(self, salt):
		"""An independent stream derived from this one and ``salt``."""
		return Rng(hash_u64(self.next_u64(), salt))


def hash_u64(*values):
	"""Deterministic 64-bit hash of integers and strings (not Python's hash())."""
	h = 0x243F6A8885A308D3
	for v in values:
		if isinstance(v, str):
			data = v.encode("utf-8")
			for b in data:
				h = ((h ^ b) * 0x100000001B3) & _MASK64
			v = len(data)
		h = (h ^ (int(v) & _MASK64)) & _MASK64
		h = (h + 0x9E3779B97F4A7C15) & _MASK64
		h = ((h ^ (h >> 30)) * 0xBF58476D1CE4E5B9) & _MASK64
		h = ((h ^ (h >> 27)) * 0x94D049BB133111EB) & _MASK64
		h ^= h >> 31
	return h
