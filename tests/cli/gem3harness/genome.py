# GEM-Mapper v3 (GEM3) test suite.
# Copyright (c) 2026 Jacob Schreiber <jacob.schreiber@umassmed.edu>
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of GEM-Mapper v3 (GEM3), distributed under the GNU General
# Public License version 3 or (at your option) any later version; see LICENSE.
"""Deterministic synthetic genomes with annotated hard regions.

``synthetic_genome(seed)`` returns a :class:`Genome`: contigs of random
sequence carrying interspersed repeat copies (exact and diverged), tandem and
homopolymer low-complexity runs, N runs (including at contig ends),
lowercase soft-masked stretches and very short contigs. Every feature is
recorded in ``Genome.features`` so a read simulator can avoid or target it.
The same seed gives the same bytes on every platform and Python version.
"""

import hashlib
import os

from .rng import Rng

COMPLEMENT = {
	"A": "T", "C": "G", "G": "C", "T": "A", "N": "N",
	"a": "t", "c": "g", "g": "c", "t": "a", "n": "n",
	"R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
	"B": "V", "V": "B", "D": "H", "H": "D",
}


def revcomp(seq):
	return "".join(COMPLEMENT.get(c, c) for c in reversed(seq))


class Feature(object):
	"""A half-open interval [start, end) on ``contig`` of a given ``kind``.

	kinds: ``repeat`` (copy of a repeat family; ``family`` and ``diverged``
	set), ``tandem`` (short-unit tandem repeat), ``homopolymer``, ``N``,
	``softmask`` (lowercase), ``short_contig`` (the whole contig).
	"""

	__slots__ = ("kind", "contig", "start", "end", "info")

	def __init__(self, kind, contig, start, end, **info):
		self.kind = kind
		self.contig = contig
		self.start = start
		self.end = end
		self.info = info

	def overlaps(self, contig, start, end, margin=0):
		return (contig == self.contig and start < self.end + margin
				and self.start - margin < end)

	def __repr__(self):
		return "Feature(%s, %s:%d-%d, %r)" % (self.kind, self.contig, self.start, self.end, self.info)


class Genome(object):
	"""Ordered contigs (name -> sequence, case preserved) plus features."""

	def __init__(self, contigs, features=()):
		self.names = [name for name, _ in contigs]
		self.seqs = dict(contigs)
		self.features = list(features)

	def __getitem__(self, name):
		return self.seqs[name]

	def __iter__(self):
		return iter(self.names)

	def __len__(self):
		return len(self.names)

	def lengths(self):
		return [(n, len(self.seqs[n])) for n in self.names]

	def total_length(self):
		return sum(len(s) for s in self.seqs.values())

	def upper(self, name):
		return self.seqs[name].upper()

	def fasta_text(self, line_width=60):
		out = []
		for name in self.names:
			seq = self.seqs[name]
			out.append(">" + name + "\n")
			for i in range(0, len(seq), line_width):
				out.append(seq[i:i + line_width] + "\n")
		return "".join(out)

	def write_fasta(self, path, line_width=60):
		with open(path, "w") as fh:
			fh.write(self.fasta_text(line_width))
		return path

	def digest(self):
		return hashlib.sha256(self.fasta_text().encode("ascii")).hexdigest()

	def hard_regions(self, contig, start, end, margin=0, kinds=None):
		"""Features overlapping [start, end) on ``contig`` (soft-masking excluded by default)."""
		if kinds is None:
			kinds = ("repeat", "tandem", "homopolymer", "N", "short_contig")
		return [f for f in self.features
				if f.kind in kinds and f.overlaps(contig, start, end, margin)]


DEFAULT_CONTIGS = (
	("chr1", 120000),
	("chr2", 80000),
	("chr3", 40000),
	("chrShort1", 400),
	("chrShort2", 90),
)


def synthetic_genome(seed=1, contigs=DEFAULT_CONTIGS, repeat_families=2,
		copies_per_family=4, repeat_length=(400, 1200), divergence=0.03,
		n_runs=3, tandems=4, homopolymers=4, softmask=6, edge_n=True):
	"""Build a :class:`Genome`.

	Contigs shorter than 1 kbp stay plain random sequence (tagged
	``short_contig``). Longer contigs get, in order: random sequence; repeat
	copies (the first copy of each family exact, the others diverged by
	``divergence`` substitutions per base); tandem repeats (unit 2-6 bp,
	60-300 bp long); homopolymers (25-80 bp); N runs (20-600 bp) and, when
	``edge_n``, an N run at the start of the first contig and at the end of
	the second; lowercase soft-masked stretches (100-800 bp, some over
	repeats).
	"""
	rng = Rng(seed)
	seqs = []
	for name, length in contigs:
		seqs.append([name, list(rng.dna(length))])
	features = []
	long_contigs = [i for i, (name, s) in enumerate(seqs) if len(s) >= 1000]
	for i, (name, s) in enumerate(seqs):
		if len(s) < 1000:
			features.append(Feature("short_contig", name, 0, len(s)))

	occupied = {}  # contig index -> list of (start, end)

	def free_slot(length, pad=50):
		for _ in range(1000):
			ci = rng.choice(long_contigs)
			L = len(seqs[ci][1])
			if L < length + 2 * pad + 200:
				continue
			start = rng.randint(pad + 100, L - length - pad - 100)
			end = start + length
			if all(end + pad <= a or b + pad <= start for a, b in occupied.get(ci, [])):
				occupied.setdefault(ci, []).append((start, end))
				return ci, start, end
		raise RuntimeError("could not place a feature of length %d" % length)

	if long_contigs:
		for fam in range(repeat_families):
			unit = rng.dna(rng.randint(repeat_length[0], repeat_length[1]))
			for copy in range(copies_per_family):
				body = list(unit)
				diverged = copy > 0 and divergence > 0
				if diverged:
					for j in range(len(body)):
						if rng.bernoulli(divergence):
							body[j] = rng.choice([b for b in "ACGT" if b != body[j]])
				if rng.bernoulli(0.5):
					body = list(revcomp("".join(body)))
				ci, start, end = free_slot(len(body))
				seqs[ci][1][start:end] = body
				features.append(Feature("repeat", seqs[ci][0], start, end,
					family=fam, copy=copy, diverged=diverged))
		for _ in range(tandems):
			unit = rng.dna(rng.randint(2, 6))
			length = rng.randint(60, 300)
			body = (unit * (length // len(unit) + 1))[:length]
			ci, start, end = free_slot(length)
			seqs[ci][1][start:end] = list(body)
			features.append(Feature("tandem", seqs[ci][0], start, end, unit=unit))
		for _ in range(homopolymers):
			base = rng.choice("ACGT")
			length = rng.randint(25, 80)
			ci, start, end = free_slot(length)
			seqs[ci][1][start:end] = [base] * length
			features.append(Feature("homopolymer", seqs[ci][0], start, end, base=base))
		for _ in range(n_runs):
			length = rng.randint(20, 600)
			ci, start, end = free_slot(length)
			seqs[ci][1][start:end] = ["N"] * length
			features.append(Feature("N", seqs[ci][0], start, end))
		if edge_n:
			first = seqs[long_contigs[0]]
			length = rng.randint(50, 300)
			first[1][0:length] = ["N"] * length
			features.append(Feature("N", first[0], 0, length))
			if len(long_contigs) > 1:
				second = seqs[long_contigs[1]]
				L = len(second[1])
				length = rng.randint(50, 300)
				second[1][L - length:L] = ["N"] * length
				features.append(Feature("N", second[0], L - length, L))
		repeats = [f for f in features if f.kind == "repeat"]
		for k in range(softmask):
			if k % 2 == 0 and repeats:
				target = rng.choice(repeats)
				ci = [x[0] for x in seqs].index(target.contig)
				start, end = target.start, target.end
			else:
				length = rng.randint(100, 800)
				ci = rng.choice(long_contigs)
				L = len(seqs[ci][1])
				start = rng.randint(400, L - length - 400)
				end = start + length
			body = seqs[ci][1]
			for j in range(start, end):
				body[j] = body[j].lower()
			features.append(Feature("softmask", seqs[ci][0], start, end))
	contigs_out = [(name, "".join(s)) for name, s in seqs]
	features.sort(key=lambda f: ([n for n, _ in contigs_out].index(f.contig), f.start, f.kind))
	return Genome(contigs_out, features)


def write_fasta_records(path, records, line_width=60):
	"""Write ``[(name, seq), ...]`` as FASTA (used for odd hand-made references)."""
	with open(path, "w") as fh:
		for name, seq in records:
			fh.write(">" + name + "\n")
			if line_width:
				for i in range(0, len(seq), line_width):
					fh.write(seq[i:i + line_width] + "\n")
			else:
				fh.write(seq + "\n")
	return os.fspath(path)
