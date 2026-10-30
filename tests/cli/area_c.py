# GEM-Mapper v3 (GEM3) test suite.
# Copyright (c) 2026 Jacob Schreiber <jacob.schreiber@umassmed.edu>
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of GEM-Mapper v3 (GEM3), distributed under the GNU General
# Public License version 3 or (at your option) any later version; see LICENSE.
"""Helpers for area C tests (paired-end, bisulfite, match selection, MAPQ).

``repeat_genome()`` is a small reference whose repeats are designed by hand,
so the expected placements of a read are known exactly.
"""

from gem3harness import Genome, Rng, revcomp
from gem3harness.sam import READ1, parse_sam

UNIT_LENGTH = 100


class RepeatGenome(object):
	"""One contig, ``chrR`` (100 kbp of random sequence), carrying:

	- ``strata``: a 100 bp unit placed 10 times, once exact and then 2, 3
	  and 4 copies with exactly 1, 2 and 3 substitutions; ``strata_copies``
	  lists (0-based position, strand, substitutions). A read equal to the
	  unit (``strata_read``) has these 10 placements and no others.
	- ``ties``: a 100 bp unit placed 6 times exactly and twice with 2
	  substitutions (``ties_copies``, ``ties_read``).
	- ``dup``: a 2,000 bp block copied exactly at ``dup_starts``.
	- ``triple``: a 300 bp unit at ``triple_copies`` (position, strand),
	  each copy between unique flanks.
	"""

	def __init__(self, genome, strata_read, strata_copies, ties_read, ties_copies, dup_starts,
			dup_length, triple, triple_copies):
		self.genome = genome
		self.strata_read = strata_read
		self.strata_copies = strata_copies
		self.ties_read = ties_read
		self.ties_copies = ties_copies
		self.dup_starts = dup_starts
		self.dup_length = dup_length
		self.triple = triple
		self.triple_copies = triple_copies

	@property
	def seq(self):
		return self.genome["chrR"]


def _place_unit(rng, bg, unit, plan, start, step):
	copies = []
	pos = start
	for subs, count in plan:
		for _ in range(count):
			body = list(unit)
			offsets = set()
			while len(offsets) < subs:
				offsets.add(rng.randint(10, UNIT_LENGTH - 11))
			for o in sorted(offsets):
				body[o] = rng.choice([b for b in "ACGT" if b != body[o]])
			strand = "-" if rng.bernoulli(0.5) else "+"
			text = "".join(body)
			bg[pos:pos + UNIT_LENGTH] = list(revcomp(text) if strand == "-" else text)
			copies.append((pos, strand, subs))
			pos += step
	return copies


def repeat_genome(seed=77):
	rng = Rng(seed)
	bg = list(rng.dna(100000))
	strata = rng.dna(UNIT_LENGTH)
	strata_copies = _place_unit(rng, bg, strata, [(0, 1), (1, 2), (2, 3), (3, 4)], 2000, 3000)
	ties = rng.dna(UNIT_LENGTH)
	ties_copies = _place_unit(rng, bg, ties, [(0, 6), (2, 2)], 40000, 2500)
	dup = rng.dna(2000)
	dup_starts = (62000, 70000)
	for s in dup_starts:
		bg[s:s + len(dup)] = list(dup)
	triple = rng.dna(300)
	triple_copies = ((80000, "+"), (86000, "-"), (92000, "+"))
	for p, st in triple_copies:
		bg[p:p + len(triple)] = list(revcomp(triple) if st == "-" else triple)
	genome = Genome([("chrR", "".join(bg))])
	return RepeatGenome(genome, strata, strata_copies, ties, ties_copies, dup_starts, len(dup),
		triple, triple_copies)


def gap_affine_score(ref, pos, ops, seq, match=1, mismatch=4, gap_open=6, gap_extend=1):
	"""Gap-affine score of an alignment (gem3-mapper's AS:i): a k-base gap
	costs gap_open + k * gap_extend; N never matches."""
	score = 0
	r, q = pos, 0
	for op, n in ops:
		if op == "M":
			for k in range(n):
				same = ref[r + k].upper() == seq[q + k].upper() and seq[q + k].upper() != "N"
				score += match if same else -mismatch
			r += n
			q += n
		elif op == "I":
			score -= gap_open + gap_extend * n
			q += n
		elif op == "D":
			score -= gap_open + gap_extend * n
			r += n
		elif op == "S":
			q += n
	return score


def event_distance(ref, pos, ops, seq):
	"""Mismatched bases plus the number of insertion and deletion runs."""
	events = 0
	r, q = pos, 0
	for op, n in ops:
		if op == "M":
			events += sum(1 for k in range(n) if ref[r + k].upper() != seq[q + k].upper()
				or seq[q + k].upper() == "N")
			r += n
			q += n
		elif op in "ID":
			events += 1
			if op == "I":
				q += n
			else:
				r += n
		elif op == "S":
			q += n
	return events


def primary_pairs(sam):
	"""qname -> (read-1 record, read-2 record), primary records only."""
	out = {}
	for rec in sam.records:
		if not rec.is_primary:
			continue
		ends = out.setdefault(rec.qname, [None, None])
		idx = 0 if rec.has(READ1) else 1
		assert ends[idx] is None, ("two primary records for", rec.qname, idx)
		ends[idx] = rec
	return dict((k, tuple(v)) for k, v in out.items())


def split_predictors(stdout):
	"""--mapq-model dump-predictors writes one line of predictors per mapped
	read to stdout next to the SAM. Returns (Sam, [predictor fields])."""
	text = stdout.decode("ascii") if isinstance(stdout, bytes) else stdout
	sam_lines, predictors = [], []
	for line in text.split("\n"):
		if not line:
			continue
		if line.startswith("@") or line.count("\t") >= 10:
			sam_lines.append(line)
		else:
			predictors.append(line.split("  "))
	return parse_sam("\n".join(sam_lines) + "\n"), predictors


PREDICTOR_FIELDS = ("tag", "delta", "wdelta", "edit", "event", "swg", "equal", "sigma", "mapq_1",
	"mapq_2", "sub_edit", "sub_event", "sub_swg", "sub_sigma", "scand", "acand", "mcs1", "mcs2",
	"len", "arlen", "mrlen", "kmer")
