# GEM-Mapper v3 (GEM3) test suite.
# Copyright (c) 2026 Jacob Schreiber <jacob.schreiber@umassmed.edu>
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of GEM-Mapper v3 (GEM3), distributed under the GNU General
# Public License version 3 or (at your option) any later version; see LICENSE.
"""Area B helpers (search and alignment): reference implementations and inputs.

Reference implementations, written from the definitions rather than from
gem3-mapper's code:

- ``gotoh_best``: the best gap-affine score of a whole read against any
  substring of a reference window (match +A, mismatch -B, a gap of k bases
  -(O + k*X); N never matches), the score gem3-mapper reports as AS.
- ``alignment_score``: the gap-affine score of a given SAM alignment.
- ``edit_distance_best``: the least edit distance of a whole read against
  any substring of a window.
- ``myers_hits`` / ``hamming_hits``: every place in a contig, on both
  strands, where a read aligns within k edits (bit-parallel Myers) or k
  mismatches (no gaps).

Inputs:

- ``planted_genome``: random contigs carrying copies of a probe read with a
  known number of edits, so the set of loci within E errors is known.
- ``chimeras``: reads whose first or last L bases come from the genome and
  whose rest matches nothing, for local alignment.
"""

from gem3harness import Genome, Rng, revcomp

INF = float("inf")


def _same(a, b):
	return a == b and a != "N"


def alignment_score(ref, pos, ops, seq, A=1, B=4, O=6, X=1):
	"""Gap-affine score of SEQ aligned at 0-based ``pos`` with CIGAR ``ops``
	(soft clips score nothing)."""
	score = 0
	r, q = pos, 0
	for op, n in ops:
		if op == "M":
			for k in range(n):
				score += A if _same(ref[r + k].upper(), seq[q + k].upper()) else -B
			r += n
			q += n
		elif op == "I":
			score -= O + X * n
			q += n
		elif op == "D":
			score -= O + X * n
			r += n
		elif op == "S":
			q += n
	return score


def gotoh_best(read, window, A=1, B=4, O=6, X=1):
	"""Best gap-affine score of all of ``read`` against any substring of
	``window`` (both ends of the window free)."""
	read = read.upper()
	window = window.upper()
	m = len(read)
	NEG = -10 ** 9
	go = O + X
	# Column-major over the window; H/F per read row, E carried along rows.
	H = [0] + [-(O + X * i) for i in range(1, m + 1)]
	F = [NEG] * (m + 1)
	best = H[m]
	for c in window:
		Hn = [0] * (m + 1)
		E = NEG
		diag = H[0]
		for i in range(1, m + 1):
			# F: gap consuming window bases (SAM D), from the previous column
			F[i] = max(H[i] - go, F[i] - X)
			# E: gap consuming read bases (SAM I), down this column
			E = max(Hn[i - 1] - go, E - X)
			s = diag + (A if _same(read[i - 1], c) else -B)
			diag = H[i]
			Hn[i] = max(s, E, F[i])
		H = Hn
		if H[m] > best:
			best = H[m]
	return best


def edit_distance_best(read, window):
	"""Least edit distance of all of ``read`` against any substring of ``window``."""
	read = read.upper()
	window = window.upper()
	m = len(read)
	prev = list(range(m + 1))
	best = prev[m]
	for c in window:
		cur = [0] * (m + 1)
		for i in range(1, m + 1):
			cur[i] = min(prev[i - 1] + (0 if _same(read[i - 1], c) else 1), prev[i] + 1, cur[i - 1] + 1)
		prev = cur
		best = min(best, cur[m])
	return best


def myers_hits(read, text, k):
	"""0-based text end positions (inclusive) where all of ``read`` aligns to
	a substring of ``text`` ending there with at most ``k`` edits, with the
	distance: [(end, distance)]. N matches nothing."""
	read = read.upper()
	m = len(read)
	mask = (1 << m) - 1
	high = 1 << (m - 1)
	peq = {}
	for i, c in enumerate(read):
		if c != "N":
			peq[c] = peq.get(c, 0) | (1 << i)
	pv, mv, score = mask, 0, m
	out = []
	for j, c in enumerate(text.upper()):
		eq = peq.get(c, 0)
		xv = eq | mv
		xh = (((eq & pv) + pv) ^ pv) | eq
		ph = mv | (~(xh | pv) & mask)
		mh = pv & xh
		if ph & high:
			score += 1
		elif mh & high:
			score -= 1
		ph = (ph << 1) & mask
		mh = (mh << 1) & mask
		pv = mh | (~(xv | ph) & mask)
		mv = ph & xv
		if score <= k:
			out.append((j, score))
	return out


def hamming_hits(read, text, k):
	"""[(start, mismatches)] for every offset where ``read`` placed without
	gaps has at most ``k`` mismatches (N matches nothing)."""
	read = read.upper()
	text = text.upper()
	m = len(read)
	out = []
	for s in range(len(text) - m + 1):
		mism = 0
		for a, b in zip(read, text[s:s + m]):
			if not _same(a, b):
				mism += 1
				if mism > k:
					break
		if mism <= k:
			out.append((s, mism))
	return out


def loci_within(genome, read, k, model="edit"):
	"""Loci where ``read`` aligns within ``k`` edits (or mismatches with
	``model='hamming'``) on either strand: a set of (contig, strand,
	end_or_start, distance) with one entry per cluster of adjacent hits,
	keeping its smallest distance. For 'edit' the coordinate is the 0-based
	end of the forward-strand alignment; for 'hamming' the 0-based start."""
	out = set()
	for name in genome:
		seq = genome[name]
		for strand, r in (("+", read), ("-", revcomp(read))):
			if model == "edit":
				hits = myers_hits(r, seq, k)
			else:
				hits = hamming_hits(r, seq, k)
			cluster = []
			for pos, d in hits + [(None, None)]:
				if cluster and (pos is None or pos - cluster[-1][0] > len(read) // 2):
					best = min(cluster, key=lambda x: (x[1], x[0]))
					out.add((name, strand, best[0], best[1]))
					cluster = []
				if pos is not None:
					cluster.append((pos, d))
	return out


def _edit(rng, s, subs, ins, dels):
	"""``s`` with spaced substitutions, insertions and deletions, none in the
	first or last 10 bases."""
	s = list(s)
	n = subs + ins + dels
	if n == 0:
		return "".join(s)
	step = (len(s) - 20) // (n + 1)
	slots = [10 + step * (i + 1) for i in range(n)]
	kinds = ["S"] * subs + ["I"] * ins + ["D"] * dels
	for pos, kind in sorted(zip(slots, kinds), reverse=True):
		if kind == "S":
			s[pos] = rng.choice([b for b in "ACGT" if b != s[pos]])
		elif kind == "I":
			s.insert(pos, rng.choice([b for b in "ACGT" if b != s[pos] and b != s[pos - 1]]))
		else:
			del s[pos]
	return "".join(s)


# (contig, 0-based start, strand, substitutions, insertions, deletions)
PLANTED_COPIES = (
	("pA", 500, "+", 0, 0, 0),
	("pA", 1500, "-", 1, 0, 0),
	("pA", 2500, "+", 2, 0, 0),
	("pA", 3500, "+", 0, 0, 1),
	("pB", 500, "-", 0, 2, 0),
	("pB", 1500, "+", 3, 0, 0),
	("pB", 2500, "-", 4, 0, 0),
	("pA", 4500, "+", 6, 0, 0),
)


class Copy(object):
	__slots__ = ("contig", "start", "strand", "subs", "indels", "length")

	def __init__(self, contig, start, strand, subs, indels, length):
		self.contig = contig
		self.start = start
		self.strand = strand
		self.subs = subs
		self.indels = indels
		self.length = length

	@property
	def edits(self):
		return self.subs + self.indels

	def __repr__(self):
		return "Copy(%s:%d%s subs=%d indels=%d)" % (self.contig, self.start, self.strand, self.subs, self.indels)


def planted_genome(seed=3, read_length=100, copies=PLANTED_COPIES):
	"""Two random contigs (pA 8 kbp, pB 6 kbp) carrying edited copies of one
	random probe. Returns (Genome, probe, [Copy])."""
	rng = Rng(seed)
	contigs = {"pA": list(rng.dna(8000)), "pB": list(rng.dna(6000))}
	probe = rng.dna(read_length)
	out = []
	for contig, start, strand, subs, ins, dels in copies:
		body = _edit(rng, probe, subs, ins, dels)
		placed = body if strand == "+" else revcomp(body)
		contigs[contig][start:start + len(placed)] = list(placed)
		out.append(Copy(contig, start, strand, subs, ins + dels, len(placed)))
	genome = Genome([(name, "".join(contigs[name])) for name in ("pA", "pB")])
	return genome, probe, out


def chimeras(genome, n, genomic_lengths, seed=1, read_length=100, contig="chr1"):
	"""Reads of ``read_length`` whose genomic part (``genomic_lengths``,
	cycled) is taken from ``contig`` and whose rest is random sequence that
	differs from the reference at every base next to the genomic part, so
	no alignment can extend into it. Half have the genomic part on the left.
	Returns [(name, seq, (contig, start, genomic_length, left))]."""
	rng = Rng(seed)
	ref = genome[contig]
	out = []
	i = 0
	while len(out) < n:
		L = genomic_lengths[i % len(genomic_lengths)]
		left = i % 2 == 0
		i += 1
		start = rng.randint(1000, len(ref) - 1000)
		part = ref[start:start + L].upper()
		if "N" in part or genome.hard_regions(contig, start - 200, start + L + 200):
			continue
		tail = read_length - L
		if left:
			# junk[j] never equals the reference base that would follow
			junk = "".join(rng.choice([b for b in "ACGT" if b != ref[start + L + j].upper()]) for j in range(tail))
			seq = part + junk
		else:
			junk = "".join(rng.choice([b for b in "ACGT" if b != ref[start - tail + j].upper()]) for j in range(tail))
			seq = junk + part
		out.append(("chim%d_%d" % (len(out), L), seq, (contig, start, L, left)))
	return out


def fastq_text(records, qual_char="I"):
	"""[(name, seq)] or [(name, seq, ...)] as FASTQ bytes with constant quality."""
	return "".join("@%s\n%s\n+\n%s\n" % (r[0], r[1], qual_char * len(r[1])) for r in records).encode("ascii")


def mapped(records):
	return [r for r in records if not r.is_unmapped]


def fill_secondary_seq(records):
	"""Secondary records (flag 256) carry SEQ '*'; give them the primary's
	sequence on their own strand so they can be checked against the
	reference. Returns the records (modified in place)."""
	primary = {}
	for r in records:
		if r.is_primary and r.seq != "*":
			primary[r.qname] = r
	for r in records:
		if r.seq == "*" and r.qname in primary:
			p = primary[r.qname]
			r.seq = p.seq if r.strand == p.strand else revcomp(p.seq)
	return records
