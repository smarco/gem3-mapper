# GEM-Mapper v3 (GEM3) test suite.
# Copyright (c) 2026 Jacob Schreiber <jacob.schreiber@umassmed.edu>
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of GEM-Mapper v3 (GEM3), distributed under the GNU General
# Public License version 3 or (at your option) any later version; see LICENSE.
"""Read simulators that record each read's true origin.

Every simulated read carries a :class:`Truth`: contig, 0-based leftmost
reference position, strand, and the true alignment (CIGAR, MD, NM in SAM
terms, i.e. on the forward reference strand). Edits are substitutions and
small insertions/deletions placed at least ``edge`` bases from either end of
the read and ``spacing`` bases from each other, so the true alignment is the
natural one; indels inside homopolymers can still be placed elsewhere by an
aligner with the same NM.

Paired reads: :func:`simulate_pe` with orientation FR, RF, FF or RR and a
normal insert-size distribution; fragments shorter than the read length read
through into adapter sequence (the truth CIGAR soft-clips the adapter part).

Bisulfite: :func:`simulate_bisulfite` converts unmethylated cytosines on the
original strand (C-to-T); in a directional library read 1 shows C-to-T and
read 2 G-to-A, in PBAT the reverse, and a non-directional library mixes both.
Methylation is a deterministic function of the site, so all reads covering a
cytosine agree.
"""

from .genome import revcomp
from .rng import Rng, hash_u64
from .sam import cigar_to_string, md_nm_from_alignment

ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"


class Truth(object):
	"""Where a read came from, in SAM conventions (forward reference strand).

	``pos`` is 0-based (SAM POS minus one). ``cigar`` is a list of (op, length)
	including soft clips for adapter bases; ``md`` and ``nm`` cover the
	aligned part. ``genomic_length`` counts read bases taken from the genome.
	"""

	__slots__ = ("contig", "pos", "strand", "cigar", "md", "nm", "ref_span",
		"genomic_length", "edits", "extra")

	def __init__(self, contig, pos, strand, cigar, md, nm, ref_span,
			genomic_length, edits, extra=None):
		self.contig = contig
		self.pos = pos
		self.strand = strand
		self.cigar = cigar
		self.md = md
		self.nm = nm
		self.ref_span = ref_span
		self.genomic_length = genomic_length
		self.edits = edits
		self.extra = extra if extra is not None else {}

	@property
	def cigar_string(self):
		return cigar_to_string(self.cigar)

	@property
	def end(self):
		return self.pos + self.ref_span

	def __repr__(self):
		return "Truth(%s:%d %s %s MD:%s NM:%d)" % (self.contig, self.pos, self.strand,
			self.cigar_string, self.md, self.nm)


class SimRead(object):
	"""A simulated read: ``name``, ``seq`` and ``qual`` as sequenced, plus ``truth``."""

	__slots__ = ("name", "seq", "qual", "truth")

	def __init__(self, name, seq, qual, truth):
		self.name = name
		self.seq = seq
		self.qual = qual
		self.truth = truth

	def __repr__(self):
		return "SimRead(%s, %d bp, %r)" % (self.name, len(self.seq), self.truth)


class SimPair(object):
	"""Two mates and the fragment they were read from (``frag_start``, ``frag_end``, 0-based half-open)."""

	__slots__ = ("r1", "r2", "contig", "frag_start", "frag_end", "orientation")

	def __init__(self, r1, r2, contig, frag_start, frag_end, orientation):
		self.r1 = r1
		self.r2 = r2
		self.contig = contig
		self.frag_start = frag_start
		self.frag_end = frag_end
		self.orientation = orientation

	@property
	def fragment_length(self):
		return self.frag_end - self.frag_start


def _qualities(rng, length, lo=20, hi=40):
	return "".join(chr(33 + rng.randint(lo, hi)) for _ in range(length))


def _place_edits(rng, length, subs, ins, dels, edge, spacing, max_indel):
	"""Pick read offsets for the edits; returns sorted [(offset, kind, size)]."""
	kinds = ["S"] * subs + ["I"] * ins + ["D"] * dels
	if not kinds:
		return []
	for _ in range(200):
		chosen = []
		ok = True
		for kind in kinds:
			size = 1 if kind == "S" else rng.randint(1, max_indel)
			for _ in range(200):
				hi = length - edge - (size if kind == "I" else 0)
				if hi < edge:
					ok = False
					break
				off = rng.randint(edge, hi)
				if all(off >= o + (s if k == "I" else 1) + spacing or off + (size if kind == "I" else 1) + spacing <= o
						for o, k, s in chosen):
					chosen.append((off, kind, size))
					break
			else:
				ok = False
			if not ok:
				break
		if ok:
			return sorted(chosen)
	raise ValueError("cannot place %d edits in a %d bp read with edge=%d spacing=%d"
		% (len(kinds), length, edge, spacing))


def _apply_edits(rng, ref, ref_start, length, edits):
	"""Build ``length`` read bases from ``ref[ref_start:]`` applying ``edits``.

	Returns (read_bases, cigar_ops, ref_span). Read bases are on the forward
	reference strand and upper-case.
	"""
	read = []
	ops = []

	def push(op, n):
		if ops and ops[-1][0] == op:
			ops[-1] = (op, ops[-1][1] + n)
		else:
			ops.append((op, n))

	by_offset = dict((off, (kind, size)) for off, kind, size in edits)
	r = ref_start
	q = 0
	while q < length:
		if q in by_offset:
			kind, size = by_offset[q]
			if kind == "S":
				orig = ref[r].upper()
				read.append(rng.choice([b for b in "ACGT" if b != orig]))
				push("M", 1)
				r += 1
				q += 1
				continue
			if kind == "I":
				prev = read[-1] if read else ""
				nxt = ref[r].upper()
				bases = []
				for _ in range(size):
					# Inserted bases differ from both neighbours so the
					# insertion cannot slide left or right.
					bases.append(rng.choice([b for b in "ACGT" if b != prev and b != nxt]))
					prev = bases[-1]
				read.extend(bases)
				push("I", size)
				q += size
				continue
			if kind == "D":
				push("D", size)
				r += size
				# fall through to the match at this read offset
		read.append(ref[r].upper())
		push("M", 1)
		r += 1
		q += 1
	return "".join(read), ops, r - ref_start


def _edit_counts(rng, length, sub_rate, ins_rate, del_rate, subs, ins, dels):
	def count(fixed, rate):
		if fixed is not None:
			return fixed
		return sum(1 for _ in range(length) if rng.bernoulli(rate)) if rate > 0 else 0
	return count(subs, sub_rate), count(ins, ins_rate), count(dels, del_rate)


def _span_for(length, edits):
	return length - sum(s for _, k, s in edits if k == "I") + sum(s for _, k, s in edits if k == "D")


def _pick_locus(rng, genome, span, unique_only, margin, contigs):
	names = contigs or [n for n in genome.names if len(genome[n]) >= span]
	weights = [len(genome[n]) for n in names]
	total = float(sum(weights))
	for _ in range(10000):
		x = rng.random() * total
		for name, w in zip(names, weights):
			if x < w:
				break
			x -= w
		L = len(genome[name])
		if L < span:
			continue
		pos = rng.randint(0, L - span)
		seg = genome[name][pos:pos + span]
		if "N" in seg or "n" in seg:
			continue
		if unique_only and genome.hard_regions(name, pos, pos + span, margin=margin):
			continue
		return name, pos
	raise RuntimeError("no locus of span %d found" % span)


def _make_truth(genome, contig, pos, strand, read_fwd, ops, span, edits, clip_left=0, clip_right=0, extra=None):
	ref = genome[contig]
	md, nm = md_nm_from_alignment(ref, pos, ops, read_fwd)
	cigar = list(ops)
	if clip_left:
		cigar = [("S", clip_left)] + cigar
	if clip_right:
		cigar = cigar + [("S", clip_right)]
	return Truth(contig, pos, strand, cigar, md, nm, span,
		sum(n for op, n in ops if op in "MI"), edits, extra)


def simulate_se(genome, n, length=100, seed=1, sub_rate=0.0, ins_rate=0.0, del_rate=0.0,
		subs=None, ins=None, dels=None, edge=8, spacing=8, max_indel=3, unique_only=True,
		margin=None, strand=None, contigs=None, name_prefix="r", qual_range=(20, 40)):
	"""Single-end reads of ``length`` bp.

	Edit counts are fixed per read when ``subs``/``ins``/``dels`` are given,
	otherwise drawn per base from the rates. ``unique_only`` keeps reads
	(plus ``margin``, default one read length) away from repeats, tandems,
	homopolymers, N runs and short contigs; reads never cover an N.
	``strand`` forces '+' or '-'; by default each read picks one at random.
	"""
	rng = Rng(hash_u64(seed, "se"))
	margin = length if margin is None else margin
	reads = []
	for i in range(n):
		ns, ni, nd = _edit_counts(rng, length, sub_rate, ins_rate, del_rate, subs, ins, dels)
		edits = _place_edits(rng, length, ns, ni, nd, edge, spacing, max_indel)
		span = _span_for(length, edits)
		contig, pos = _pick_locus(rng, genome, span, unique_only, margin, contigs)
		read_fwd, ops, span2 = _apply_edits(rng, genome[contig], pos, length, edits)
		assert span2 == span
		st = strand or ("+" if rng.bernoulli(0.5) else "-")
		qual_fwd = _qualities(rng, length, *qual_range)
		truth = _make_truth(genome, contig, pos, st, read_fwd, ops, span, edits)
		if st == "+":
			seq, qual = read_fwd, qual_fwd
		else:
			seq, qual = revcomp(read_fwd), qual_fwd[::-1]
		reads.append(SimRead("%s%d" % (name_prefix, i), seq, qual, truth))
	return reads


def _mate(rng, genome, contig, start, end, at_left, forward, length, edits_spec, edge, spacing,
		max_indel, adapter, qual_range):
	"""One mate of a fragment [start, end): reads ``length`` bases from the
	fragment end it starts at, on the strand given by ``forward``."""
	frag_len = end - start
	genomic = min(length, frag_len)
	ns, ni, nd = edits_spec(genomic)
	edits = _place_edits(rng, genomic, ns, ni, nd, edge, spacing, max_indel) if genomic > 2 * edge else []
	span = _span_for(genomic, edits)
	# Deletions lengthen the reference span; keep it inside the fragment.
	if span > frag_len:
		edits = []
		span = genomic
	pos = start if at_left else end - span
	read_fwd, ops, span2 = _apply_edits(rng, genome[contig], pos, genomic, edits)
	assert span2 == span
	tail = length - genomic
	adapter_bases = ""
	if tail > 0:
		adapter_bases = (adapter * (tail // len(adapter) + 1))[:tail]
	qual = _qualities(rng, length, *qual_range)
	if forward:
		seq = read_fwd + adapter_bases
		truth = _make_truth(genome, contig, pos, "+", read_fwd, ops, span, edits,
			clip_right=tail, extra={"adapter": tail})
	else:
		seq = revcomp(read_fwd) + adapter_bases
		truth = _make_truth(genome, contig, pos, "-", read_fwd, ops, span, edits,
			clip_left=tail, extra={"adapter": tail})
	return seq, qual, truth


def simulate_pe(genome, n, length=100, insert_mean=300, insert_sd=30, insert_min=None,
		insert_max=None, orientation="FR", seed=1, sub_rate=0.0, ins_rate=0.0, del_rate=0.0,
		subs=None, ins=None, dels=None, edge=8, spacing=8, max_indel=3, unique_only=True,
		margin=None, contigs=None, name_prefix="p", qual_range=(20, 40),
		adapters=(ADAPTER_R1, ADAPTER_R2)):
	"""Paired reads; returns a list of :class:`SimPair`.

	``orientation`` describes the mates as they map, leftmost mate first:
	``FR`` (Illumina paired-end), ``RF`` (mate-pair), ``FF`` or ``RR``. Half of
	the fragments come from each strand, so read 1 is the leftmost mate for
	half of the pairs. Insert sizes are normal(``insert_mean``,
	``insert_sd``) clipped to [``insert_min``, ``insert_max``] (default
	[1, mean + 6 sd]); a fragment shorter than ``length`` reads through into
	the adapter given for that mate.
	"""
	if orientation not in ("FR", "RF", "FF", "RR"):
		raise ValueError("orientation must be FR, RF, FF or RR")
	rng = Rng(hash_u64(seed, "pe", orientation))
	insert_min = 1 if insert_min is None else insert_min
	insert_max = int(insert_mean + 6 * insert_sd) if insert_max is None else insert_max
	margin = length if margin is None else margin
	pairs = []

	def edits_spec(genomic):
		return _edit_counts(rng, genomic, sub_rate, ins_rate, del_rate, subs, ins, dels)

	left_fwd = orientation[0] == "F"
	right_fwd = orientation[1] == "F"
	for i in range(n):
		frag = int(round(rng.normal(insert_mean, insert_sd)))
		frag = max(insert_min, min(insert_max, frag))
		contig, start = _pick_locus(rng, genome, max(frag, length), unique_only, margin, contigs)
		end = start + frag
		# Left mate reads from the fragment start, right mate from its end.
		r1_left = rng.bernoulli(0.5)
		left = _mate(rng, genome, contig, start, end, True, left_fwd, length, edits_spec, edge,
			spacing, max_indel, adapters[0 if r1_left else 1], qual_range)
		right = _mate(rng, genome, contig, start, end, False, right_fwd, length, edits_spec, edge,
			spacing, max_indel, adapters[1 if r1_left else 0], qual_range)
		m1, m2 = (left, right) if r1_left else (right, left)
		name = "%s%d" % (name_prefix, i)
		pairs.append(SimPair(SimRead(name, m1[0], m1[1], m1[2]), SimRead(name, m2[0], m2[1], m2[2]),
			contig, start, end, orientation))
	return pairs


def methylated(contig, pos, strand_of_c, context, seed, cpg_rate, chh_rate):
	"""Whether the cytosine at ``pos`` on ``strand_of_c`` is methylated.

	A deterministic function of the site and ``seed``: every read covering the
	site sees the same state. ``context`` is 'CG' or 'CH'.
	"""
	rate = cpg_rate if context == "CG" else chh_rate
	if rate <= 0.0:
		return False
	if rate >= 1.0:
		return True
	return (hash_u64(seed, contig, pos, 1 if strand_of_c == "+" else 2) >> 11) * (1.0 / (1 << 53)) < rate


def bisulfite_convert(genome, contig, start, end, strand, seed=1, cpg_rate=0.8, chh_rate=0.01):
	"""Bisulfite-convert the original ``strand`` of ``genome[contig][start:end]``.

	Returns (converted, calls): ``converted`` is the original strand read
	5'->3' after unmethylated C became T (so for strand '-' it is the reverse
	complement of the reference, converted). ``calls`` maps each reference
	position of a cytosine on that strand to True (methylated, kept C) or
	False (converted).
	"""
	ref = genome[contig]
	calls = {}
	if strand == "+":
		out = []
		for p in range(start, end):
			b = ref[p].upper()
			if b == "C":
				ctx = "CG" if p + 1 < len(ref) and ref[p + 1].upper() == "G" else "CH"
				m = methylated(contig, p, "+", ctx, seed, cpg_rate, chh_rate)
				calls[p] = m
				out.append("C" if m else "T")
			else:
				out.append(b)
		return "".join(out), calls
	out = []
	for p in range(end - 1, start - 1, -1):
		b = ref[p].upper()
		if b == "G":
			ctx = "CG" if p - 1 >= 0 and ref[p - 1].upper() == "C" else "CH"
			m = methylated(contig, p, "-", ctx, seed, cpg_rate, chh_rate)
			calls[p] = m
			out.append("C" if m else "T")
		else:
			out.append(revcomp(b))
	return "".join(out), calls


def simulate_bisulfite(genome, n, length=100, paired=True, library="directional",
		insert_mean=250, insert_sd=25, seed=1, cpg_rate=0.8, chh_rate=0.01,
		unique_only=True, margin=None, contigs=None, name_prefix="b", qual_range=(20, 40)):
	"""Bisulfite reads without sequencing errors.

	``library``: ``directional`` (read 1 from the converted original top or
	bottom strand, so read 1 shows C-to-T and read 2 G-to-A), ``pbat``
	(read 1 from the strand complementary to the converted one: read 1 G-to-A,
	read 2 C-to-T) or ``non-directional`` (each fragment picks one of the two
	at random). Returns :class:`SimPair` objects (``paired=True``) or
	:class:`SimRead` objects (read 1 only). Each truth's ``extra`` holds
	``origin`` (OT, OB, CTOT or CTOB), ``conversion`` ('CT' or 'GA', as seen in
	the read) and ``calls`` (reference position -> methylated).
	"""
	if library not in ("directional", "pbat", "non-directional"):
		raise ValueError("library must be directional, pbat or non-directional")
	rng = Rng(hash_u64(seed, "bs", library))
	margin = length if margin is None else margin
	out = []
	for i in range(n):
		frag = max(length if not paired else 1, int(round(rng.normal(insert_mean, insert_sd))))
		if not paired:
			frag = max(frag, length)
		contig, start = _pick_locus(rng, genome, max(frag, length), unique_only, margin, contigs)
		end = start + frag
		strand = "+" if rng.bernoulli(0.5) else "-"
		conv, calls = bisulfite_convert(genome, contig, start, end, strand, seed, cpg_rate, chh_rate)
		pbat = library == "pbat" or (library == "non-directional" and rng.bernoulli(0.5))
		# The converted original strand, 5'->3', is ``conv``; its complement
		# (CTOT/CTOB) is revcomp(conv).
		g1 = min(length, frag)
		top_name = "OT" if strand == "+" else "OB"
		comp_name = "CTOT" if strand == "+" else "CTOB"
		reads_conv = conv[:g1]                      # from the converted strand
		reads_comp = revcomp(conv)[:g1]             # from its complement
		ra = (reads_conv, top_name, "CT")
		rb = (reads_comp, comp_name, "GA")
		first, second = (rb, ra) if pbat else (ra, rb)
		name = "%s%d" % (name_prefix, i)
		mates = []
		for mate_idx, (bases, origin, conversion) in enumerate((first, second)):
			if mate_idx == 1 and not paired:
				break
			# Which reference interval and strand does this mate align to?
			from_conv = origin in ("OT", "OB")
			if strand == "+":
				fwd = from_conv
				pos = start if from_conv else end - g1
			else:
				fwd = not from_conv
				pos = end - g1 if from_conv else start
			# Forward-strand read sequence (as SAM SEQ would show it).
			read_fwd = bases if fwd else revcomp(bases)
			ops = [("M", g1)]
			tail = length - g1
			adapter = ADAPTER_R1 if mate_idx == 0 else ADAPTER_R2
			adapter_bases = (adapter * (tail // len(adapter) + 1))[:tail] if tail > 0 else ""
			seq = bases + adapter_bases
			ref_seg = genome[contig][pos:pos + g1].upper()
			# Bisulfite mismatches are expected; MD/NM here are against the
			# unconverted reference (what SAM MD/NM describe).
			md, nm = md_nm_from_alignment(genome[contig], pos, ops, read_fwd)
			cigar = list(ops)
			if tail:
				cigar = cigar + [("S", tail)] if fwd else [("S", tail)] + cigar
			mate_calls = dict((p, m) for p, m in calls.items() if pos <= p < pos + g1)
			truth = Truth(contig, pos, "+" if fwd else "-", cigar, md, nm, g1, g1, [],
				{"origin": origin, "conversion": conversion, "calls": mate_calls,
				 "adapter": tail, "ref": ref_seg})
			mates.append(SimRead(name, seq, _qualities(rng, length, *qual_range), truth))
		if paired:
			out.append(SimPair(mates[0], mates[1], contig, start, end, "FR"))
		else:
			out.append(mates[0])
	return out
