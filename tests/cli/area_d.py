# GEM-Mapper v3 (GEM3) test suite.
# Copyright (c) 2026 Jacob Schreiber <jacob.schreiber@umassmed.edu>
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of GEM-Mapper v3 (GEM3), distributed under the GNU General
# Public License version 3 or (at your option) any later version; see LICENSE.
"""Helpers for the area D tests (SAM/MAP output, report file, CLI, threads).

Expected values are derived here from the SAM specification (SAMv1 1.4),
from the reference sequence and from the reads that went in, never from
gem-mapper's own writers:

- ``se_record_problems`` / ``pair_problems``: FLAG bits, mate fields and
  TLEN as SAMv1 defines them, unmapped-record placement as SAMv1
  recommends, SEQ/QUAL orientation, NM/MD and the gap-affine AS recomputed
  from the reference.
- ``map_cigar``: the GEM MAP-format CIGAR of an alignment, built from the
  SAM CIGAR, the read and the reference.
- ``expected_report_se`` / ``expected_report_pe``: the ``--report-file``
  JSON recomputed from the input reads and the SAM records.

Datasets (``Dataset``) write their inputs once per pytest session into the
session's gem3 work directory and memoize gem-mapper runs, so the area D
modules share runs instead of repeating them.
"""

import hashlib
import json
import os

from gem3harness import (Rng, check_record_against_reference, default_threads, parse_sam, revcomp, simulate_bisulfite, simulate_pe, simulate_se,
	synthetic_genome, write_reads, write_split)
from gem3harness.sam import (MATE_REVERSE, MATE_UNMAPPED, PAIRED, PROPER_PAIR, READ1, READ2, REVERSE,
	SECONDARY, SUPPLEMENTARY, UNMAPPED)

# Default --gap-affine-penalties (match, mismatch, gap open, gap extension).
MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND = 1, 4, 6, 1

# Contig names gem-mapper's report classifies as controls (report_stats.h).
SEQUENCING_CONTROL = "NC_001422.1"
UNDERCONVERSION_CONTROL = "NC_001416.1"
OVERCONVERSION_CONTROL = "NC_001604.1"


class Read(object):
	"""A read for ``write_reads``: ``name``, ``seq``, ``qual``."""

	__slots__ = ("name", "seq", "qual", "truth")

	def __init__(self, name, seq, qual, truth=None):
		self.name = name
		self.seq = seq
		self.qual = qual
		self.truth = truth


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

def _shuffle(rng, items):
	items = list(items)
	for i in range(len(items) - 1, 0, -1):
		j = rng.randint(0, i)
		items[i], items[j] = items[j], items[i]
	return items


def _qual(rng, n):
	return "".join(chr(33 + rng.randint(20, 40)) for _ in range(n))


def repeat_reads(genome, n, length, seed, prefix):
	"""Reads cut from the exact copy (copy 0) of each repeat family, half of
	them reverse-complemented: each has several equally good placements."""
	rng = Rng(seed)
	copies = [f for f in genome.features if f.kind == "repeat" and f.info["copy"] == 0]
	out = []
	for i in range(n):
		f = copies[i % len(copies)]
		start = rng.randint(f.start, f.end - length)
		seq = genome[f.contig][start:start + length].upper()
		if i % 2:
			seq = revcomp(seq)
		out.append(Read("%s%d" % (prefix, i), seq, _qual(rng, length)))
	return out


def se_mix(genome, seed=101):
	"""Single-end reads covering the writer's cases: unique reads with
	substitutions and indels (both strands), multi-mapping repeat reads,
	random (unmappable) reads, reads with N bases and 60 bp reads."""
	rng = Rng(seed)
	unique = simulate_se(genome, 110, 100, seed=seed, sub_rate=0.02, ins_rate=0.003,
		del_rate=0.003, name_prefix="u")
	short = simulate_se(genome, 10, 60, seed=seed + 1, subs=1, name_prefix="s")
	withn = simulate_se(genome, 10, 100, seed=seed + 2, subs=1, name_prefix="n")
	ns = []
	for r in withn:
		seq = list(r.seq)
		for k in (17, 50, 81):
			seq[k] = "N"
		ns.append(Read(r.name, "".join(seq), r.qual, r.truth))
	rep = repeat_reads(genome, 24, 100, seed + 3, "m")
	rnd = [Read("x%d" % i, rng.dna(100), _qual(rng, 100)) for i in range(12)]
	return _shuffle(rng, list(unique) + list(short) + ns + rep + rnd)


def pe_mix(genome, seed=202):
	"""Pairs covering the paired writer's cases: proper FR pairs, pairs
	from a repeat copy (alternative pairings), pairs with one random mate
	(unpaired output, one end unmapped) and pairs with both mates random."""
	rng = Rng(seed)
	pairs = simulate_pe(genome, 100, 100, insert_mean=320, insert_sd=40, insert_min=220,
		seed=seed, sub_rate=0.01, ins_rate=0.002, del_rate=0.002, name_prefix="p")
	copies = [f for f in genome.features if f.kind == "repeat" and f.info["copy"] == 0
		and f.end - f.start >= 500]
	reps = []
	for i in range(12):
		f = copies[i % len(copies)]
		frag = rng.randint(250, min(400, f.end - f.start))
		start = rng.randint(f.start, f.end - frag)
		left = genome[f.contig][start:start + 100].upper()
		right = revcomp(genome[f.contig][start + frag - 100:start + frag].upper())
		m1, m2 = (left, right) if i % 2 else (right, left)
		reps.append((Read("q%d" % i, m1, _qual(rng, 100)), Read("q%d" % i, m2, _qual(rng, 100))))
	out = [(p.r1, p.r2) for p in pairs] + reps
	result = []
	for i, (a, b) in enumerate(out):
		kind = i % 8
		if kind == 3:
			b = Read(b.name, rng.dna(100), b.qual)
		elif kind == 5:
			a = Read(a.name, rng.dna(100), a.qual)
		elif kind == 7 and i % 16 == 7:
			a = Read(a.name, rng.dna(100), a.qual)
			b = Read(b.name, rng.dna(100), b.qual)
		result.append((a, b))
	return _shuffle(rng, result)


def pe_tie_pairs(genome, seed=77, n=150):
	"""Pairs from fragments overlapping repeats, tandem repeats and
	homopolymers (placements that tie), mixed with unique pairs. Their output
	depends on the template-length estimate, which gem-mapper restarts at
	every input buffer (README, Determinism)."""
	rng = Rng(seed)
	feats = [f for f in genome.features if f.kind in ("tandem", "homopolymer", "repeat")]
	out = []
	for i in range(n):
		f = feats[i % len(feats)]
		frag = rng.randint(200, 450)
		start = rng.randint(max(0, f.start - frag + 60), min(len(genome[f.contig]) - frag, f.end - 60))
		s = genome[f.contig][start:start + frag].upper()
		if "N" in s:
			continue
		a, b = s[:100], revcomp(s[-100:])
		if i % 2:
			a, b = b, a
		out.append((Read("t%d" % i, a, _qual(rng, 100)), Read("t%d" % i, b, _qual(rng, 100))))
	uniq = simulate_pe(genome, 200, 100, seed=seed + 1, insert_mean=300, insert_sd=60, name_prefix="w")
	return _shuffle(rng, out + [(p.r1, p.r2) for p in uniq])


def control_genome(seed=31):
	"""A small genome whose contigs carry the report's control names."""
	return synthetic_genome(seed, contigs=(("chrA", 40000), (SEQUENCING_CONTROL, 5000),
		(UNDERCONVERSION_CONTROL, 6000), (OVERCONVERSION_CONTROL, 6000)),
		repeat_families=0, tandems=0, homopolymers=0, n_runs=0, softmask=0, edge_n=False)


# ---------------------------------------------------------------------------
# Datasets with memoized runs
# ---------------------------------------------------------------------------

_DATASETS = {}


class Dataset(object):
	"""Inputs written once per session plus memoized gem-mapper runs.

	``reads`` (single-end) or ``pairs`` (list of (read1, read2)) are the
	inputs in file order; ``paths`` are the FASTQ files. ``run(args,
	threads)`` returns a RunResult (cached on the arguments)."""

	def __init__(self, gem3, index, name, reads=None, pairs=None):
		self.gem3 = gem3
		self.index = index
		self.reads = reads
		self.pairs = pairs
		d = os.path.join(gem3.work_dir, "area_d")
		if not os.path.isdir(d):
			os.makedirs(d)
		self.dir = d
		if reads is not None:
			self.paths = [write_reads(os.path.join(d, name + ".fq"), reads)]
		else:
			p1 = os.path.join(d, name + "_1.fq")
			p2 = os.path.join(d, name + "_2.fq")
			write_split(p1, p2, [_Pair(a, b) for a, b in pairs])
			self.paths = [p1, p2]
		self._runs = {}

	@property
	def paired(self):
		return self.pairs is not None

	def run(self, args=(), threads=None, inputs=None, check=True, **kw):
		threads = default_threads() if threads is None else threads
		inputs = self.paths if inputs is None else inputs
		key = (tuple(args), threads, tuple(inputs), tuple(sorted(kw.items())), check)
		if key not in self._runs:
			self._runs[key] = self.gem3.run_map(self.index, args=list(args), inputs=inputs,
				threads=threads, check=check, **kw)
		return self._runs[key]


class _Pair(object):
	__slots__ = ("r1", "r2")

	def __init__(self, r1, r2):
		self.r1 = r1
		self.r2 = r2


def dataset(gem3, index, name, make):
	"""The session's Dataset ``name`` for ``index`` (``make()`` -> ('se', reads) or ('pe', pairs))."""
	key = (gem3.work_dir, index, name)
	if key not in _DATASETS:
		kind, items = make()
		if kind == "se":
			_DATASETS[key] = Dataset(gem3, index, name, reads=items)
		else:
			_DATASETS[key] = Dataset(gem3, index, name, pairs=items)
	return _DATASETS[key]


def se_dataset(gem3, genome, index):
	return dataset(gem3, index, "se_mix", lambda: ("se", se_mix(genome)))


def pe_dataset(gem3, genome, index):
	return dataset(gem3, index, "pe_mix", lambda: ("pe", pe_mix(genome)))


def bs_dataset(gem3, genome, index, name="bs_pairs", seed=303, n=120):
	def make():
		bs = simulate_bisulfite(genome, n, 100, library="directional", seed=seed,
			insert_mean=260, insert_sd=30)
		return "pe", [(p.r1, p.r2) for p in bs]
	return dataset(gem3, index, name, make)


# ---------------------------------------------------------------------------
# SAM checks (SAMv1)
# ---------------------------------------------------------------------------

def gap_affine_score(ref, pos, ops, seq):
	"""Gap-affine score with the default penalties; a k-base gap costs O + k*X."""
	score = 0
	r, q = pos, 0
	for op, n in ops:
		if op == "M":
			for k in range(n):
				same = ref[r + k].upper() == seq[q + k].upper() and seq[q + k].upper() != "N"
				score += MATCH if same else -MISMATCH
			r += n
			q += n
		elif op in ("I", "D"):
			score -= GAP_OPEN + GAP_EXTEND * n
			if op == "I":
				q += n
			else:
				r += n
		elif op == "S":
			q += n
	return score


def forward_seq(read_seq, strand):
	return read_seq if strand == "+" else revcomp(read_seq)


def mismatches(ref, pos, ops, seq, conversion=None):
	"""Substituted bases of an alignment. With ``conversion`` 'C' (C-to-T)
	or 'G' (G-to-A) both sequences are converted first, as gem-mapper
	aligns bisulfite reads to the converted reference."""
	conv = {"C": ("C", "T"), "G": ("G", "A")}.get(conversion)
	n = 0
	r, q = pos, 0
	for op, k in ops:
		if op == "M":
			for i in range(k):
				rb, qb = ref[r + i].upper(), seq[q + i].upper()
				if conv:
					rb = conv[1] if rb == conv[0] else rb
					qb = conv[1] if qb == conv[0] else qb
				if rb != qb or rb == "N":
					n += 1
			r += k
			q += k
		elif op in ("I", "S"):
			q += k
		elif op == "D":
			r += k
	return n


def alignment_problems(rec, seq, genome):
	"""CIGAR/NM/MD problems of a mapped record whose forward SEQ is ``seq``."""
	clone = _Clone(rec, seq)
	return check_record_against_reference(clone, genome[rec.rname])


class _Clone(object):
	"""A SamRecord stand-in with SEQ replaced (secondary records print '*')."""

	def __init__(self, rec, seq):
		self.cigar_ops = rec.cigar_ops
		self.cigar = rec.cigar
		self.seq = seq
		self.tags = rec.tags
		self.ref_start = rec.ref_start
		self.ref_end = rec.ref_end


def unmapped_problems(rec, read):
	p = []
	if (rec.rname, rec.pos, rec.mapq, rec.cigar) != ("*", 0, 0, "*") and not (rec.flag & PAIRED):
		p.append("unmapped SE record has RNAME/POS/MAPQ/CIGAR %s/%d/%d/%s" % (rec.rname, rec.pos, rec.mapq, rec.cigar))
	if rec.seq != read.seq or rec.qual != read.qual:
		p.append("unmapped record SEQ/QUAL not the read as given")
	for t in ("NM", "MD", "AS", "XA", "XS"):
		if t in rec.tags:
			p.append("unmapped record carries %s" % t)
	return p


def se_record_problems(rec, read, genome, primary=True):
	"""SAMv1 checks of a single-end record against its input read."""
	p = []
	allowed = UNMAPPED | REVERSE | SECONDARY
	if rec.flag & ~allowed:
		p.append("flag %d has bits outside 0x4|0x10|0x100" % rec.flag)
	if rec.qname != read.name:
		p.append("QNAME %s != %s" % (rec.qname, read.name))
	if (rec.rnext, rec.pnext, rec.tlen) != ("*", 0, 0):
		p.append("single-end RNEXT/PNEXT/TLEN %s/%d/%d" % (rec.rnext, rec.pnext, rec.tlen))
	if rec.flag & UNMAPPED:
		if rec.flag & (REVERSE | SECONDARY):
			p.append("unmapped record with 0x10/0x100")
		return p + unmapped_problems(rec, read)
	if bool(rec.flag & SECONDARY) == primary:
		p.append("0x100 is %s on a %s record" % (bool(rec.flag & SECONDARY), "primary" if primary else "secondary"))
	seq = forward_seq(read.seq, rec.strand)
	if primary:
		if rec.seq != seq:
			p.append("SEQ is not the read on the forward strand")
		if rec.qual != (read.qual if rec.strand == "+" else read.qual[::-1]):
			p.append("QUAL is not the read's qualities on the forward strand")
	elif (rec.seq, rec.qual) != ("*", "*"):
		p.append("secondary record prints SEQ/QUAL")
	p += alignment_problems(rec, seq, genome)
	if not p:
		score = gap_affine_score(genome[rec.rname], rec.ref_start, rec.cigar_ops, seq)
		if rec.tags.get("AS") != score:
			p.append("AS:%s but the gap-affine score is %d" % (rec.tags.get("AS"), score))
	return p


def template_length(a, b):
	"""SAMv1 TLEN of two mapped records on one contig: from the leftmost
	mapped base to the rightmost, signed + on the leftmost record."""
	left = min(a.ref_start, b.ref_start)
	right = max(a.ref_end, b.ref_end)
	return right - left


def pair_problems(r1, r2, read1, read2, genome):
	"""SAMv1 checks of a primary pair (records of end 1 and end 2)."""
	p = []
	for rec, read, mate, first in ((r1, read1, r2, True), (r2, read2, r1, False)):
		tag = "end%d" % (1 if first else 2)
		if rec.qname != read.name:
			p.append("%s QNAME %s != %s" % (tag, rec.qname, read.name))
		if not rec.flag & PAIRED:
			p.append("%s lacks 0x1" % tag)
		if bool(rec.flag & READ1) != first or bool(rec.flag & READ2) == first:
			p.append("%s has wrong 0x40/0x80 (flag %d)" % (tag, rec.flag))
		if rec.flag & (SECONDARY | SUPPLEMENTARY):
			p.append("%s primary record has 0x100/0x800" % tag)
		if bool(rec.flag & MATE_UNMAPPED) != bool(mate.flag & UNMAPPED):
			p.append("%s 0x8 disagrees with the mate's 0x4" % tag)
		if not mate.flag & UNMAPPED and bool(rec.flag & MATE_REVERSE) != bool(mate.flag & REVERSE):
			p.append("%s 0x20 disagrees with the mate's 0x10" % tag)
		if rec.flag & PROPER_PAIR and (rec.flag & UNMAPPED or mate.flag & UNMAPPED):
			p.append("%s has 0x2 with an end unmapped" % tag)
		if bool(rec.flag & PROPER_PAIR) != bool(mate.flag & PROPER_PAIR):
			p.append("%s 0x2 differs between the mates" % tag)
		if rec.flag & UNMAPPED:
			if rec.cigar != "*" or rec.mapq != 0:
				p.append("%s unmapped with CIGAR/MAPQ" % tag)
			if rec.seq != read.seq or rec.qual != read.qual:
				p.append("%s unmapped record SEQ/QUAL not the read as given" % tag)
			if mate.flag & UNMAPPED:
				if (rec.rname, rec.pos, rec.rnext, rec.pnext, rec.tlen) != ("*", 0, "*", 0, 0):
					p.append("%s both unmapped but placed" % tag)
			else:
				# SAMv1 recommended practice: an unmapped mate takes the mapped
				# mate's RNAME and POS.
				if (rec.rname, rec.pos) != (mate.rname, mate.pos):
					p.append("%s unmapped not placed at its mate" % tag)
				if rec.rnext not in ("=", mate.rname) or rec.pnext != mate.pos or rec.tlen != 0:
					p.append("%s unmapped: RNEXT/PNEXT/TLEN %s/%d/%d" % (tag, rec.rnext, rec.pnext, rec.tlen))
			continue
		seq = forward_seq(read.seq, rec.strand)
		if rec.seq != seq or rec.qual != (read.qual if rec.strand == "+" else read.qual[::-1]):
			p.append("%s SEQ/QUAL not on the forward strand" % tag)
		probs = alignment_problems(rec, seq, genome)
		p += ["%s %s" % (tag, x) for x in probs]
		if not probs:
			score = gap_affine_score(genome[rec.rname], rec.ref_start, rec.cigar_ops, seq)
			if rec.tags.get("AS") != score:
				p.append("%s AS:%s but gap-affine score %d" % (tag, rec.tags.get("AS"), score))
		if mate.flag & UNMAPPED:
			if (rec.rnext, rec.pnext, rec.tlen) != ("=", rec.pos, 0):
				p.append("%s with mate unmapped: RNEXT/PNEXT/TLEN %s/%d/%d" % (tag, rec.rnext, rec.pnext, rec.tlen))
			continue
		rnext = "=" if mate.rname == rec.rname else mate.rname
		if (rec.rnext, rec.pnext) != (rnext, mate.pos):
			p.append("%s RNEXT/PNEXT %s/%d, mate at %s/%d" % (tag, rec.rnext, rec.pnext, mate.rname, mate.pos))
		if rec.flag & PROPER_PAIR:
			if mate.rname != rec.rname:
				p.append("%s proper pair across contigs" % tag)
			else:
				tl = template_length(rec, mate)
				if abs(rec.tlen) != tl:
					p.append("%s |TLEN| %d, SAMv1 span %d" % (tag, abs(rec.tlen), tl))
				if rec.pos < mate.pos and rec.tlen <= 0 or rec.pos > mate.pos and rec.tlen >= 0:
					p.append("%s TLEN sign %d (pos %d, mate %d)" % (tag, rec.tlen, rec.pos, mate.pos))
			if rec.mapq != mate.mapq:
				p.append("%s proper pair MAPQs differ (template MAPQ)" % tag)
		elif rec.tlen != 0:
			p.append("%s not a proper pair but TLEN %d" % (tag, rec.tlen))
	return p


def primaries_by_name(sam):
	out = {}
	for rec in sam.records:
		if rec.is_primary:
			out.setdefault(rec.qname, []).append(rec)
	return out


def xa_entries(rec):
	"""XA:Z entries as (rname, strand, pos, cigar, distance)."""
	out = []
	if "XA" not in rec.tags:
		return out
	body = rec.tags["XA"]
	assert body.endswith(";"), body
	for e in body[:-1].split(";"):
		rname, spos, cigar, dist = e.split(",")
		assert spos[0] in "+-", e
		out.append((rname, spos[0], int(spos[1:]), cigar, int(dist)))
	return out


def soft_clipped(ops):
	return sum(n for op, n in ops if op == "S")


# ---------------------------------------------------------------------------
# MAP format
# ---------------------------------------------------------------------------

def map_cigar(rname, strand, pos, ops, read_seq, genome):
	"""GEM MAP CIGAR (format 2) of an alignment, in read orientation.

	Matches print as run lengths, a substitution as the reference base (its
	complement on the reverse strand), a SAM deletion as '>k+', a SAM
	insertion as '>k-' and a soft clip as '(k)'. ``pos`` is 1-based;
	``read_seq`` is the read as sequenced."""
	ref = genome[rname]
	seq = forward_seq(read_seq, strand)
	elems = []
	r, q = pos - 1, 0
	for op, n in ops:
		if op == "M":
			for k in range(n):
				rb, qb = ref[r + k].upper(), seq[q + k].upper()
				if rb == qb and rb != "N":
					if elems and isinstance(elems[-1], int):
						elems[-1] += 1
					else:
						elems.append(1)
				else:
					elems.append(("X", rb))
			r += n
			q += n
		elif op == "I":
			elems.append(("I", n))
			q += n
		elif op == "D":
			elems.append(("D", n))
			r += n
		elif op == "S":
			elems.append(("S", n))
			q += n
		else:
			raise ValueError(op)
	comp = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
	if strand == "-":
		elems = [e if isinstance(e, int) else (("X", comp[e[1]]) if e[0] == "X" else e) for e in reversed(elems)]
	out = []
	for e in elems:
		if isinstance(e, int):
			out.append(str(e))
		elif e[0] == "X":
			out.append(e[1])
		elif e[0] == "I":
			out.append(">%d-" % e[1])
		elif e[0] == "D":
			out.append(">%d+" % e[1])
		else:
			out.append("(%d)" % e[1])
	return "".join(out)


def parse_map_line(line):
	"""A MAP line -> dict(tag, read, qual, counters, maps); ``maps`` is the
	raw column ('-' when unmapped)."""
	f = line.split("\t")
	if len(f) == 5:
		return dict(tag=f[0], read=f[1], qual=f[2], counters=f[3], maps=f[4])
	if len(f) == 4:
		return dict(tag=f[0], read=f[1], qual=None, counters=f[2], maps=f[3])
	raise ValueError("MAP line with %d fields: %r" % (len(f), line))


def parse_se_maps(column):
	"""'chr:+:12:100:::60,...' -> [(rname, strand, pos, cigar, mapq)]."""
	if column == "-":
		return []
	out = []
	for m in column.split(","):
		head, mapq = m.rsplit(":::", 1)
		rname, strand, pos, cigar = head.split(":", 3)
		out.append((rname, strand, int(pos), cigar, int(mapq)))
	return out


def parse_pe_maps(column):
	"""'c:s:p:cig::c:s:p:cig:::mapq,...' -> [((end1), (end2), mapq)]."""
	if column == "-":
		return []
	out = []
	for m in column.split(","):
		body, mapq = m.rsplit(":::", 1)
		e1, e2 = body.split("::")
		a = e1.split(":", 3)
		b = e2.split(":", 3)
		out.append(((a[0], a[1], int(a[2]), a[3]), (b[0], b[1], int(b[2]), b[3]), int(mapq)))
	return out


def counters_ok(counters):
	"""MAP counters column: ':'-separated counts, at most one '+' (the
	maximum complete stratum)."""
	if counters.count("+") > 1:
		return False
	parts = counters.replace("+", ":").split(":")
	return all(x.isdigit() for x in parts)


def counter_values(counters):
	return [int(x) for x in counters.replace("+", ":").split(":")]


# ---------------------------------------------------------------------------
# --report-file
# ---------------------------------------------------------------------------

READ_TYPES = ("General", "SequencingControl", "UnderConversionControl", "OverConversionControl")
CONVERSION_TYPES = ("GeneralC2T", "GeneralG2A", "UnderConversionControlC2T",
	"UnderConversionControlG2A", "OverConversionControlC2T", "OverConversionControlG2A")


def read_type(rname):
	"""Report classification of an alignment by contig name."""
	for i, name in enumerate((SEQUENCING_CONTROL, UNDERCONVERSION_CONTROL, OVERCONVERSION_CONTROL)):
		if name in rname:
			return i + 1
	return 0


def _bases(seq):
	c = {"A": 0, "C": 0, "G": 0, "T": 0}
	for b in seq:
		if b in c:
			c[b] += 1
	return c


class _Acc(object):
	def __init__(self, ends):
		self.ends = ends
		self.reads = [[0] * ends for _ in READ_TYPES]
		self.unmapped = [0] * ends
		self.bs = {"C2T": [0] * ends, "G2A": [0] * ends}
		self.bases = {"Overall": [_bases("") for _ in range(ends)]}
		self.mapq = [0] * 256
		self.readlen = [{} for _ in range(ends)]
		self.mism = [{} for _ in range(ends)]
		self.tlen = {}
		self.correct = 0

	def add_bases(self, kind, end, seq):
		if kind not in self.bases:
			self.bases[kind] = [_bases("") for _ in range(self.ends)]
		for b, n in _bases(seq).items():
			self.bases[kind][end][b] += n

	def count(self, d, key):
		d[key] = d.get(key, 0) + 1

	def mapped_end(self, end, rec, read, genome, bisulfite):
		seq = forward_seq(read.seq, rec.strand)
		conv = rec.tags["XB"] if bisulfite else None
		self.count(self.mism[end], mismatches(genome[rec.rname], rec.ref_start, rec.cigar_ops, seq, conv))
		if bisulfite:
			self.bs["C2T" if conv == "C" else "G2A"][end] += 1

	def as_json(self, mapper_type, read_group=None):
		e = self.ends
		out = {"MapperType": mapper_type}
		if read_group is not None:
			out["ReadGroup"] = read_group
		reads = {}
		for i, name in enumerate(READ_TYPES):
			if i == 0 or any(self.reads[i]):
				reads[name] = self.reads[i]
		reads["Unmapped"] = self.unmapped
		out["Reads"] = reads
		if any(sum(v) for v in self.bs.values()):
			out["NumReadsBS"] = self.bs
		if e == 2:
			out["CorrectPairs"] = self.correct
		bc = {}
		for kind in ("Overall",) + CONVERSION_TYPES:
			if kind in self.bases:
				bc[kind] = {b: [self.bases[kind][k][b] for k in range(e)] for b in "ACGT"}
		out["BaseCounts"] = bc
		top = max([i for i in range(256) if self.mapq[i]] or [0])
		out["HistMapq"] = self.mapq[:top + 1]
		out["HistReadLen"] = [{str(k): v for k, v in sorted(d.items())} for d in self.readlen]
		out["HistMismatch"] = [{str(k): v for k, v in sorted(d.items())} for d in self.mism]
		if e == 2:
			out["HistTemplateLen"] = {str(k): v for k, v in sorted(self.tlen.items())}
		return out


def strip_report_n(report):
	"""The report without its 'N' base counts (known bug: one extra N per
	read, see known-bugs/D-report-basecounts-extra-n.md)."""
	r = json.loads(json.dumps(report))
	for kind in r.get("BaseCounts", {}).values():
		kind.pop("N", None)
	return r


def expected_report_se(reads, sam, genome, read_group=None):
	"""``--report-file`` contents for single-end reads, from the input reads
	and the primary SAM records (N base counts left out)."""
	acc = _Acc(1)
	prim = primaries_by_name(sam)
	for read in reads:
		rec = prim[read.name][0]
		acc.add_bases("Overall", 0, read.seq)
		acc.count(acc.readlen[0], len(read.seq))
		if rec.flag & UNMAPPED:
			acc.unmapped[0] += 1
			continue
		acc.mapped_end(0, rec, read, genome, False)
		acc.mapq[rec.mapq] += 1
		acc.reads[read_type(rec.rname)][0] += 1
	return acc.as_json("Single", read_group)


def expected_report_pe(pairs, sam, genome, bisulfite=False, discordant=()):
	"""``--report-file`` contents for pairs (N base counts left out).

	A pair printed as a proper pair (0x2) is a paired map: both ends count
	as mapped under the first end's read type, its template MAPQ goes to
	HistMapq, and, when concordant, it counts in CorrectPairs and |TLEN| in
	HistTemplateLen. ``discordant`` names the pairs the input was built to
	pair as discordant (SAM does not mark them). Otherwise each end counts
	on its own and HistMapq gets one 0. Conversion base counts (bisulfite)
	come from paired maps with MAPQ > 0.
	"""
	acc = _Acc(2)
	prim = primaries_by_name(sam)
	cnv = {1: ("GeneralC2T", None, "UnderConversionControlC2T", "OverConversionControlC2T"),
		2: ("GeneralG2A", None, "UnderConversionControlG2A", "OverConversionControlG2A")}
	for a, b in pairs:
		recs = prim[a.name]
		r1 = [r for r in recs if r.flag & READ1][0]
		r2 = [r for r in recs if r.flag & READ2][0]
		for end, read in ((0, a), (1, b)):
			acc.add_bases("Overall", end, read.seq)
			acc.count(acc.readlen[end], len(read.seq))
		if r1.flag & PROPER_PAIR:
			rt = read_type(r1.rname)
			if a.name not in discordant:
				acc.correct += 1
				acc.count(acc.tlen, abs(r1.tlen))
			acc.mapq[r1.mapq] += 1
			for end, rec, read in ((0, r1, a), (1, r2, b)):
				acc.reads[rt][end] += 1
				acc.mapped_end(end, rec, read, genome, bisulfite)
				if bisulfite and r1.mapq > 0:
					kind = cnv[1 if rec.tags["XB"] == "C" else 2][rt]
					if kind:
						acc.add_bases(kind, end, read.seq)
			continue
		acc.mapq[0] += 1
		for end, rec, read in ((0, r1, a), (1, r2, b)):
			if rec.flag & UNMAPPED:
				acc.unmapped[end] += 1
			else:
				acc.reads[read_type(rec.rname)][end] += 1
				acc.mapped_end(end, rec, read, genome, bisulfite)
	return acc.as_json("Paired")


def digest(data):
	return hashlib.sha256(data).hexdigest()


def body_lines(data):
	"""SAM records (no header) as a list of str."""
	text = data.decode("ascii") if isinstance(data, bytes) else data
	return [l for l in text.split("\n") if l and not l.startswith("@")]
