# GEM-Mapper v3 (GEM3) test suite.
# Copyright (c) 2026 Jacob Schreiber <jacob.schreiber@umassmed.edu>
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of GEM-Mapper v3 (GEM3), distributed under the GNU General
# Public License version 3 or (at your option) any later version; see LICENSE.
"""Pilot end-to-end tests: simulated reads land where they came from.

Expected values come from the simulator's truth and from recomputing NM, MD
and the gap-affine score from the reference (default --gap-affine-penalties
1,4,6,1: match +1, mismatch -4, a gap of length k costs 6 + k).
"""

import pytest

from gem3harness import (bisulfite_nm, check_record_against_reference, cigar_is_valid,
	revcomp, simulate_bisulfite, simulate_pe, simulate_se, write_reads, write_split)
from gem3harness.sam import PAIRED, PROPER_PAIR, READ1, READ2, MATE_REVERSE, REVERSE


def gap_affine_score(ref, pos, ops, seq, match=1, mismatch=4, gap_open=6, gap_extend=1):
	"""Gap-affine score of an alignment, as gem3-mapper's AS:i reports it."""
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


def primary_by_name(sam):
	out = {}
	for rec in sam.records:
		if rec.is_primary:
			out.setdefault(rec.qname, []).append(rec)
	return out


def _forward_seq(read):
	return read.seq if read.truth.strand == "+" else revcomp(read.seq)


@pytest.mark.parametrize("edits", [
	{"subs": 0},
	{"subs": 2, "ins": 1, "dels": 1},
	{"sub_rate": 0.02, "ins_rate": 0.003, "del_rate": 0.003},
], ids=["exact", "2sub-1ins-1del", "rates"])
def test_se_unique_reads_map_to_truth(gem3, genome, genome_index, tmp_path, threads, edits):
	reads = simulate_se(genome, 300, 100, seed=21, **edits)
	fq = write_reads(tmp_path / "se.fq", reads)
	res = gem3.run_map(genome_index, inputs=[fq], threads=threads)
	sam = res.sam()
	assert sam.sq == genome.lengths()
	by_name = primary_by_name(sam)
	assert sorted(by_name) == sorted(r.name for r in reads)
	for read in reads:
		recs = by_name[read.name]
		assert len(recs) == 1, recs
		rec, t = recs[0], read.truth
		assert not rec.is_unmapped, (read, rec)
		assert (rec.rname, rec.pos, rec.strand) == (t.contig, t.pos + 1, t.strand), (read, rec)
		# SEQ is the read on the forward reference strand; QUAL follows it.
		assert rec.seq == _forward_seq(read)
		assert rec.qual == (read.qual if t.strand == "+" else read.qual[::-1])
		assert cigar_is_valid(rec.cigar_ops, 100), rec
		assert check_record_against_reference(rec, genome[rec.rname]) == [], rec
		assert rec.tags["NM"] <= t.nm, (read, rec)
		# AS is the alignment's gap-affine score, and the aligner finds an
		# alignment at least as good as the simulated one.
		score = gap_affine_score(genome[rec.rname], rec.ref_start, rec.cigar_ops, rec.seq)
		assert rec.tags["AS"] == score, rec
		assert score >= gap_affine_score(genome[t.contig], t.pos, t.cigar, _forward_seq(read)), (read, rec)
		assert 0 <= rec.mapq <= 255


def test_pe_fr_pairs_map_to_truth(gem3, genome, genome_index, tmp_path, threads):
	pairs = simulate_pe(genome, 200, 100, insert_mean=300, insert_sd=30, seed=31,
		sub_rate=0.01, ins_rate=0.002, del_rate=0.002)
	fq1, fq2 = write_split(tmp_path / "r1.fq", tmp_path / "r2.fq", pairs)
	res = gem3.run_map(genome_index, inputs=[fq1, fq2], threads=threads)
	sam = res.sam()
	by_name = primary_by_name(sam)
	assert len(by_name) == len(pairs)
	for p in pairs:
		recs = by_name[p.r1.name]
		assert len(recs) == 2, recs
		first = [r for r in recs if r.has(READ1)]
		second = [r for r in recs if r.has(READ2)]
		assert len(first) == 1 and len(second) == 1, recs
		a, b = first[0], second[0]
		for rec, read, mate in ((a, p.r1, b), (b, p.r2, a)):
			t = read.truth
			assert rec.has(PAIRED) and rec.has(PROPER_PAIR), rec
			assert (rec.rname, rec.pos, rec.strand) == (t.contig, t.pos + 1, t.strand), (read, rec)
			assert check_record_against_reference(rec, genome[rec.rname]) == [], rec
			# Mate fields describe the other record.
			assert rec.rnext == "=" and rec.pnext == mate.pos
			assert rec.has(MATE_REVERSE) == mate.has(REVERSE)
		# TLEN: fragment length, positive on the leftmost mate (SAMv1 1.4).
		frag = p.frag_end - p.frag_start
		left, right = (a, b) if a.pos < b.pos else (b, a)
		assert (left.tlen, right.tlen) == (frag, -frag), (p.frag_start, p.frag_end, a, b)


def test_bisulfite_directional_pairs_map_to_truth(gem3, genome, genome_bs_index, tmp_path, threads):
	"""Directional bisulfite pairs on a -b index with the default
	--bisulfite-read (inferred): every pair maps to its origin, XB:A tells
	which strand was converted (C: C-to-T on the forward reference, G: G-to-A),
	and NM counts no conversion as an edit."""
	pairs = simulate_bisulfite(genome, 150, 100, library="directional", seed=41,
		cpg_rate=0.8, chh_rate=0.01)
	fq1, fq2 = write_split(tmp_path / "b1.fq", tmp_path / "b2.fq", pairs)
	res = gem3.run_map(genome_bs_index, inputs=[fq1, fq2], threads=threads)
	sam = res.sam()
	# The -b index adds no contigs to the header.
	assert sam.sq == genome.lengths()
	by_name = primary_by_name(sam)
	for p in pairs:
		recs = by_name[p.r1.name]
		assert len(recs) == 2
		for rec in recs:
			read = p.r1 if rec.has(READ1) else p.r2
			t = read.truth
			assert (rec.rname, rec.pos, rec.strand) == (t.contig, t.pos + 1, t.strand), (read, rec)
			assert rec.has(PROPER_PAIR)
			# SEQ is the read as sequenced (unconverted), on the forward strand.
			assert rec.seq == _forward_seq(read)
			converted = "C" if t.extra["origin"] in ("OT", "CTOT") else "G"
			assert rec.tags["XB"] == converted, (t.extra["origin"], rec)
			assert rec.tags["NM"] == bisulfite_nm(genome[rec.rname], rec.ref_start, rec.cigar_ops,
				rec.seq, converted) == 0, rec


# sha256 of the SAM output without @PG/@CO, from base 6e155e3 (see README.md,
# "Determinism", for the builds and thread counts it was checked on).
PINNED = {
	"se": "c97ed817dd095a3f6c0d72626716a7ee6f44d884a1c7387e8f3702fad4f1003d",
	"pe": "637d510193941a930463c86993244fb4b25e78a22744d4303682e6b16d1e57b9",
	"bisulfite": "205d4916404f66661369a003f40bb177b9d6923f79a89b2e8c8f0cb609424e2b",
}


@pytest.mark.parametrize("case", sorted(PINNED))
def test_pilot_output_bytes_are_pinned(gem3, genome, genome_index, genome_bs_index, tmp_path, threads, case):
	"""Pins the current output: the SAM bytes (minus the variable @PG/@CO
	lines) for fixed simulated inputs, including repeat and N-adjacent reads."""
	from gem3harness import sam_digest
	if case == "se":
		reads = simulate_se(genome, 400, 100, seed=51, sub_rate=0.01, ins_rate=0.002,
			del_rate=0.002, unique_only=False)
		fq = write_reads(tmp_path / "se.fq", reads)
		res = gem3.run_map(genome_index, inputs=[fq], threads=threads)
	elif case == "pe":
		pairs = simulate_pe(genome, 300, 100, insert_mean=250, insert_sd=60, seed=52,
			sub_rate=0.01, ins_rate=0.002, del_rate=0.002, unique_only=False)
		fq1, fq2 = write_split(tmp_path / "1.fq", tmp_path / "2.fq", pairs)
		res = gem3.run_map(genome_index, inputs=[fq1, fq2], threads=threads)
	else:
		pairs = simulate_bisulfite(genome, 200, 100, seed=53, unique_only=False)
		fq1, fq2 = write_split(tmp_path / "1.fq", tmp_path / "2.fq", pairs)
		res = gem3.run_map(genome_bs_index, inputs=[fq1, fq2], threads=threads)
	assert sam_digest(res.stdout) == PINNED[case]
