# GEM-Mapper v3 (GEM3) test suite.
# Copyright (c) 2026 Jacob Schreiber <jacob.schreiber@umassmed.edu>
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of GEM-Mapper v3 (GEM3), distributed under the GNU General
# Public License version 3 or (at your option) any later version; see LICENSE.
"""Bisulfite mapping on a -b index: --bisulfite-read modes for single-end
reads, directional pairs, XB:A, NM without conversions, SEQ as sequenced,
and the control-sequence options counted in --report-file.

Expected values come from the simulator's truth (origin strand OT, OB,
CTOT, CTOB) and from ``bisulfite_nm``. XB:A is 'C' when the read aligns to
the C-to-T converted forward strand (OT and CTOT reads), 'G' otherwise.

Not tested (known-bugs/): single-end '2' and 'interleaved' act as '1';
paired-end 'non-stranded' is never applied.
"""

import json

import pytest

from gem3harness import (Genome, Rng, SimRead, bisulfite_nm, revcomp, sam_digest, simulate_bisulfite,
	write_reads, write_split)
from gem3harness.sam import PROPER_PAIR, READ1
from area_c import primary_pairs


def converted_strand(read):
	return "C" if read.truth.extra["origin"] in ("OT", "CTOT") else "G"


def forward_seq(read):
	return read.seq if read.truth.strand == "+" else revcomp(read.seq)


def check_bisulfite_record(rec, read, genome):
	t = read.truth
	assert (rec.rname, rec.pos, rec.strand) == (t.contig, t.pos + 1, t.strand), (read, rec)
	assert rec.seq == forward_seq(read), rec      # the read as sequenced, not converted
	assert rec.tags["XB"] == converted_strand(read), (t.extra["origin"], rec)
	assert "MD" not in rec.tags
	assert rec.tags["NM"] == bisulfite_nm(genome[rec.rname], rec.ref_start, rec.cigar_ops, rec.seq,
		rec.tags["XB"]), rec


def with_substitutions(reads, rng, k):
	"""Copies of the reads with k A->C substitutions each, 10+ bp apart; an
	A->C change is never a bisulfite conversion in either direction."""
	out = []
	for r in reads:
		seq = list(r.seq)
		offsets = [o for o in range(10, len(seq) - 10, 12) if seq[o] == "A"]
		chosen = offsets[:k]
		for o in chosen:
			seq[o] = "C"
		out.append((SimRead(r.name, "".join(seq), r.qual, r.truth), len(chosen)))
	return out


@pytest.fixture(scope="module")
def directional(genome):
	return simulate_bisulfite(genome, 120, 100, library="directional", seed=101, cpg_rate=0.8, chh_rate=0.02)


@pytest.mark.parametrize("mode", ["inferred", "1", "non-stranded"])
def test_se_read1(gem3, genome, genome_bs_index, directional, tmp_path, threads, mode):
	"""Read 1 of directional pairs (C-to-T), named with /1: every mode that
	converts C-to-T places each read at its origin, with the same records."""
	reads = [p.r1 for p in directional]
	fq = write_reads(tmp_path / "r1.fq", reads, suffix="/1")
	sam = gem3.run_map(genome_bs_index, inputs=[fq], args=["--bisulfite-read", mode], threads=threads).sam()
	assert len(sam.records) == len(reads)
	for rec, read in zip(sam.records, reads):
		check_bisulfite_record(rec, read, genome)
		assert rec.tags["NM"] == 0
	ref = gem3.run_map(genome_bs_index, inputs=[fq], threads=threads).sam()
	assert [r.line for r in sam.records] == [r.line for r in ref.records]


@pytest.mark.parametrize("named,mode", [(True, "inferred"), (False, "non-stranded"), (True, "non-stranded")])
def test_se_read2(gem3, genome, genome_bs_index, directional, tmp_path, threads, named, mode):
	"""Read 2 of directional pairs (G-to-A): converted G-to-A when named /2
	(inferred) or, in non-stranded mode, because it has more C than G."""
	reads = [p.r2 for p in directional]
	fq = write_reads(tmp_path / "r2.fq", reads, suffix="/2" if named else "")
	sam = gem3.run_map(genome_bs_index, inputs=[fq], args=["--bisulfite-read", mode], threads=threads).sam()
	for rec, read in zip(sam.records, reads):
		check_bisulfite_record(rec, read, genome)
		assert rec.tags["NM"] == 0


def test_se_mixed_reads_non_stranded(gem3, genome, genome_bs_index, tmp_path, threads):
	"""Unnamed reads from both ends of directional and PBAT libraries, with
	two A->C substitutions each: non-stranded mode picks each read's
	conversion; NM counts the substitutions and no conversion."""
	rng = Rng(102)
	d = simulate_bisulfite(genome, 50, 100, library="directional", seed=103)
	b = simulate_bisulfite(genome, 50, 100, library="pbat", seed=104)
	reads = []
	for i, p in enumerate(d + b):
		for j, r in enumerate((p.r1, p.r2)):
			reads.append(SimRead("m%d_%d" % (i, j), r.seq, r.qual, r.truth))
	subst = with_substitutions(reads, rng, 2)
	fq = write_reads(tmp_path / "m.fq", [r for r, _ in subst])
	sam = gem3.run_map(genome_bs_index, inputs=[fq], args=["--bisulfite-read", "non-stranded"],
		threads=threads).sam()
	assert len(sam.records) == len(reads)
	for rec, (read, k) in zip(sam.records, subst):
		check_bisulfite_record(rec, read, genome)
		assert rec.tags["NM"] == k, (rec, k)


def test_pe_directional_pairs(gem3, genome, genome_bs_index, tmp_path, threads):
	"""Directional pairs with two A->C substitutions per mate: proper pairs
	at their origin, TLEN = fragment length, both mates on the same
	converted strand (XB), NM = the substitutions."""
	rng = Rng(105)
	pairs = simulate_bisulfite(genome, 80, 100, library="directional", seed=106, insert_mean=300, insert_sd=30)
	pairs = [p for p in pairs if p.fragment_length > 100]
	r1 = with_substitutions([p.r1 for p in pairs], rng, 2)
	r2 = with_substitutions([p.r2 for p in pairs], rng, 2)
	for p, (a, _), (b, _) in zip(pairs, r1, r2):
		p.r1, p.r2 = a, b
	fq1, fq2 = write_split(tmp_path / "1.fq", tmp_path / "2.fq", pairs)
	by = primary_pairs(gem3.run_map(genome_bs_index, inputs=[fq1, fq2], threads=threads).sam())
	for p, (_, k1), (_, k2) in zip(pairs, r1, r2):
		a, b = by[p.r1.name]
		check_bisulfite_record(a, p.r1, genome)
		check_bisulfite_record(b, p.r2, genome)
		assert (a.tags["NM"], b.tags["NM"]) == (k1, k2)
		assert a.has(PROPER_PAIR) and b.has(PROPER_PAIR) and a.has(READ1)
		assert a.tags["XB"] == b.tags["XB"]
		left, right = (a, b) if a.pos < b.pos else (b, a)
		assert (left.tlen, right.tlen) == (p.fragment_length, -p.fragment_length)


# Control sequences -------------------------------------------------------------------

def control_genome():
	rng = Rng(107)
	return Genome([("chrA", rng.dna(30000)), ("NC_001422.1", rng.dna(5000)), ("NC_001416.1", rng.dna(5000)),
		("NC_001604.1", rng.dna(5000)), ("ctlPhi", rng.dna(5000)), ("ctlLambda", rng.dna(5000)),
		("ctlT7", rng.dna(5000))])


@pytest.mark.parametrize("names", ["default", "options"])
def test_control_sequences_in_report(gem3, tmp_path, threads, names):
	"""Reads on the sequencing, under-conversion and over-conversion control
	sequences are counted apart in the report (Reads: SequencingControl,
	UnderConversionControl, OverConversionControl, General). The defaults
	are phiX (NC_001422.1), lambda (NC_001416.1) and T7 (NC_001604.1); the
	options name other sequences."""
	g = control_genome()
	index = gem3.index(g, bisulfite=True)
	plan = (("chrA", 6), ("NC_001422.1", 3), ("NC_001416.1", 2), ("NC_001604.1", 4), ("ctlPhi", 5),
		("ctlLambda", 1), ("ctlT7", 2))
	reads = []
	for k, (contig, n) in enumerate(plan):
		for i in range(n):
			reads += simulate_bisulfite(g, 1, 100, paired=False, seed=1000 * k + i, contigs=[contig],
				name_prefix="%s_%d_" % (contig, i))
	fq = write_reads(tmp_path / "c.fq", reads, suffix="/1")
	args = ["-i", "c.fq", "--report-file", "report.json"]
	if names == "options":
		args += ["--control_sequence", "ctlPhi", "--underconversion_sequence", "ctlLambda",
			"--overconversion_sequence", "ctlT7"]
		want = {"SequencingControl": 5, "UnderConversionControl": 1, "OverConversionControl": 2,
			"General": 6 + 3 + 2 + 4}
	else:
		want = {"SequencingControl": 3, "UnderConversionControl": 2, "OverConversionControl": 4,
			"General": 6 + 5 + 1 + 2}
	with open(fq, "rb") as fh:
		data = fh.read()
	res = gem3.run_map(index, args=args, files={"c.fq": data}, threads=threads)
	sam = res.sam()
	for rec, read in zip(sam.records, reads):
		assert (rec.rname, rec.pos) == (read.truth.contig, read.truth.pos + 1), rec
	report = json.loads(res.files["report.json"].decode("ascii"))
	got = dict((k, v) for k, v in report["Reads"].items() if k != "Unmapped")
	assert got == dict((k, [v]) for k, v in want.items())
	assert report["Reads"]["Unmapped"] == [0]


# Pinned bytes ----------------------------------------------------------------------------

PINNED = {
	"se-inferred": "f065cdfd599f1f7bb8791e236253582ece14169ec7abbe4bd100340c02d2a6fb",
	"se-non-stranded": "39622a8e586eb1c00534923580a0841b8b9b6ab597d28bec1eadabac71184134",
	"pe-classify": "59ebb1fe1582d7e2f4344b011b4aa85b0ba31315310cee4dda275f66afcb0c7d",
	"pe-M1": "790173f686cea0618abfcc84a362f55664ff364a9b828c3092033896e9e9b175",
}


@pytest.mark.parametrize("case", sorted(PINNED))
def test_bisulfite_output_bytes_are_pinned(gem3, genome, genome_bs_index, tmp_path, threads, case):
	"""Pins the current bisulfite output on fixed inputs, repeats included."""
	if case.startswith("se-"):
		d = simulate_bisulfite(genome, 120, 100, library="non-directional", seed=108, unique_only=False)
		reads = [SimRead("s%d_%d" % (i, j), r.seq, r.qual, r.truth)
			for i, p in enumerate(d) for j, r in enumerate((p.r1, p.r2))]
		if case == "se-inferred":
			d = simulate_bisulfite(genome, 120, 100, library="directional", seed=110, unique_only=False)
			fq1 = write_reads(tmp_path / "a.fq", [p.r1 for p in d], suffix="/1")
			fq2 = write_reads(tmp_path / "b.fq", [p.r2 for p in d], suffix="/2")
			data = b""
			for path in (fq1, fq2):
				with open(path, "rb") as fh:
					data += fh.read()
			res = gem3.run_map(genome_bs_index, args=["-i", "r.fq"], files={"r.fq": data}, threads=threads)
		else:
			fq = write_reads(tmp_path / "s.fq", reads)
			res = gem3.run_map(genome_bs_index, inputs=[fq], args=["--bisulfite-read", "non-stranded"],
				threads=threads)
	else:
		pairs = simulate_bisulfite(genome, 150, 100, library="directional", seed=109, unique_only=False)
		fq1, fq2 = write_split(tmp_path / "1.fq", tmp_path / "2.fq", pairs)
		args = ["--mapq-model", "classify"] if case == "pe-classify" else ["-M", 1]
		res = gem3.run_map(genome_bs_index, inputs=[fq1, fq2], args=args, threads=threads)
	assert sam_digest(res.stdout) == PINNED[case]
