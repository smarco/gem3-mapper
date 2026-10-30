# GEM-Mapper v3 (GEM3) test suite.
# Copyright (c) 2026 Jacob Schreiber <jacob.schreiber@umassmed.edu>
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of GEM-Mapper v3 (GEM3), distributed under the GNU General
# Public License version 3 or (at your option) any later version; see LICENSE.
"""Self-tests of the harness: generators, simulator truth, SAM helpers, build."""

import bz2
import gzip
import os

import pytest

from gem3harness import (Rng, bisulfite_nm, cigar_is_valid, cigar_query_length, cigar_ref_span,
	md_equal, md_is_wellformed, md_nm_from_alignment, parse_cigar, parse_sam, revcomp,
	simulate_bisulfite, simulate_pe, simulate_se, synthetic_genome, write_interleaved,
	write_reads, write_split)
from gem3harness.build import REPO_ROOT, TOOLS
from gem3harness.sam import md_ref_span


# ---------------------------------------------------------------------------
# Random numbers and the synthetic genome
# ---------------------------------------------------------------------------

def test_rng_matches_splitmix64_reference():
	"""SplitMix64 seeded with 0 starts with the published reference outputs."""
	r = Rng(0)
	assert [r.next_u64() for _ in range(3)] == [0xE220A8397B1DCDAF, 0x6E789E6AA1B965F4,
		0x06C45D188009454F]


def test_synthetic_genome_is_pinned():
	"""Pins the bytes of synthetic_genome(1): fixtures and any pinned SAM depend on them."""
	g = synthetic_genome(1)
	assert g.digest() == "7e6cb58b6321f330b8c2f6f54a272e007fa7f180eead61888cc40404cc5c9d47"
	assert synthetic_genome(1).fasta_text() == g.fasta_text()
	assert synthetic_genome(2).digest() != g.digest()


def test_synthetic_genome_features():
	g = synthetic_genome(1)
	kinds = set(f.kind for f in g.features)
	assert kinds == {"repeat", "tandem", "homopolymer", "N", "softmask", "short_contig"}
	for f in g.features:
		body = g[f.contig][f.start:f.end]
		if f.kind == "N":
			assert set(body.upper()) == {"N"}
		elif f.kind == "homopolymer":
			assert set(body.upper()) == {f.info["base"]}
		elif f.kind == "softmask":
			assert body == body.lower()
		elif f.kind == "tandem":
			unit = f.info["unit"]
			assert body.upper() == (unit * len(body))[:len(body)]
	# Copy 0 of each repeat family is exact; the others diverged a little.
	fams = {}
	for f in g.features:
		if f.kind == "repeat":
			fams.setdefault(f.info["family"], []).append(f)
	for copies in fams.values():
		assert len(copies) >= 2
		seqs = [g[f.contig][f.start:f.end].upper() for f in copies]
		lens = set(len(s) for s in seqs)
		assert len(lens) == 1
	# An N run sits at the start of chr1 and the end of chr2; two contigs are short.
	assert g["chr1"][0] == "N" and g["chr2"][-1] == "N"
	assert sorted(L for _, L in g.lengths())[:2] == [90, 400]
	assert set("".join(g.seqs.values())) <= set("ACGTNacgtn")


# ---------------------------------------------------------------------------
# CIGAR / MD / NM helpers against hand-worked examples
# ---------------------------------------------------------------------------

REF = "ACGTACGTAC"


@pytest.mark.parametrize("pos,cigar,seq,md,nm", [
	(0, "10M", "ACGTACGTAC", "10", 0),
	(0, "10M", "ACGAACGTAC", "3T6", 1),
	(0, "10M", "TTGTACGTAC", "0A0C8", 2),        # adjacent mismatches need the 0
	(0, "3M1I6M", "ACGGTACGTA", "9", 1),          # insertion adds to NM only
	(0, "4M2D4M", "ACGTGTAC", "4^AC4", 2),
	(0, "4M2D4M", "ACGTCTAC", "4^AC0G3", 3),      # mismatch right after a deletion
	(2, "2S8M", "GGGTACGTAC", "8", 0),            # soft clips are skipped
	(0, "10M", "ACNTACGTAC", "2G7", 1),           # a read N is a mismatch
	(1, "3=1X2=", "CGTTCG", "3A2", 1),
])
def test_md_nm_hand_examples(pos, cigar, seq, md, nm):
	assert md_nm_from_alignment(REF, pos, parse_cigar(cigar), seq) == (md, nm)
	assert md_is_wellformed(md)
	assert md_ref_span(md) == cigar_ref_span(parse_cigar(cigar))


def test_md_nm_sam_spec_example():
	"""Read r001/1 of the SAM specification's section 1.1 example
	(8M2I4M1D3M at position 7); MD and NM worked out by hand."""
	ref = "AGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT"
	md, nm = md_nm_from_alignment(ref, 6, parse_cigar("8M2I4M1D3M"), "TTAGATAAAGGATACTG")
	assert (md, nm) == ("12^G3", 3)


def test_md_nm_lowercase_and_reference_n():
	assert md_nm_from_alignment("acgtACGTac", 0, [("M", 10)], "ACGTACGTAC") == ("10", 0)
	# A reference N never matches, not even a read N.
	assert md_nm_from_alignment("ACNTA", 0, [("M", 5)], "ACNTA") == ("2N2", 1)


def test_bisulfite_nm():
	# C-to-T over a reference C is free with conversion C, not with G.
	assert bisulfite_nm("ACGTC", 0, [("M", 5)], "ATGTT", "C") == 0
	assert bisulfite_nm("ACGTC", 0, [("M", 5)], "ATGTT", "G") == 2
	assert bisulfite_nm("ACGTC", 0, [("M", 5)], "ACATC", "G") == 0
	assert bisulfite_nm("ACGTC", 0, [("M", 5)], "ACATC", "C") == 1
	assert bisulfite_nm("ACGTC", 0, [("M", 2), ("D", 1), ("M", 2)], "ATTT", "C") == 1


@pytest.mark.parametrize("cigar,length,valid", [
	("10M", 10, True), ("5S5M5S", 15, True), ("5H5S5M", 10, True), ("3S2M1I4M3H", 10, True),
	("5M5S5M", 15, False), ("5M5H5M", 10, False), ("0M10M", 10, False), ("10M", 11, False),
	("5S5I", 10, False), ("4M2D4M", 8, True),
])
def test_cigar_validity(cigar, length, valid):
	assert cigar_is_valid(parse_cigar(cigar), length) is valid


def test_cigar_lengths_and_md_equal():
	ops = parse_cigar("3S4M2I3M1D5M2S")
	assert cigar_query_length(ops) == 19
	assert cigar_ref_span(ops) == 13
	assert md_equal("10A5", "10A5") and not md_equal("10A5", "10C5")
	with pytest.raises(ValueError):
		parse_cigar("10M3")


def test_sam_parser_types():
	text = ("@HD\tVN:1.4\n@SQ\tSN:c1\tLN:100\n"
		"r1\t16\tc1\t5\t60\t3S7M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tNM:i:1\tMD:Z:2A4\t"
		"XB:A:C\tZF:f:1.5\tZB:B:i,1,-2,3\tZH:H:1AE3\n")
	sam = parse_sam(text)
	assert sam.sq == [("c1", 100)]
	r = sam.records[0]
	assert (r.qname, r.flag, r.rname, r.pos, r.mapq, r.cigar) == ("r1", 16, "c1", 5, 60, "3S7M")
	assert r.is_reverse and r.is_primary and not r.is_unmapped
	assert (r.ref_start, r.ref_end) == (4, 11)
	assert r.tags == {"NM": 1, "MD": "2A4", "XB": "C", "ZF": 1.5, "ZB": ("i", [1, -2, 3]),
		"ZH": "1AE3"}
	assert r.tag_types["NM"] == "i"


# ---------------------------------------------------------------------------
# Simulator truth
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def g():
	return synthetic_genome(1)


def _forward(read):
	return read.seq if read.truth.strand == "+" else revcomp(read.seq)


def test_se_truth_is_self_consistent(g):
	reads = simulate_se(g, 300, 100, seed=5, subs=2, ins=1, dels=1)
	assert len(set(r.name for r in reads)) == 300
	for r in reads:
		t = r.truth
		fwd = _forward(r)
		assert len(r.seq) == len(r.qual) == 100
		assert cigar_query_length(t.cigar) == 100
		assert cigar_ref_span(t.cigar) == t.ref_span
		# Exactly the requested edits: 2 substitutions, 1 insertion, 1 deletion.
		kinds = sorted(k for _, k, _ in t.edits)
		assert kinds == ["D", "I", "S", "S"]
		assert t.nm == sum(s for _, _, s in t.edits)
		assert md_nm_from_alignment(g[t.contig], t.pos, t.cigar, fwd) == (t.md, t.nm)
		# Unique reads stay out of repeats, N runs and short contigs.
		assert not g.hard_regions(t.contig, t.pos, t.end)
		assert "N" not in g[t.contig][t.pos:t.end].upper()


def test_se_reads_without_edits_are_reference_substrings(g):
	for r in simulate_se(g, 100, 75, seed=9):
		t = r.truth
		assert _forward(r) == g[t.contig][t.pos:t.pos + 75].upper()
		assert (t.cigar_string, t.md, t.nm) == ("75M", "75", 0)


def test_se_is_deterministic(g):
	a = simulate_se(g, 50, 80, seed=3, sub_rate=0.02)
	b = simulate_se(g, 50, 80, seed=3, sub_rate=0.02)
	assert [(r.seq, r.qual, r.truth.pos) for r in a] == [(r.seq, r.qual, r.truth.pos) for r in b]


@pytest.mark.parametrize("orientation", ["FR", "RF", "FF", "RR"])
def test_pe_truth_orientation_and_fragment(g, orientation):
	pairs = simulate_pe(g, 200, 100, insert_mean=300, insert_sd=40, orientation=orientation,
		seed=7, subs=1)
	r1_left = 0
	for p in pairs:
		left, right = sorted([p.r1, p.r2], key=lambda r: (r.truth.pos, r.truth is p.r2.truth))
		assert left.truth.pos == p.frag_start
		assert right.truth.end == p.frag_end
		assert (left.truth.strand, right.truth.strand) == (
			"+" if orientation[0] == "F" else "-", "+" if orientation[1] == "F" else "-")
		r1_left += left is p.r1
		for r in (p.r1, p.r2):
			t = r.truth
			assert md_nm_from_alignment(g[t.contig], t.pos, [op for op in t.cigar if op[0] != "S"],
				_forward(r) if not t.extra.get("adapter") else _forward_genomic(r)) == (t.md, t.nm)
	# Both strands of the fragment are sampled.
	assert 60 < r1_left < 140


def _forward_genomic(read):
	t = read.truth
	n = t.genomic_length
	fwd = read.seq[:n]
	return fwd if t.strand == "+" else revcomp(fwd)


def test_pe_adapter_read_through(g):
	pairs = simulate_pe(g, 100, 100, insert_mean=60, insert_sd=15, insert_min=30, seed=11)
	short = [p for p in pairs if p.fragment_length < 100]
	assert short
	for p in short:
		tail = 100 - p.fragment_length
		for r in (p.r1, p.r2):
			t = r.truth
			assert t.extra["adapter"] == tail
			assert t.genomic_length == p.fragment_length
			assert cigar_query_length(t.cigar) == 100
			# The adapter follows the genomic bases in the read as sequenced.
			assert r.seq[p.fragment_length:].startswith("AGATCGGAAGAGC"[:tail])
			assert _forward_genomic(r) == g[t.contig][t.pos:t.pos + p.fragment_length].upper()


@pytest.mark.parametrize("library", ["directional", "pbat", "non-directional"])
def test_bisulfite_truth(g, library):
	pairs = simulate_bisulfite(g, 150, 100, library=library, seed=13, cpg_rate=0.75, chh_rate=0.02)
	first_conv = set()
	for p in pairs:
		assert p.r1.truth.extra["conversion"] != p.r2.truth.extra["conversion"]
		first_conv.add(p.r1.truth.extra["conversion"])
		for r in (p.r1, p.r2):
			t = r.truth
			n = t.genomic_length
			fwd = r.seq[:n] if t.strand == "+" else revcomp(r.seq[:n])
			ref = g[t.contig][t.pos:t.pos + n].upper()
			# Converted strand in reference coordinates: C->T if the original
			# top strand was converted, G->A if the bottom one.
			top = t.extra["origin"] in ("OT", "CTOT")
			conv = "C" if top else "G"
			assert bisulfite_nm(g[t.contig], t.pos, [("M", n)], fwd, conv) == 0
			for i, (rb, qb) in enumerate(zip(ref, fwd)):
				if rb == conv:
					assert qb in (conv, "T" if conv == "C" else "A")
					assert (qb == conv) == t.extra["calls"][t.pos + i]
				else:
					assert qb == rb
	expected = {"directional": {"CT"}, "pbat": {"GA"}, "non-directional": {"CT", "GA"}}[library]
	assert first_conv == expected


def test_bisulfite_methylation_rates(g):
	pairs = simulate_bisulfite(g, 400, 100, seed=17, cpg_rate=0.8, chh_rate=0.0)
	cg = []
	for p in pairs:
		t = p.r1.truth
		top = t.extra["origin"] in ("OT", "CTOT")
		for pos, m in t.extra["calls"].items():
			ref = g[t.contig]
			is_cg = (ref[pos + 1:pos + 2].upper() == "G") if top else (ref[pos - 1:pos].upper() == "C")
			if is_cg:
				cg.append(m)
			else:
				assert not m
	assert 0.7 < sum(cg) / float(len(cg)) < 0.9


# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------

def test_writers_roundtrip(tmp_path, g):
	reads = simulate_se(g, 5, 30, seed=1)
	plain = write_reads(tmp_path / "r.fq", reads)
	gz = write_reads(tmp_path / "r.fq.gz", reads, compression="gz")
	bz = write_reads(tmp_path / "r.fq.bz2", reads, compression="bz2")
	data = open(plain, "rb").read()
	assert gzip.open(gz).read() == data
	assert bz2.open(bz).read() == data
	assert data.count(b"\n") == 20 and data.startswith(b"@r0\n")
	# gzip bytes do not depend on the time or file name.
	gz2 = write_reads(tmp_path / "other.fq.gz", reads, compression="gz")
	assert open(gz, "rb").read() == open(gz2, "rb").read()
	fa = write_reads(tmp_path / "r.fa", reads, fmt="fasta")
	assert open(fa).read().startswith(">r0\n" + reads[0].seq + "\n")
	pairs = simulate_pe(g, 3, 30, insert_mean=100, insert_sd=5, seed=1)
	il = open(write_interleaved(tmp_path / "i.fq", pairs)).read().split("\n")
	assert il[0] == "@p0/1" and il[4] == "@p0/2" and il[8] == "@p1/1"
	a, b = write_split(tmp_path / "1.fq", tmp_path / "2.fq", pairs)
	assert open(a).read().split("\n")[1] == pairs[0].r1.seq
	assert open(b).read().split("\n")[1] == pairs[0].r2.seq


# ---------------------------------------------------------------------------
# The build under test
# ---------------------------------------------------------------------------

def test_build_leaves_source_tree_untouched(gem3_build):
	"""The out-of-tree build adds, removes and modifies nothing in the source tree."""
	assert gem3_build.tree_changes == []
	if gem3_build.root is not None and gem3_build.info.get("target") != "external":
		root = os.path.realpath(gem3_build.root)
		assert not root.startswith(os.path.realpath(REPO_ROOT) + os.sep)


def test_build_has_tools(gem3_build):
	for tool in TOOLS[:2]:
		assert os.access(gem3_build.tool(tool), os.X_OK)


def test_index_is_built_in_its_own_directory(gem3, genome):
	"""gem-indexer's <name>.sa.tmp goes to its working directory; the harness
	runs it in the index's directory and the temporary file is gone after."""
	idx = gem3.index(genome)
	d = os.path.dirname(idx)
	assert os.path.exists(idx)
	assert not [f for f in os.listdir(d) if f.endswith(".sa.tmp")]
	assert not [f for f in os.listdir(os.getcwd()) if f.endswith(".sa.tmp")]


def test_strip_variable_header():
	from gem3harness import sam_digest, strip_variable_header
	a = b"@HD\tVN:1.4\n@PG\tID:GEM\tCL:x -t 1\n@CO\tTM:1\nr\t4\t*\t0\t0\t*\t*\t0\t0\tA\tI\n"
	b = b"@HD\tVN:1.4\n@PG\tID:GEM\tCL:x -t 8\n@CO\tTM:2\nr\t4\t*\t0\t0\t*\t*\t0\t0\tA\tI\n"
	assert strip_variable_header(a) == b"@HD\tVN:1.4\nr\t4\t*\t0\t0\t*\t*\t0\t0\tA\tI\n"
	assert sam_digest(a) == sam_digest(b)
	assert sam_digest(a) != sam_digest(a.replace(b"\tA\t", b"\tC\t"))
