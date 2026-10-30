# GEM-Mapper v3 (GEM3) test suite.
# Copyright (c) 2026 Jacob Schreiber <jacob.schreiber@umassmed.edu>
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of GEM-Mapper v3 (GEM3), distributed under the GNU General
# Public License version 3 or (at your option) any later version; see LICENSE.
"""SAM writer (src/io/output_sam.c): header, core fields, flags, mate
fields, optional tags, compact XA against --sam-compact=false, read groups,
FASTA input and --sam-gem-compatible.

Expected values come from the SAM specification (SAMv1), the reference and
the input reads (see area_d.py). Known bugs left out: XA's last field counts
soft-clipped bases (D-xa-distance-counts-clipping.md), RNEXT of an unmapped
mate names the contig instead of '=' (D-rnext-unmapped-mate.md; accepted
either way here),
paired-end --sam-compact=false can crash (D-pe-sam-compact-false-crash.md).
"""

import os
import pwd
import socket
import time

import pytest

import area_d
from area_d import (body_lines, pair_problems, primaries_by_name, se_record_problems, xa_entries)
from gem3harness import parse_sam, revcomp, sam_digest, write_reads
from gem3harness.sam import READ1, READ2, SECONDARY, UNMAPPED


@pytest.fixture(scope="module")
def se(gem3, genome, genome_index):
	return area_d.se_dataset(gem3, genome, genome_index)


@pytest.fixture(scope="module")
def pe(gem3, genome, genome_index):
	return area_d.pe_dataset(gem3, genome, genome_index)


def version(gem3):
	from gem3harness import run_process
	res = run_process([gem3.mapper, "--version"])
	assert res.returncode == 0
	return res.stderr_text.strip()


# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------

def test_header_lines(gem3, genome, se, tmp_path):
	"""@HD, one @SQ per contig in FASTA order with its length, @PG with the
	version and the command line, @CO with time, working directory, host and
	user (output_sam_print_header, system_get_cwd/hostname/user_name)."""
	t0 = time.time()
	res = gem3.run_map(se.index, inputs=se.paths, threads=1, cwd=str(tmp_path))
	t1 = time.time()
	sam = res.sam()
	h = sam.header_lines
	assert h[0] == "@HD\tVN:1.4\tSO:unsorted"
	sq = h[1:1 + len(genome.lengths())]
	assert sq == ["@SQ\tSN:%s\tLN:%d" % (n, L) for n, L in genome.lengths()]
	assert [l[:3] for l in h] == ["@HD"] + ["@SQ"] * len(sq) + ["@PG", "@CO"]
	pg = sam.header("@PG")[0]
	assert (pg["ID"], pg["PN"], pg["VN"]) == ("GEM", "gem-mapper", version(gem3))
	assert pg["CL"] == " ".join(res.argv)
	co = sam.header("@CO")[0]
	assert co["WD"] == os.path.realpath(str(tmp_path))
	assert co["HN"] == socket.gethostname()
	assert co["UN"] == pwd.getpwuid(os.geteuid()).pw_name
	# TM:YYYY/M/D HH:MM:SS CET, local time at the run
	stamp = co["TM"]
	assert stamp.endswith(" CET")
	got = time.mktime(time.strptime(stamp[:-4], "%Y/%m/%d %H:%M:%S"))
	assert int(t0) - 1 <= got <= t1 + 1, (stamp, t0, t1)


def test_bisulfite_index_header_lists_plain_contigs(gem3, genome, genome_bs_index, tmp_path):
	"""A bisulfite index stores converted copies of every contig; the header
	lists each original contig once."""
	reads = area_d.se_mix(genome)[:5]
	fq = write_reads(tmp_path / "r.fq", reads)
	sam = gem3.run_map(genome_bs_index, inputs=[fq]).sam()
	assert sam.sq == genome.lengths()


def test_command_line_escaped_in_pg(gem3, se):
	"""Tabs and newlines in arguments are written as \\t and \\n in @PG CL; a
	real tab in -r reaches the @RG line as a tab."""
	res = se.run(["-r", "@RG\tID:tab\tSM:x"], threads=1)
	head = res.stdout.decode().split("\n")
	pg = [l for l in head if l.startswith("@PG")][0]
	assert pg.count("\t") == 4  # ID, PN, VN, CL: no raw tab from the argument
	assert "-r @RG\\tID:tab\\tSM:x" in pg
	assert "@RG\tID:tab\tSM:x" in head


# ---------------------------------------------------------------------------
# Single-end records
# ---------------------------------------------------------------------------

def test_se_records_follow_spec(genome, se):
	"""One primary record per read in input order; FLAG, RNAME..TLEN,
	SEQ/QUAL orientation, NM/MD/AS of every primary against the reference;
	unmapped records as SAMv1 describes them."""
	sam = se.run().sam()
	prim = [r for r in sam.records if r.is_primary]
	assert [r.qname for r in prim] == [r.name for r in se.reads]
	for rec, read in zip(prim, se.reads):
		assert se_record_problems(rec, read, genome) == [], (rec, read.seq)
	kinds = {}
	for rec, read in zip(prim, se.reads):
		kinds.setdefault(read.name[0], set()).add(bool(rec.flag & UNMAPPED))
	# The input covers mapped (unique, repeat, short, N) and unmapped reads.
	assert kinds["u"] == {False} and kinds["m"] == {False} and True in kinds["x"]
	assert any("XA" in r.tags for r in prim)


def test_se_tag_order(se):
	"""Optional fields of a mapped record are NM, MD, AS then XA (compact)."""
	for line in body_lines(se.run().stdout):
		f = line.split("\t")
		tags = [x[:2] for x in f[11:]]
		if f[1] == "4":
			assert tags == []
		else:
			assert tags in (["NM", "MD", "AS"], ["NM", "MD", "AS", "XA"]), line


def test_se_compact_xa_lists_the_secondary_records(genome, se):
	"""Compact XA:Z of a primary lists the alignments that
	--sam-compact=false prints as secondary records, in the same order; the
	primary records are the same apart from XA. Secondary records follow
	SAMv1 (0x100, SEQ/QUAL '*') and their NM/MD/AS match the reference."""
	compact = se.run().sam()
	full = se.run(["--sam-compact=false"]).sam()
	byname = full.by_name()
	reads = {r.name: r for r in se.reads}
	n_alt = 0
	for rec in compact.records:
		assert rec.is_primary
		recs = byname[rec.qname]
		assert recs[0].line == "\t".join(rec.line.split("\t")[:11] + [f for f in rec.line.split("\t")[11:]
			if not f.startswith("XA:")])
		secondary = recs[1:]
		assert all(r.flag & SECONDARY for r in secondary)
		got = [(r.rname, r.strand, r.pos, r.cigar) for r in secondary]
		xa = xa_entries(rec)
		assert [e[:4] for e in xa] == got, rec.qname
		for e, r in zip(xa, secondary):
			assert se_record_problems(r, reads[rec.qname], genome, primary=False) == [], r
			# XA's last field is NM for unclipped alignments (clipped: known bug).
			if area_d.soft_clipped(r.cigar_ops) == 0:
				assert e[4] == r.tags["NM"], (e, r)
			n_alt += 1
	assert n_alt >= 20


def test_se_max_reported_matches_limits_records(se):
	"""-M n prints at most n records per read with --sam-compact=false and at
	most n-1 XA entries in compact mode; the first n alignments are the same
	as with the default -M 5."""
	full = se.run(["--sam-compact=false"]).sam().by_name()
	two = se.run(["--sam-compact=false", "-M", "2"]).sam().by_name()
	compact = se.run(["-M", "2"]).sam()
	n_one = 0
	for rec in compact.records:
		xa = xa_entries(rec)
		assert len(xa) <= 1
		# A single alternative is still listed in XA
		assert [e[:4] for e in xa] == [(r.rname, r.strand, r.pos, r.cigar) for r in two[rec.qname][1:]]
		n_one += len(xa) == 1
	assert n_one >= 10
	for name, recs in two.items():
		assert 1 <= len(recs) <= 2
		assert [(r.rname, r.pos, r.cigar) for r in recs] == [(r.rname, r.pos, r.cigar) for r in full[name][:len(recs)]]


def test_fasta_input_prints_star_qual(gem3, se, tmp_path):
	"""FASTA reads print QUAL '*' (SAMv1) and otherwise the same records as
	the FASTQ reads mapped with -q ignore."""
	fa = write_reads(tmp_path / "r.fa", se.reads, fmt="fasta")
	a = gem3.run_map(se.index, inputs=[fa]).sam()
	b = se.run(["-q", "ignore"]).sam()
	assert len(a.records) == len(b.records)
	for x, y in zip(a.records, b.records):
		fx, fy = x.line.split("\t"), y.line.split("\t")
		assert fx[10] == "*"
		assert fx[:10] + fx[11:] == fy[:10] + fy[11:]


def test_qname_drops_comment_and_mate_suffix(gem3, genome, se, tmp_path):
	"""QNAME is the read name up to the first space, without /1 or /2."""
	reads = se.reads[:6]
	text = []
	for i, r in enumerate(reads):
		name = [r.name + " 1:N:0:ACGT", r.name + "/1", r.name + "/2 comment", r.name + "\tx",
			r.name, r.name + " a b c"][i]
		text.append("@%s\n%s\n+\n%s\n" % (name, r.seq, r.qual))
	p = tmp_path / "names.fq"
	p.write_text("".join(text))
	sam = gem3.run_map(se.index, inputs=[str(p)]).sam()
	assert [r.qname for r in sam.records if r.is_primary] == [r.name for r in reads]


# ---------------------------------------------------------------------------
# Read groups
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("header,rg_line,rg_id", [
	("@RG\\tID:grp1\\tSM:s1\\tPL:ILLUMINA", "@RG\tID:grp1\tSM:s1\tPL:ILLUMINA", "grp1"),
	("@RG\\tSM:s2\\tID:lane.7\\n", "@RG\tSM:s2\tID:lane.7", "lane.7"),
], ids=["escaped-tabs", "id-not-first-trailing-newline"])
def test_read_group(genome, se, pe, header, rg_line, rg_id):
	"""-r: the header gets the @RG line (escapes \\t/\\n translated, a trailing
	newline dropped) between @SQ and @PG, and every record, single or paired,
	RG:Z:<ID> as its first optional field; nothing else changes."""
	for ds in (se, pe):
		res = ds.run(["-r", header])
		sam = res.sam()
		types = [l[:3] for l in sam.header_lines]
		assert types.index("@RG") == types.index("@PG") - 1
		assert [l for l in sam.header_lines if l.startswith("@RG")] == [rg_line]
		plain = body_lines(ds.run().stdout)
		got = body_lines(res.stdout)
		assert len(plain) == len(got)
		for a, b in zip(plain, got):
			f = b.split("\t")
			assert f[11] == "RG:Z:" + rg_id
			assert "\t".join(f[:11] + f[12:]) == a


# ---------------------------------------------------------------------------
# Paired-end records
# ---------------------------------------------------------------------------

def test_pe_records_follow_spec(genome, pe):
	"""Two primary records per pair (end 1 first), in input order; FLAG
	bits, RNEXT/PNEXT, TLEN (leftmost to rightmost mapped base, + on the
	leftmost), the placement of an unmapped mate and NM/MD/AS, from SAMv1."""
	sam = pe.run().sam()
	prim = [r for r in sam.records if r.is_primary]
	assert len(prim) == 2 * len(pe.pairs)
	seen = {"proper": 0, "one": 0, "none": 0}
	for k, (a, b) in enumerate(pe.pairs):
		r1, r2 = prim[2 * k], prim[2 * k + 1]
		assert r1.flag & READ1 and r2.flag & READ2
		assert pair_problems(r1, r2, a, b, genome) == [], (r1, r2)
		um = bool(r1.flag & UNMAPPED) + bool(r2.flag & UNMAPPED)
		seen["proper" if r1.flag & 2 else ("none" if um == 2 else "one")] += 1
	assert seen["proper"] > 70 and seen["one"] > 10 and seen["none"] >= 3, seen


def test_pe_compact_xa_per_end(genome, pe):
	"""Paired XA:Z. In a proper pair, entry i of end 1 and entry i of end 2
	are the two ends of alternative pairing i: equally many, on one contig.
	A mapped end of an unpaired pair lists its own alternatives. Every
	entry's CIGAR matches the reference and the read."""
	sam = pe.run().sam()
	prim = [r for r in sam.records if r.is_primary]
	pairs = {a.name: (a, b) for a, b in pe.pairs}
	n_paired = n_single = 0
	for k in range(0, len(prim), 2):
		r1, r2 = prim[k], prim[k + 1]
		x1, x2 = xa_entries(r1), xa_entries(r2)
		a, b = pairs[r1.qname]
		if r1.flag & 2:
			assert len(x1) == len(x2), (r1, r2)
			assert [e[0] for e in x1] == [e[0] for e in x2]
			n_paired += len(x1)
		else:
			n_single += len(x1) + len(x2)
		for xs, read in ((x1, a), (x2, b)):
			for e in xs:
				rec = _record_from_xa(e)
				assert area_d.alignment_problems(rec, area_d.forward_seq(read.seq, e[1]), genome) == [], e
	assert n_paired >= 5 and n_single >= 2, (n_paired, n_single)


class _XaRecord(object):
	def __init__(self, e):
		from gem3harness import parse_cigar, cigar_ref_span
		self.rname = e[0]
		self.cigar = e[3]
		self.cigar_ops = parse_cigar(e[3])
		self.ref_start = e[2] - 1
		self.ref_end = self.ref_start + cigar_ref_span(self.cigar_ops)
		self.tags = {}


def _record_from_xa(e):
	return _XaRecord(e)


def test_bisulfite_records_carry_xb_not_md(genome, gem3, genome_bs_index):
	"""Bisulfite output: mapped records carry NM, XB (C or G) and AS, no MD."""
	bs = area_d.bs_dataset(gem3, genome, genome_bs_index, n=60)
	for line in body_lines(bs.run().stdout):
		f = line.split("\t")
		tags = [x[:4] for x in f[11:]]
		if int(f[1]) & UNMAPPED:
			continue
		assert tags[:3] == ["NM:i", "XB:A", "AS:i"], line
		assert f[12] in ("XB:A:C", "XB:A:G")


# ---------------------------------------------------------------------------
# --sam-gem-compatible
# ---------------------------------------------------------------------------

def test_gem_compatible_adds_x4_x5(se, pe):
	"""--sam-gem-compatible prints the --sam-compact=false records with X4:Z
	(MAP counters) and X5:Z (MAP matches) appended to each primary record
	and no XA; unmapped primaries get X5:Z:-."""
	full = body_lines(se.run(["--sam-compact=false"]).stdout)
	gem = body_lines(se.run(["--sam-gem-compatible"]).stdout)
	assert len(full) == len(gem)
	for a, b in zip(full, gem):
		fa, fb = a.split("\t"), b.split("\t")
		if int(fa[1]) & SECONDARY:
			assert a == b
			continue
		assert [x[:5] for x in fb[-2:]] == ["X4:Z:", "X5:Z:"], b
		assert "\t".join(fb[:-2]) == a
		if int(fa[1]) & UNMAPPED:
			assert fb[-1] == "X5:Z:-"
	# Paired: the compact records plus X4/X5 on every primary, no XA.
	for line in body_lines(pe.run(["--sam-gem-compatible"]).stdout):
		f = line.split("\t")
		assert not any(x.startswith("XA:") for x in f[11:])
		if not int(f[1]) & SECONDARY:
			assert [x[:5] for x in f[-2:]] == ["X4:Z:", "X5:Z:"], line


def test_sam_compact_explicit_true_is_default(se):
	"""--sam-compact=true and --sam-gem-compatible=false are the defaults."""
	d = sam_digest(se.run().stdout)
	assert sam_digest(se.run(["--sam-compact=true"]).stdout) == d
	assert sam_digest(se.run(["--sam-compact"]).stdout) == d
	assert sam_digest(se.run(["--sam-gem-compatible=false"]).stdout) == d


# ---------------------------------------------------------------------------
# Byte pins
# ---------------------------------------------------------------------------

# sha256 of the SAM output without @PG/@CO on base 6e155e3 (pins the current
# output; the speed search gates on byte identity).
PINNED = {
	"se-default": "2711545c5118e477b70783ef4ae33238bba80a11344eaae108dae0b8849c6472",
	"se-noncompact": "d7fee76642ccba7fb2c052310c7a9f63b13fe0776d3154761f948586e6eba834",
	"se-gemcompat": "ddbb39486c4cef394a15832e8a877cc91377fdf14cb95fa6489dfaa7faeb27dc",
	"se-readgroup": "4716771e0b2fca89071ebd66000c951cd1973249189cf6c6f7613759de9b02f1",
	"pe-default": "c7c82ae5087589bfaf11f5ef0811791278c9663a62164399b41b1ca12fe72cfb",
	"pe-gemcompat": "61c7e6364e91feeccd2fa0f6376e87cacdf69b3da411e74dd211459868674757",
}

PIN_ARGS = {
	"se-default": [],
	"se-noncompact": ["--sam-compact=false"],
	"se-gemcompat": ["--sam-gem-compatible"],
	"se-readgroup": ["-r", "@RG\\tID:grp1\\tSM:s1"],
	"pe-default": [],
	"pe-gemcompat": ["--sam-gem-compatible"],
}


@pytest.mark.parametrize("case", sorted(PINNED))
def test_sam_output_bytes_are_pinned(se, pe, case):
	"""Pins the current output of the SAM writer's modes on the mixed inputs."""
	data = (se if case.startswith("se") else pe).run(PIN_ARGS[case])
	assert sam_digest(data.stdout) == PINNED[case]
