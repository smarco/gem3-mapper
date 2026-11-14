# GEM-Mapper v3 (GEM3) test suite.
# Copyright (c) 2026 Jacob Schreiber <jacob.schreiber@umassmed.edu>
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of GEM-Mapper v3 (GEM3), distributed under the GNU General
# Public License version 3 or (at your option) any later version; see LICENSE.
"""Paired-end mapping: pair orientation, template length (-l/-L, TLEN),
discordant pairs (--discordant-pair-search), mates in repeats, ties,
unpaired and unmapped mates, interleaved input (-p), -M and the MAPQ models.

Expected values come from the simulator's truth (fragment, positions,
strands), the designed repeats of ``area_c.repeat_genome()`` and the SAM
specification: TLEN is the span from the leftmost to the rightmost mapped
base, positive on the leftmost mate. Defaults from ``--help``: concordant
orientation FR, discordant RF, concordant layouts separate and overlap,
discordant contain, -L 10000, discordant search if-no-concordant.

Not tested (known-bugs/C-tlen-sam-spec.md, C-pair-layout-end2-contains-end1.md,
C-mapq-depends-on-max-reported-matches.md, D-pe-sam-compact-false-crash.md):
TLEN of mates at the same position, of overlapping RF pairs and of unpaired
mates; pairs whose read 2 contains read 1; MAPQ under -M; uncompacted output
of pairs whose mates have different numbers of alignments.
"""

import pytest

from gem3harness import (Genome, Rng, check_record_against_reference, sam_digest, simulate_pe, write_interleaved,
	write_split, revcomp)
from gem3harness.sam import (MATE_REVERSE, MATE_UNMAPPED, PAIRED, PROPER_PAIR, READ1, READ2, REVERSE,
	SECONDARY, UNMAPPED)
from area_c import PREDICTOR_FIELDS, primary_pairs, repeat_genome, split_predictors


def fastq(records):
	return "".join("@%s\n%s\n+\n%s\n" % (name, seq, "I" * len(seq)) for name, seq in records)


def map_pairs(gem3, index, pairs, args, threads, tmp_path):
	fq1, fq2 = write_split(tmp_path / "1.fq", tmp_path / "2.fq", pairs)
	return gem3.run_map(index, inputs=[fq1, fq2], args=list(args), threads=threads)


def map_records(gem3, index, recs1, recs2, args, threads, tmp_path):
	(tmp_path / "1.fq").write_text(fastq(recs1))
	(tmp_path / "2.fq").write_text(fastq(recs2))
	return gem3.run_map(index, inputs=[str(tmp_path / "1.fq"), str(tmp_path / "2.fq")], args=list(args),
		threads=threads)


def at_truth(rec, read):
	t = read.truth
	return (not rec.is_unmapped and (rec.rname, rec.pos, rec.strand) == (t.contig, t.pos + 1, t.strand))


def check_mates(rec, mate):
	"""Flag and mate fields of two mapped mates (SAMv1 1.4)."""
	assert rec.has(PAIRED) and mate.has(PAIRED)
	assert rec.has(READ1) != rec.has(READ2) and rec.has(READ1) == mate.has(READ2)
	assert not rec.has(UNMAPPED) and not rec.has(MATE_UNMAPPED)
	assert rec.has(MATE_REVERSE) == mate.has(REVERSE)
	assert rec.rnext == ("=" if rec.rname == mate.rname else mate.rname)
	assert rec.pnext == mate.pos


def check_proper_pair(rec1, rec2, pair):
	"""Both mates at their origin, flagged as a proper pair with SAM TLEN."""
	for rec, read, mate in ((rec1, pair.r1, rec2), (rec2, pair.r2, rec1)):
		assert at_truth(rec, read), (read, rec)
		assert rec.has(PROPER_PAIR), rec
		check_mates(rec, mate)
	span = max(rec1.ref_end, rec2.ref_end) - min(rec1.ref_start, rec2.ref_start)
	assert span == pair.fragment_length
	left, right = (rec1, rec2) if rec1.pos < rec2.pos else (rec2, rec1)
	assert (left.tlen, right.tlen) == (span, -span), (rec1, rec2)
	assert rec1.mapq == rec2.mapq


def check_unpaired(rec1, rec2, pair):
	"""Both mates mapped at their origin but not as a pair."""
	for rec, read, mate in ((rec1, pair.r1, rec2), (rec2, pair.r2, rec1)):
		assert at_truth(rec, read), (read, rec)
		assert not rec.has(PROPER_PAIR), rec
		check_mates(rec, mate)


# Orientation and discordant search ---------------------------------------------

_SIM = {}


def sim_pairs(genome, orientation, n=80, **kw):
	"""Unique-region pairs, fragments 210-390 bp (100 bp mates never touch)."""
	key = (orientation, n, tuple(sorted(kw.items())))
	if key not in _SIM:
		args = dict(insert_mean=300, insert_sd=30, insert_min=210, insert_max=390, seed=81,
			orientation=orientation, sub_rate=0.01)
		args.update(kw)
		_SIM[key] = simulate_pe(genome, n, 100, **args)
	return _SIM[key]


@pytest.mark.parametrize("mode", ["if-no-concordant", "always", "never"])
@pytest.mark.parametrize("orientation", ["FR", "RF", "FF", "RR"])
def test_orientation_and_discordant_search(gem3, genome, genome_index, tmp_path, threads, orientation, mode):
	"""FR pairs are concordant: proper pairs with TLEN = fragment length.
	RF pairs are discordant: reported as pairs with MAPQ 0 and TLEN = fragment
	length unless --discordant-pair-search never, which leaves them unpaired.
	FF and RR pairs are never paired. Every mate is placed at its origin."""
	pairs = sim_pairs(genome, orientation)
	sam = map_pairs(gem3, genome_index, pairs, ["--discordant-pair-search", mode], threads, tmp_path).sam()
	by = primary_pairs(sam)
	assert len(by) == len(pairs)
	for p in pairs:
		r1, r2 = by[p.r1.name]
		if orientation == "FR":
			check_proper_pair(r1, r2, p)
			assert r1.mapq > 3
		elif orientation == "RF" and mode != "never":
			for rec, read, mate in ((r1, p.r1, r2), (r2, p.r2, r1)):
				assert at_truth(rec, read), (read, rec)
				check_mates(rec, mate)
				assert rec.mapq == 0, rec
			left, right = (r1, r2) if r1.pos < r2.pos else (r2, r1)
			assert (left.tlen, right.tlen) == (p.fragment_length, -p.fragment_length)
		else:
			check_unpaired(r1, r2, p)
		for rec in (r1, r2):
			assert check_record_against_reference(rec, genome[rec.rname]) == [], rec


# Template length ------------------------------------------------------------------

@pytest.mark.parametrize("limits", [("-L", 400), ("-l", 300), ("both", (260, 450))])
def test_template_length_limits(gem3, genome, genome_index, tmp_path, threads, limits):
	"""-l/-L bound the template length of a pair (inclusive): pairs outside
	are not paired at all, not even as discordant pairs."""
	pairs = sim_pairs(genome, "FR", n=120, insert_mean=350, insert_sd=110, insert_min=201, insert_max=600, seed=82)
	opt, value = limits
	if opt == "both":
		lo, hi = value
		args = ["-l", lo, "-L", hi]
	elif opt == "-L":
		lo, hi = 0, value
		args = ["-L", value]
	else:
		lo, hi = value, 10000
		args = ["-l", value]
	for mode in ("if-no-concordant", "always"):
		sam = map_pairs(gem3, genome_index, pairs, args + ["--discordant-pair-search", mode], threads, tmp_path).sam()
		by = primary_pairs(sam)
		inside = 0
		for p in pairs:
			r1, r2 = by[p.r1.name]
			if lo <= p.fragment_length <= hi:
				check_proper_pair(r1, r2, p)
				inside += 1
			else:
				check_unpaired(r1, r2, p)
		assert 0 < inside < len(pairs)


def test_template_length_limit_boundaries(gem3, genome, genome_index, tmp_path, threads):
	"""Fragments of exactly 300 bp: -L 300 and -l 300 keep them, -L 299 and
	-l 301 do not."""
	pairs = sim_pairs(genome, "FR", n=20, insert_sd=0, insert_min=300, insert_max=300, seed=83)
	for args, paired in ((["-L", 300], True), (["-L", 299], False), (["-l", 300], True), (["-l", 301], False)):
		by = primary_pairs(map_pairs(gem3, genome_index, pairs, args, threads, tmp_path).sam())
		for p in pairs:
			r1, r2 = by[p.r1.name]
			(check_proper_pair if paired else check_unpaired)(r1, r2, p)


def test_long_templates_are_discordant(gem3, genome, genome_index, tmp_path, threads):
	"""Pairs far beyond the expected template length (3,000 bp against an
	estimate starting at 0-1,000 bp and learnt from 300 bp pairs) but within
	-L are discordant: paired with MAPQ 0 by default, unpaired with
	--discordant-pair-search never or with -L below their length."""
	normal = sim_pairs(genome, "FR", n=100, insert_sd=20, seed=84)
	long_ = sim_pairs(genome, "FR", n=15, insert_mean=3000, insert_sd=0, insert_min=3000, insert_max=3000,
		seed=85, name_prefix="long")
	pairs = normal + long_
	by = primary_pairs(map_pairs(gem3, genome_index, pairs, [], threads, tmp_path).sam())
	for p in normal:
		check_proper_pair(*by[p.r1.name], pair=p)
	for p in long_:
		r1, r2 = by[p.r1.name]
		for rec, read, mate in ((r1, p.r1, r2), (r2, p.r2, r1)):
			assert at_truth(rec, read) and rec.mapq == 0, rec
			check_mates(rec, mate)
		left, right = (r1, r2) if r1.pos < r2.pos else (r2, r1)
		assert (left.tlen, right.tlen) == (3000, -3000)
	for args in (["--discordant-pair-search", "never"], ["-L", 2999, "--discordant-pair-search", "always"]):
		by = primary_pairs(map_pairs(gem3, genome_index, pairs, args, threads, tmp_path).sam())
		for p in normal:
			check_proper_pair(*by[p.r1.name], pair=p)
		for p in long_:
			check_unpaired(*by[p.r1.name], pair=p)


def test_discordant_search_only_without_concordant(gem3, tmp_path, threads):
	"""A pair with a concordant placement and a discordant one (mate 2 also
	matches 8,000 bp away, beyond the expected template length): every mode
	reports the concordant pair as primary; only 'always' keeps the
	discordant pair, whose mate-2 placement then appears in mate 2's XA.
	(Uncompacted output is not used: known-bugs/D-pe-sam-compact-false-crash.md.)"""
	rg = repeat_genome()
	index = gem3.index(rg.genome)
	seq = rg.seq
	s = rg.dup_starts[0] + 50          # mate 2 in the duplicated block
	other = rg.dup_starts[1] + 50      # its second placement
	m1 = seq[s - 300:s - 200]          # mate 1 in unique sequence before the block
	m2 = revcomp(seq[s:s + 100])
	for mode in ("if-no-concordant", "never", "always"):
		sam = map_records(gem3, index, [("p", m1)], [("p", m2)],
			["--discordant-pair-search", mode, "-M", "all"], threads, tmp_path).sam()
		r1, r2 = primary_pairs(sam)["p"]
		assert len(sam.records) == 2
		assert (r1.pos, r1.tlen, r2.pos, r2.tlen) == (s - 300 + 1, 400, s + 1, -400)
		assert r1.has(PROPER_PAIR) and r2.has(PROPER_PAIR)
		if mode == "always":
			assert r2.tags["XA"] == "chrR,-%d,100M,0;" % (other + 1)
		else:
			assert "XA" not in r2.tags and "XA" not in r1.tags


# Repeats ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def rg():
	return repeat_genome()


@pytest.fixture(scope="module")
def rg_index(gem3, rg):
	return gem3.index(rg.genome)


def triple_pairs(rg):
	"""Mate 1 in unique sequence next to a copy of the 300 bp repeat, mate 2
	inside that copy (three equally good placements on its own)."""
	seq = rg.seq
	recs1, recs2, truth = [], [], []
	for k, (p, st) in enumerate(rg.triple_copies):
		for j, off in enumerate((0, 60, 150)):
			if st == "+":
				a, b = p - 250 + off, p + 100 + off      # mate 1 left (+), mate 2 inside (-)
				recs1.append(("t%d_%d" % (k, j), seq[a:a + 100]))
				recs2.append(("t%d_%d" % (k, j), revcomp(seq[b:b + 100])))
			else:
				a, b = p + 300 + 150 + off, p + off       # mate 2 inside (+), mate 1 right (-)
				recs1.append(("t%d_%d" % (k, j), revcomp(seq[a:a + 100])))
				recs2.append(("t%d_%d" % (k, j), seq[b:b + 100]))
			truth.append(("t%d_%d" % (k, j), a, b))
	return recs1, recs2, truth


def test_repeat_mate_placed_next_to_its_partner(gem3, rg, rg_index, tmp_path, threads):
	"""A mate with three identical placements is placed at the copy next to
	its uniquely placed partner: a proper pair with the true template length
	and a MAPQ above a tie's."""
	recs1, recs2, truth = triple_pairs(rg)
	by = primary_pairs(map_records(gem3, rg_index, recs1, recs2, [], threads, tmp_path).sam())
	for name, a, b in truth:
		r1, r2 = by[name]
		assert (r1.pos - 1, r2.pos - 1) == (a, b), (r1, r2)
		assert r1.has(PROPER_PAIR) and r2.has(PROPER_PAIR)
		check_mates(r1, r2)
		span = max(a, b) + 100 - min(a, b)
		left, right = (r1, r2) if r1.pos < r2.pos else (r2, r1)
		assert (left.tlen, right.tlen) == (span, -span)
		assert r1.mapq == r2.mapq and r1.mapq > 3
		assert r1.tags["NM"] == r2.tags["NM"] == 0
	# On its own, mate 2 has three placements
	(tmp_path / "se.fq").write_text(fastq(recs2[:1]))
	se = gem3.run_map(rg_index, inputs=[str(tmp_path / "se.fq")], threads=threads).sam()
	assert len(se.records[0].tags["XA"].rstrip(";").split(";")) == 2


def dup_pairs(rg):
	"""Pairs inside the 2,000 bp block duplicated exactly at two places."""
	seq = rg.seq
	d0, d1 = rg.dup_starts
	recs1, recs2, offsets = [], [], []
	for k, (off, frag, r1_left) in enumerate(((100, 300, True), (500, 350, False), (1200, 420, True))):
		left = seq[d0 + off:d0 + off + 100]
		right = revcomp(seq[d0 + off + frag - 100:d0 + off + frag])
		recs1.append(("d%d" % k, left if r1_left else right))
		recs2.append(("d%d" % k, right if r1_left else left))
		offsets.append(("d%d" % k, off, frag))
	return recs1, recs2, offsets


def test_pair_tie_between_two_copies(gem3, rg, rg_index, tmp_path, threads):
	"""A pair inside an exactly duplicated block has two equally good
	placements: the primary is at either copy with MAPQ <= 3, each mate's
	XA names its placement in the other copy, and -M all with uncompacted
	output reports both pairs (the second as secondary)."""
	recs1, recs2, offsets = dup_pairs(rg)
	d0, d1 = rg.dup_starts
	sam = map_records(gem3, rg_index, recs1, recs2, [], threads, tmp_path).sam()
	by = primary_pairs(sam)
	for name, off, frag in offsets:
		r1, r2 = by[name]
		base = min(r1.pos, r2.pos) - 1 - off
		assert base in (d0, d1)
		other = d1 if base == d0 else d0
		for rec in (r1, r2):
			assert rec.has(PROPER_PAIR) and rec.mapq <= 3, rec
			xa = rec.tags["XA"].rstrip(";").split(";")
			assert xa == ["chrR,%s%d,100M,0" % ("-" if rec.is_reverse else "+", rec.pos - base + other)], rec
		left, right = (r1, r2) if r1.pos < r2.pos else (r2, r1)
		assert (left.tlen, right.tlen) == (frag, -frag)
	full = map_records(gem3, rg_index, recs1, recs2, ["-M", "all", "--sam-compact=false"], threads, tmp_path).sam()
	for name, off, frag in offsets:
		recs = [r for r in full.records if r.qname == name]
		assert len(recs) == 4
		assert [r.has(SECONDARY) for r in recs] == [False, False, True, True]
		bases = sorted(set(min(a.pos, b.pos) - 1 - off for a, b in (recs[:2], recs[2:])))
		assert bases == [d0, d1]
		for rec in recs[2:]:
			assert rec.seq == "*" and rec.mapq == 0
	one = map_records(gem3, rg_index, recs1, recs2, ["-M", 1], threads, tmp_path).sam()
	assert len(one.records) == 2 * len(offsets)
	assert all("XA" not in r.tags for r in one.records)


def test_mate_rescue_in_high_copy_repeat(gem3, tmp_path, threads):
	"""Mate 2 lies in a 100 bp unit repeated 300 times (more placements than
	the 100 matches a paired-end search keeps per end); mate 1 is unique.
	The mate is found next to its partner by extending mate 1's alignment
	(mate rescue): a proper pair with the true template length."""
	rng = Rng(96)
	bg = list(rng.dna(200000))
	unit = rng.dna(100)
	copies = [1000 + 600 * k for k in range(300)]
	for p in copies:
		bg[p:p + 100] = list(unit)
	g = Genome([("chrU", "".join(bg))])
	index = gem3.index(g)
	seq = g["chrU"]
	chosen = [copies[k] for k in (3, 50, 120, 250, 299)]
	recs1 = [("q%d" % j, seq[p - 250:p - 150]) for j, p in enumerate(chosen)]
	recs2 = [("q%d" % j, revcomp(seq[p:p + 100])) for j, p in enumerate(chosen)]
	by = primary_pairs(map_records(gem3, index, recs1, recs2, [], threads, tmp_path).sam())
	for j, p in enumerate(chosen):
		r1, r2 = by["q%d" % j]
		assert (r1.pos, r2.pos) == (p - 250 + 1, p + 1), (r1, r2)
		assert (r1.tlen, r2.tlen) == (350, -350)
		assert r1.has(PROPER_PAIR) and r2.has(PROPER_PAIR)
		check_mates(r1, r2)


def tandem_array_pairs():
	"""A 60 kbp contig with 25 tandem copies of a 60 bp unit at 30,000; 1,000
	unique pairs (so the template-length estimate is reliable) and then 60
	pairs whose mates both fall inside the array, where many pairs tie on
	score and template length but differ in layout and position."""
	rng = Rng(606)
	bg = rng.dna(60000)
	array = rng.dna(60) * 25
	seq = bg[:30000] + array + bg[30000 + len(array):]
	g = Genome([("chrT", seq)])
	pairs = simulate_pe(g, 1000, 100, seed=607, insert_mean=300, insert_sd=30, name_prefix="u")
	recs1 = [(p.r1.name, p.r1.seq) for p in pairs]
	recs2 = [(p.r2.name, p.r2.seq) for p in pairs]
	for i in range(60):
		frag = rng.randint(180, 420)
		start = rng.randint(30000, 30000 + len(array) - frag)
		a, b = seq[start:start + 100], revcomp(seq[start + frag - 100:start + frag])
		if i % 2:
			a, b = b, a
		recs1.append(("tr%d" % i, a))
		recs2.append(("tr%d" % i, b))
	return g, recs1, recs2


# sha256 of the SAM output without @PG/@CO on base 6e155e3, the same from the
# five builds README.md lists at -t 1 and 4 (pins the current output). Tied
# pairs are ordered by the order in which they are enumerated
# (known-bugs/r2-pe-tied-pair-order.md); "Mall" changes if the layout
# tie-break of paired_matches_cmp_swg_score changes.
TANDEM_PINNED = {
	"default": "4c73924b4d1d8321e1e7fb055461b9fdf5a1b494c737d1f8c2361ca8488a8a3b",
	"Mall": "110dc9747907afd1652c801c04860e3d3297873a57916c8efa59e8cd569108e3",
}


@pytest.mark.parametrize("case", sorted(TANDEM_PINNED))
def test_tandem_array_pairs_are_pinned(gem3, tmp_path, threads, case):
	"""Pins the current output of pairs inside a tandem array (tied pairs)."""
	g, recs1, recs2 = tandem_array_pairs()
	args = {"default": [], "Mall": ["-M", "all"]}[case]
	res = map_records(gem3, gem3.index(g), recs1, recs2, args, threads, tmp_path)
	assert sam_digest(res.stdout) == TANDEM_PINNED[case]


# Unpaired and unmapped mates ------------------------------------------------------

def test_unmapped_and_cross_contig_mates(gem3, genome, genome_index, tmp_path, threads):
	"""A mate that cannot map leaves its partner mapped alone (mate-unmapped
	flag; the unmapped record takes the partner's RNAME and POS, as SAM
	recommends). Mates on different contigs are reported unpaired with the
	other contig in RNEXT and TLEN 0."""
	pairs = sim_pairs(genome, "FR", n=12, seed=86)
	rng = Rng(87)
	noise = [rng.dna(100) for _ in range(4)]
	recs1 = [(p.r1.name, p.r1.seq) for p in pairs]
	recs2 = [(p.r2.name, p.r2.seq) for p in pairs]
	for k in range(4):
		recs2[k] = (recs2[k][0], noise[k])
	# Mate 2 of pairs 4-7 replaced by mate 2 of a pair on another contig
	for k in range(4, 8):
		donor = next(q for q in pairs[8:] if q.contig != pairs[k].contig)
		recs2[k] = (recs2[k][0], donor.r2.seq)
	by = primary_pairs(map_records(gem3, genome_index, recs1, recs2, [], threads, tmp_path).sam())
	for k, p in enumerate(pairs):
		r1, r2 = by[p.r1.name]
		if k < 4:
			assert at_truth(r1, p.r1)
			assert r1.has(MATE_UNMAPPED) and not r1.has(PROPER_PAIR) and not r1.has(UNMAPPED)
			assert r2.has(UNMAPPED) and r2.cigar == "*" and r2.mapq == 0
			assert (r2.rname, r2.pos) == (r1.rname, r1.pos)
			assert (r1.rnext, r1.pnext) == ("=", r1.pos)
			assert r1.tlen == r2.tlen == 0
		elif k < 8:
			donor = next(q for q in pairs[8:] if q.contig != p.contig)
			assert at_truth(r1, p.r1) and at_truth(r2, donor.r2)
			assert not r1.has(PROPER_PAIR) and not r2.has(PROPER_PAIR)
			check_mates(r1, r2)
			assert r1.rnext == r2.rname != r1.rname
			assert r1.tlen == r2.tlen == 0
		else:
			check_proper_pair(r1, r2, p)


def test_interleaved_input_equals_split_input(gem3, genome, genome_index, tmp_path, threads):
	"""-p with interleaved FASTQ gives the same SAM as -1/-2."""
	pairs = simulate_pe(genome, 150, 100, insert_mean=250, insert_sd=60, seed=88, sub_rate=0.01,
		unique_only=False)
	split = map_pairs(gem3, genome_index, pairs, [], threads, tmp_path)
	il = write_interleaved(tmp_path / "il.fq", pairs)
	inter = gem3.run_map(genome_index, inputs=[il], args=["-p"], threads=threads)
	assert sam_digest(inter.stdout) == sam_digest(split.stdout)


# MAPQ models ------------------------------------------------------------------------

@pytest.fixture(scope="module")
def mixed_pairs(genome):
	"""Pairs from everywhere (repeats included) with sequencing errors."""
	return simulate_pe(genome, 250, 100, insert_mean=260, insert_sd=50, seed=89, sub_rate=0.012,
		ins_rate=0.002, del_rate=0.002, unique_only=False)


def _without_mapq(rec):
	f = rec.line.split("\t")
	return f[:4] + f[5:]


def test_pe_mapq_model_none_and_classify(gem3, genome_index, mixed_pairs, tmp_path, threads):
	"""--mapq-model none and classify change only MAPQ. none: 0 everywhere
	(pins the current value). classify: paired records take the paired
	model's classes (discordant 0, local 1, ties 2 or 10-40, 50-80, 100-160,
	161); unpaired records keep single-end MAPQs."""
	gem = map_pairs(gem3, genome_index, mixed_pairs, [], threads, tmp_path).sam()
	none = map_pairs(gem3, genome_index, mixed_pairs, ["--mapq-model", "none"], threads, tmp_path).sam()
	cls = map_pairs(gem3, genome_index, mixed_pairs, ["--mapq-model", "classify"], threads, tmp_path).sam()
	assert [_without_mapq(r) for r in gem.records] == [_without_mapq(r) for r in none.records]
	assert [_without_mapq(r) for r in gem.records] == [_without_mapq(r) for r in cls.records]
	assert set(r.mapq for r in none.records) == {0}
	allowed = {0, 1, 2} | set(range(10, 41)) | set(range(50, 81)) | set(range(100, 162))
	for g, c in zip(gem.records, cls.records):
		if c.has(PROPER_PAIR) and c.is_primary:
			assert c.mapq in allowed, c
			if g.mapq == 0:
				assert c.mapq in (0, 1, 2) or 10 <= c.mapq <= 40, (g, c)
		elif not c.has(PROPER_PAIR):
			assert c.mapq == g.mapq
	assert any(c.mapq == 161 for c in cls.records)


def test_pe_dump_predictors(gem3, genome, genome_index, mixed_pairs, tmp_path, threads):
	"""--mapq-model dump-predictors prints one line per paired read with
	the read-1 name as tag. Recomputed from the two records: edit =
	(L - NM1 - NM2) / L and swg = (AS1 + AS2) / L with L the sum of the mate
	lengths, and mcs2 0."""
	res = map_pairs(gem3, genome_index, mixed_pairs, ["--mapq-model", "dump-predictors"], threads, tmp_path)
	sam, predictors = split_predictors(res.stdout)
	by = primary_pairs(sam)
	paired = dict((name, recs) for name, recs in by.items() if recs[0].has(PROPER_PAIR))
	tags = {}
	for fields in predictors:
		assert len(fields) == len(PREDICTOR_FIELDS), fields
		d = dict(zip(PREDICTOR_FIELDS, fields))
		assert d["tag"].endswith("/1")
		tags[d["tag"][:-2]] = d
	assert sorted(tags) == sorted(paired)
	for name, (r1, r2) in paired.items():
		d = tags[name]
		L = float(len(r1.seq) + len(r2.seq))
		assert abs(float(d["edit"]) - (L - r1.tags["NM"] - r2.tags["NM"]) / L) < 1e-6, (r1, r2, d)
		assert abs(float(d["swg"]) - (r1.tags["AS"] + r2.tags["AS"]) / L) < 1e-6, (r1, r2, d)
		assert int(d["mcs2"]) == 0


# Pinned bytes -------------------------------------------------------------------------

PINNED = {
	"default": "2e0f71689aecb4029a9067edaf33c9a05a572a831150c0c7d5866bac02fea55a",
	"always": "50fb11dfa817bc2ecaaa63e0085158d8233dcfec1551e2de55b8ee6a8a977be6",
	"never": "c8ab755d4d49152f7b1b9bd5fca129b81d9666cecffe8bd5a08b205f89e91121",
	"L400": "97b01f48695ad0cddbe3e11f5e3a7fe2c1968009e771315fbfeddfe32c1beaa0",
	"l260": "bc962ee55755f56d35e5bfc4b9eb1642ad44192e4ed30152927ea18152e2a30b",
	"M1": "d7e6effc663faf1bc39569291bad8b8dd00ce7d0e9874ce4ab0d54d7a0cca4f2",
	"Mall": "2e0f71689aecb4029a9067edaf33c9a05a572a831150c0c7d5866bac02fea55a",
	"classify": "c326193414bf80185ff86db781b43c35fc3d0e3170d32bd5a9b6111335d17e16",
	"rg-default": "4f756bbaf7ff61412e697610fd330d18bf2b16647d4687d62602e8c588020f98",
	"rg-Mall-full": "7eee790c99dab9696e5b5cbb649db5af5cc8bfe19ac3b5d3e17ae589bc997bc8",
	"rg-always": "6e8ca826a443b55b3c04a49cd81460f171e5a692baee069cd6d64ce2754c1c50",
}
PIN_ARGS = {
	"default": [],
	"always": ["--discordant-pair-search", "always"],
	"never": ["--discordant-pair-search", "never"],
	"L400": ["-L", 400],
	"l260": ["-l", 260],
	"M1": ["-M", 1],
	"Mall": ["-M", "all"],
	"classify": ["--mapq-model", "classify"],
	"rg-default": [],
	# Every end here has as many alignments as its template has pairs, so
	# XS reads no stale slot (README "Determinism", D-pe-sam-compact-false-crash.md).
	"rg-Mall-full": ["-M", "all", "--sam-compact=false"],
	"rg-always": ["--discordant-pair-search", "always", "-M", "all"],
}


@pytest.fixture(scope="module")
def pin_pairs(genome):
	"""Pairs from everywhere with errors, plus RF, FF, long-template pairs and
	pairs with a random mate 2, in one input."""
	pairs = simulate_pe(genome, 200, 100, insert_mean=260, insert_sd=60, seed=90, sub_rate=0.012,
		ins_rate=0.002, del_rate=0.002, unique_only=False)
	pairs += simulate_pe(genome, 20, 100, insert_mean=400, insert_sd=60, seed=91, orientation="RF",
		name_prefix="rf", sub_rate=0.01)
	pairs += simulate_pe(genome, 10, 100, insert_mean=300, insert_sd=30, seed=92, orientation="FF",
		name_prefix="ff")
	pairs += simulate_pe(genome, 10, 100, insert_mean=2500, insert_sd=400, seed=93, name_prefix="long")
	rng = Rng(94)
	extra = simulate_pe(genome, 10, 100, insert_mean=300, insert_sd=30, seed=95, name_prefix="noise")
	for p in extra:
		p.r2.seq = rng.dna(100)
	return pairs + extra


@pytest.mark.parametrize("case", sorted(PINNED))
def test_paired_output_bytes_are_pinned(gem3, genome_index, pin_pairs, rg, rg_index, tmp_path, threads, case):
	"""Pins the current output of paired-end mapping options on fixed inputs."""
	if case.startswith("rg-"):
		t1, t2, _ = triple_pairs(rg)
		d1, d2, _ = dup_pairs(rg)
		res = map_records(gem3, rg_index, t1 + d1, t2 + d2, PIN_ARGS[case], threads, tmp_path)
	else:
		res = map_pairs(gem3, genome_index, pin_pairs, PIN_ARGS[case], threads, tmp_path)
	assert sam_digest(res.stdout) == PINNED[case]
