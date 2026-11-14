# GEM-Mapper v3 (GEM3) test suite.
# Copyright (c) 2026 Jacob Schreiber <jacob.schreiber@umassmed.edu>
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of GEM-Mapper v3 (GEM3), distributed under the GNU General
# Public License version 3 or (at your option) any later version; see LICENSE.
"""Mapping modes, the complete-search error (-E) and the alignment models.

Expected values come from the simulator's truth, from recomputing CIGAR, NM,
MD and the gap-affine score against the reference, and from a brute-force
search of a small planted genome (area_b.loci_within: bit-parallel Myers for
edit distance, a gapless scan for Hamming distance) for every locus within
E errors of a probe read.
"""

import pytest

from gem3harness import (check_record_against_reference, cigar_is_valid, revcomp, simulate_se,
	write_reads)
import area_b

MODES = ["fast", "sensitive", "hybrid", "customed", "complete", "complete-pure",
	"complete-partition", "complete-partition-filter"]

# One gem-mapper run is ~20 ms; the planted genome is 14 kbp.
NO_LOCAL = ["--alignment-local", "never"]
ALL_MATCHES = ["-M", "all", "--sam-compact=false"]


def _forward_seq(read):
	return read.seq if read.truth.strand == "+" else revcomp(read.seq)


@pytest.fixture(scope="module")
def unique_reads(genome):
	"""Unique 100 bp reads with 0-3 edits (subs and 1-3 bp indels at least 8 bp
	from the ends), so every mode's search depth (default -E 0.04, 4 errors)
	covers them."""
	reads = simulate_se(genome, 60, 100, seed=61, subs=2, ins=1, max_indel=1)
	reads += simulate_se(genome, 60, 100, seed=62, subs=1, dels=1, max_indel=2)
	reads += simulate_se(genome, 30, 100, seed=63)
	return reads


_gotoh_cache = {}


def _best_score(seq, ref, start, end, pen=(1, 4, 6, 1)):
	key = (seq, id(ref), start, end, pen)
	if key not in _gotoh_cache:
		_gotoh_cache[key] = area_b.gotoh_best(seq, ref[max(0, start - 12):end + 12], *pen)
	return _gotoh_cache[key]


@pytest.mark.parametrize("mode", MODES)
def test_mode_maps_unique_reads_to_truth(gem3, genome, genome_index, tmp_path, threads, unique_reads, mode):
	"""Every mode places unique reads with at most 3 edits at their origin;
	the alignment is consistent with the reference and AS is the optimal
	gap-affine score (1,4,6,1) of the read in a window around it."""
	fq = write_reads(tmp_path / "r.fq", unique_reads)
	res = gem3.run_map(genome_index, inputs=[fq], args=["--mapping-mode", mode], threads=threads)
	sam = res.sam()
	assert [r.qname for r in sam.records] == [r.name for r in unique_reads]
	for read, rec in zip(unique_reads, sam.records):
		t = read.truth
		assert not rec.is_unmapped, (mode, read, rec)
		assert (rec.rname, rec.pos, rec.strand) == (t.contig, t.pos + 1, t.strand), (mode, read, rec)
		assert rec.seq == _forward_seq(read)
		assert cigar_is_valid(rec.cigar_ops, 100)
		assert check_record_against_reference(rec, genome[rec.rname]) == [], rec
		assert rec.tags["NM"] <= t.nm, (read, rec)
		ref = genome[rec.rname]
		score = area_b.alignment_score(ref, rec.ref_start, rec.cigar_ops, rec.seq)
		assert rec.tags["AS"] == score, rec
		assert score == _best_score(rec.seq, ref, rec.ref_start, rec.ref_end), rec


@pytest.mark.parametrize("a,b", [("sensitive", "hybrid"), ("complete", "customed")])
def test_mode_aliases_are_identical(gem3, genome, genome_index, tmp_path, threads, a, b):
	"""--help lists 'sensitive|hybrid' and 'customed|complete' as one mode each."""
	reads = simulate_se(genome, 200, 100, seed=64, sub_rate=0.02, ins_rate=0.004, del_rate=0.004,
		unique_only=False, spacing=4)
	fq = write_reads(tmp_path / "r.fq", reads)
	from gem3harness import strip_variable_header
	out = [strip_variable_header(gem3.run_map(genome_index, inputs=[fq], args=["--mapping-mode", m],
		threads=threads).stdout) for m in (a, b)]
	assert out[0] == out[1]


# --------------------------------------------------------------------------
# Planted genome: copies of one probe at known edit distances
# --------------------------------------------------------------------------

@pytest.fixture(scope="module")
def planted(gem3):
	genome, probe, copies = area_b.planted_genome()
	index = gem3.index(genome)
	edit = area_b.loci_within(genome, probe, 7, "edit")
	ham = area_b.loci_within(genome, probe, 7, "hamming")
	return genome, probe, copies, index, {"edit": edit, "hamming": ham}


def test_planted_genome_matches_its_design(planted):
	"""The brute-force search finds exactly the planted copies (the random
	background has nothing within 7 errors), each at its planted distance."""
	genome, probe, copies, index, loci = planted
	edit = sorted((c, s, d) for c, s, e, d in loci["edit"])
	assert edit == sorted((c.contig, c.strand, c.edits) for c in copies)
	ham = sorted((c, s, start, d) for c, s, start, d in loci["hamming"])
	assert ham == sorted((c.contig, c.strand, c.start, c.subs) for c in copies if c.indels == 0)


def _reported(res, model):
	out = set()
	for rec in area_b.mapped(res.sam().records):
		coord = rec.ref_end - 1 if model == "edit" else rec.ref_start
		out.add((rec.rname, rec.strand, coord, rec.tags["NM"]))
	return out


def _expected(loci, model, max_error):
	return set(x for x in loci[model] if x[3] <= max_error)


def _map_probe(gem3, index, probe, args, threads):
	data = area_b.fastq_text([("probe", probe)])
	return gem3.run_map(index, args=list(args) + ["-i", "probe.fq"], files={"probe.fq": data},
		threads=threads, timeout=60)


@pytest.mark.parametrize("model", ["edit", "hamming"])
@pytest.mark.parametrize("mode", ["complete-brute-force", "complete-partition"])
def test_neighborhood_search_reports_exactly_loci_within_E(gem3, planted, threads, mode, model):
	"""The pure neighborhood searches report every locus within -E errors
	(edit distance, or mismatches under --alignment-model hamming) and
	nothing else, each with NM equal to its distance; -E given as a count or
	as a fraction of the read length is the same search."""
	genome, probe, copies, index, loci = planted
	am = ["--alignment-model", "hamming"] if model == "hamming" else []
	for E in (0, 1, 2, 3, 4, 6):
		res = _map_probe(gem3, index, probe, ["--mapping-mode", mode, "-E", str(E)] + am
			+ ALL_MATCHES + NO_LOCAL, threads)
		assert _reported(res, model) == _expected(loci, model, E), (mode, model, E)
		for rec in area_b.mapped(area_b.fill_secondary_seq(res.sam().records)):
			assert check_record_against_reference(rec, genome[rec.rname]) == [], rec
			if model == "hamming":
				assert rec.cigar == "100M", rec
		if E in (2, 4):
			frac = _map_probe(gem3, index, probe, ["--mapping-mode", mode, "-E", "0.0%d" % E] + am
				+ ALL_MATCHES + NO_LOCAL, threads)
			assert _reported(frac, model) == _reported(res, model)


@pytest.mark.parametrize("model", ["edit", "hamming"])
@pytest.mark.parametrize("mode", ["fast", "sensitive", "complete", "complete-pure", "complete-partition-filter"])
def test_search_finds_every_locus_within_E(gem3, planted, threads, mode, model):
	"""The filtering and hybrid modes report at least every locus within -E
	errors (they may report more, up to --alignment-max-error), and every
	reported alignment is consistent with the reference."""
	genome, probe, copies, index, loci = planted
	am = ["--alignment-model", "hamming"] if model == "hamming" else []
	for E in (0, 2, 4):
		res = _map_probe(gem3, index, probe, ["--mapping-mode", mode, "-E", str(E)] + am
			+ ALL_MATCHES + NO_LOCAL, threads)
		got = _reported(res, model)
		assert _expected(loci, model, E) <= got, (mode, model, E, got)
		for rec in area_b.mapped(area_b.fill_secondary_seq(res.sam().records)):
			assert check_record_against_reference(rec, genome[rec.rname]) == [], rec
		# The primary is the exact copy
		primary = [r for r in res.sam().records if r.is_primary]
		assert len(primary) == 1 and primary[0].tags["NM"] == 0, primary


def test_brute_force_matches_partition_on_reads_with_edits(gem3, planted, threads):
	"""complete-brute-force and complete-partition are two complete searches
	of the same space: the same records on reads carrying edits, except for
	MAPQ (the MAPQ model reads search metrics that differ between them)."""
	genome, probe, copies, index, loci = planted
	rng = area_b.Rng(5)
	reads = []
	for i in range(6):
		s = list(probe)
		for j in range(i % 3):
			q = 20 + 25 * j
			s[q] = rng.choice([b for b in "ACGT" if b != s[q]])
		reads.append(("q%d" % i, "".join(s) if i < 3 else revcomp("".join(s))))
	files = {"q.fq": area_b.fastq_text(reads)}
	outs = []
	for mode in ("complete-brute-force", "complete-partition"):
		res = gem3.run_map(index, args=["-i", "q.fq", "--mapping-mode", mode, "-E", "3"] + ALL_MATCHES
			+ NO_LOCAL, files=files, threads=threads, timeout=60)
		lines = [l.split("\t") for l in res.stdout_text.split("\n") if l and not l.startswith("@")]
		outs.append([f[:4] + f[5:] for f in lines])
	assert len(outs[0]) >= 6 * 4
	assert outs[0] == outs[1]


# --------------------------------------------------------------------------
# exact mode and alignment models
# --------------------------------------------------------------------------

def test_exact_mode_maps_only_exact_reads(gem3, genome, genome_index, tmp_path, threads):
	"""--mapping-mode exact: reads identical to the reference map with 100M
	and NM 0 at their origin; a read with one substitution is unmapped."""
	exact = simulate_se(genome, 40, 100, seed=65)
	one = simulate_se(genome, 40, 100, seed=66, subs=1, name_prefix="s")
	fq = write_reads(tmp_path / "r.fq", exact + one)
	sam = gem3.run_map(genome_index, inputs=[fq], args=["--mapping-mode", "exact"], threads=threads).sam()
	by = sam.by_name()
	for read in exact:
		(rec,) = by[read.name]
		assert (rec.rname, rec.pos, rec.strand, rec.cigar, rec.tags["NM"]) == (
			read.truth.contig, read.truth.pos + 1, read.truth.strand, "100M", 0), rec
	for read in one:
		(rec,) = by[read.name]
		assert rec.is_unmapped, rec


def test_hamming_model_aligns_without_gaps(gem3, genome, genome_index, tmp_path, threads):
	"""--alignment-model hamming (no local fallback): substitution-only reads
	map at their origin with 100M and NM = mismatches; reads with an indel
	are either unmapped or placed with 100M and NM = the mismatches of that
	gapless placement."""
	subs = simulate_se(genome, 50, 100, seed=67, subs=3)
	indel = simulate_se(genome, 30, 100, seed=68, subs=1, ins=1, name_prefix="i")
	fq = write_reads(tmp_path / "r.fq", subs + indel)
	sam = gem3.run_map(genome_index, inputs=[fq], args=["--alignment-model", "hamming"] + NO_LOCAL,
		threads=threads).sam()
	by = sam.by_name()
	for read in subs:
		(rec,) = by[read.name]
		assert (rec.rname, rec.pos, rec.strand, rec.cigar, rec.tags["NM"]) == (
			read.truth.contig, read.truth.pos + 1, read.truth.strand, "100M", 3), rec
	for rec in area_b.mapped(sam.records):
		assert rec.cigar == "100M", rec
		assert check_record_against_reference(rec, genome[rec.rname]) == [], rec
		# AS is the mismatch count under this model
		assert rec.tags["AS"] == rec.tags["NM"], rec


def test_edit_model_reports_least_edit_distance(gem3, genome, genome_index, tmp_path, threads):
	"""--alignment-model edit: reads map at their origin; NM (and AS) equal
	the least edit distance of the read against the reference around it."""
	reads = simulate_se(genome, 50, 100, seed=69, subs=2, ins=1, dels=1, max_indel=3)
	fq = write_reads(tmp_path / "r.fq", reads)
	sam = gem3.run_map(genome_index, inputs=[fq], args=["--alignment-model", "edit"] + NO_LOCAL,
		threads=threads).sam()
	for read, rec in zip(reads, sam.records):
		t = read.truth
		assert (rec.rname, rec.strand) == (t.contig, t.strand), (read, rec)
		assert abs(rec.pos - (t.pos + 1)) <= 3, (read, rec)
		assert check_record_against_reference(rec, genome[rec.rname]) == [], rec
		ref = genome[rec.rname]
		best = area_b.edit_distance_best(rec.seq, ref[max(0, t.pos - 10):t.end + 10])
		assert rec.tags["NM"] == best == rec.tags["AS"], (read, rec, best)


def test_pseudoalignment_places_reads_near_origin(gem3, genome, genome_index, tmp_path, threads):
	"""--alignment-model none|pseudoalignment: reads are placed on their
	strand within the error budget of their origin, with a CIGAR covering
	the read. (Its NM and MD do not describe the reference:
	known-bugs/B-pseudoalignment-nm-md.md.)"""
	reads = simulate_se(genome, 40, 100, seed=70, subs=2, dels=1, max_indel=3)
	fq = write_reads(tmp_path / "r.fq", reads)
	for model in ("none", "pseudoalignment"):
		sam = gem3.run_map(genome_index, inputs=[fq], args=["--alignment-model", model] + NO_LOCAL,
			threads=threads).sam()
		for read, rec in zip(reads, sam.records):
			t = read.truth
			assert (rec.rname, rec.strand) == (t.contig, t.strand), (read, rec)
			assert abs(rec.pos - (t.pos + 1)) <= 12, (read, rec)
			assert cigar_is_valid(rec.cigar_ops, 100), rec
