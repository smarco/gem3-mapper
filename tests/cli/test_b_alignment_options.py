# GEM-Mapper v3 (GEM3) test suite.
# Copyright (c) 2026 Jacob Schreiber <jacob.schreiber@umassmed.edu>
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of GEM-Mapper v3 (GEM3), distributed under the GNU General
# Public License version 3 or (at your option) any later version; see LICENSE.
"""Alignment thresholds, local alignment, gap-affine scoring, banding,
end curation, candidate verification, the q-gram filter and -c.

Expected values are derived from the option's meaning in --help and the
read's construction: a 100 bp read with k spaced substitutions has k edits,
95..100-k identity and gap-affine score A*(100-k) - B*k; a value below 1 is
a fraction of the read length, 1 or more an absolute number. AS is checked
against a brute-force Gotoh alignment (area_b.gotoh_best) with the same
penalties. Options whose behaviour contradicts --help are known bugs and not
tested here: -B and --gap-affine-penalties
(known-bugs/B-mismatch-penalty-sign.md), --alignment-curation
(B-alignment-curation-ignored.md), the --help defaults of
--alignment-global-min-score and --alignment-curation-min-end-context
(B-help-defaults.md), -c other than 'correct'/'none' (B-check-alignments-crash.md)
and --qgram-filter with 0 tiles (B-qgram-filter-zero-tiles.md).
"""

import pytest

from gem3harness import (check_record_against_reference, simulate_se, strip_variable_header,
	write_reads)
import area_b

NO_LOCAL = ["--alignment-local", "never"]


@pytest.fixture(scope="module")
def sub_reads(genome):
	"""{k: reads with exactly k substitutions}, k = 0..9, at least 6 bp apart
	so every read keeps an exact seed of 10 bp or more."""
	return dict((k, simulate_se(genome, 8, 100, seed=200 + k, subs=k, spacing=6, edge=6,
		name_prefix="k%d_" % k)) for k in range(10))


def _mapped_names(res):
	return set(r.qname for r in res.sam().records if not r.is_unmapped)


def _map(gem3, index, reads, tmp_path, args, threads, name="r.fq"):
	fq = write_reads(tmp_path / name, reads)
	return gem3.run_map(index, inputs=[fq], args=args, threads=threads)


@pytest.mark.parametrize("model", ["gap-affine", "hamming"])
@pytest.mark.parametrize("value,nominal", [("0", 0), ("3", 3), ("5", 5), ("8", 8), ("0.05", 5), ("0.08", 8)])
def test_max_error_bounds_mapped_reads(gem3, genome, genome_index, tmp_path, threads, sub_reads, value, nominal, model):
	"""-e: a read with k substitutions maps (at its origin, NM k) exactly
	when k <= the error budget, under the gap-affine and the Hamming model;
	a fraction is of the read length."""
	reads = [r for k in sorted(sub_reads) for r in sub_reads[k]]
	res = _map(gem3, genome_index, reads, tmp_path, ["-e", value, "--alignment-model", model] + NO_LOCAL, threads)
	by = res.sam().by_name()
	for k in sorted(sub_reads):
		for read in sub_reads[k]:
			(rec,) = by[read.name]
			if k <= nominal:
				assert (rec.rname, rec.pos, rec.strand, rec.tags.get("NM")) == (
					read.truth.contig, read.truth.pos + 1, read.truth.strand, k), (value, read, rec)
			else:
				assert rec.is_unmapped, (value, read, rec)


def test_max_error_count_and_fraction_agree(gem3, genome, genome_index, tmp_path, threads, sub_reads):
	"""-e 5 and -e 0.05 are the same budget on 100 bp reads: same bytes."""
	reads = [r for k in sorted(sub_reads) for r in sub_reads[k]]
	outs = [strip_variable_header(_map(gem3, genome_index, reads, tmp_path, ["-e", v], threads).stdout)
		for v in ("5", "0.05")]
	assert outs[0] == outs[1]


@pytest.mark.parametrize("args,max_k", [
	(["--alignment-global-min-identity", "0.95"], 5),   # 100-k >= 95
	(["--alignment-global-min-identity", "0.9"], 9),    # 100-k >= 90 (k <= 9 here)
	(["--alignment-global-min-identity", "93"], 7),     # absolute: 100-k >= 93
	(["--alignment-global-min-score", "0.9"], 2),       # 100-5k >= 90
	(["--alignment-global-min-score", "0.7"], 6),       # 100-5k >= 70
	(["--alignment-global-min-score", "80"], 4),        # absolute: 100-5k >= 80
	(["--alignment-global-min-score", "0.9", "-A", "2"], 3),  # 2(100-k)-4k >= 180
	(["--alignment-global-min-score", "0.9", "-A", "3"], 4),  # 3(100-k)-4k >= 270
])
def test_global_identity_and_score_thresholds(gem3, genome, genome_index, tmp_path, threads, sub_reads, args, max_k):
	"""--alignment-global-min-identity (matching bases) and
	--alignment-global-min-score (gap-affine score; a fraction is of the
	read length times the match score): reads below either are not
	reported as global alignments (unmapped without local alignment)."""
	reads = [r for k in sorted(sub_reads) for r in sub_reads[k]]
	res = _map(gem3, genome_index, reads, tmp_path, args + ["-e", "12"] + NO_LOCAL, threads)
	got = _mapped_names(res)
	want = set(r.name for k in sorted(sub_reads) if k <= max_k for r in sub_reads[k])
	assert got == want, (args, sorted(got ^ want))


# --------------------------------------------------------------------------
# Local alignment
# --------------------------------------------------------------------------

@pytest.fixture(scope="module")
def chimera_reads(genome):
	return area_b.chimeras(genome, 32, [30, 45, 55, 65], seed=9)


def _check_local(rec, info, genome):
	contig, start, L, left = info
	clip = 100 - L
	assert (rec.rname, rec.strand) == (contig, "+"), rec
	assert rec.cigar == ("%dM%dS" % (L, clip) if left else "%dS%dM" % (clip, L)), rec
	assert rec.pos == start + 1, rec
	assert rec.tags["NM"] == 0 and rec.tags["AS"] == L, rec
	assert check_record_against_reference(rec, genome[rec.rname]) == [], rec


@pytest.mark.parametrize("args,min_genomic", [
	([], 40),                                                  # defaults: identity 40, score 20
	(["--alignment-local", "if-unmapped"], 40),
	(["--alignment-local-min-identity", "50"], 50),
	(["--alignment-local-min-identity", "0.5"], 50),
	(["--alignment-local-min-identity", "60"], 60),
	(["--alignment-local-min-score", "50"], 50),
	(["--alignment-local-min-score", "0.6"], 60),
	(["--alignment-local", "never"], 101),
])
def test_local_alignment_of_chimeras(gem3, genome, genome_index, tmp_path, threads, chimera_reads, args, min_genomic):
	"""Reads whose first or last L bases come from the genome (L <= 65, so
	no global alignment passes the 80% identity default): with local alignment
	they map at the genomic part with LM/S CIGARs, AS = L and NM 0 when L
	meets --alignment-local-min-identity (matching bases) and
	--alignment-local-min-score; otherwise, or with --alignment-local never,
	they are unmapped."""
	files = {"c.fq": area_b.fastq_text(chimera_reads)}
	res = gem3.run_map(genome_index, args=["-i", "c.fq"] + args, files=files, threads=threads)
	by = res.sam().by_name()
	for name, seq, info in chimera_reads:
		(rec,) = by[name]
		if info[2] >= min_genomic:
			_check_local(rec, info, genome)
		else:
			assert rec.is_unmapped, (args, info, rec)


# --------------------------------------------------------------------------
# Gap-affine scoring and the aligner
# --------------------------------------------------------------------------

@pytest.fixture(scope="module")
def indel_reads(genome):
	reads = simulate_se(genome, 25, 100, seed=210, subs=2, ins=1, dels=1, max_indel=4)
	reads += simulate_se(genome, 15, 100, seed=211, subs=1, dels=1, max_indel=4, edge=12, name_prefix="d")
	return reads


@pytest.mark.parametrize("args,pen", [
	([], (1, 4, 6, 1)),
	(["-A", "2"], (2, 4, 6, 1)),
	(["--matching-score", "3"], (3, 4, 6, 1)),
	(["-O", "2"], (1, 4, 2, 1)),
	(["--gap-open-penalty", "15"], (1, 4, 15, 1)),
	(["-X", "3"], (1, 4, 6, 3)),
	(["--gap-extension-penalty", "2", "-O", "0"], (1, 4, 0, 2)),
	(["--alignment-force-full-swg"], (1, 4, 6, 1)),
	(["--alignment-force-full-swg=true", "-O", "10"], (1, 4, 10, 1)),
])
def test_alignment_is_optimal_under_penalties(gem3, genome, genome_index, tmp_path, threads, indel_reads, args, pen):
	"""-A/-O/-X (and --alignment-force-full-swg): every read maps at its
	origin; AS is the gap-affine score of the reported alignment under the
	given penalties and equals the best score any alignment of the read
	achieves in a window around it (brute-force Gotoh)."""
	res = _map(gem3, genome_index, indel_reads, tmp_path, args + NO_LOCAL, threads)
	for read, rec in zip(indel_reads, res.sam().records):
		t = read.truth
		assert not rec.is_unmapped, (args, read)
		assert (rec.rname, rec.strand) == (t.contig, t.strand), (args, read, rec)
		assert abs(rec.pos - (t.pos + 1)) <= 4, (args, read, rec)
		ref = genome[rec.rname]
		assert check_record_against_reference(rec, ref) == [], rec
		score = area_b.alignment_score(ref, rec.ref_start, rec.cigar_ops, rec.seq, *pen)
		assert rec.tags["AS"] == score, (args, rec, score)
		best = area_b.gotoh_best(rec.seq, ref[max(0, t.pos - 15):t.end + 15], *pen)
		assert score == best, (args, read, rec, best)


def _deletion_reads(genome, sizes, seed):
	rng = area_b.Rng(seed)
	ref = genome["chr1"]
	out = []
	while len(out) < 2 * len(sizes):
		n = sizes[len(out) % len(sizes)]
		p = rng.randint(1000, len(ref) - 1000)
		if genome.hard_regions("chr1", p - 50, p + 150 + n):
			continue
		seq = (ref[p:p + 50] + ref[p + 50 + n:p + 100 + n]).upper()
		if "N" in seq:
			continue
		out.append(("del%d_%d" % (n, len(out)), seq, (p, n)))
	return out


@pytest.mark.parametrize("value,band", [("0.1", 10), ("0.2", 20), (None, 20), ("30", 30)])
def test_max_bandwidth_limits_gap_length(gem3, genome, genome_index, tmp_path, threads, value, band):
	"""--alignment-max-bandwidth: a read missing n reference bases in its
	middle aligns end to end with one n-base deletion when n is within the
	band; otherwise it gets no such alignment (here: a local alignment of one
	half)."""
	reads = _deletion_reads(genome, [8, 15, 25], seed=12)
	args = ["--alignment-max-bandwidth", value] if value else []
	res = gem3.run_map(genome_index, args=["-i", "d.fq"] + args, files={"d.fq": area_b.fastq_text(reads)},
		threads=threads)
	by = res.sam().by_name()
	for name, seq, (p, n) in reads:
		(rec,) = by[name]
		ops = rec.cigar_ops
		gapped = [op for op in ops if op[0] == "D"] == [("D", n)] and "S" not in rec.cigar
		if n <= band:
			assert gapped and rec.pos == p + 1 and rec.tags["AS"] == 100 - 6 - n, (value, n, rec)
			assert check_record_against_reference(rec, genome["chr1"]) == [], rec
		else:
			assert not gapped, (value, n, rec)


def _end_mismatch_reads(genome, offsets, seed):
	rng = area_b.Rng(seed)
	ref = genome["chr1"]
	out = []
	for off in offsets:
		for side in "LR":
			while True:
				p = rng.randint(1000, len(ref) - 1000)
				if not genome.hard_regions("chr1", p - 20, p + 120):
					break
			s = list(ref[p:p + 100].upper())
			q = off if side == "L" else 99 - off
			s[q] = rng.choice([b for b in "ACGT" if b != s[q]])
			out.append(("end%s%d" % (side, off), "".join(s), (p, side, off)))
	return out


@pytest.mark.parametrize("value,context", [(None, 0), ("3", 3), ("6", 6), ("0.05", 5)])
def test_curation_clips_mismatches_near_ends(gem3, genome, genome_index, tmp_path, threads, value, context):
	"""--alignment-curation-min-end-context c: a mismatch fewer than c bases
	from a read end is soft-clipped with the bases beyond it (AS counts only
	the aligned matches); one at c bases or more stays a mismatch (AS
	100 - 5). The --help default is not the code's
	(known-bugs/B-help-defaults.md); the code's is 0."""
	reads = _end_mismatch_reads(genome, [1, 2, 3, 5, 8], seed=13)
	args = ["--alignment-curation-min-end-context", value] if value else []
	res = gem3.run_map(genome_index, args=["-i", "e.fq"] + args, files={"e.fq": area_b.fastq_text(reads)},
		threads=threads)
	by = res.sam().by_name()
	for name, seq, (p, side, off) in reads:
		(rec,) = by[name]
		if off < context:
			clip = off + 1
			want = ("%dS%dM" % (clip, 100 - clip)) if side == "L" else ("%dM%dS" % (100 - clip, clip))
			pos = p + 1 + (clip if side == "L" else 0)
			assert (rec.cigar, rec.pos, rec.tags["AS"], rec.tags["NM"]) == (want, pos, 100 - clip, 0), (value, name, rec)
		else:
			assert (rec.cigar, rec.pos, rec.tags["AS"], rec.tags["NM"]) == ("100M", p + 1, 95, 1), (value, name, rec)
		assert check_record_against_reference(rec, genome["chr1"]) == [], rec


# --------------------------------------------------------------------------
# Candidate verification, q-gram filter, check
# --------------------------------------------------------------------------

@pytest.fixture(scope="module")
def mixed_fastq(genome, chimera_reads):
	"""Unique and repeat reads with substitutions and indels, chimeras,
	N-rich reads and random reads."""
	from gem3harness import fastq_record
	rng = area_b.Rng(17)
	reads = simulate_se(genome, 150, 100, seed=220, sub_rate=0.02, ins_rate=0.004, del_rate=0.004,
		unique_only=False, spacing=4)
	text = "".join(fastq_record(r.name, r.seq, r.qual) for r in reads)
	text += area_b.fastq_text(chimera_reads[:12]).decode()
	for i in range(8):
		s = list(genome["chr2"][3000 + 5000 * i:3100 + 5000 * i].upper())
		for k in range(3, 100, 19):
			s[k] = "N"
		text += fastq_record("nrich%d" % i, "".join(s), "I" * 100)
	for i in range(6):
		text += fastq_record("rand%d" % i, rng.dna(100), "I" * 100)
	return text.encode("ascii")


def test_candidate_verification_strategies(gem3, genome, genome_index, threads, mixed_fastq):
	"""--candidate-verification BPM and chained (the default) verify the
	same candidates: identical output, and every mapped record is consistent
	with the reference."""
	outs = {}
	for v in (None, "BPM", "chained"):
		args = ["--candidate-verification", v] if v else []
		res = gem3.run_map(genome_index, args=["-i", "m.fq"] + args, files={"m.fq": mixed_fastq}, threads=threads)
		outs[v] = strip_variable_header(res.stdout)
		for rec in area_b.mapped(res.sam().records):
			assert check_record_against_reference(rec, genome[rec.rname]) == [], (v, rec)
	assert outs["BPM"] == outs["chained"] == outs[None]


def test_qgram_filter_never_changes_output(gem3, genome, genome_index, threads, mixed_fastq):
	"""--qgram-filter <tiles>,<k>: the k-mer filter only drops candidates
	that cannot be within the error budget, so any valid setting leaves the
	output as it is. A k outside 3..7 is replaced by 7 with a warning."""
	base = strip_variable_header(gem3.run_map(genome_index, args=["-i", "m.fq"],
		files={"m.fq": mixed_fastq}, threads=threads).stdout)
	for q in ("1,4", "2,6", "3,5", "4,7", "5,3", "6,6"):
		res = gem3.run_map(genome_index, args=["-i", "m.fq", "--qgram-filter", q], files={"m.fq": mixed_fastq},
			threads=threads)
		assert strip_variable_header(res.stdout) == base, q
	res9 = gem3.run_map(genome_index, args=["-i", "m.fq", "--qgram-filter", "2,9"],
		files={"m.fq": mixed_fastq}, threads=threads)
	assert "Invalid proposed k-mer length. Set to k=7" in res9.stderr_text
	assert strip_variable_header(res9.stdout) == base


@pytest.mark.parametrize("value", ["correct", "none"])
def test_check_alignments_correct(gem3, genome, genome_index, threads, mixed_fastq, value):
	"""-c correct re-checks every reported CIGAR against the reference and
	prints a 'Match check failed' block for a wrong one: none here, and the
	output is unchanged. -c none checks nothing."""
	base = gem3.run_map(genome_index, args=["-i", "m.fq"], files={"m.fq": mixed_fastq}, threads=threads)
	res = gem3.run_map(genome_index, args=["-i", "m.fq", "-c", value], files={"m.fq": mixed_fastq},
		threads=threads)
	assert "Match check failed" not in res.stderr_text
	assert strip_variable_header(res.stdout) == strip_variable_header(base.stdout)


def test_check_alignments_correct_paired(gem3, genome, genome_index, tmp_path, threads):
	"""-c correct on paired-end input checks both ends' alignments: no
	failure report, same output as without -c."""
	from gem3harness import simulate_pe, write_split
	pairs = simulate_pe(genome, 60, 100, insert_mean=300, insert_sd=30, seed=230, sub_rate=0.02,
		ins_rate=0.004, del_rate=0.004, unique_only=False, spacing=4)
	fq1, fq2 = write_split(tmp_path / "1.fq", tmp_path / "2.fq", pairs)
	base = gem3.run_map(genome_index, inputs=[fq1, fq2], threads=threads)
	res = gem3.run_map(genome_index, inputs=[fq1, fq2], args=["-c", "correct"], threads=threads)
	assert "Match check failed" not in res.stderr_text
	assert strip_variable_header(res.stdout) == strip_variable_header(base.stdout)
	assert len([r for r in res.sam().records if not r.is_unmapped]) >= 100
