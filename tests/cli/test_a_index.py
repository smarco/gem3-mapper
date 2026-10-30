# GEM-Mapper v3 (GEM3) test suite.
# Copyright (c) 2026 Jacob Schreiber <jacob.schreiber@umassmed.edu>
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of GEM-Mapper v3 (GEM3), distributed under the GNU General
# Public License version 3 or (at your option) any later version; see LICENSE.
"""Area A: loading the GEM index (-I) and using it to place reads.

Expected values come from the FASTA: the indexed text and locator
intervals rebuilt by ``area_a.indexed_text``, k-mer suffix-array intervals
counted by brute force, and reads cut from known positions (contig ends,
next to N runs, short contigs, IUPAC codes in the reference).
"""

import os
import re

import pytest

from gem3harness import (Rng, check_record_against_reference, revcomp, sam_digest, simulate_pe,
	simulate_se, synthetic_genome, write_fasta_records, write_reads, write_split)
from area_a import assert_maps_to_truth, indexed_text, kmer_intervals, primary

TINY = [
	# 7 N (kept), lowercase bases, an IUPAC code, a 70 N run (removed) and
	# a contig that starts with 55 N (removed)
	("c1", "ACGTACGTTAGCATCGATCGACTAGCTAGCTAGCATCGATCAGCTACGATCGACTGACTGATCGATGCATGCATCGATCGATGCATGC"
		"NNNNNNNacgtagctagctagctgactgatcgatcgtagctagctagcatcgRTTGCA"),
	("c2", "TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA" + "N" * 70 + "ACGATCGATCGATCGTAGCTAGCTAGTCGATCGATCGA"),
	("c3", "N" * 55 + "GATTACAGATTACACCGGTTAACGGTTCAGT"),
]


def _dev_output(gem3, tmp_path, records, extra_index_args=()):
	fa = write_fasta_records(str(tmp_path / "tiny.fa"), records)
	index = gem3.build_index(fa, extra_args=extra_index_args)
	fq = write_reads(tmp_path / "r.fq", [])
	res = gem3.run_map(index, inputs=[fq], args=["--verbose=dev"], threads=1, quiet=False)
	return res.stderr_text


def test_verbose_dev_prints_the_loaded_index(gem3, tmp_path):
	"""--verbose=dev prints the archive as loaded: the locator intervals,
	the text lengths, the sampling and the rank table's k-mer intervals,
	which must equal the ones counted on the text rebuilt from the FASTA."""
	err = _dev_output(gem3, tmp_path, TINY)
	full, fwd_length, intervals = indexed_text(TINY)
	assert "Index.Type GEM.FM-Index.DNA.Classic.full" in err
	assert re.search(r"Ns\.Threshold 50\n", err)
	assert re.search(r"Archive\.Forward\.Length +%d\n" % fwd_length, err)
	assert re.search(r"Archive\.RL\.Text +No\n", err)
	# Locator: forward intervals, uncalled N runs and reverse-complement intervals
	printed = re.findall(r"\[\d+\] (\w+)\tIdx=\[(\d+),(\d+)\)\tRL\.Idx=\[\d+,\d+\)\tText=([^:\t]+):([+-]):\[(\d+),(\d+)\)", err)
	got = [(k, int(b), int(e), name, st, int(tb), int(te)) for k, b, e, name, st, tb, te in printed]
	assert got == intervals
	assert re.search(r"Intervals\.Size \d+ MB \(%d intervals" % len(intervals), err)
	# One SA sample every 4 text positions (gem-indexer's CPU default)
	assert re.search(r"Text\.Sampling\.Rate 2\^2 \(4\)", err)
	assert re.search(r"PackedArray\.Length %d\n" % ((len(full) + 3) // 4), err)
	# Rank table: every k-mer of length 1-3 with its suffix-array interval
	content = re.findall(r"\[([ACGT]+)\] (\d+) (\d+)\n", err)
	assert len(content) == 4 + 16 + 64
	intervals_by_k = {k: kmer_intervals(full, k) for k in (1, 2, 3)}
	for kmer, lo, hi in content:
		assert (int(lo), int(hi)) == intervals_by_k[len(kmer)](kmer), kmer


def test_verbose_dev_run_length_index(gem3, tmp_path):
	"""A run-length index (gem-indexer -r) is reported as such, with the RL
	sampling of text positions."""
	err = _dev_output(gem3, tmp_path, TINY, extra_index_args=("-r",))
	assert "Index.Type GEM.FM-Index.DNA.RunLength.full" in err
	assert re.search(r"Archive\.RL\.Text +Yes\n", err)
	assert "[GEM]>Sampled-RL" in err


@pytest.fixture(scope="module")
def mixed_reads(genome):
	"""Single-end reads with errors, including repeats, tandems and
	homopolymers (unique_only=False), and pairs."""
	se = simulate_se(genome, 400, 100, seed=91, sub_rate=0.01, ins_rate=0.002, del_rate=0.002, unique_only=False)
	pe = simulate_pe(genome, 150, 100, insert_mean=300, insert_sd=40, seed=92, sub_rate=0.01, unique_only=False)
	return se, pe


@pytest.mark.parametrize("sampling", [("-s", "8"), ("-s", "32"), ("-S", "8"), ("-S", "4", "-s", "16")],
	ids=["text8", "text32", "sa8", "sa4-text16"])
def test_index_sampling_does_not_change_output(gem3, genome, genome_index, mixed_reads, tmp_path, threads, sampling):
	"""Suffix-array positions are decoded from samples taken every s text
	positions (-s) or every S suffix-array positions (-S); the decoded
	positions, and so the SAM output, do not depend on the sampling."""
	se, pe = mixed_reads
	fa = os.path.join(gem3.work_dir, "genome-%s.fa" % genome.digest()[:16])
	index = gem3.build_index(fa, extra_args=sampling, threads=4)
	fq = write_reads(tmp_path / "se.fq", se)
	fq1, fq2 = write_split(tmp_path / "1.fq", tmp_path / "2.fq", pe)
	for inputs in ([fq], [fq1, fq2]):
		base = gem3.run_map(genome_index, inputs=inputs, threads=threads)
		other = gem3.run_map(index, inputs=inputs, threads=threads)
		assert sam_digest(other.stdout) == sam_digest(base.stdout)


# sha256 of the SAM output (without @PG/@CO) on the run-length index, base 6e155e3
RL_PINNED = "85979448324f4edebd0669fb57ec37ec037be25035ca6c593da48086331c31d5"


def test_run_length_index(gem3, genome, mixed_reads, tmp_path, threads):
	"""Reads map to their origin on a run-length index; the output is
	pinned (pins the current output), and paired-end mode is refused."""
	se, pe = mixed_reads
	fa = os.path.join(gem3.work_dir, "genome-%s.fa" % genome.digest()[:16])
	index = gem3.build_index(fa, extra_args=("-r",), threads=4)
	unique = simulate_se(genome, 200, 100, seed=93, subs=2)
	fq = write_reads(tmp_path / "u.fq", unique)
	assert_maps_to_truth(gem3.run_map(index, inputs=[fq], threads=threads).sam(), unique, genome)
	fq = write_reads(tmp_path / "se.fq", se)
	res = gem3.run_map(index, inputs=[fq], threads=threads)
	for rec in res.sam().records:
		if not rec.is_unmapped:
			assert check_record_against_reference(rec, genome[rec.rname]) == [], rec
	assert sam_digest(res.stdout) == RL_PINNED
	fq1, fq2 = write_split(tmp_path / "1.fq", tmp_path / "2.fq", pe[:10])
	res = gem3.run_map(index, inputs=[fq1, fq2], threads=threads, check=False)
	assert res.returncode != 0 and not res.crashed
	assert "Archive RL-text not supported for paired-end mode" in res.stderr_text


def _edge_reads(genome, rng):
	"""Exact reads at the first and last bases of every contig, against both
	sides of every N run, covering whole short contigs; both strands."""
	out = []

	def add(contig, start, length, label):
		seq = genome[contig][start:start + length].upper()
		if len(seq) != length or "N" in seq:
			return
		for strand in "+-":
			name = "%s_%s_%d_%s" % (label, contig, start, "f" if strand == "+" else "r")
			out.append((name, seq if strand == "+" else revcomp(seq), contig, start, strand))

	for contig in genome:
		L = len(genome[contig])
		seq = genome[contig].upper()
		first = next((i for i, c in enumerate(seq) if c != "N"), None)
		last = max((i for i, c in enumerate(seq) if c != "N"), default=None)
		for length in (100, 60, 36):
			if first is not None:
				add(contig, first, length, "start%d" % length)
			if last is not None:
				add(contig, last + 1 - length, length, "end%d" % length)
		if L <= 100:
			add(contig, 0, L, "whole")
	for f in genome.features:
		if f.kind != "N":
			continue
		for length in (100, 50):
			add(f.contig, f.start - length, length, "beforeN%d" % length)
			add(f.contig, f.end, length, "afterN%d" % length)
	return out


def test_reads_at_contig_ends_and_n_runs(gem3, genome, genome_index, tmp_path, threads):
	"""Positions next to the ends of the indexed intervals (contig ends,
	removed N runs, the 90 bp contig) translate back to the right contig,
	position and strand."""
	reads = _edge_reads(genome, Rng(5))
	assert len(reads) > 40
	fq = tmp_path / "edge.fq"
	with open(str(fq), "w") as fh:
		for name, seq, _, _, _ in reads:
			fh.write("@%s\n%s\n+\n%s\n" % (name, seq, "I" * len(seq)))
	sam = gem3.run_map(genome_index, inputs=[str(fq)], threads=threads).sam()
	by_name = primary(sam)
	for name, seq, contig, start, strand in reads:
		rec = by_name[name]
		assert not rec.is_unmapped, name
		assert check_record_against_reference(rec, genome[rec.rname]) == [], rec
		# Reads from the exact repeat copies or tandems may have equal placements
		if genome.hard_regions(contig, start, start + len(seq), kinds=("repeat", "tandem", "homopolymer")):
			continue
		assert (rec.rname, rec.pos, rec.strand, rec.cigar) == (contig, start + 1, strand, "%dM" % len(seq)), (name, rec)
		assert rec.tags["NM"] == 0


def test_reads_do_not_map_across_contig_boundaries(gem3, genome, genome_index, tmp_path, threads):
	"""A read made of the end of one contig and the start of the next must
	not be reported as one alignment across the separator."""
	names = genome.names
	reads = []
	for a, b in zip(names, names[1:]):
		left = genome[a][-50:].upper()
		right = genome[b][:50].upper()
		if "N" in left or "N" in right:
			continue
		reads.append(("junction_%s_%s" % (a, b), left + right))
	assert reads
	fq = tmp_path / "j.fq"
	with open(str(fq), "w") as fh:
		for name, seq in reads:
			fh.write("@%s\n%s\n+\n%s\n" % (name, seq, "I" * len(seq)))
	sam = gem3.run_map(genome_index, inputs=[str(fq)], threads=threads).sam()
	for rec in sam.records:
		if rec.is_unmapped:
			continue
		assert rec.ref_end <= len(genome[rec.rname]), rec
		assert check_record_against_reference(rec, genome[rec.rname]) == [], rec


def test_reference_case_and_iupac_codes(gem3, tmp_path, threads):
	"""Lowercase reference bases are indexed as bases (MD reports them in
	upper case); IUPAC codes are indexed as N, so a read over one has a
	mismatch there."""
	rng = Rng(77)
	body = rng.dna(3000)
	ref = body[:1000] + body[1000:1200].lower() + body[1200:1500] + "R" + body[1501:]
	fa = write_fasta_records(str(tmp_path / "iupac.fa"), [("chrI", ref)])
	index = gem3.build_index(fa)
	reads = [
		("soft", body[1050:1150]),               # inside the lowercase stretch
		("edge", body[1150:1250]),               # half lowercase
		("iupac", body[1450:1550]),              # over the R at 1500
		("iupac_rc", revcomp(body[1450:1550])),
	]
	fq = tmp_path / "r.fq"
	with open(str(fq), "w") as fh:
		for name, seq in reads:
			fh.write("@%s\n%s\n+\n%s\n" % (name, seq, "I" * len(seq)))
	by_name = primary(gem3.run_map(index, inputs=[str(fq)], threads=threads).sam())
	for name, start, nm in (("soft", 1050, 0), ("edge", 1150, 0), ("iupac", 1450, 1), ("iupac_rc", 1450, 1)):
		rec = by_name[name]
		assert (rec.rname, rec.pos) == ("chrI", start + 1), rec
		assert rec.tags["NM"] == nm, rec
		# The index holds N where the FASTA has R
		assert check_record_against_reference(rec, ref.replace("R", "N")) == [], rec
	assert by_name["iupac"].tags["MD"] == "50N49"


@pytest.mark.parametrize("case", ["missing", "empty", "not_gem", "fasta"])
def test_index_errors(gem3, genome, tmp_path, case):
	"""A missing, empty or foreign index file stops gem-mapper with a
	message, without crashing."""
	fq = write_reads(tmp_path / "r.fq", simulate_se(genome, 5, 100, seed=1))
	if case == "missing":
		index, message = str(tmp_path / "nope.gem"), "Could not stat file"
	elif case == "empty":
		index, message = str(tmp_path / "empty.gem"), "Mapping zero-bytes File"
		open(index, "wb").close()
	elif case == "not_gem":
		index, message = str(tmp_path / "junk.gem"), "Wrong GEM-Index Model"
		with open(index, "wb") as fh:
			fh.write(b"\x01" * 4096)
	else:
		index, message = genome.write_fasta(str(tmp_path / "g.fa")), "Wrong GEM-Index Model"
	res = gem3.run_map(index, inputs=[fq], check=False, threads=1)
	assert res.returncode == 1 and not res.crashed, res.describe()
	assert message in res.stderr_text
	assert res.stdout == b""
