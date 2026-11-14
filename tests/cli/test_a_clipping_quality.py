# GEM-Mapper v3 (GEM3) test suite.
# Copyright (c) 2026 Jacob Schreiber <jacob.schreiber@umassmed.edu>
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of GEM-Mapper v3 (GEM3), distributed under the GNU General
# Public License version 3 or (at your option) any later version; see LICENSE.
"""Area A: read preparation, --clipping and -q/--quality-format.

Clipping removes read ends before the search (uncalled: N; masked:
lowercase; fixed,L,R: L and R bases); they appear as soft clips (S) in the
CIGAR, on the CIGAR's right for reverse-strand records, and POS moves to
the first aligned base. Quality strings are decoded with the chosen ASCII
offset, so the same phred values in two encodings must map identically.
"""

import pytest

from gem3harness import Rng, check_record_against_reference, revcomp, sam_digest, simulate_se, write_reads
from area_a import primary, records_without_qual


def _unique_starts(genome, rng, n, length=100):
	"""Positions on chr1/chr2 of reads away from repeats and N runs."""
	out = []
	while len(out) < n:
		contig = rng.choice(["chr1", "chr2"])
		start = rng.randint(2000, len(genome[contig]) - 2000)
		seq = genome[contig][start:start + length]
		if "N" in seq.upper() or genome.hard_regions(contig, start, start + length, margin=length):
			continue
		out.append((contig, start))
	return out


def _run(gem3, index, reads, args, threads):
	text = "".join("@%s\n%s\n+\n%s\n" % (name, seq, "I" * len(seq)) for name, seq in reads)
	return gem3.run_map(index, args=["-i", "r.fq"] + list(args), files={"r.fq": text.encode()}, threads=threads)


def _expected_cigar(left, right, length, reverse):
	"""Soft clips of a read's first ``left`` and last ``right`` bases, in
	reference orientation."""
	if reverse:
		left, right = right, left
	ops = []
	if left:
		ops.append("%dS" % left)
	ops.append("%dM" % (length - left - right))
	if right:
		ops.append("%dS" % right)
	return "".join(ops)


@pytest.fixture(scope="module")
def starts(genome):
	return _unique_starts(genome, Rng(141), 24)


@pytest.mark.parametrize("mode", ["default", "uncalled", "bare"])
def test_clipping_uncalled(gem3, genome, genome_index, starts, threads, mode):
	"""Default (and --clipping, --clipping=uncalled): N at either end is
	clipped; an N inside the read is not."""
	reads, expected = [], {}
	for i, (contig, start) in enumerate(starts):
		ref = genome[contig][start:start + 100].upper()
		left, right = [(3, 0), (0, 5), (2, 2), (0, 0)][i % 4]
		seq = "N" * left + ref[left:100 - right] + "N" * right
		if i % 8 == 7:
			seq = seq[:50] + "N" + seq[51:]  # internal N: a mismatch, not a clip
		reverse = i % 2 == 1
		name = "u%d" % i
		reads.append((name, revcomp(seq) if reverse else seq))
		# A reverse-strand record's left/right clips swap
		read_left, read_right = (right, left) if reverse else (left, right)
		expected[name] = (contig, start + left + 1, _expected_cigar(read_left, read_right, 100, reverse), "-" if reverse else "+")
	args = {"default": [], "uncalled": ["--clipping=uncalled"], "bare": ["--clipping"]}[mode]
	by_name = primary(_run(gem3, genome_index, reads, args, threads).sam())
	for name, (contig, pos, cigar, strand) in expected.items():
		rec = by_name[name]
		assert (rec.rname, rec.pos, rec.cigar, rec.strand) == (contig, pos, cigar, strand), rec


def test_clipping_masked(gem3, genome, genome_index, starts, threads):
	"""--clipping=masked: lowercase bases at either end are clipped;
	lowercase bases inside the read are kept."""
	reads, expected = [], {}
	for i, (contig, start) in enumerate(starts):
		ref = genome[contig][start:start + 100].upper()
		left, right = [(5, 6), (0, 4), (7, 0), (0, 0)][i % 4]
		seq = ref[:left].lower() + ref[left:100 - right] + ref[100 - right:].lower()
		seq = seq[:40] + seq[40:45].lower() + seq[45:]  # internal lowercase, not clipped
		reverse = i % 2 == 1
		name = "m%d" % i
		# revcomp keeps the case; in reference orientation the clips stay where they are in seq
		reads.append((name, revcomp(seq) if reverse else seq))
		expected[name] = (contig, start + left + 1, _expected_cigar(left, right, 100, False), "-" if reverse else "+")
	by_name = primary(_run(gem3, genome_index, reads, ["--clipping=masked"], threads).sam())
	for name, (contig, pos, cigar, strand) in expected.items():
		rec = by_name[name]
		assert (rec.rname, rec.pos, rec.cigar, rec.strand) == (contig, pos, cigar, strand), rec
		assert rec.tags["NM"] == 0


@pytest.mark.parametrize("left,right", [(5, 3), (0, 10), (12, 0), (0, 0)])
def test_clipping_fixed(gem3, genome, genome_index, starts, threads, left, right):
	"""--clipping=fixed,L,R removes the read's first L and last R bases."""
	reads, expected = [], {}
	for i, (contig, start) in enumerate(starts):
		ref = genome[contig][start:start + 100].upper()
		reverse = i % 2 == 1
		name = "f%d" % i
		reads.append((name, revcomp(ref) if reverse else ref))
		# In reference orientation a reverse read's first bases are on the right
		ref_left = right if reverse else left
		expected[name] = (contig, start + ref_left + 1, _expected_cigar(left, right, 100, reverse), "-" if reverse else "+")
	by_name = primary(_run(gem3, genome_index, reads, ["--clipping=fixed,%d,%d" % (left, right)], threads).sam())
	for name, (contig, pos, cigar, strand) in expected.items():
		rec = by_name[name]
		assert (rec.rname, rec.pos, rec.cigar, rec.strand) == (contig, pos, cigar, strand), rec


def test_clipping_none(gem3, genome, genome_index, starts, threads):
	"""--clipping=none: reads without N at their ends are unchanged by
	uncalled clipping, so both give the same records; reads with terminal N
	still map to their origin with alignments consistent with the reference."""
	clean = [("c%d" % i, genome[c][s:s + 100].upper()) for i, (c, s) in enumerate(starts)]
	none = _run(gem3, genome_index, clean, ["--clipping=none"], threads)
	default = _run(gem3, genome_index, clean, [], threads)
	assert sam_digest(none.stdout) == sam_digest(default.stdout)
	reads = [("n%d" % i, "NNN" + genome[c][s + 3:s + 97].upper() + "NNN") for i, (c, s) in enumerate(starts)]
	by_name = primary(_run(gem3, genome_index, reads, ["--clipping=none"], threads).sam())
	for (name, _), (contig, start) in zip(reads, starts):
		rec = by_name[name]
		assert rec.rname == contig and start <= rec.pos - 1 <= start + 3, rec
		assert check_record_against_reference(rec, genome[contig]) == [], rec


def test_clipping_fixed_longer_than_read(gem3, genome, genome_index, starts, threads):
	"""Clipping the whole read leaves nothing to map: the read is unmapped."""
	contig, start = starts[0]
	reads = [("all", genome[contig][start:start + 100].upper())]
	sam = _run(gem3, genome_index, reads, ["--clipping=fixed,200,0"], threads).sam()
	assert [r.is_unmapped for r in sam.records] == [True]


@pytest.mark.parametrize("value", ["bogus", "fixed,5", "masked,1"])
def test_clipping_rejects_bad_values(gem3, genome_index, threads, value):
	res = gem3.run_map(genome_index, args=["-i", "r.fq", "--clipping=" + value], files={"r.fq": b""}, threads=threads, check=False)
	assert res.returncode == 1 and not res.crashed
	assert "Option '--clipping' must be" in res.stderr_text


@pytest.fixture(scope="module")
def quality_reads(genome):
	"""Reads from anywhere with errors and qualities spread over 2-40, so
	that qualities at mismatches change MAPQ."""
	return simulate_se(genome, 300, 100, seed=81, sub_rate=0.03, ins_rate=0.003, del_rate=0.003,
		unique_only=False, qual_range=(2, 40))


class _Read(object):
	def __init__(self, name, seq, qual):
		self.name, self.seq, self.qual = name, seq, qual


def test_quality_offsets(gem3, genome_index, quality_reads, tmp_path, threads):
	"""Phred+33 qualities read with offset-33 (the default) and the same
	phred values written as phred+64 read with offset-64 give identical
	records apart from QUAL; -q ignore gives the records of the same reads
	given as FASTA (no qualities); qualities do change some MAPQ."""
	fq33 = write_reads(tmp_path / "q33.fq", quality_reads)
	fq64 = write_reads(tmp_path / "q64.fq", [_Read(r.name, r.seq, "".join(chr(ord(c) + 31) for c in r.qual))
		for r in quality_reads])
	fa = write_reads(tmp_path / "q.fa", quality_reads, fmt="fasta")
	default = gem3.run_map(genome_index, inputs=[fq33], threads=threads).stdout
	q33 = gem3.run_map(genome_index, inputs=[fq33], args=["-q", "offset-33"], threads=threads).stdout
	q64 = gem3.run_map(genome_index, inputs=[fq64], args=["--quality-format", "offset-64"], threads=threads).stdout
	ignore = gem3.run_map(genome_index, inputs=[fq33], args=["-q", "ignore"], threads=threads).stdout
	fasta = gem3.run_map(genome_index, inputs=[fa], threads=threads).stdout
	assert sam_digest(default) == sam_digest(q33)
	assert records_without_qual(q33) == records_without_qual(q64)
	assert records_without_qual(ignore) == records_without_qual(fasta)
	assert records_without_qual(q33) != records_without_qual(ignore)


def test_quality_format_rejects_bad_values(gem3, genome_index, threads):
	res = gem3.run_map(genome_index, args=["-i", "r.fq", "-q", "phred"], files={"r.fq": b""}, threads=threads, check=False)
	assert res.returncode == 1 and not res.crashed
	assert "Option '-q|--quality-format' must be 'ignore', 'offset-33' or 'offset-64'" in res.stderr_text
