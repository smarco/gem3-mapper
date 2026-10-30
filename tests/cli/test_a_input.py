# GEM-Mapper v3 (GEM3) test suite.
# Copyright (c) 2026 Jacob Schreiber <jacob.schreiber@umassmed.edu>
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of GEM-Mapper v3 (GEM3), distributed under the GNU General
# Public License version 3 or (at your option) any later version; see LICENSE.
"""Area A: read input (-i, -1/-2, -p input parsing, -z, -j, pipes,
--input-model) as seen in gem-mapper's output.

The same reads given in every supported way (FASTQ or FASTA; plain, gzip,
bzip2; file, standard input or pipe; LF or CRLF; blocks and buffers of
other sizes) must give the same SAM records. Read names, sequences and
qualities in the SAM are checked against the FASTQ records and the SAM
specification; malformed input must stop gem-mapper with a message naming
the line.
"""

import re

import pytest

from gem3harness import (compress_bytes, fastq_record, revcomp, sam_digest, simulate_pe,
	simulate_se, write_interleaved, write_reads, write_split)
from area_a import assert_maps_to_truth, primary, records_without_qual


@pytest.fixture(scope="module")
def se_reads(genome):
	"""Unique reads with substitutions and indels (truth known), plus reads
	from repeats and low-complexity regions."""
	unique = simulate_se(genome, 250, 100, seed=101, sub_rate=0.01, ins_rate=0.002, del_rate=0.002)
	other = simulate_se(genome, 150, 75, seed=102, sub_rate=0.01, unique_only=False, name_prefix="m")
	return unique, unique + other


@pytest.fixture(scope="module")
def se_fastq_bytes(se_reads):
	_, reads = se_reads
	return "".join(fastq_record(r.name, r.seq, r.qual) for r in reads).encode("ascii")


@pytest.fixture(scope="module")
def se_reference_sam(gem3, genome, genome_index, se_reads, se_fastq_bytes, tmp_path_factory):
	"""Output for the plain FASTQ file, checked against the simulator's truth."""
	unique, reads = se_reads
	d = tmp_path_factory.mktemp("se")
	path = d / "reads.fq"
	path.write_bytes(se_fastq_bytes)
	res = gem3.run_map(genome_index, inputs=[str(path)], threads=4)
	sam = res.sam()
	assert_maps_to_truth(type(sam)(sam.header_lines, [r for r in sam.records if r.qname.startswith("r")]),
		unique, genome)
	by_name = primary(sam)
	for r in reads:
		rec = by_name[r.name]
		# SEQ and QUAL are the read's, reverse-complemented and reversed on the reverse strand
		assert rec.seq == (revcomp(r.seq) if rec.is_reverse else r.seq), rec
		assert rec.qual == (r.qual[::-1] if rec.is_reverse else r.qual), rec
	return res.stdout


def _variants(data):
	"""(label, file bytes, how to pass them, extra arguments)."""
	crlf = data.replace(b"\n", b"\r\n")
	return [
		("gzip_file", compress_bytes(data, "gz"), "file", []),
		("bzip2_file", compress_bytes(data, "bz2"), "file", []),
		("crlf", crlf, "file", []),
		("crlf_gzip", compress_bytes(crlf, "gz"), "file", []),
		("no_final_newline", data[:-1], "file", []),
		("two_gzip_members", compress_bytes(data[:len(data) // 2], "gz") + compress_bytes(data[len(data) // 2:], "gz"), "file", []),
		("stdin", data, "stdin", []),
		("stdin_gzip", compress_bytes(data, "gz"), "stdin", ["-z"]),
		("stdin_bzip2", compress_bytes(data, "bz2"), "stdin", ["-j"]),
		("pipe", data, "pipe", []),
		("pipe_gzip", compress_bytes(data, "gz"), "pipe_gzip", []),
		# -z/-j are consulted only for standard input: files are detected by content
		("gzip_file_with_j", compress_bytes(data, "gz"), "file", ["-j"]),
	]


@pytest.mark.parametrize("variant", [v[0] for v in _variants(b"")])
def test_input_forms_give_identical_output(gem3, genome_index, se_fastq_bytes, se_reference_sam, tmp_path, threads, variant):
	"""Compressed, CRLF, unterminated, standard-input and piped copies of
	the FASTQ map exactly like the plain file."""
	label, data, how, extra = [v for v in _variants(se_fastq_bytes) if v[0] == variant][0]
	if how == "file":
		path = tmp_path / ("reads.fq.gz" if "gzip" in label else "reads.fq")
		path.write_bytes(data)
		res = gem3.run_map(genome_index, inputs=[str(path)], args=extra, threads=threads)
	elif how == "stdin":
		res = gem3.run_map(genome_index, args=extra, stdin=data, threads=threads)
	else:
		path = tmp_path / "reads.bin"
		path.write_bytes(data)
		command = ("cat %s |" if how == "pipe" else "gzip -dc %s|") % path
		res = gem3.run_map(genome_index, args=["-i", command] + extra, threads=threads)
	assert sam_digest(res.stdout) == sam_digest(se_reference_sam)


def test_fasta_input(gem3, genome_index, se_reads, se_reference_sam, tmp_path, threads):
	"""FASTA reads give the FASTQ output with QUAL '*' (no qualities: the
	same as -q ignore on the FASTQ)."""
	_, reads = se_reads
	fa = write_reads(tmp_path / "r.fa", reads, fmt="fasta")
	res = gem3.run_map(genome_index, inputs=[fa], threads=threads)
	assert all(rec.qual == "*" for rec in res.sam().records)
	fq = write_reads(tmp_path / "r.fq", reads)
	ignored = gem3.run_map(genome_index, inputs=[fq], args=["-q", "ignore"], threads=threads)
	assert records_without_qual(res.stdout) == records_without_qual(ignored.stdout)


@pytest.mark.parametrize("model", ["64K,3,16K", "128K,2,8K", "256K,4,32K", "1M,2,64K", "32M,4,4M"])
def test_input_model_does_not_change_single_end_output(gem3, genome_index, se_fastq_bytes, se_reference_sam, tmp_path, threads, model):
	"""--input-model block_size,num_blocks,buffer_size changes how the input
	is cut into buffers, not the single-end output. (Blocks of 4K or less
	lose lines, always with a few hundred bytes and now and then at -t 4:
	known-bugs/A-input-model-small-blocks-drop-lines.md.)"""
	path = tmp_path / "reads.fq"
	path.write_bytes(se_fastq_bytes)
	res = gem3.run_map(genome_index, inputs=[str(path)], args=["--input-model", model], threads=threads)
	assert sam_digest(res.stdout) == sam_digest(se_reference_sam)


def test_read_names(gem3, genome, genome_index, tmp_path, threads):
	"""QNAME is the header up to the first space or tab, without a /1, /2 or
	/3 suffix; a CASAVA 1.8 comment is dropped; a comment field ending in
	/1-/3 is joined to the name with '_' (SRA style)."""
	seq = genome["chr1"][20000:20100].upper()
	headers = [
		("plain", "plain"), ("end1/1", "end1"), ("end2/2", "end2"), ("end3/3", "end3"), ("x/4", "x/4"),
		("cas 1:N:0:ACGTAC", "cas"), ("cas2 2:Y:18:ATCACG", "cas2"),
		("tab\tsome comment", "tab"), ("spc  two  spaces", "spc"),
		("SRR384920.1 HWI-ST382_0049:1:1:1217:1879/1", "SRR384920.1_HWI-ST382_0049:1:1:1217:1879"),
	]
	text = "".join("@%s\n%s\n+\n%s\n" % (h, seq, "I" * 100) for h, _ in headers)
	res = gem3.run_map(genome_index, args=["-i", "r.fq"], files={"r.fq": text.encode()}, threads=threads)
	names = [rec.qname for rec in res.sam().records]
	assert names == [q for _, q in headers]
	assert all(rec.pos == 20001 and rec.rname == "chr1" for rec in res.sam().records)


def test_read_characters(gem3, genome, genome_index, tmp_path, threads):
	"""Read bases keep their case; IUPAC codes and other characters inside a
	read become N in SEQ (the others with a warning naming the line) and
	count as mismatches; reads with lowercase bases map like uppercase ones."""
	seq = genome["chr2"][30000:30100].upper()
	reads = [
		("upper", seq), ("lower", seq.lower()), ("mixed", seq[:40].lower() + seq[40:]),
		("iupac", seq[:30] + "R" + seq[31:60] + "y" + seq[61:]),
		("dot", seq[:50] + "." + seq[51:]),
	]
	text = "".join("@%s\n%s\n+\n%s\n" % (n, s, "I" * 100) for n, s in reads)
	res = gem3.run_map(genome_index, args=["-i", "r.fq"], files={"r.fq": text.encode()}, threads=threads)
	by_name = primary(res.sam())
	for name, s in reads:
		rec = by_name[name]
		assert (rec.rname, rec.pos, rec.strand) == ("chr2", 30001, "+"), rec
		expected = "".join(c if c in "ACGTNacgtn" else "N" for c in s)
		assert rec.seq == expected, rec
	assert by_name["iupac"].tags["NM"] == 2 and by_name["iupac"].tags["MD"] == "30%s29%s39" % (seq[30], seq[60])
	assert by_name["dot"].tags["NM"] == 1
	# The warning names the file (line numbers come from an uninitialized
	# counter: known-bugs/A-input-line-counter-uninitialized.md)
	assert re.search(r"Parsing FASTA/FASTQ error\(r\.fq:\d+\)\. Input read sequence contains bad character", res.stderr_text)


@pytest.mark.parametrize("case", ["no_header", "late_no_header", "no_plus", "short_qualities", "low_quality"])
def test_malformed_fastq_is_rejected(gem3, genome, genome_index, threads, case):
	"""Malformed records stop gem-mapper (exit 1, no crash) with the parse
	error. (The line number in the message comes from an uninitialized
	counter, known-bugs/A-input-line-counter-uninitialized.md: release builds
	print the 1-based line given here, make debug one less; not checked.)"""
	seq = genome["chr1"][5000:5060].upper()
	good = "@ok\n%s\n+\n%s\n" % (seq, "I" * 60)
	text, line, message = {
		"no_header": ("ok\n%s\n+\n%s\n" % (seq, "I" * 60), 1, "Beginning Symbol ('>' or '@') not found"),
		"late_no_header": (good * 2 + "bad\n%s\n+\n%s\n" % (seq, "I" * 60), 9, "Beginning Symbol ('>' or '@') not found"),
		"no_plus": (good + "@x\n%s\n-\n%s\n" % (seq, "I" * 60), 7, "Wrong separator symbol ('+')"),
		"short_qualities": (good + "@x\n%s\n+\n%s\n" % (seq, "I" * 57), 8, "Different read and qualities lengths"),
		"low_quality": (good + "@x\n%s\n+\n%s\n" % (seq, "I" * 20 + "\x1f" + "I" * 39), 8,
			"Input read qualities contains invalid quality value"),
	}[case]
	res = gem3.run_map(genome_index, args=["-i", "r.fq"], files={"r.fq": text.encode()}, threads=threads, check=False)
	assert res.returncode == 1 and not res.crashed, res.describe()
	assert re.search(r"Parsing FASTA/FASTQ error\(r\.fq:\d+\)\. " + re.escape(message), res.stderr_text), res.stderr_text[:500]
	if case == "short_qualities":
		# The error report shows the record being parsed
		assert "\n%s\n+\n%s" % (seq, "I" * 57) in res.stderr_text, res.stderr_text


def test_missing_input_and_failing_pipe(gem3, genome_index, threads):
	"""An input file that does not exist, or a pipe whose command fails,
	stops gem-mapper with exit 1 and a message."""
	res = gem3.run_map(genome_index, args=["-i", "nonexistent.fq"], threads=threads, check=False)
	assert res.returncode == 1 and not res.crashed
	assert "Could not open file 'nonexistent.fq'" in res.stderr_text
	res = gem3.run_map(genome_index, args=["-i", "cat nonexistent.fq |"], threads=threads, check=False)
	assert res.returncode == 1 and not res.crashed
	assert "Could not close file 'cat nonexistent.fq '" in res.stderr_text


@pytest.mark.parametrize("how", ["file", "stdin"])
def test_empty_input(gem3, genome, genome_index, threads, how):
	"""No reads: a header with one @SQ per contig and no records."""
	if how == "file":
		res = gem3.run_map(genome_index, args=["-i", "e.fq"], files={"e.fq": b""}, threads=threads)
	else:
		res = gem3.run_map(genome_index, stdin=b"", threads=threads)
	sam = res.sam()
	assert sam.records == [] and sam.sq == genome.lengths()


@pytest.mark.parametrize("bad", ["1M,2", "1X,2,1M", "1M,2,4Q"])
def test_input_model_rejects_bad_values(gem3, genome_index, threads, bad):
	"""A wrong number of fields or a size that does not parse is an error.
	(num_blocks 0, or text read as 0, hangs: known-bugs/A-input-model-zero-blocks-hang.md.)"""
	res = gem3.run_map(genome_index, args=["-i", "r.fq", "--input-model", bad], files={"r.fq": b""}, threads=threads, check=False)
	assert res.returncode == 1 and not res.crashed
	assert "--input-model" in res.stderr_text


@pytest.fixture(scope="module")
def pairs(genome):
	return simulate_pe(genome, 150, 100, insert_mean=300, insert_sd=30, seed=111, sub_rate=0.01)


@pytest.fixture(scope="module")
def pe_reference_sam(gem3, genome_index, pairs, tmp_path_factory):
	d = tmp_path_factory.mktemp("pe")
	fq1, fq2 = write_split(d / "1.fq", d / "2.fq", pairs)
	return gem3.run_map(genome_index, inputs=[fq1, fq2], threads=4).stdout


@pytest.mark.parametrize("form", ["split_gzip", "split_bzip2", "split_crlf", "interleaved", "interleaved_gzip",
	"interleaved_stdin", "split_casava"])
def test_paired_input_forms_give_identical_output(gem3, genome_index, pairs, pe_reference_sam, tmp_path, threads, form):
	"""Pairs read from two files (-1/-2, any compression, CRLF, CASAVA ends
	instead of /1 /2) or one interleaved file (-p, file or standard input)
	give the same SAM; each mate is named without its /1 or /2."""
	if form.startswith("split"):
		comp = {"split_gzip": "gz", "split_bzip2": "bz2"}.get(form)
		suffixes = (" 1:N:0:ACGT", " 2:N:0:ACGT") if form == "split_casava" else ("/1", "/2")
		fq1, fq2 = write_split(tmp_path / "1.fq", tmp_path / "2.fq", pairs, compression=comp, suffixes=suffixes)
		if form == "split_crlf":
			for p in (fq1, fq2):
				with open(p, "rb") as fh:
					data = fh.read()
				with open(p, "wb") as fh:
					fh.write(data.replace(b"\n", b"\r\n"))
		res = gem3.run_map(genome_index, inputs=[fq1, fq2], threads=threads)
	elif form == "interleaved_stdin":
		path = write_interleaved(tmp_path / "i.fq", pairs)
		with open(path, "rb") as fh:
			res = gem3.run_map(genome_index, args=["-p"], stdin=fh.read(), threads=threads)
	else:
		comp = "gz" if form == "interleaved_gzip" else None
		path = write_interleaved(tmp_path / "i.fq", pairs, compression=comp)
		res = gem3.run_map(genome_index, inputs=[path], args=["-p"], threads=threads)
	assert sam_digest(res.stdout) == sam_digest(pe_reference_sam)
	names = [r.qname for r in res.sam().records]
	assert names == [n for p in pairs for n in (p.r1.name, p.r1.name)]


# Unequal split files and odd interleaved files report the error and then
# hang (known-bugs/D-unsynced-pairs-hang.md): not tested.
