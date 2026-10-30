# GEM-Mapper v3 (GEM3) test suite.
# Copyright (c) 2026 Jacob Schreiber <jacob.schreiber@umassmed.edu>
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of GEM-Mapper v3 (GEM3), distributed under the GNU General
# Public License version 3 or (at your option) any later version; see LICENSE.
"""Output destinations, compression and the threads that share them:
-o, --gzip-output, --bzip-output, --output-model, standard input, -t, and
the ordering of output blocks (src/io/output_file.c, buffered_output_file.c,
src/mapper/mapper.c, mapper_io.c).

The reference for every comparison is the plain single-threaded run: the
other modes must write the same SAM bytes (the @PG line, which holds the
command line, is left out). README "Determinism": the SAM is identical at
every -t and keeps the input order; paired output depends on the input
buffer size, so paired comparisons keep --input-model fixed.
"""

import bz2
import gzip
import zlib

import pytest

import area_d
from area_d import body_lines
from gem3harness import sam_digest, strip_variable_header

# Small input buffers so a few hundred reads span many input blocks, which
# the mapper threads take in turn and whose output must come back in order.
SMALL_INPUT = ["--input-model", "1M,c,4K"]


@pytest.fixture(scope="module")
def se(gem3, genome, genome_index):
	return area_d.se_dataset(gem3, genome, genome_index)


@pytest.fixture(scope="module")
def pe(gem3, genome, genome_index):
	return area_d.pe_dataset(gem3, genome, genome_index)


@pytest.fixture(scope="module")
def bs(gem3, genome, genome_bs_index):
	return area_d.bs_dataset(gem3, genome, genome_bs_index)


def plain(data):
	return strip_variable_header(data)


def test_output_file_matches_stdout(gem3, se, tmp_path):
	"""-o <file> writes the stdout bytes to the file and nothing to stdout."""
	res = gem3.run_map(se.index, inputs=se.paths, args=["-o", "out.sam"], cwd=str(tmp_path))
	assert res.stdout == b""
	assert plain((tmp_path / "out.sam").read_bytes()) == plain(se.run().stdout)


@pytest.mark.parametrize("mode", ["gzip", "bzip"])
@pytest.mark.parametrize("to_file", [False, True], ids=["stdout", "file"])
def test_compressed_output(gem3, se, tmp_path, mode, to_file):
	"""--gzip-output / --bzip-output write one valid gzip or bzip2 stream
	(to stdout or the -o file) that decompresses to the plain SAM."""
	args = ["--%s-output" % mode] + (["-o", "out.sam.x"] if to_file else []) + SMALL_INPUT
	res = gem3.run_map(se.index, inputs=se.paths, args=args, cwd=str(tmp_path))
	data = (tmp_path / "out.sam.x").read_bytes() if to_file else res.stdout
	if to_file:
		assert res.stdout == b""
	if mode == "gzip":
		assert data[:3] == b"\x1f\x8b\x08"
		d = zlib.decompressobj(16 + zlib.MAX_WBITS)
		text = d.decompress(data) + d.flush()
		assert d.eof and d.unused_data == b""  # a single, complete member
		assert gzip.decompress(data) == text
	else:
		assert data[:3] == b"BZh"
		text = bz2.decompress(data)
	assert plain(text) == plain(se.run().stdout)


def test_last_compression_option_wins(gem3, se, tmp_path):
	"""--gzip-output and --bzip-output together: the later one applies."""
	res = gem3.run_map(se.index, inputs=se.paths, args=["--gzip-output", "--bzip-output"])
	assert plain(bz2.decompress(res.stdout)) == plain(se.run().stdout)
	res = gem3.run_map(se.index, inputs=se.paths, args=["--bzip-output", "--gzip-output"])
	assert plain(gzip.decompress(res.stdout)) == plain(se.run().stdout)


def test_standard_input(gem3, se):
	"""Without -i the reads come from standard input."""
	with open(se.paths[0], "rb") as fh:
		data = fh.read()
	res = gem3.run_map(se.index, stdin=data)
	assert plain(res.stdout) == plain(se.run().stdout)


@pytest.mark.parametrize("model", ["4M,5c", "64K,3", "4K,2", "1K,1", "512,4"])
def test_output_model_does_not_change_output(se, pe, model):
	"""--output-model <buffer_size>,<num_buffers>: any buffer size that holds
	a record and any number of buffers, including fewer buffers than
	threads, give the same bytes at 1 and 4 threads (records larger than a
	buffer are split across buffer extensions)."""
	for ds in (se, pe):
		ref = plain(ds.run(SMALL_INPUT, threads=1).stdout)
		for t in (1, 4):
			got = ds.run(SMALL_INPUT + ["--output-model", model], threads=t).stdout
			assert plain(got) == ref, (model, t, ds.paired)


@pytest.mark.parametrize("kind", ["se", "pe", "bs"])
def test_thread_count_does_not_change_output(se, pe, bs, kind):
	"""-t 1, 2, 3, 4 and 8 give the same SAM bytes over many input blocks,
	and records follow the input order."""
	ds = {"se": se, "pe": pe, "bs": bs}[kind]
	ref = ds.run(SMALL_INPUT, threads=1)
	for t in (2, 3, 4, 8):
		assert plain(ds.run(SMALL_INPUT, threads=t).stdout) == plain(ref.stdout), t
	names = [l.split("\t", 1)[0] for l in body_lines(ref.stdout)]
	order = []
	for n in names:
		if not order or order[-1] != n:
			order.append(n)
	want = [r.name for r in ds.reads] if ds.reads is not None else [a.name for a, b in ds.pairs]
	assert order == want


def test_single_end_output_ignores_input_buffer_size(se):
	"""Single-end output does not depend on the input buffer size."""
	ref = plain(se.run().stdout)
	assert plain(se.run(SMALL_INPUT).stdout) == ref
	assert plain(se.run(["--input-model", "64K,c,2K"]).stdout) == ref


def test_threads_read_whole_inputs(gem3, se):
	"""-t larger than the number of input blocks maps every read once."""
	res = se.run(threads=8)
	names = [r.qname for r in res.sam().records if r.is_primary]
	assert names == [r.name for r in se.reads]


# sha256 of the multi-block SAM output on base 6e155e3 (pins the current
# output; checked at -t 1 and 4).
PINNED = {
	"se": "2711545c5118e477b70783ef4ae33238bba80a11344eaae108dae0b8849c6472",
	"pe": "c7c82ae5087589bfaf11f5ef0811791278c9663a62164399b41b1ca12fe72cfb",
	"bs": "497814ce34d7c053cf037f01192d304b6ae7ff2400fb4507e6b1bd3fe99d200c",
}


@pytest.mark.parametrize("kind", sorted(PINNED))
def test_multiblock_output_bytes_are_pinned(se, pe, bs, threads, kind):
	"""Pins the current output over many input blocks at the suite's -t."""
	ds = {"se": se, "pe": pe, "bs": bs}[kind]
	assert sam_digest(ds.run(SMALL_INPUT, threads=threads).stdout) == PINNED[kind]
