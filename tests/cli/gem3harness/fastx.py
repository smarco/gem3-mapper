# GEM-Mapper v3 (GEM3) test suite.
# Copyright (c) 2026 Jacob Schreiber <jacob.schreiber@umassmed.edu>
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of GEM-Mapper v3 (GEM3), distributed under the GNU General
# Public License version 3 or (at your option) any later version; see LICENSE.
"""FASTQ/FASTA writers: plain, gzip or bzip2; single, interleaved or split.

``compression`` is None, 'gz' or 'bz2'. Output bytes are deterministic (the
gzip header carries no file name and mtime 0).
"""

import bz2
import gzip
import io
import os


def _open_write(path, compression):
	path = os.fspath(path)
	if compression is None:
		return open(path, "wb")
	if compression == "gz":
		raw = open(path, "wb")
		return _ClosingGzip(raw)
	if compression == "bz2":
		return bz2.BZ2File(path, "wb")
	raise ValueError("compression must be None, 'gz' or 'bz2'")


class _ClosingGzip(object):
	def __init__(self, raw):
		self.raw = raw
		self.gz = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)

	def write(self, data):
		return self.gz.write(data)

	def close(self):
		self.gz.close()
		self.raw.close()

	def __enter__(self):
		return self

	def __exit__(self, *exc):
		self.close()


def fastq_record(name, seq, qual, suffix=""):
	if len(seq) != len(qual):
		raise ValueError("read %s: %d bases but %d qualities" % (name, len(seq), len(qual)))
	return "@%s%s\n%s\n+\n%s\n" % (name, suffix, seq, qual)


def fasta_record(name, seq, suffix=""):
	return ">%s%s\n%s\n" % (name, suffix, seq)


def _records_text(reads, fmt, suffix=""):
	parts = []
	for r in reads:
		if fmt == "fastq":
			parts.append(fastq_record(r.name, r.seq, r.qual, suffix))
		elif fmt == "fasta":
			parts.append(fasta_record(r.name, r.seq, suffix))
		else:
			raise ValueError("fmt must be 'fastq' or 'fasta'")
	return "".join(parts)


def write_reads(path, reads, fmt="fastq", compression=None, suffix=""):
	"""Write :class:`SimRead`-like objects (``name``, ``seq``, ``qual``)."""
	data = _records_text(reads, fmt, suffix).encode("ascii")
	with _open_write(path, compression) as fh:
		fh.write(data)
	return os.fspath(path)


def write_interleaved(path, pairs, fmt="fastq", compression=None, suffixes=("/1", "/2")):
	"""Mates alternate: read 1 then read 2 of each pair, names suffixed by ``suffixes``."""
	parts = []
	for p in pairs:
		parts.append(_records_text([p.r1], fmt, suffixes[0]))
		parts.append(_records_text([p.r2], fmt, suffixes[1]))
	with _open_write(path, compression) as fh:
		fh.write("".join(parts).encode("ascii"))
	return os.fspath(path)


def write_split(path1, path2, pairs, fmt="fastq", compression=None, suffixes=("/1", "/2")):
	"""Read 1 of every pair to ``path1`` and read 2 to ``path2``."""
	write_reads(path1, [p.r1 for p in pairs], fmt, compression, suffixes[0])
	write_reads(path2, [p.r2 for p in pairs], fmt, compression, suffixes[1])
	return os.fspath(path1), os.fspath(path2)


def compress_bytes(data, compression):
	"""In-memory compression with the same deterministic settings."""
	if compression is None:
		return data
	if compression == "gz":
		buf = io.BytesIO()
		with gzip.GzipFile(filename="", mode="wb", fileobj=buf, mtime=0) as gz:
			gz.write(data)
		return buf.getvalue()
	if compression == "bz2":
		return bz2.compress(data)
	raise ValueError("compression must be None, 'gz' or 'bz2'")
