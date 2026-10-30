# GEM-Mapper v3 (GEM3) test suite.
# Copyright (c) 2026 Jacob Schreiber <jacob.schreiber@umassmed.edu>
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of GEM-Mapper v3 (GEM3), distributed under the GNU General
# Public License version 3 or (at your option) any later version; see LICENSE.
"""Helpers for area A's tests (index load and input).

``indexed_text`` rebuilds, from FASTA records, the text gem-indexer indexes
and the locator intervals gem-mapper prints with ``--verbose=dev``;
``kmer_intervals`` gives suffix-array intervals of k-mers over that text by
brute force. ``primary`` and ``assert_maps_to_truth`` check simulated reads.
"""

import bisect

from gem3harness import check_record_against_reference, revcomp

# Encoded order of the index alphabet: A < C < G < T < N < separator
_RANK = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4, "|": 5}
NS_THRESHOLD = 50  # gem-indexer's default --strip-unknown-bases-threshold


def indexed_text(records, ns_threshold=NS_THRESHOLD):
	"""(full_text, forward_length, intervals) for FASTA ``records`` [(name, seq)].

	Rules (gem-indexer, regular index): bases are upper-cased; N and IUPAC
	codes are N; a run of N shorter than ``ns_threshold`` stays in the text,
	a longer one is dropped, ends the current interval and becomes an
	``uncalled`` interval of index length 0; every non-empty interval is
	followed by one separator '|'. The reverse-complement intervals follow
	in reverse order, then a final separator. ``intervals`` holds
	(type, idx_begin, idx_end, name, strand, text_begin, text_end).
	"""
	fwd = []
	intervals = []
	for name, seq in records:
		seq = "".join(c if c.upper() in "ACGT" else "N" for c in seq).upper()
		i = 0
		start_text = 0
		start_idx = len(fwd)
		length = 0
		while i < len(seq):
			if seq[i] == "N":
				j = i
				while j < len(seq) and seq[j] == "N":
					j += 1
				if j - i < ns_threshold:
					fwd.extend("N" * (j - i))
					length += j - i
				else:
					if length > 0:
						intervals.append(("chromosome", start_idx, start_idx + length, name, "+",
							start_text, start_text + length))
						fwd.append("|")
					intervals.append(("uncalled", len(fwd), len(fwd), name, "+", i, j))
					start_idx = len(fwd)
					start_text = j
					length = 0
				i = j
			else:
				fwd.append(seq[i])
				length += 1
				i += 1
		if length > 0:
			intervals.append(("chromosome", start_idx, start_idx + length, name, "+",
				start_text, start_text + length))
			fwd.append("|")
	forward = "".join(fwd)
	pos = len(forward)
	rc_intervals = []
	for (kind, b, e, name, strand, tb, te) in reversed(intervals):
		rc_intervals.append((kind, pos, pos + (e - b), name, "-", tb, te))
		pos += e - b
		if e > b:
			pos += 1  # separator
	full = forward + revcomp(forward[:-1]) + "||"
	return full, len(forward), intervals + rc_intervals


def kmer_intervals(text, k):
	"""A function kmer -> (lo, hi): the interval of an ACGT k-mer among the
	rotations of the cyclic ``text`` sorted by _RANK (for an absent k-mer,
	lo == hi is where it would be inserted)."""
	n = len(text)
	keys = []
	for i in range(n):
		window = (text[i:i + k] if i + k <= n else text[i:] + text[:i + k - n])
		keys.append(tuple(_RANK[c] for c in window))
	keys.sort()

	def interval(kmer):
		key = tuple(_RANK[c] for c in kmer)
		return bisect.bisect_left(keys, key), bisect.bisect_right(keys, key)
	return interval


def primary(sam):
	"""qname -> the primary record (one per read in single-end output)."""
	out = {}
	for rec in sam.records:
		if rec.is_primary:
			assert rec.qname not in out, rec
			out[rec.qname] = rec
	return out


def forward_seq(read):
	return read.seq if read.truth.strand == "+" else revcomp(read.seq)


def assert_maps_to_truth(sam, reads, genome, allow_unmapped=0):
	"""Every read has a primary record at its true contig, position and
	strand (up to ``allow_unmapped`` unmapped), with SEQ the read on the
	forward strand and CIGAR/NM/MD consistent with the reference."""
	by_name = primary(sam)
	assert sorted(by_name) == sorted(r.name for r in reads)
	unmapped = 0
	for read in reads:
		rec = by_name[read.name]
		if rec.is_unmapped:
			unmapped += 1
			continue
		t = read.truth
		assert (rec.rname, rec.pos, rec.strand) == (t.contig, t.pos + 1, t.strand), (read, rec)
		assert rec.seq == forward_seq(read), rec
		assert check_record_against_reference(rec, genome[rec.rname]) == [], rec
	assert unmapped <= allow_unmapped, "%d reads unmapped" % unmapped


def records_without_qual(data):
	"""SAM records as tuples without QUAL (column 11), for comparing runs
	that differ only in how qualities were given."""
	out = []
	for line in (data.decode("ascii") if isinstance(data, bytes) else data).split("\n"):
		if line and not line.startswith("@"):
			f = line.split("\t")
			out.append(tuple(f[:10] + f[11:]))
	return out
