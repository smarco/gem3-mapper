# GEM-Mapper v3 (GEM3) test suite.
# Copyright (c) 2026 Jacob Schreiber <jacob.schreiber@umassmed.edu>
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of GEM-Mapper v3 (GEM3), distributed under the GNU General
# Public License version 3 or (at your option) any later version; see LICENSE.
"""MAP writer (-F MAP, src/io/output_map.c) and the MAP fields of
--sam-gem-compatible (X4/X5).

Each MAP match is checked against the SAM record of the same alignment: its
contig, strand and 1-based position, and its CIGAR rebuilt in GEM's MAP
syntax from the SAM CIGAR, the read and the reference (area_d.map_cigar).
--map-format 1 is a known bug (D-map-format-1-not-implemented.md); names
with comments are left out (the MAP tag column joins them without spaces).
"""

import pytest

import area_d
from area_d import (body_lines, counter_values, counters_ok, map_cigar, parse_map_line,
	parse_pe_maps, parse_se_maps, xa_entries)
from gem3harness import parse_cigar, sam_digest
from gem3harness.sam import READ1, READ2, SECONDARY, UNMAPPED


@pytest.fixture(scope="module")
def se(gem3, genome, genome_index):
	return area_d.se_dataset(gem3, genome, genome_index)


@pytest.fixture(scope="module")
def pe(gem3, genome, genome_index):
	return area_d.pe_dataset(gem3, genome, genome_index)


def map_lines(res):
	text = res.stdout.decode("ascii")
	assert text.endswith("\n")
	return text[:-1].split("\n")


def test_se_map_matches_sam(genome, se):
	"""-F MAP, single end: one line per read in input order (name, read and
	qualities as given, counters, matches); the matches are the
	--sam-compact=false records in order, with the record's MAPQ."""
	lines = map_lines(se.run(["-F", "MAP"]))
	full = se.run(["--sam-compact=false"]).sam().by_name()
	assert len(lines) == len(se.reads)
	n_multi = 0
	for line, read in zip(lines, se.reads):
		m = parse_map_line(line)
		assert (m["tag"], m["read"], m["qual"]) == (read.name, read.seq, read.qual)
		assert counters_ok(m["counters"]), m["counters"]
		recs = full[read.name]
		maps = parse_se_maps(m["maps"])
		if recs[0].flag & UNMAPPED:
			assert maps == []
			continue
		want = [(r.rname, r.strand, r.pos,
			map_cigar(r.rname, r.strand, r.pos, r.cigar_ops, read.seq, genome), r.mapq) for r in recs]
		assert maps == want, (line, [r.line[:80] for r in recs])
		assert sum(counter_values(m["counters"])) >= len(maps)
		n_multi += len(maps) > 1
	assert n_multi >= 10


def test_pe_map_matches_sam(genome, pe):
	"""-F MAP, paired end: a proper pair prints one line (name without /1,
	'read1 read2', 'qual1 qual2', counters, pairings); pairing 0 is the
	primary pair with its template MAPQ and pairings 1.. are the XA entries
	of the two ends. Other pairs print one single-end line per end, named
	<name>/1 and <name>/2."""
	lines = map_lines(pe.run(["-F", "MAP"]))
	prim = [r for r in pe.run().sam().records if r.is_primary]
	k = 0
	n_proper = n_single = 0
	for i, (a, b) in enumerate(pe.pairs):
		r1, r2 = prim[2 * i], prim[2 * i + 1]
		if r1.flag & 2:
			m = parse_map_line(lines[k])
			k += 1
			assert (m["tag"], m["read"], m["qual"]) == (a.name, a.seq + " " + b.seq, a.qual + " " + b.qual)
			assert counters_ok(m["counters"])
			maps = parse_pe_maps(m["maps"])
			want = [((r1.rname, r1.strand, r1.pos, map_cigar(r1.rname, r1.strand, r1.pos, r1.cigar_ops, a.seq, genome)),
				(r2.rname, r2.strand, r2.pos, map_cigar(r2.rname, r2.strand, r2.pos, r2.cigar_ops, b.seq, genome)))]
			for e1, e2 in zip(xa_entries(r1), xa_entries(r2)):
				want.append(tuple((e[0], e[1], e[2], map_cigar(e[0], e[1], e[2], parse_cigar(e[3]), read.seq, genome))
					for e, read in ((e1, a), (e2, b))))
			assert [x[:2] for x in maps] == want, lines[k - 1]
			assert maps[0][2] == r1.mapq == r2.mapq
			n_proper += 1
			continue
		for rec, read, suffix in ((r1, a, "/1"), (r2, b, "/2")):
			m = parse_map_line(lines[k])
			k += 1
			assert (m["tag"], m["read"], m["qual"]) == (read.name + suffix, read.seq, read.qual)
			maps = parse_se_maps(m["maps"])
			if rec.flag & UNMAPPED:
				assert maps == []
				continue
			want = [(rec.rname, rec.strand, rec.pos, map_cigar(rec.rname, rec.strand, rec.pos, rec.cigar_ops, read.seq, genome))]
			for e in xa_entries(rec):
				want.append((e[0], e[1], e[2], map_cigar(e[0], e[1], e[2], parse_cigar(e[3]), read.seq, genome)))
			assert [x[:4] for x in maps] == want, lines[k - 1]
			assert maps[0][4] == rec.mapq
			n_single += 1
	assert k == len(lines)
	assert n_proper > 50 and n_single > 10


def test_map_format_3_differs_only_in_homopolymer_marks(se, pe):
	"""--map-format 3 prints indels flagged as homopolymer with '&' where
	format 2 prints '>'; nothing else differs. --map-format 2 is the default.
	(Only the run-length text path flags homopolymer indels; on these inputs
	none is flagged, so the two formats must print the same bytes.)"""
	for ds in (se, pe):
		v2 = ds.run(["-F", "MAP"]).stdout
		assert ds.run(["-F", "MAP", "--map-format", "2"]).stdout == v2
		v3 = ds.run(["-F", "MAP", "--map-format", "3"]).stdout
		assert len(v3) == len(v2)
		assert v3.replace(b"&", b">") == v2


def test_map_output_is_case_insensitive_option(se):
	"""-F takes MAP/SAM in any case."""
	assert se.run(["-F", "map"]).stdout == se.run(["-F", "MAP"]).stdout
	assert sam_digest(se.run(["-F", "sam"]).stdout) == sam_digest(se.run().stdout)


def test_gem_compatible_x4_x5_equal_map_columns(se, pe):
	"""--sam-gem-compatible: on a single-end primary X4:Z is the MAP counters
	column and X5:Z the MAP matches column. On paired output X4:Z of each end
	is the MAP line's counters column and X5:Z lists that end's alignment in
	each pairing (contig, strand, position, CIGAR) in order."""
	lines = map_lines(se.run(["-F", "MAP"]))
	recs = [r for r in se.run(["--sam-gem-compatible"]).sam().records if not r.flag & SECONDARY]
	assert len(lines) == len(recs)
	for line, rec in zip(lines, recs):
		m = parse_map_line(line)
		assert (rec.tags["X4"], rec.tags["X5"]) == (m["counters"], m["maps"]), (line, rec)
	lines = map_lines(pe.run(["-F", "MAP"]))
	prim = [r for r in pe.run(["--sam-gem-compatible"]).sam().records if not r.flag & SECONDARY]
	k = 0
	for i in range(0, len(prim), 2):
		r1, r2 = prim[i], prim[i + 1]
		if r1.flag & 2:
			m = parse_map_line(lines[k])
			k += 1
			maps = parse_pe_maps(m["maps"])
			for rec, end in ((r1, 0), (r2, 1)):
				assert rec.tags["X4"] == m["counters"]
				x5 = [x[:4] for x in parse_se_maps(rec.tags["X5"])]
				assert x5 == [p[end] for p in maps], (rec, line)
		else:
			for rec in (r1, r2):
				m = parse_map_line(lines[k])
				k += 1
				assert (rec.tags["X4"], rec.tags["X5"]) == (m["counters"], m["maps"])
	assert k == len(lines)


# sha256 of the MAP output on base 6e155e3 (pins the current output).
PINNED = {
	"se-map2": "e4c0f0e830218c4f66af609dad913040fcf971aee1484a309bdbdde9d361d013",
	"pe-map2": "3ad6c9ba26d702fa9846847d79ada7885e98c2b4e81d4e625935c8d16dfc6d89",
}


@pytest.mark.parametrize("case", sorted(PINNED))
def test_map_output_bytes_are_pinned(se, pe, case):
	"""Pins the current MAP output (no header, so the whole file)."""
	ds = se if case.startswith("se") else pe
	assert area_d.digest(ds.run(["-F", "MAP"]).stdout) == PINNED[case]
