# GEM-Mapper v3 (GEM3) test suite.
# Copyright (c) 2026 Jacob Schreiber <jacob.schreiber@umassmed.edu>
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of GEM-Mapper v3 (GEM3), distributed under the GNU General
# Public License version 3 or (at your option) any later version; see LICENSE.
"""--report-file (src/stats/report_stats.c, merged across threads in
src/mapper/mapper.c).

The expected JSON is recomputed from the input reads and the SAM records of
the same run (area_d.expected_report_se/pe). The genome carries contigs
named like gem-mapper's sequencing, under- and over-conversion controls so
that every read class appears. N base counts are left out: the report
counts one N too many per read (known-bugs/D-report-basecounts-extra-n.md).
"""

import json

import pytest

import area_d
from area_d import expected_report_pe, expected_report_se, strip_report_n
from gem3harness import Rng, simulate_bisulfite, simulate_pe, simulate_se

REPORT = ["--report-file", "report.json"]


@pytest.fixture(scope="module")
def cgenome():
	return area_d.control_genome()


@pytest.fixture(scope="module")
def cindex(gem3, cgenome):
	return gem3.index(cgenome)


@pytest.fixture(scope="module")
def cbs_index(gem3, cgenome):
	return gem3.index(cgenome, bisulfite=True)


def _se_reads(genome):
	rng = Rng(9)
	reads = simulate_se(genome, 150, 100, seed=91, sub_rate=0.02, ins_rate=0.003, del_rate=0.003,
		name_prefix="c")
	reads += simulate_se(genome, 20, 70, seed=92, subs=1, name_prefix="k")
	rnd = [area_d.Read("z%d" % i, rng.dna(90), "I" * 90) for i in range(10)]
	return area_d._shuffle(rng, list(reads) + rnd)


def _pairs(genome):
	rng = Rng(10)
	pairs = simulate_pe(genome, 100, 100, insert_mean=300, insert_sd=30, insert_min=220, seed=93,
		sub_rate=0.01, name_prefix="c")
	out = []
	for i, p in enumerate(pairs):
		a, b = p.r1, p.r2
		if i % 9 == 4:
			b = area_d.Read(b.name, rng.dna(100), b.qual)
		elif i % 9 == 6:
			a = area_d.Read(a.name, rng.dna(100), a.qual)
		elif i % 18 == 8:
			a = area_d.Read(a.name, rng.dna(100), a.qual)
			b = area_d.Read(b.name, rng.dna(100), b.qual)
		out.append((a, b))
	return out


@pytest.fixture(scope="module")
def se(gem3, cgenome, cindex):
	return area_d.dataset(gem3, cindex, "report_se", lambda: ("se", _se_reads(cgenome)))


@pytest.fixture(scope="module")
def pe(gem3, cgenome, cindex):
	return area_d.dataset(gem3, cindex, "report_pe", lambda: ("pe", _pairs(cgenome)))


@pytest.fixture(scope="module")
def bs(gem3, cgenome, cbs_index):
	return area_d.bs_dataset(gem3, cgenome, cbs_index, name="report_bs", seed=94, n=150)


def _bs_pairs_with_long_templates(genome):
	"""Directional bisulfite pairs, then pairs with 3,000 bp templates: far
	beyond the learnt template length, so paired as discordant with MAPQ 0."""
	pairs = simulate_bisulfite(genome, 150, 100, library="directional", seed=94, insert_mean=260,
		insert_sd=30)
	pairs += simulate_bisulfite(genome, 12, 100, library="directional", seed=95, insert_mean=3000,
		insert_sd=0, contigs=["chrA"], name_prefix="long")
	return [(p.r1, p.r2) for p in pairs]


@pytest.fixture(scope="module")
def bs_long(gem3, cgenome, cbs_index):
	return area_d.dataset(gem3, cbs_index, "report_bs_long",
		lambda: ("pe", _bs_pairs_with_long_templates(cgenome)))


def report(res):
	return json.loads(res.files["report.json"].decode("ascii"))


def test_report_single_end(cgenome, se):
	"""Single end: MapperType, mapped reads per class (general and the three
	controls) and unmapped reads, base counts, MAPQ, read-length and
	mismatch histograms of the primary alignments."""
	res = se.run(REPORT)
	got = report(res)
	want = expected_report_se(se.reads, res.sam(), cgenome)
	assert strip_report_n(got) == want
	assert list(got) == list(want)
	assert set(got["Reads"]) == {"General", "SequencingControl", "UnderConversionControl",
		"OverConversionControl", "Unmapped"}
	assert got["Reads"]["Unmapped"][0] >= 10


def test_report_paired_end(cgenome, pe):
	"""Paired end: the same per end, plus CorrectPairs and the
	template-length histogram of the proper pairs."""
	res = pe.run(REPORT)
	got = report(res)
	want = expected_report_pe(pe.pairs, res.sam(), cgenome)
	assert strip_report_n(got) == want
	assert list(got) == list(want)
	assert got["CorrectPairs"] > 50 and min(got["Reads"]["Unmapped"]) > 5


def test_report_bisulfite_pairs(cgenome, bs, bs_long):
	"""Bisulfite pairs: NumReadsBS (ends aligned C-to-T or G-to-A) and the
	base counts of reads in pairs with MAPQ > 0 by conversion and class;
	pairs with MAPQ 0 (here discordant long templates) are left out of the
	conversion counts."""
	for ds in (bs, bs_long):
		res = ds.run(REPORT)
		got = report(res)
		long_ = set(a.name for a, b in ds.pairs if a.name.startswith("long"))
		want = expected_report_pe(ds.pairs, res.sam(), cgenome, bisulfite=True, discordant=long_)
		assert strip_report_n(got) == want
		assert "NumReadsBS" in got and "UnderConversionControlC2T" in got["BaseCounts"]
	prim = [r for r in res.sam().records if r.is_primary and r.qname.startswith("long")]
	assert len(prim) == 24 and all(r.mapq == 0 and not r.is_unmapped for r in prim)


def test_report_read_group(gem3, genome, genome_index):
	"""-r: the report repeats the @RG line with tabs written as \\t."""
	ds = area_d.se_dataset(gem3, genome, genome_index)
	res = ds.run(["-r", "@RG\\tID:g7\\tSM:s"] + REPORT)
	got = report(res)
	assert list(got)[:2] == ["MapperType", "ReadGroup"]
	assert got["ReadGroup"] == "@RG\tID:g7\tSM:s"
	assert strip_report_n(got) == expected_report_se(ds.reads, res.sam(), genome,
		read_group="@RG\tID:g7\tSM:s")


@pytest.mark.parametrize("kind", ["se", "pe", "bs"])
def test_report_same_at_any_thread_count(se, pe, bs, kind):
	"""Per-thread statistics merge to the same report at -t 1, 3 and 8 over
	many input blocks (the whole file, byte for byte)."""
	ds = {"se": se, "pe": pe, "bs": bs}[kind]
	args = REPORT + ["--input-model", "1M,c,4K"]
	ref = ds.run(args, threads=1).files["report.json"]
	for t in (3, 8):
		assert ds.run(args, threads=t).files["report.json"] == ref, t


def test_report_json_layout(se):
	"""The report is JSON with tab indentation and HistMapq printed eight
	values per line (pins the current layout)."""
	text = se.run(REPORT).files["report.json"].decode("ascii")
	assert text.startswith("{\n\t\"MapperType\": \"Single\",\n")
	assert text.endswith("}\n")
	lines = text.split("\n")
	i = lines.index("\t\"HistMapq\": [")
	row = lines[i + 1]
	assert row.startswith("\t\t") and row.count(",") == 8
