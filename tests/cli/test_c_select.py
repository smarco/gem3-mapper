# GEM-Mapper v3 (GEM3) test suite.
# Copyright (c) 2026 Jacob Schreiber <jacob.schreiber@umassmed.edu>
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of GEM-Mapper v3 (GEM3), distributed under the GNU General
# Public License version 3 or (at your option) any later version; see LICENSE.
"""Single-end match selection (-M, -m) and MAPQ models (--mapq-model).

Expected placements come from ``area_c.repeat_genome()``, whose repeat
copies and their substitutions are designed by hand, and from the
simulator's truth. MAPQ expectations follow the SAM specification (a read
with two equally good placements is wrong at least half the time, so MAPQ
<= 3) or the model classes in src/archive/score/.

Not tested (known-bugs/C-mapq-depends-on-max-reported-matches.md): MAPQ under
-M below the number of equally good placements; the byte pins record it.
"""

import pytest

from gem3harness import Rng, SimRead, md_equal, md_nm_from_alignment, sam_digest, simulate_se, write_reads
from gem3harness.sam import SECONDARY
from area_c import PREDICTOR_FIELDS, event_distance, gap_affine_score, repeat_genome, split_predictors
from gem3harness.genome import revcomp


@pytest.fixture(scope="module")
def rg():
	return repeat_genome()


@pytest.fixture(scope="module")
def rg_index(gem3, rg):
	return gem3.index(rg.genome)


def fastq(records):
	return "".join("@%s\n%s\n+\n%s\n" % (name, seq, "I" * len(seq)) for name, seq in records)


def map_records(gem3, index, records, args, threads, tmp_path, name="r.fq"):
	path = tmp_path / name
	path.write_text(fastq(records))
	return gem3.run_map(index, inputs=[str(path)], args=list(args), threads=threads)


def placements(sam, qname):
	return [(r.pos - 1, r.strand, r.tags["NM"]) for r in sam.records if r.qname == qname]


def designed(copies, read_strand):
	"""(position, strand as the read maps, substitutions) of each copy."""
	flip = {"+": "-", "-": "+"}
	return sorted((p, s if read_strand == "+" else flip[s], d) for p, s, d in copies)


def xa_entries(rec):
	out = []
	for item in rec.tags.get("XA", "").split(";"):
		if item:
			contig, spos, cigar, nm = item.split(",")
			out.append((contig, int(spos[1:]) - 1, spos[0], cigar, int(nm)))
	return out


# -M ------------------------------------------------------------------------

def test_max_reported_all_finds_every_designed_copy(gem3, rg, rg_index, tmp_path, threads):
	"""-M all with uncompacted SAM: one record per designed copy of the unit
	(1 exact, 2/3/4 with 1/2/3 substitutions), primary first and exact,
	NM equal to the substitutions, CIGAR, NM and MD consistent with the
	reference; the same for the reverse-complemented read."""
	reads = [("fwd", rg.strata_read), ("rev", revcomp(rg.strata_read))]
	sam = map_records(gem3, rg_index, reads, ["-M", "all", "--sam-compact=false"], threads, tmp_path).sam()
	for name, strand in (("fwd", "+"), ("rev", "-")):
		recs = [r for r in sam.records if r.qname == name]
		assert sorted(placements(sam, name)) == designed(rg.strata_copies, strand)
		assert not recs[0].has(SECONDARY) and recs[0].tags["NM"] == 0
		assert all(r.has(SECONDARY) for r in recs[1:])
		assert [r.tags["NM"] for r in recs] == sorted(r.tags["NM"] for r in recs)
		read = dict(reads)[name]
		for r in recs:
			# Secondary records carry SEQ '*': use the read on the record's strand
			seq = r.seq if r.seq != "*" else (read if r.strand == "+" else revcomp(read))
			assert r.cigar == "100M"
			md, nm = md_nm_from_alignment(rg.seq, r.ref_start, r.cigar_ops, seq)
			assert r.tags["NM"] == nm and md_equal(r.tags["MD"], md), r
			assert r.tags["AS"] == gap_affine_score(rg.seq, r.ref_start, r.cigar_ops, seq)


@pytest.mark.parametrize("m", [1, 2, 3, 5, 7])
def test_max_reported_matches_keeps_the_best(gem3, rg, rg_index, tmp_path, threads, m):
	"""-M m reports the m best placements (fewest substitutions), as
	separate records (--sam-compact=false) or as the primary record's XA
	list; both forms name the same placements."""
	reads = [("fwd", rg.strata_read), ("rev", revcomp(rg.strata_read))]
	full = map_records(gem3, rg_index, reads, ["-M", m, "--sam-compact=false"], threads, tmp_path).sam()
	compact = map_records(gem3, rg_index, reads, ["-M", m], threads, tmp_path).sam()
	best = sorted(d for _, _, d in rg.strata_copies)[:m]
	for name, strand in (("fwd", "+"), ("rev", "-")):
		got = placements(full, name)
		assert len(got) == m
		assert sorted(d for _, _, d in got) == best
		assert set(got) <= set(designed(rg.strata_copies, strand))
		recs = [r for r in compact.records if r.qname == name]
		assert len(recs) == 1
		alt = [(p, st, nm) for _, p, st, _, nm in xa_entries(recs[0])]
		assert len(alt) == m - 1
		assert sorted([(recs[0].pos - 1, recs[0].strand, recs[0].tags["NM"])] + alt) == sorted(got)
		if m == 1:
			assert "XA" not in recs[0].tags


def test_max_reported_matches_default_is_five(gem3, rg, rg_index, tmp_path, threads):
	"""Without -M, five placements are reported (--help: default=5)."""
	sam = map_records(gem3, rg_index, [("fwd", rg.strata_read)], ["--sam-compact=false"], threads,
		tmp_path).sam()
	assert len(sam.records) == 5


# -m ------------------------------------------------------------------------

def test_min_reported_strata_reports_whole_strata(gem3, rg, rg_index, tmp_path, threads):
	"""-m k reports every placement of the first k strata (edit distances),
	even beyond the default of five matches: the tie unit has six exact
	copies; without -m only five of them are reported."""
	sam = map_records(gem3, rg_index, [("t", rg.ties_read)], ["--sam-compact=false"], threads, tmp_path).sam()
	exact = set(p for p in designed(rg.ties_copies, "+") if p[2] == 0)
	got = placements(sam, "t")
	assert len(got) == 5 and set(got) <= exact
	for k in (1, 2):
		sam = map_records(gem3, rg_index, [("t", rg.ties_read)], ["-m", k, "--sam-compact=false"], threads,
			tmp_path).sam()
		assert exact <= set(placements(sam, "t"))
	# Strata of the strata unit: 0 (1 copy), 1 (2), 2 (3), 3 (4)
	for k, max_subs in ((1, 0), (2, 1), (3, 2), (4, 3)):
		sam = map_records(gem3, rg_index, [("s", rg.strata_read)], ["-m", k, "--sam-compact=false"], threads,
			tmp_path).sam()
		want = set(p for p in designed(rg.strata_copies, "+") if p[2] <= max_subs)
		got = set(placements(sam, "s"))
		assert want <= got, (k, sorted(want - got))
		assert got <= set(designed(rg.strata_copies, "+"))


def test_min_reported_strata_fraction_of_read_length(gem3, rg, rg_index, tmp_path, threads):
	"""-m below 1 is a fraction of the read length, rounded down: 0.02 of
	100 bp is 2 strata; 0.005 is 0, the same as not giving -m."""
	reads = [("s", rg.strata_read), ("t", rg.ties_read), ("t_rc", revcomp(rg.ties_read))]

	def digest(args):
		return sam_digest(map_records(gem3, rg_index, reads, args + ["--sam-compact=false"], threads,
			tmp_path).stdout)
	assert digest(["-m", "0.02"]) == digest(["-m", "2"])
	assert digest(["-m", "0.005"]) == digest([])
	assert digest([]) != digest(["-m", "2"])


@pytest.mark.parametrize("args,message", [
	(["-M", "0"], "'--max-reported-matches' must be greater than zero"),
	(["-m", "1", "-M", "3"], "cannot be used at the same time"),
	(["-M", "3", "-m", "1"], "cannot be used at the same time"),
	(["--mapq-model", "best"], "Option '--mapq-model' must be in"),
	(["--discordant-pair-search", "sometimes"], "must be 'always'|'if-no-concordant'|'never'"),
	(["--bisulfite-read", "3"], "must be '1'|'2'|'interleaved'|'inferred'|'non-stranded'"),
])
def test_option_value_errors(gem3, genome_index, tmp_path, args, message):
	"""Invalid values of the reporting, MAPQ, discordant-search and bisulfite
	options are rejected with a message naming the valid ones."""
	fq = tmp_path / "r.fq"
	fq.write_text(fastq([("r", "ACGT" * 25)]))
	res = gem3.run_map(genome_index, inputs=[str(fq)], args=args, threads=1, check=False)
	assert res.returncode != 0 and not res.crashed
	assert message in res.stderr_text
	assert b"\t" not in res.stdout


# MAPQ ------------------------------------------------------------------------

def test_mapq_ties_and_unique_reads(gem3, genome, genome_index, rg, rg_index, tmp_path, threads):
	"""A read with two or more equally good placements gets MAPQ <= 3 (SAM:
	-10 log10 Pr{wrong} with Pr >= 1/2); exact reads from unique sequence
	get a high MAPQ."""
	ties = [("ties", rg.ties_read), ("dup", rg.seq[rg.dup_starts[0] + 500:rg.dup_starts[0] + 600]),
		("dup_rc", revcomp(rg.seq[rg.dup_starts[1] + 900:rg.dup_starts[1] + 1000]))]
	sam = map_records(gem3, rg_index, ties, [], threads, tmp_path).sam()
	for rec in sam.records:
		assert not rec.is_unmapped and rec.tags["NM"] == 0
		assert rec.mapq <= 3, rec
	reads = simulate_se(genome, 200, 100, seed=71)
	fq = write_reads(tmp_path / "u.fq", reads)
	sam = gem3.run_map(genome_index, inputs=[fq], threads=threads).sam()
	for rec in sam.records:
		assert rec.tags["NM"] == 0 and "XA" not in rec.tags
		assert rec.mapq >= 20, rec


@pytest.fixture(scope="module")
def mixed_se(genome):
	"""Reads from everywhere with sequencing errors, plus reads taken from
	every repeat copy, tandem repeat and homopolymer of the genome."""
	reads = simulate_se(genome, 300, 100, seed=73, sub_rate=0.015, ins_rate=0.002, del_rate=0.002,
		unique_only=False)
	rng = Rng(74)
	for i, f in enumerate(f for f in genome.features if f.kind in ("repeat", "tandem", "homopolymer")):
		for j in range(3):
			start = max(0, min(len(genome[f.contig]) - 100, f.start - 50 + rng.randint(0, max(1, f.end - f.start))))
			seq = genome[f.contig][start:start + 100].upper()
			if "N" in seq:
				continue
			if rng.bernoulli(0.5):
				seq = revcomp(seq)
			reads.append(SimRead("f%d_%d" % (i, j), seq, "I" * len(seq), None))
	# Reads with 14 substitutions (beyond the 12% global alignment limit):
	# local alignments or unmapped
	for r in simulate_se(genome, 12, 100, seed=75, subs=14, spacing=2, edge=2):
		reads.append(SimRead("loc" + r.name, r.seq, r.qual, r.truth))
	return reads


def _without_mapq(rec):
	f = rec.line.split("\t")
	return f[:4] + f[5:]


def test_mapq_model_none(gem3, genome_index, mixed_se, tmp_path, threads):
	"""--mapq-model none changes only MAPQ, which is 0 everywhere (pins the
	current value)."""
	fq = write_reads(tmp_path / "m.fq", mixed_se)
	gem = gem3.run_map(genome_index, inputs=[fq], threads=threads).sam()
	none = gem3.run_map(genome_index, inputs=[fq], args=["--mapq-model", "none"], threads=threads).sam()
	assert [_without_mapq(r) for r in gem.records] == [_without_mapq(r) for r in none.records]
	assert set(r.mapq for r in none.records) == {0}
	assert len(set(r.mapq for r in gem.records)) > 2


def test_mapq_model_classify_classes(gem3, genome_index, mixed_se, tmp_path, threads):
	"""--mapq-model classify changes only MAPQ; the values are the model's
	classes: unique reads 240-250 (well separated) or 110-170, reads with
	other placements 220-230 or 10-70, ties 0 or 2, local alignments 1."""
	fq = write_reads(tmp_path / "m.fq", mixed_se)
	gem = gem3.run_map(genome_index, inputs=[fq], threads=threads).sam()
	cls = gem3.run_map(genome_index, inputs=[fq], args=["--mapq-model", "classify"], threads=threads).sam()
	assert [_without_mapq(r) for r in gem.records] == [_without_mapq(r) for r in cls.records]
	unique = set(range(110, 171)) | set(range(240, 251)) | {1}
	multi = set(range(10, 71)) | set(range(220, 231)) | {0, 1, 2}
	for g, c in zip(gem.records, cls.records):
		if c.is_unmapped:
			assert c.mapq == 0
		elif "XA" in c.tags:
			assert c.mapq in multi, c
		else:
			assert c.mapq in unique | multi, c
	assert any(c.mapq >= 240 for c in cls.records)


def test_mapq_model_dump_predictors(gem3, genome, genome_index, mixed_se, tmp_path, threads):
	"""--mapq-model dump-predictors scores like the gem model and prints one
	line of predictors per mapped read. The predictors that the SAM record
	determines are recomputed: edit = (L - NM) / L and event = (L - events) / L
	for unclipped alignments, swg = AS / L (match score 1), template sigma 1
	and both mate MAPQs 0 for single-end reads, delta -1 exactly for reads
	with one placement (-M all)."""
	fq = write_reads(tmp_path / "m.fq", mixed_se)
	args = ["-M", "all"]
	gem = gem3.run_map(genome_index, inputs=[fq], args=args, threads=threads).sam()
	res = gem3.run_map(genome_index, inputs=[fq], args=args + ["--mapq-model", "dump-predictors"],
		threads=threads)
	sam, predictors = split_predictors(res.stdout)
	assert [r.line for r in sam.records] == [r.line for r in gem.records]
	by_tag = {}
	for fields in predictors:
		assert len(fields) == len(PREDICTOR_FIELDS), fields
		d = dict(zip(PREDICTOR_FIELDS, fields))
		assert d["tag"] not in by_tag
		by_tag[d["tag"]] = d
	mapped = [r for r in sam.records if not r.is_unmapped]
	assert sorted(by_tag) == sorted(r.qname for r in mapped)
	for rec in mapped:
		d = by_tag[rec.qname]
		L = float(len(rec.seq))
		if "S" not in rec.cigar:  # local alignments also count their clipped bases
			assert abs(float(d["edit"]) - (L - rec.tags["NM"]) / L) < 1e-6, (rec, d)
			ev = event_distance(genome[rec.rname], rec.ref_start, rec.cigar_ops, rec.seq)
			assert abs(float(d["event"]) - (L - ev) / L) < 1e-6, (rec, d)
		assert abs(float(d["swg"]) - rec.tags["AS"] / L) < 1e-6, (rec, d)
		assert float(d["sigma"]) == 1.0 and float(d["mapq_1"]) == 0.0 and float(d["mapq_2"]) == 0.0
		assert int(d["mcs2"]) == 0
		assert (int(d["delta"]) == -1) == ("XA" not in rec.tags), (rec, d)


# Pinned bytes ------------------------------------------------------------------

# sha256 of the SAM output without @PG/@CO from base 6e155e3 (gcc 13 release;
# checked on the other builds listed in README.md at -t 1 and 4).
PINNED = {
	"M1": "209e2a26604941ed435736e2e093cfdf64b6ae1fe54baf9c1781bfa7e04475c6",
	"M3": "cc92053a7688356872512b43fa8dab90c5cda5b6ccf602616f3f878ba7cb7fe1",
	"Mall-full": "1d420f73973d96077ad1e74380e7277929260b3ef8d1fb823ae20d2a85f9304a",
	"m1": "f993f13749fa0a7c6da3562d45fcfe6f0ffc08756eefdbfad9a26a46dd657bb0",
	"m0.03": "6e0f90b34ea7a46996cfb0d2eed433e8c88344a01645b7a58323fa8a051a8f7b",
	"classify": "9db220135536445efaebbc6adf9c505fd34144375b334c1177ca036c94059960",
	"rg-Mall-full": "8a9c12819a1c14a8ba03ac1374e7b032183c92839d05c282415ce02ee5e875d9",
	"rg-m1-full": "8a9c12819a1c14a8ba03ac1374e7b032183c92839d05c282415ce02ee5e875d9",
	"rg-M2": "a7cdb6647e85ff4a2e2c9d1e7880e4cb7a73c67bc027559246d4b579f6e1d5b0",
}
PIN_ARGS = {
	"M1": ["-M", "1"],
	"M3": ["-M", "3"],
	"Mall-full": ["-M", "all", "--sam-compact=false"],
	"m1": ["-m", "1"],
	"m0.03": ["-m", "0.03", "--sam-compact=false"],
	"classify": ["--mapq-model", "classify"],
	"rg-Mall-full": ["-M", "all", "--sam-compact=false"],
	"rg-m1-full": ["-m", "1", "--sam-compact=false"],
	"rg-M2": ["-M", "2"],
}


@pytest.mark.parametrize("case", sorted(PINNED))
def test_select_output_bytes_are_pinned(gem3, genome_index, mixed_se, rg, rg_index, tmp_path, threads, case):
	"""Pins the current output of match selection and MAPQ options on fixed
	inputs (reads from repeats included)."""
	if case.startswith("rg-"):
		reads = [("s", rg.strata_read), ("s_rc", revcomp(rg.strata_read)), ("t", rg.ties_read),
			("t_rc", revcomp(rg.ties_read)), ("d", rg.seq[rg.dup_starts[0] + 10:rg.dup_starts[0] + 110]),
			("x", rg.triple[50:150])]
		res = map_records(gem3, rg_index, reads, PIN_ARGS[case], threads, tmp_path)
	else:
		fq = write_reads(tmp_path / "m.fq", mixed_se)
		res = gem3.run_map(genome_index, inputs=[fq], args=PIN_ARGS[case], threads=threads)
	assert sam_digest(res.stdout) == PINNED[case]
