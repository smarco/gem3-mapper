# GEM-Mapper v3 (GEM3) test suite.
# Copyright (c) 2026 Jacob Schreiber <jacob.schreiber@umassmed.edu>
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of GEM-Mapper v3 (GEM3), distributed under the GNU General
# Public License version 3 or (at your option) any later version; see LICENSE.
"""Area A: seeding (--candidate-generation, --candidate-generation-adaptive).

Each region-profile strategy must find the simulated reads at their origin
with alignments consistent with the reference; the output of each is also
pinned, since the seeding decides which candidates are verified. Strategies
and values that crash or misplace reads on base 6e155e3 are left out (see
known-bugs/A-candidate-generation-strategies.md and
candidate-generation-adaptive.md).
"""

import re

import pytest

from gem3harness import check_record_against_reference, sam_digest, simulate_se, write_reads
from area_a import assert_maps_to_truth, digest_without_mapq, primary

# name -> (arguments, reads with indels for the pinned run?, guaranteed
# substitutions): a read with at most that many substitutions leaves one of
# the strategy's disjoint exact-search regions error-free (pigeonhole), so it
# must be found; None for the multi-stage default and adaptive profiles,
# which must find every simulated read, indels included.
STRATEGIES = {
	"default": ([], True, None),
	"fixed": (["--candidate-generation", "fixed,20,20,0,100"], True, 4),                # 5 regions
	"fixed_overlapping": (["--candidate-generation", "fixed,12,6,0,50"], True, 7),      # 8 disjoint of 15
	"fixed_gapped": (["--candidate-generation", "fixed,20,25,0,100"], True, 3),         # 4 regions
	"CKS": (["--candidate-generation", "CKS,3,20,100"], True, 2),                        # best 3 of 5
	"OPS": (["--candidate-generation", "OPS,4,20,100"], True, 3),                        # 4 windows
	"factors": (["--candidate-generation", "factors,4,0,100"], True, 3),                 # 4 x 25 bp
	"MEM": (["--candidate-generation", "MEM,100"], False, 6),
	"SMEM": (["--candidate-generation", "SMEM,100"], False, 6),
	"adaptive_20_4_1": (["--candidate-generation-adaptive", "20,4,1"], True, None),
	"adaptive_5_2_2": (["--candidate-generation-adaptive", "5,2,2"], True, None),
	"adaptive_100_8_1": (["--candidate-generation-adaptive", "100,8,1"], True, None),
}

# sha256 of the SAM output (without @PG/@CO) of each strategy on
# _pinned_reads(), base 6e155e3 (pins the current output)
PINNED = {
	"default": "cf5eb18d5f165645977127d80158e031001de33a8c80d333d416f005920eb1e7",
	"fixed": "5438322922895fd88b2312d0c940d8cc09e45ceb7d5a153b102458c5b8fdd0bb",
	"fixed_overlapping": "a3710b3a616bd4eab71e891647105fbfbd5471b80c4012ca1a61b50bf94799b9",
	"fixed_gapped": "4fbfc4860ee47e685fb92b0c75d5ef7ffd8ba339a1c57b4230697e99bc899d9a",
	"CKS": "4a91d2f0a76966724575d5c37297147567c4b6f5beb242a199576737e68f303a",
	"OPS": "20f73474cbf5821628fbc06c2df6d0576483683aa81c59e9678e01b5f2aff2a2",
	"factors": "7694d511e18c78f43e62049ae93ea5c8975989135561e2f91a4dd7691bae1b2c",
	"MEM": "430fa4cda0f31e9d59472ce099be9a324c17b74c8c9769aa4eba075e1a79c4f2",  # digest_without_mapq
	"SMEM": "70ff9d9b3b2f04c07435ac5886c0ab51ed98dbe3a96ee5873e08d2d987de6f61",
	"adaptive_20_4_1": "cf5eb18d5f165645977127d80158e031001de33a8c80d333d416f005920eb1e7",
	"adaptive_5_2_2": "847b983264d90fe3c19035b5b0708ccf1bcdd184fb054e874b240d50a8d20d37",
	"adaptive_100_8_1": "afeba0fc2648ab439c0cbdf737adae066c30e6ca48e46bc9d9715f3d11251a21",
}


@pytest.fixture(scope="module")
def truth_reads(genome):
	"""Unique reads: with substitutions and indels, and with exactly 0-6
	substitutions (40 of each)."""
	with_indels = simulate_se(genome, 250, 100, seed=121, sub_rate=0.02, ins_rate=0.003, del_rate=0.003)
	by_subs = []
	for k in range(7):
		by_subs += simulate_se(genome, 40, 100, seed=130 + k, subs=k, name_prefix="s%d_" % k)
	return with_indels, by_subs


def _pinned_reads(genome, indels):
	# Reads from anywhere, repeats and low-complexity regions included
	if indels:
		return simulate_se(genome, 300, 100, seed=123, sub_rate=0.015, ins_rate=0.002, del_rate=0.002, unique_only=False)
	return simulate_se(genome, 300, 100, seed=124, sub_rate=0.015, unique_only=False)


@pytest.mark.parametrize("name", sorted(STRATEGIES))
def test_strategy_maps_reads_to_origin(gem3, genome, genome_index, truth_reads, tmp_path, threads, name):
	"""Reads land at their origin (all of them, or all within the strategy's
	guaranteed number of substitutions); every alignment reported is
	consistent with the reference."""
	args, _, guaranteed = STRATEGIES[name]
	with_indels, by_subs = truth_reads
	if guaranteed is None:
		fq = write_reads(tmp_path / "r.fq", with_indels)
		assert_maps_to_truth(gem3.run_map(genome_index, inputs=[fq], args=args, threads=threads).sam(),
			with_indels, genome)
	fq = write_reads(tmp_path / "s.fq", by_subs)
	sam = gem3.run_map(genome_index, inputs=[fq], args=args, threads=threads).sam()
	by_name = primary(sam)
	for read in by_subs:
		rec = by_name[read.name]
		if rec.is_unmapped:
			assert guaranteed is not None and len(read.truth.edits) > guaranteed, (read, rec)
			continue
		assert check_record_against_reference(rec, genome[rec.rname]) == [], rec
		t = read.truth
		assert (rec.rname, rec.pos, rec.strand) == (t.contig, t.pos + 1, t.strand), (read, rec)


@pytest.mark.parametrize("name", sorted(STRATEGIES))
def test_strategy_output_is_pinned(gem3, genome, genome_index, tmp_path, threads, name):
	"""Pins the current output of each strategy on reads from repeats and
	low-complexity regions too; every mapped record stays consistent with
	the reference."""
	args, indels, _ = STRATEGIES[name]
	fq = write_reads(tmp_path / "r.fq", _pinned_reads(genome, indels))
	res = gem3.run_map(genome_index, inputs=[fq], args=args, threads=threads)
	for rec in res.sam().records:
		if not rec.is_unmapped:
			assert check_record_against_reference(rec, genome[rec.rname]) == [], rec
	if name == "MEM":
		# One read's MAPQ differs between -Ofast and -O3 builds with the same
		# predictors (floating point, src/archive/score/); MAPQ is left out
		assert digest_without_mapq(res.stdout) == PINNED[name]
	else:
		assert sam_digest(res.stdout) == PINNED[name]


@pytest.mark.parametrize("group", [
	# On 100 bp reads each of these is the five regions [0,20), ..., [80,100), all searched
	["fixed,20,20,0,100", "CKS,5,20,100", "OPS,5,20,100", "factors,5,0,100"],
	# The adaptive profile's documented defaults are 20,4,1
	[None, "adaptive:20,4,1"],
], ids=["five-20bp-regions", "adaptive-defaults"])
def test_equivalent_settings_give_identical_output(gem3, genome, genome_index, tmp_path, threads, group):
	"""Settings that define the same regions give byte-identical output."""
	fq = write_reads(tmp_path / "r.fq", _pinned_reads(genome, True))
	digests = set()
	for value in group:
		if value is None:
			args = []
		elif value.startswith("adaptive:"):
			args = ["--candidate-generation-adaptive", value.split(":")[1]]
		else:
			args = ["--candidate-generation", value]
		digests.add(sam_digest(gem3.run_map(genome_index, inputs=[fq], args=args, threads=threads).stdout))
	assert len(digests) == 1


def test_test_strategy_prints_every_profile(gem3, genome, genome_index, tmp_path):
	"""'test' computes the profile of every strategy for each read, prints
	them to stderr and keeps no region, so no read is mapped."""
	reads = simulate_se(genome, 6, 100, seed=125, subs=1)
	fq = write_reads(tmp_path / "r.fq", reads)
	res = gem3.run_map(genome_index, inputs=[fq], args=["--candidate-generation", "test,5,20,20,0,100"], threads=1)
	assert all(rec.is_unmapped for rec in res.sam().records)
	labels = re.findall(r"\[GEM\]>Region\.Profile \((\w+)\)", res.stderr_text)
	order = ["Fixed", "CKS", "OPS", "Factors", "MEM", "SMEM", "Adaptive", "OPP"]
	assert labels[:8] == order
	for label in order:
		assert labels.count(label) >= len(reads), label


@pytest.mark.parametrize("value,message", [
	("bogus", "Invalid <strategy> from option '--candidate-generation'"),
	("fixed,20,20", "Invalid usage '--candidate-generation fixed,"),
	("CKS,5,20", "Invalid usage '--candidate-generation CKS,"),
	("OPS,5", "Invalid usage '--candidate-generation OPS,"),
	("factors,5", "Invalid usage '--candidate-generation factors,"),
	("MEM", "Invalid usage '--candidate-generation MEM,"),
	("OPP", "Invalid usage '--candidate-generation OPP,"),
	("test,1,2", "Invalid usage '--candidate-generation test,"),
])
def test_candidate_generation_rejects_bad_values(gem3, genome_index, value, message):
	res = gem3.run_map(genome_index, args=["-i", "r.fq", "--candidate-generation", value],
		files={"r.fq": b""}, threads=1, check=False)
	assert res.returncode == 1 and not res.crashed, res.describe()
	assert message in res.stderr_text


@pytest.mark.parametrize("value", ["20,4", "20"])
def test_candidate_generation_adaptive_rejects_bad_values(gem3, genome_index, value):
	res = gem3.run_map(genome_index, args=["-i", "r.fq", "--candidate-generation-adaptive", value],
		files={"r.fq": b""}, threads=1, check=False)
	assert res.returncode == 1 and not res.crashed, res.describe()
	assert "--candidate-generation-adaptive" in res.stderr_text
