# GEM-Mapper v3 (GEM3) test suite.
# Copyright (c) 2026 Jacob Schreiber <jacob.schreiber@umassmed.edu>
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of GEM-Mapper v3 (GEM3), distributed under the GNU General
# Public License version 3 or (at your option) any later version; see LICENSE.
"""Pins the current output of the search and alignment options.

Each case maps one fixed single-end input with one option set and compares
sha256 of the SAM output without @PG/@CO (sam_digest) with the value from
base 6e155e3. The inputs mix what the search and aligner treat differently:
unique and repeat reads with substitutions and indels, reads with a large
deletion, end mismatches, chimeras (local alignment), N-rich, random, 36 bp
and 250 bp reads. These tests check no property; the other test_b_* files
do. A failure here means the bytes changed, which the speed search's gate
forbids; inspect the difference before updating a digest.

complete-partition and complete-partition-filter compute MAPQ from match
metrics that mode never sets, so their MAPQ depends on what the memory held
(known-bugs/r1-complete-partition-mapq-uninitialized.md): those two pins
leave MAPQ out.
"""

import pytest

from gem3harness import fastq_record, sam_digest, simulate_se
from area_a import digest_without_mapq
import area_b


def _pin_input(genome):
	rng = area_b.Rng(41)
	parts = []
	reads = simulate_se(genome, 160, 100, seed=301, sub_rate=0.02, ins_rate=0.004, del_rate=0.004,
		unique_only=False, spacing=4)
	reads += simulate_se(genome, 40, 100, seed=302, subs=6, spacing=5, name_prefix="s")
	reads += simulate_se(genome, 40, 36, seed=303, sub_rate=0.02, unique_only=False, edge=4, spacing=4,
		name_prefix="short")
	# long9, and reads missing 12 or 18 bases, get MAPQ 59 from -Ofast builds
	# and 0 from -O3/-Og builds; they are left out so the pins hold on every
	# build the README lists.
	reads += [r for r in simulate_se(genome, 30, 250, seed=304, sub_rate=0.02, ins_rate=0.004,
		del_rate=0.004, unique_only=False, name_prefix="long") if r.name != "long9"]
	parts.extend(fastq_record(r.name, r.seq, r.qual) for r in reads)
	for name, seq, info in area_b.chimeras(genome, 16, [35, 45, 60], seed=42):
		parts.append(fastq_record(name, seq, "I" * len(seq)))
	ref = genome["chr1"]
	for i, n in enumerate((6, 24)):
		p = 20000 + 7000 * i
		parts.append(fastq_record("del%d" % n, (ref[p:p + 50] + ref[p + 50 + n:p + 100 + n]).upper(), "I" * 100))
	for i in range(10):
		p = 41000 + 3000 * i
		s = list(ref[p:p + 100].upper())
		q = (i % 5) if i < 5 else 99 - (i % 5)
		s[q] = "A" if s[q] != "A" else "C"
		parts.append(fastq_record("end%d" % i, "".join(s), "I" * 100))
	for i in range(8):
		s = list(genome["chr2"][3000 + 5000 * i:3100 + 5000 * i].upper())
		for k in range(2 + i % 3, 100, 17):
			s[k] = "N"
		parts.append(fastq_record("nrich%d" % i, "".join(s), "I" * 100))
	for i in range(6):
		parts.append(fastq_record("rand%d" % i, rng.dna(100), "I" * 100))
	return "".join(parts).encode("ascii")


@pytest.fixture(scope="module")
def pin_fastq(genome):
	return _pin_input(genome)


# Base 6e155e3; the same at -t 1 and 4 from gcc 13 release, -O3 and debug,
# gcc 11 (-O3) and clang 18 builds.
PINNED = {
	"bandwidth": "8e29bb0ffc2aca5f511e4e53cc92028e653bccdf20e0473c405fe92229b2bcd6",
	"bpm": "5741d7eff2060689f5b710144d927dc71f5882b06a9af7f5e5a544a7ba60bdee",
	"complete": "52e78074c0ea8592d2b4af6549d90187d09f5115310436a4be6330bf41cbde02",
	"complete-E2-s2": "52e78074c0ea8592d2b4af6549d90187d09f5115310436a4be6330bf41cbde02",
	"complete-E8-s0": "959cebb117453cdea308081979b62b97a122066406c7c2dd6816aa8b7ac3bba1",
	"complete-all": "85bf14aa049a7e48af7b49bc036130db13dc448a3a8d742aeb32b7e23a143c67",
	"complete-partition": "59f1cd5382993f861d2ab82cc0791025b93a1cc69d77811eef46094e329e2d98",  # without MAPQ
	"complete-partition-filter": "8221e714e8bef8dcb592512db38cd6cf52a00946147291b41844ac77412224e9",  # without MAPQ
	"complete-pure": "52e78074c0ea8592d2b4af6549d90187d09f5115310436a4be6330bf41cbde02",
	"default": "5741d7eff2060689f5b710144d927dc71f5882b06a9af7f5e5a544a7ba60bdee",
	"e20": "2d64ae2effdc7c6b84fde75e67e75a6a6454eb0f601832cc8a369df829fbeab2",
	"e5": "cd33bc20eee18cab6cf90d78d2817d642987d4e1b306f91c6173dfe7d53d65dd",
	"edit": "34db65d1341b8f1b728696747ccd84a5a10020f8ee062548f68a5e775c6e0fc5",
	"end-context": "339929286741445dc97bbf944eb1c2d9e5384d52b2b11110e55629c81a52a00b",
	"exact": "dfcefeb3a1fb694a26a3be167d5df12dc88ff634ef633a46d7629d342a5a8fff",
	"full-swg": "b0869a0268d0b03282deea21cf3ac498e8e2883190a9c01b65a34c8b1d13f25c",
	"global-thresholds": "6b6ff5e280bee8349bf67f928a98cb90b9cdc2f8f5441661b88cb202fb7ea342",
	"hamming": "3dd3550399f0824ff5a91374a64da7344be3f700e9af8bd9fc3ff64ca110a451",
	"hamming-complete": "e3493ac88ac89da04404843d8d29a45ceac4219572595a7b27d91ab19ae75c62",
	"local-never": "69f55a8fd44b6850b30e144cdc3078581cd7b94ae12082dbcaf5cb1b073b9123",
	"local-thresholds": "14c63fd66a742a7b28cbe45e70575ddf9d3735ba610658a5c275768b6fb0825e",
	"penalties": "a54b252fa765e32698e31b44ea0745d47c3fc58093dd2678d90f88673fcec718",
	"qgram": "5741d7eff2060689f5b710144d927dc71f5882b06a9af7f5e5a544a7ba60bdee",
	"scaffolding-3": "b71ecc1730474ef272c39cfff53d8404f05e1b5455d40957c2a8085390a52df1",
	"scaffolding-30": "a5edc777e39d31563d48b2c3b0940f4acff231a19d0a3a143aa98b47b268a433",
	"sensitive": "5c6944ce288d70db4ebf8bc352c01f7a2dd9b06a4b59ad7029b44fd0c75fcf5d",
}

CONFIGS = {
	"default": [],
	"sensitive": ["--mapping-mode", "sensitive"],
	"complete": ["--mapping-mode", "complete"],
	"complete-E8-s0": ["--mapping-mode", "complete", "-E", "0.08", "-s", "0"],
	"complete-E2-s2": ["--mapping-mode", "complete", "-E", "2", "-s", "2"],
	"complete-all": ["--mapping-mode", "complete", "-M", "all", "--sam-compact=false"],
	"complete-pure": ["--mapping-mode", "complete-pure"],
	"complete-partition": ["--mapping-mode", "complete-partition"],
	"complete-partition-filter": ["--mapping-mode", "complete-partition-filter", "-E", "0.03"],
	"exact": ["--mapping-mode", "exact"],
	"e5": ["-e", "0.05"],
	"e20": ["-e", "0.2"],
	"hamming": ["--alignment-model", "hamming"],
	"hamming-complete": ["--alignment-model", "hamming", "--mapping-mode", "complete"],
	"edit": ["--alignment-model", "edit"],
	"local-never": ["--alignment-local", "never"],
	"local-thresholds": ["--alignment-local-min-identity", "60", "--alignment-local-min-score", "30"],
	"global-thresholds": ["--alignment-global-min-identity", "0.9", "--alignment-global-min-score", "0.6"],
	"full-swg": ["--alignment-force-full-swg"],
	"bandwidth": ["--alignment-max-bandwidth", "0.05"],
	"scaffolding-30": ["--alignment-scaffolding-min-matching_length", "30"],
	"scaffolding-3": ["--alignment-scaffolding-min-matching_length", "3"],
	"end-context": ["--alignment-curation-min-end-context", "4"],
	"bpm": ["--candidate-verification", "BPM"],
	"qgram": ["--qgram-filter", "4,5"],
	"penalties": ["-A", "2", "-O", "8", "-X", "2"],
}

# Pins compared without MAPQ (see the module docstring)
WITHOUT_MAPQ = ("complete-partition", "complete-partition-filter")


@pytest.mark.parametrize("config", sorted(CONFIGS))
def test_b_output_bytes_are_pinned(gem3, genome_index, threads, pin_fastq, config):
	"""Pins the current output (see the module docstring)."""
	res = gem3.run_map(genome_index, args=["-i", "p.fq"] + CONFIGS[config], files={"p.fq": pin_fastq},
		threads=threads, timeout=120)
	digest = digest_without_mapq(res.stdout) if config in WITHOUT_MAPQ else sam_digest(res.stdout)
	assert digest == PINNED[config], config
