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
"""

import pytest

from gem3harness import fastq_record, sam_digest, simulate_se
import area_b


def _pin_input(genome):
	rng = area_b.Rng(41)
	parts = []
	reads = simulate_se(genome, 160, 100, seed=301, sub_rate=0.02, ins_rate=0.004, del_rate=0.004,
		unique_only=False, spacing=4)
	reads += simulate_se(genome, 40, 100, seed=302, subs=6, spacing=5, name_prefix="s")
	reads += simulate_se(genome, 40, 36, seed=303, sub_rate=0.02, unique_only=False, edge=4, spacing=4,
		name_prefix="short")
	reads += simulate_se(genome, 30, 250, seed=304, sub_rate=0.02, ins_rate=0.004, del_rate=0.004,
		unique_only=False, name_prefix="long")
	parts.extend(fastq_record(r.name, r.seq, r.qual) for r in reads)
	for name, seq, info in area_b.chimeras(genome, 16, [35, 45, 60], seed=42):
		parts.append(fastq_record(name, seq, "I" * len(seq)))
	ref = genome["chr1"]
	for i, n in enumerate((6, 12, 18, 24)):
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


# Base 6e155e3, gcc 13 release, -t 1 and 4.
PINNED = {
	"bandwidth": "9248b084566fad9efa7004075b1710fc13eb71d7d3dd9e3729710ca7c3c4ddfd",
	"bpm": "aea7cf011196d060c600725f77e1156d8d7c8281776ec482420fdf119d96f620",
	"complete": "b2c64254c5314de47c1b0bcc32106ab82e32984a3fbbf36921714e3dc2b271a9",
	"complete-E2-s2": "b2c64254c5314de47c1b0bcc32106ab82e32984a3fbbf36921714e3dc2b271a9",
	"complete-E8-s0": "187f73426ec8bfcf1b25f7c96503898420fed120f06b1c8b6ea667ebc8a378b1",
	"complete-all": "4bdb2511f4ae915abcbca3510c1a6d78a08cf53adbb73157c1efa83b04850d33",
	"complete-partition": "a835851718328424cf64e18438056b5e9d8c9c377f798150a76359568e6932f9",
	"complete-partition-filter": "ed6b1b39be61990e3aa8e7e65e8f8ed7ea4e84d04d97b942fb4b4810fbca7395",
	"complete-pure": "b2c64254c5314de47c1b0bcc32106ab82e32984a3fbbf36921714e3dc2b271a9",
	"default": "aea7cf011196d060c600725f77e1156d8d7c8281776ec482420fdf119d96f620",
	"e20": "2da43ab5799d57f841029349f12e410684e3812fc700b60a20ea9fa33e1011ac",
	"e5": "896a3080bca5bf08f7424552f3a1b76f80fa6baa9ab4bf0128c48916c831b1f6",
	"edit": "0e3b0e0a453bef163b0c3a334199340a9883702a95cfacc7b380ddc434f60461",
	"end-context": "e5222cf9534fae7ecd4a06dfd74127f7a88e60a5608ec60b897abd735be4c02e",
	"exact": "76f2aef10e56d6099b9f5c17ac6e2efda4b2d0bb94717d2e0371e3f7d3e645ad",
	"full-swg": "c9e091804f26e132c608fbbc20901f07aae0be8ba584d593dff06b1e567a90d4",
	"global-thresholds": "178aa8b09a6f7864a3a41a32a5b26d777f81d5d4a0635cdb9dee3ba9317896c8",
	"hamming": "0a17359c4076b7828d4d9910c725dde8eff932e8a1537856ef0b76db7c9deda8",
	"hamming-complete": "7d89bc7d0ca0f523a739f4055b3de033a5eb08ac864abd1febb7d6df8bba56af",
	"local-never": "fa47e1d39818e37745e9f3d3bf66124570449d5c14e056b7b6605d7824afebd6",
	"local-thresholds": "b3df2b37ec8ef0dc25d63de97dfa1bc31aaf51eeb68c250b53844ce70c30cd3f",
	"penalties": "b28e9cca6dde8ccfd8c0aaa88257467b8034af62dc8ce8d7f86cc45618e4dfbd",
	"qgram": "aea7cf011196d060c600725f77e1156d8d7c8281776ec482420fdf119d96f620",
	"scaffolding-3": "7805a9910c1da8bf4ac042748e4e0769156cf56560e0e46d803dd0f15f7bafa6",
	"scaffolding-30": "89f60b2d38c077df34d18d662fce8f06896ecb60f3077e47a7f70beb14ed8f06",
	"sensitive": "2460eb826d9736c50f1f01eb2e0123231a44f8aeac191892a795dd054829b1d9",
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
	"bandwidth": ["--alignment-max-bandwidth", "0.1"],
	"scaffolding-30": ["--alignment-scaffolding-min-matching_length", "30"],
	"scaffolding-3": ["--alignment-scaffolding-min-matching_length", "3"],
	"end-context": ["--alignment-curation-min-end-context", "4"],
	"bpm": ["--candidate-verification", "BPM"],
	"qgram": ["--qgram-filter", "4,5"],
	"penalties": ["-A", "2", "-O", "8", "-X", "2"],
}


@pytest.mark.parametrize("config", sorted(CONFIGS))
def test_b_output_bytes_are_pinned(gem3, genome_index, threads, pin_fastq, config):
	"""Pins the current output (see the module docstring)."""
	res = gem3.run_map(genome_index, args=["-i", "p.fq"] + CONFIGS[config], files={"p.fq": pin_fastq},
		threads=threads, timeout=120)
	assert sam_digest(res.stdout) == PINNED[config], config
