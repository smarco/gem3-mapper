# GEM-Mapper v3 (GEM3) test suite.
# Copyright (c) 2026 Jacob Schreiber <jacob.schreiber@umassmed.edu>
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of GEM-Mapper v3 (GEM3), distributed under the GNU General
# Public License version 3 or (at your option) any later version; see LICENSE.
"""Session fixtures shared by the C and CLI layers (see README.md)."""

import os
import sys

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "cli"))

import pytest  # noqa: E402

from gem3harness import Gem3, gem3_from_env, synthetic_genome, default_threads  # noqa: E402
from gem3harness.build import tree_snapshot  # noqa: E402


@pytest.fixture(scope="session")
def gem3_build():
	"""The build under test (GEM3_BIN_DIR, or an out-of-tree build of this tree).

	``tree_changes`` lists source-tree paths the build added, removed or
	modified (it must stay empty)."""
	before = tree_snapshot()
	build = gem3_from_env()
	after = tree_snapshot()
	build.tree_changes = sorted(
		[("added", p) for p in set(after) - set(before)]
		+ [("removed", p) for p in set(before) - set(after)]
		+ [("modified", p) for p in set(before) & set(after) if before[p] != after[p]])
	return build


@pytest.fixture(scope="session")
def gem3(gem3_build, tmp_path_factory):
	"""A :class:`gem3harness.Gem3` bound to the build; indexes are cached for the session."""
	return Gem3(gem3_build, str(tmp_path_factory.mktemp("gem3")))


@pytest.fixture(scope="session")
def genome():
	"""The default synthetic genome (seed 1, ~240 kbp, five contigs)."""
	return synthetic_genome(seed=1)


@pytest.fixture(scope="session")
def genome_index(gem3, genome):
	return gem3.index(genome)


@pytest.fixture(scope="session")
def genome_bs_index(gem3, genome):
	return gem3.index(genome, bisulfite=True)


@pytest.fixture
def threads():
	"""-t for tests that do not pin it (GEM3_TEST_THREADS, default 1)."""
	return default_threads()


def pytest_report_header(config):
	return ["gem3: GEM3_BIN_DIR=%s GEM3_TEST_THREADS=%s GEM3_CC=%s GEM3_MAKE_TARGET=%s GEM3_MAKE_ARGS=%s"
		% (os.environ.get("GEM3_BIN_DIR", "-"), default_threads(), os.environ.get("GEM3_CC", "-"),
			os.environ.get("GEM3_MAKE_TARGET", "release"), os.environ.get("GEM3_MAKE_ARGS", "-"))]
