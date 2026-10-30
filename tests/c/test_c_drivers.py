# GEM-Mapper v3 (GEM3) test suite.
# Copyright (c) 2026 Jacob Schreiber <jacob.schreiber@umassmed.edu>
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of GEM-Mapper v3 (GEM3), distributed under the GNU General
# Public License version 3 or (at your option) any later version; see LICENSE.
"""Build and run every C driver tests/c/test_*.c (one pytest item each)."""

import os

import pytest

from gem3harness.build import build_c_driver, c_driver_sources
from gem3harness.runner import run_process

DRIVER_TIMEOUT = 120


@pytest.mark.parametrize("source", c_driver_sources(), ids=lambda p: os.path.basename(p)[:-2])
def test_c_driver(gem3_build, source):
	if not gem3_build.has_objects:
		pytest.skip("GEM3_BIN_DIR without GEM3_BUILD_DIR: no objects to link C drivers against")
	exe = build_c_driver(gem3_build, source)
	res = run_process([exe], timeout=DRIVER_TIMEOUT, collect=False)
	out = res.stdout_text
	assert res.returncode == 0 and not res.timed_out, (
		"%s exited %s%s\n%s\nstderr:\n%s" % (os.path.basename(source), res.returncode,
			" (timed out)" if res.timed_out else "", out[-8000:], res.stderr_text[-2000:]))
	assert "not ok" not in out, out[-8000:]
