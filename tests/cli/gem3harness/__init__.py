# GEM-Mapper v3 (GEM3) test suite.
# Copyright (c) 2026 Jacob Schreiber <jacob.schreiber@umassmed.edu>
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of GEM-Mapper v3 (GEM3), distributed under the GNU General
# Public License version 3 or (at your option) any later version; see LICENSE.
"""Helpers for the gem3-mapper test suite (standard library only)."""

from .genome import Feature, Genome, revcomp, synthetic_genome, write_fasta_records
from .rng import Rng, hash_u64
from .sam import (parse_sam, parse_cigar, cigar_to_string, cigar_query_length, cigar_ref_span,
	cigar_is_valid, md_nm_from_alignment, bisulfite_nm, md_is_wellformed, md_equal, check_record_against_reference,
	strip_variable_header, sam_digest, Sam, SamRecord)
from .simulate import (Truth, SimRead, SimPair, simulate_se, simulate_pe, simulate_bisulfite,
	bisulfite_convert, ADAPTER_R1, ADAPTER_R2)
from .fastx import write_reads, write_interleaved, write_split, fastq_record, fasta_record, compress_bytes
from .runner import Gem3, RunResult, run_process, default_threads
from .build import build_gem3, gem3_from_env, build_c_driver, c_driver_sources, BuildError
