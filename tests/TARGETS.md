# gem3-mapper test targets

Everything the new suite has to cover, listed before the tests were written.
**Done** is ticked (`x`) and **Test** names the covering tests when a test
exercises the row through its inputs and outputs. An unticked row says why
in **Test**: `indirect` (the function runs under the CLI tests, which check
the output it feeds, not its own results), `partial` (what is and is not
checked), or `not covered` (a known bug, unreachable, or no observable
behaviour). The ticks were made from the area writers' coverage maps
(2026-10-30); long entries are cut at 260 characters.

## How this list was made

- Functions: every `src/**/*.c` was compiled with `gcc -O0 -std=c99 -DNDEBUG
  -fcallgraph-info` (release configuration: no `GEM_DEBUG`, no `GEM_PROFILE`,
  `--disable-cuda`), and the call graph was walked from `main` in
  `src/tools/gem-mapper.c`. The two calls the graph cannot see were added by
  hand: the mapper thread entry points `mapper_se_thread` / `mapper_pe_thread`
  (passed to `pthread_create` in `src/mapper/mapper.c`) and the `qsort`
  comparators `match_trace_cmp_genomic_position` and `region_sort_cmp_value`.
  1022 functions are reachable; static inline functions from headers appear
  under the header that defines them. gem-indexer, gem-retriever and
  gem-constructor are not walked (index building is out of scope).
- Options: the `gem_mapper_options[]` table and the `gem_mapper_parse_arguments_*`
  switch cases in `src/tools/interface/mapper_arguments.c` (all visibility
  levels: user, advanced, developer; `gem-mapper --help=developer` prints them).

## Areas

Four writers work in parallel. Each owns the source files listed under its
area below and writes only its own test files:

| Area | Scope | Test files it owns | Area helpers |
|---|---|---|---|
| A | Index load and input: `src/archive/` (archive, text, locator, sampled SA), `src/fm_index/`, `src/text/` (DNA encoding, sequences, qualities), FASTA/FASTQ parsing and buffered/sliced input with gzip and bzip2 (`src/io/input_*`, `buffered_input_file.c`, `src/system/fm.c`), index memory mapping (`src/system/mm.c`), the region profile (`src/filtering/region_profile/`, `approximate_search_region_profile.c`) | `tests/cli/test_a_*.py`, `tests/c/test_a_*.c` | `tests/cli/area_a.py` |
| B | Search and alignment: `src/approximate_search/`, `src/neighborhood_search/`, `src/filtering/candidates/` and `region/`, `src/align/` (BPM, SWG, k-mer filter), `src/matches/align/` and `scaffold/`, SE search in `src/archive/search/` | `tests/cli/test_b_*.py`, `tests/c/test_b_*.c` | `tests/cli/area_b.py` |
| C | Paired-end, bisulfite, matches and MAPQ: `src/archive/search/archive_search_pe*`, `archive_select*`, `src/archive/score/`, `src/matches/` (classification, MAPQ predictors, CIGAR, paired matches), `src/text/sequence_bisulfite.c` | `tests/cli/test_c_*.py`, `tests/c/test_c_*.c` | `tests/cli/area_c.py` |
| D | Output, CLI, threads and runtime: SAM/MAP writers and output buffers (`src/io/output_*`, `buffered_output_file.c`), `src/mapper/` (threads, I/O handler, statistics), `src/stats/`, `src/system/` (allocators, errors, threads), `src/utils/`, `src/profiler/`, `src/tools/` (CLI) | `tests/cli/test_d_*.py`, `tests/c/test_d_*.c` | `tests/cli/area_d.py` |

Notes on the split, from reading the source:

- `src/search_pipeline/` is reached only from the CUDA mapper
  (`mapper_cuda_se_run`/`mapper_cuda_pe_run`); with `--disable-cuda` the CPU
  mapper (`src/mapper/mapper.c`) reads, searches and writes each read itself, one
  thread per input buffer. D therefore tests thread-count invariance and output
  order through `src/mapper/` and the buffered input/output files, not the
  search pipeline. The CUDA-only rows are listed at the end and need no tests.
- `--candidate-generation` and `--candidate-generation-adaptive` select the
  region-profile strategy, so they belong to A (region profile), not B.
- `-m`/`-M` are implemented by match selection (`archive_select.c`), which sits
  with matches and MAPQ in C.
- `src/text/sequence_bisulfite.c` (read conversion) is bisulfite, so C; the
  rest of `src/text/` is A.
- Paired-end input parsing (`-1/-2`, interleaved `-p`) is A (parsers); what is
  done with the pairs is C.

Shared harness files (`tests/conftest.py`, `tests/cli/gem3harness/`,
`tests/c/gem3_test.h`) are frozen while the writers work: put area-specific
helpers in the area's helper module, and report a needed harness change instead
of making it.

## Options

| Done | Option | Short | Argument | Visibility | Area | Notes | Test |
|---|---|---|---|---|---|---|---|
| x | `--index` | `-I` | required | user | A | memory-mapped vs loaded index | every CLI test; cli/test_a_index.py::test_index_errors (missing, empty, foreign file). gem-mapper always loads the index into memory (mm_bulk_mload_file); there is no memory-mapped mode to compare |
| x | `--input` | `-i` | required | user | A | FASTA/FASTQ, stdin when absent; a trailing '\|' opens a pipe | cli/test_a_input.py (files, gzip, bzip2, CRLF, stdin when absent, 'command \|' pipes, missing file, failing pipe, empty input) |
| x | `--i1` | `-1` | required | user | A | split paired input, end 1 (sets paired mode) | cli/test_a_input.py::test_paired_input_forms_give_identical_output (plain, gzip, bzip2, CRLF, CASAVA ends); unequal files hang (known-bugs/D-unsynced-pairs-hang.md) |
| x | `--i2` | `-2` | required | user | A | split paired input, end 2 | as --i1 |
| x | `--gzip-input` | `-z` | none | user | A | only consulted for stdin | cli/test_a_input.py stdin_gzip; gzip_file_with_j shows files are detected by content |
| x | `--bzip-input` | `-j` | none | user | A | only consulted for stdin | cli/test_a_input.py stdin_bzip2, gzip_file_with_j |
|  | `--input-model` |  | required | developer | A | block_size,num_blocks,buffer_size; the buffer size changes PE output (README); known bug: PE hangs when num_blocks < threads | partial: cli/test_a_input.py: 5 models with 64K-32M blocks give identical SE output; 1-4K blocks refilled at -t 1 (test_small_reused_blocks_single_thread, test_interleaved_pairs_with_small_buffers); bad values rejected. Not tested: small blocks at -t > 1 (known bugs r1-input-model-reader-count-race, A-input-model-small-blocks-drop-lines), one block (hang) |
| x | `--clipping` |  | optional (=) | advanced | A | none, uncalled, masked, fixed,<l>,<r> | cli/test_a_clipping_quality.py (default, bare, uncalled, masked, fixed with 4 values, clip longer than read, none: equal to uncalled on reads without terminal N, consistent alignments with terminal N; bad values) |
| x | `--quality-format` | `-q` | required | advanced | A | ignore, offset-33, offset-64 | cli/test_a_clipping_quality.py::test_quality_offsets (offset-33 = default, offset-64 = same phred values in +64, ignore = FASTA), bad value rejected |
|  | `--candidate-generation` |  | required | developer | A | region-profile strategy: fixed, CKS, OPS, factors, MEM, SMEM, OPP, test; known bug: 'adaptive' and 'adaptive-limited' rejected | partial: cli/test_a_candidate_generation.py: fixed, CKS, OPS, factors, MEM, SMEM, test (truth, pins, equivalences, bad values); not tested: OPP and fixed/factors with region_error > 0 or step 0 (known-bugs/A-candidate-generation-strategies.md), adaptive and... |
| x | `--candidate-generation-adaptive` |  | required | developer | A | <app_threshold>,<app_steps>,<app_dec> | cli/test_a_candidate_generation.py (20,4,1 = default; 5,2,2; 100,8,1; bad values) |
| x | `--mapping-mode` |  | required | user | B | fast, sensitive\|hybrid, customed\|complete, exact, complete-pure, complete-brute-force, complete-partition, complete-partition-filter | test_mode_maps_unique_reads_to_truth (8 values), test_mode_aliases_are_identical, planted tests (complete-brute-force, complete-partition exact sets; fast/sensitive/complete/complete-pure/complete-partition-filter supersets), test_exact_mode_maps_only_exact... |
| x | `--complete-search-error` | `-E` | required | advanced | B |  | test_neighborhood_search_reports_exactly_loci_within_E (E 0-6, count and fraction), pins |
|  | `--complete-strata-after-best` | `-s` | required | advanced | B |  | partial: pins (complete-E8-s0, complete-E2-s2) only; no input found where -s changes the output (filtering finds every planted copy before the strata limit applies) |
| x | `--alignment-max-error` | `-e` | required | user | B | number or fraction | test_max_error_bounds_mapped_reads (gap-affine and Hamming, count and fraction), test_max_error_count_and_fraction_agree, pins |
| x | `--alignment-max-bandwidth` |  | required | advanced | B |  | test_max_bandwidth_limits_gap_length, pins |
| x | `--alignment-global-min-identity` |  | required | user | B |  | test_global_identity_and_score_thresholds |
| x | `--alignment-global-min-score` |  | required | user | B |  | test_global_identity_and_score_thresholds (help default: known-bugs/B-help-defaults.md) |
| x | `--alignment-local` |  | required | user | B | never, if-unmapped | test_local_alignment_of_chimeras (never, if-unmapped) |
| x | `--alignment-local-min-identity` |  | required | user | B |  | test_local_alignment_of_chimeras |
| x | `--alignment-local-min-score` |  | required | user | B |  | test_local_alignment_of_chimeras |
| x | `--alignment-force-full-swg` |  | optional (=) | advanced | B |  | test_alignment_is_optimal_under_penalties, pins |
|  | `--alignment-scaffolding-min-matching_length` |  | required | developer | B |  | partial: pins only (scaffolding-3, scaffolding-30); semantics not derived |
|  | `--alignment-curation` |  | optional (=) | advanced | B |  | not covered: known-bugs/B-alignment-curation-ignored.md (no effect) |
| x | `--alignment-curation-min-end-context` |  | required | developer | B |  | test_curation_clips_mismatches_near_ends (help default: B-help-defaults.md) |
| x | `--candidate-verification` |  | required | developer | B | BPM, chained | test_candidate_verification_strategies, pins |
| x | `--qgram-filter` |  | required | developer | B | <kmer_tiles>,<qgram_length> | test_qgram_filter_never_changes_output, C test_b_kmer_filter; 0 tiles: known-bugs/B-qgram-filter-zero-tiles.md |
| x | `--alignment-model` |  | required | advanced | B | none\|pseudoalignment, hamming, edit\|levenshtein, gap-affine | test_hamming_model_aligns_without_gaps, test_edit_model_reports_least_edit_distance, test_pseudoalignment_places_reads_near_origin (NM/MD: known-bugs/B-pseudoalignment-nm-md.md), planted Hamming tests, pins |
|  | `--gap-affine-penalties` |  | required | user | B | A,B,O,X | not covered: known-bugs/B-mismatch-penalty-sign.md |
| x | `--matching-score` | `-A` | required | advanced | B |  | test_alignment_is_optimal_under_penalties, test_global_identity_and_score_thresholds |
|  | `--mismatch-penalty` | `-B` | required | advanced | B |  | not covered: known-bugs/B-mismatch-penalty-sign.md |
| x | `--gap-open-penalty` | `-O` | required | advanced | B |  | test_alignment_is_optimal_under_penalties |
| x | `--gap-extension-penalty` | `-X` | required | advanced | B |  | test_alignment_is_optimal_under_penalties |
| x | `--paired-end-alignment` | `-p` | none | user | C | interleaved paired input | test_c_paired.py::test_interleaved_input_equals_split_input |
| x | `--min-template-length` | `-l` | required | user | C |  | test_c_paired.py::test_template_length_limits, test_template_length_limit_boundaries, pin l260 |
| x | `--max-template-length` | `-L` | required | user | C |  | same, test_long_templates_are_discordant, pin L400 |
| x | `--discordant-pair-search` |  | required | user | C | always, if-no-concordant, never | test_c_paired.py::test_orientation_and_discordant_search, test_discordant_search_only_without_concordant, test_long_templates_are_discordant, pins always/never; invalid value in test_c_select.py::test_option_value_errors |
|  | `--pair-orientation` |  | required | advanced | C | known bug: hangs on every valid value | not covered: , known bug (pair-orientation-layout-hang.md) |
|  | `--discordant-pair-orientation` |  | required | advanced | C | known bug: hangs | not covered: , known bug (same) |
|  | `--pair-layout` |  | required | advanced | C | known bug: hangs | not covered: , known bug (same) |
|  | `--discordant-pair-layout` |  | required | advanced | C | known bug: hangs | not covered: , known bug (same) |
|  | `--pe-template-length` |  | required | advanced | C | known bug: 3 values rejected, 2 values crash | not covered: , known bug (pe-template-length.md) |
|  | `--bisulfite-read` |  | required | user | C | inferred, 1, 2, interleaved, non-stranded; known bugs: SE '2' and 'interleaved' act as '1', PE 'non-stranded' ignored | partial: inferred, 1, non-stranded single-end in test_c_bisulfite.py, invalid value in test_c_select.py::test_option_value_errors; '2'/'interleaved' single-end and paired 'non-stranded' are known bugs |
| x | `--underconversion_sequence` |  | required | user | C |  | test_c_bisulfite.py::test_control_sequences_in_report |
| x | `--overconversion_sequence` |  | required | user | C |  | same |
| x | `--control_sequence` |  | required | user | C |  | same, and the three defaults |
| x | `--mapq-model` |  | required | advanced | C | none, gem, classify, dump-predictors | none, gem, classify, dump-predictors: test_c_select.py and test_c_paired.py MAPQ tests; invalid value in test_option_value_errors; paired dump-predictors MAPQ is a known bug, C-pe-dump-predictors-uses-classify.md |
| x | `--min-reported-strata` | `-m` | required | advanced | C | cannot combine with -M | test_c_select.py::test_min_reported_strata_reports_whole_strata, test_min_reported_strata_fraction_of_read_length, pins m1/m0.03/rg-m1-full; -m with -M rejected in test_option_value_errors |
| x | `--max-reported-matches` | `-M` | required | user | C | number or 'all'; 0 rejected | test_c_select.py -M tests incl. all and default 5, 0 rejected in test_option_value_errors; test_c_paired.py -M 1/all; MAPQ under -M is a known bug, C-mapq-depends-on-max-reported-matches.md |
| x | `--output` | `-o` | required | user | D |  | out test_output_file_matches_stdout, test_compressed_output[file]; cli errors input/index/end is output, output-dir |
| x | `--gzip-output` |  | none | user | D |  | out test_compressed_output, test_last_compression_option_wins; Cout |
| x | `--bzip-output` |  | none | user | D |  | same |
| x | `--output-model` |  | required | developer | D | buffer_size,num_buffers | out test_output_model_does_not_change_output: 4M,5c 64K,3 4K,2 1K,1 512,4 at -t 1/4, SE/PE; cli errors); zero/non-numeric count known bug (D-zero-or-non-numeric-counts.md |
| x | `--report-file` |  | required | user | D | JSON mapping statistics | rep, SE/PE/bisulfite, thread counts, -r |
| x | `--output-format` | `-F` | required | user | D | MAP, SAM | map; case-insensitive; cli error |
| x | `--sam-compact` |  | optional (=) | user | D | known bug (upstream #37): a space-separated value is ignored | sam: =false single end, =true/bare default); paired =false crashes (known bug D-pe-sam-compact-false-crash.md); space-separated value ignored (known bug optional-argument-ignored.md |
| x | `--sam-read-group-header` | `-r` | required | user | D |  | sam test_read_group SE/PE, @PG escaping; rep ReadGroup; cli errors and 255-char ID |
| x | `--sam-gem-compatible` |  | optional (=) | developer | D |  | sam, map |
|  | `--map-format` |  | required | developer | D | 1, 2, 3 | partial: 2 and 3 (map); 1 is a known bug (D-map-format-1-not-implemented.md); 4 rejected (cli) |
| x | `--threads` | `-t` | required | user | D | number or <x>c (cores) | (out thread-count tests 1-8; Cfmt '<x>c' forms); 0/non-number known bug |
| x | `--max-memory` |  | required | developer | D |  | cli test_system_options_do_not_change_output; error |
| x | `--tmp-folder` |  | required | developer | D |  | cli test_system_options_do_not_change_output |
|  | `--check-alignments` | `-c` | required | developer | B | correct, first-best\|best, all-best, complete, none (archive_check.c) | partial: correct, none (SE and PE); first-best/best/all-best/complete: known-bugs/B-check-alignments-crash.md |
| x | `--profile` |  | optional (=) | developer | D | sum, min, max, mean, sample | cli test_profile x5; error |
| x | `--verbose` | `-v` | optional (=) | user | D | quiet, user, dev\|debug | cli test_verbose_user_log, test_verbose_quiet_and_dev; bare -v; error); space form known bug (#37 |
| x | `--progress-step` |  | required | developer | D |  | cli test_progress_step |
| x | `--version` |  | none | user | D |  | cli test_version, equals @PG VN |
| x | `--help` | `-h` | optional (=) | user | D | user, advanced, developer | cli: three levels against the option table, spellings, defaults text, invalid level |
| x | `(no option)` |  |  |  | D | argument errors: no -I, same file twice, -m with -M, unknown option, positional arguments | cli test_argument_errors: no -I, same files, -m with -M, -M 0, unknown option, CUDA options, unreachable --candidate-drop-off/--pe-extension; no arguments prints usage). Stray positional arguments are a known bug (optional-argument-ignored.md |

`optional (=)`: the value must be attached (`--sam-compact=false`, `-vquiet`);
getopt leaves a space-separated value as an ignored positional argument.

Unreachable from the command line: the parser handles `--candidate-drop-off`
(case 316) and `--pe-extension` (case 405), but neither is in
`gem_mapper_options[]`, so getopt rejects them. The CUDA options (`--gpu`,
`--gpu-devices`, `--gpu-buffers-model`) exist only in `HAVE_CUDA` builds; a
`--disable-cuda` build rejects them as unrecognized (area D tests that).

## Area A: Index load and input (252 functions in 30 files)

### `src/approximate_search/approximate_search_region_profile.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `approximate_search_region_profile_stats` | 59 | partial: runs after every profile in cli/test_a_candidate_generation.py; it only adds GEM_PROFILE counters (no output in release builds; counters reach --profile, area D) |
| x | `approximate_search_region_profile` | 109 | cli/test_a_candidate_generation.py (fixed, CKS, OPS, factors, MEM, SMEM, adaptive via --candidate-generation-adaptive, test branches; truth + pins); OPP branch not tested (known-bugs/A-candidate-generation-strategies.md); adaptive_limited branch unreachable... |

### `src/archive/archive.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `archive_read_mem` | 43 | c/test_a_fm_index.c index_* (check_text/check_bwt/check_decode/check_bsearch/check_rank_mtable/check_2query/check_locator); cli/test_a_index.py::test_index_errors (model check); every CLI test (index load/teardown) |
| x | `archive_read` | 65 | c/test_a_fm_index.c index_* (check_text/check_bwt/check_decode/check_bsearch/check_rank_mtable/check_2query/check_locator); cli/test_a_index.py::test_index_errors (model check); every CLI test (index load/teardown) |
| x | `archive_delete` | 71 | c/test_a_fm_index.c index_* (check_text/check_bwt/check_decode/check_bsearch/check_rank_mtable/check_2query/check_locator); cli/test_a_index.py::test_index_errors (model check); every CLI test (index load/teardown) |
| x | `archive_print` | 92 | cli/test_a_index.py::test_verbose_dev_prints_the_loaded_index, cli/test_a_index.py::test_verbose_dev_run_length_index |

### `src/archive/archive_text.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `archive_text_read_mem` | 64 | c/test_a_fm_index.c index_* (check_text/check_bwt/check_decode/check_bsearch/check_rank_mtable/check_2query/check_locator) (check_text: stored text = FASTA-derived text) |
| x | `archive_text_delete` | 81 | c/test_a_fm_index.c index_* (check_text/check_bwt/check_decode/check_bsearch/check_rank_mtable/check_2query/check_locator) |
| x | `archive_text_get_size` | 91 | cli/test_a_index.py::test_verbose_dev_prints_the_loaded_index (Archive.Text.Size printed; value not checked) |
| x | `archive_text_get_position_strand` | 97 | c/test_a_fm_index.c check_text |
| x | `archive_text_get_projection` | 107 | c/test_a_fm_index.c check_text (reverse-half retrieval vs naive) |
| x | `archive_text_retrieve` | 116 | c/test_a_fm_index.c check_text (4000 ranges, both strands, rc) and index_run_length (RL text); CLI NM/MD checks |
| x | `archive_text_print` | 186 | cli/test_a_index.py::test_verbose_dev_prints_the_loaded_index (RL.Text, Forward.Length checked) |

### `src/archive/archive_text_rl.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `archive_text_rl_encode` | 33 | c/test_a_fm_index.c::rl_encode (hand-worked + 500 random vs naive), index_run_length |
| x | `archive_text_rl_position_translate` | 67 | c/test_a_fm_index.c::index_run_length (every RL position within the forward text vs naive) |
| x | `archive_text_rl_get_run_length` | 95 | c/test_a_fm_index.c::rl_encode |
| x | `archive_text_rl_get_decoded_offset_exl` | 105 | c/test_a_fm_index.c::rl_encode |

### `src/archive/locator.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `locator_intervals_init` | 45 | c/test_a_fm_index.c index_* (check_text/check_bwt/check_decode/check_bsearch/check_rank_mtable/check_2query/check_locator); cli/test_a_index.py::test_verbose_dev_prints_the_loaded_index (intervals equal the FASTA-derived ones) |
| x | `locator_read_mem` | 61 | c/test_a_fm_index.c index_* (check_text/check_bwt/check_decode/check_bsearch/check_rank_mtable/check_2query/check_locator); cli/test_a_index.py::test_verbose_dev_prints_the_loaded_index (intervals equal the FASTA-derived ones) |
| x | `locator_delete` | 85 | c/test_a_fm_index.c index_* (check_text/check_bwt/check_decode/check_bsearch/check_rank_mtable/check_2query/check_locator); cli/test_a_index.py::test_verbose_dev_prints_the_loaded_index (intervals equal the FASTA-derived ones) |
| x | `locator_get_size` | 96 | c/test_a_fm_index.c check_locator (non-zero only); printed by cli/test_a_index.py::test_verbose_dev_prints_the_loaded_index |
| x | `locator_interval_get_tag` | 112 | c/test_a_fm_index.c check_locator, index_run_length |
| x | `locator_lookup_interval_index_4way` | 150 | c/test_a_fm_index.c check_locator (every index position) |
| x | `locator_lookup_interval` | 187 | c/test_a_fm_index.c check_locator |
| x | `locator_lookup_rl_interval_index` | 195 | c/test_a_fm_index.c::index_run_length |
| x | `locator_lookup_rl_interval` | 216 | c/test_a_fm_index.c::index_run_length |
| x | `locator_map` | 225 | c/test_a_fm_index.c check_locator (every forward and reverse position vs FASTA coordinates); cli/test_a_index.py::test_reads_at_contig_ends_and_n_runs |
| x | `locator_inverse_map` | 245 | c/test_a_fm_index.c check_locator (forward strand; gem-mapper only calls it forward, from archive_check) |
| x | `locator_inverse_map_position` | 271 | c/test_a_fm_index.c check_locator |
| x | `locator_interval_print` | 283 | cli/test_a_index.py::test_verbose_dev_prints_the_loaded_index (every printed interval and the interval count checked) |
| x | `locator_print_summary` | 315 | cli/test_a_index.py::test_verbose_dev_prints_the_loaded_index (every printed interval and the interval count checked) |
| x | `locator_print` | 327 | cli/test_a_index.py::test_verbose_dev_prints_the_loaded_index (every printed interval and the interval count checked) |

### `src/archive/sampled_rl.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `sampled_rl_read_mem` | 56 | c/test_a_fm_index.c::index_run_length |
| x | `sampled_rl_delete` | 78 | c/test_a_fm_index.c::index_run_length |
|  | `sampled_rl_get_size` | 87 | partial: cli/test_a_index.py::test_verbose_dev_run_length_index (printed, value not checked) |
| x | `sampled_rl_get_sample` | 97 | c/test_a_fm_index.c::index_run_length (every 100th RL position) |
| x | `sampled_rl_print` | 104 | cli/test_a_index.py::test_verbose_dev_run_length_index (section present) |

### `src/filtering/region_profile/region_profile.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `region_profile_init` | 39 | c/test_a_fm_index.c::region_profile_strategies, c/test_a_fm_index.c::region_profile_helpers |
| x | `region_profile_model_init` | 48 | c/test_a_fm_index.c::region_profile_strategies (adaptive cases start from it); cli/test_a_candidate_generation.py default == adaptive 20,4,1 |
| x | `region_profile_clear` | 63 | c/test_a_fm_index.c::region_profile_strategies, c/test_a_fm_index.c::region_profile_helpers |
| x | `region_profile_inject_mm` | 72 | c/test_a_fm_index.c::region_profile_strategies, c/test_a_fm_index.c::region_profile_helpers |
| x | `region_profile_allocate_regions` | 86 | c/test_a_fm_index.c::region_profile_strategies, c/test_a_fm_index.c::region_profile_helpers |
| x | `region_profile_has_exact_matches` | 103 | c/test_a_fm_index.c::region_profile_helpers |
| x | `region_profile_query_character` | 114 | c/test_a_fm_index.c::region_profile_strategies (rp_query_character vs naive, mtable and BWT paths) |
| x | `region_profile_query_regions` | 179 | c/test_a_fm_index.c::region_profile_strategies (fixed/CKS/factors intervals vs naive) |
| x | `region_profile_extend_last_region` | 190 | c/test_a_fm_index.c::region_profile_strategies (adaptive with allow_zero_regions; intervals vs naive) |
| x | `region_profile_fill_gaps_add` | 215 | c/test_a_fm_index.c::region_profile_helpers (hand-worked) |
| x | `region_profile_fill_gaps` | 226 | c/test_a_fm_index.c::region_profile_helpers (hand-worked) |
| x | `region_profile_merge_small_regions` | 276 | c/test_a_fm_index.c::region_profile_helpers (hand-worked) |
| x | `region_profile_compute_error_limits` | 346 | c/test_a_fm_index.c::region_profile_helpers (hand-worked) |
| x | `region_profile_compute_kmer_frequency` | 363 | c/test_a_fm_index.c::region_profile_strategies (rp_kmer_frequency vs naive mean log4 counts); feeds MAPQ, pinned via cli/test_a_candidate_generation.py |
| x | `region_search_cmp_position` | 446 | c/test_a_fm_index.c::region_profile_helpers |
| x | `region_profile_sort_by_position` | 453 | c/test_a_fm_index.c::region_profile_helpers |
| x | `region_profile_print_region_pretty` | 483 | cli/test_a_candidate_generation.py::test_test_strategy_prints_every_profile (labels only) |
| x | `region_profile_print_pretty` | 524 | cli/test_a_candidate_generation.py::test_test_strategy_prints_every_profile (labels only) |

### `src/filtering/region_profile/region_profile_adaptive.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `region_profile_generator_save_cut_point` | 54 | c/test_a_fm_index.c::region_profile_strategies (rp_adaptive: 4 models x zero regions on/off; intervals, disjointness, threshold) and cli/test_a_candidate_generation.py (default + 3 adaptive models: truth and pins) |
| x | `region_profile_generator_restart` | 59 | c/test_a_fm_index.c::region_profile_strategies (rp_adaptive: 4 models x zero regions on/off; intervals, disjointness, threshold) and cli/test_a_candidate_generation.py (default + 3 adaptive models: truth and pins) |
| x | `region_profile_generator_init` | 73 | c/test_a_fm_index.c::region_profile_strategies (rp_adaptive: 4 models x zero regions on/off; intervals, disjointness, threshold) and cli/test_a_candidate_generation.py (default + 3 adaptive models: truth and pins) |
| x | `region_profile_generator_close_region` | 96 | c/test_a_fm_index.c::region_profile_strategies (rp_adaptive: 4 models x zero regions on/off; intervals, disjointness, threshold) and cli/test_a_candidate_generation.py (default + 3 adaptive models: truth and pins) |
| x | `region_profile_generator_close_profile` | 110 | c/test_a_fm_index.c::region_profile_strategies (rp_adaptive: 4 models x zero regions on/off; intervals, disjointness, threshold) and cli/test_a_candidate_generation.py (default + 3 adaptive models: truth and pins) |
| x | `region_profile_generator_add_character` | 140 | c/test_a_fm_index.c::region_profile_strategies (rp_adaptive: 4 models x zero regions on/off; intervals, disjointness, threshold) and cli/test_a_candidate_generation.py (default + 3 adaptive models: truth and pins) |
| x | `region_profile_generator_disallow_character` | 192 | c/test_a_fm_index.c::region_profile_strategies (rp_adaptive: 4 models x zero regions on/off; intervals, disjointness, threshold) and cli/test_a_candidate_generation.py (default + 3 adaptive models: truth and pins) |
| x | `region_profile_generate_adaptive` | 248 | c/test_a_fm_index.c::region_profile_strategies (rp_adaptive: 4 models x zero regions on/off; intervals, disjointness, threshold) and cli/test_a_candidate_generation.py (default + 3 adaptive models: truth and pins) |
|  | `region_profile_generate_adaptive_limited` | 304 | partial: c/test_a_fm_index.c::region_profile_strategies (rp_adaptive: intervals, disjointness, length/threshold rule); not reachable from the CLI (known-bugs/candidate-generation-adaptive.md) |

### `src/filtering/region_profile/region_profile_fixed.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `region_profile_generate_fixed` | 97 | c/test_a_fm_index.c::region_profile_strategies (rp_fixed: regions, intervals, schedule) and cli/test_a_candidate_generation.py (fixed, overlapping, gapped) |
| x | `region_profile_generate_cks` | 134 | c/test_a_fm_index.c::region_profile_strategies (rp_cks) and cli/test_a_candidate_generation.py (CKS) |
| x | `region_profile_generate_factors` | 174 | c/test_a_fm_index.c::region_profile_strategies (rp_factors) and cli/test_a_candidate_generation.py (factors) |

### `src/filtering/region_profile/region_profile_mem.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `region_profile_compute_mem` | 39 | c/test_a_fm_index.c::region_profile_strategies (rp_mem: one left-maximal region per position) and cli/test_a_candidate_generation.py (MEM, substitution-only reads; pin without MAPQ) |
| x | `region_profile_generate_mem` | 69 | c/test_a_fm_index.c::region_profile_strategies (rp_mem: one left-maximal region per position) and cli/test_a_candidate_generation.py (MEM, substitution-only reads; pin without MAPQ) |
| x | `smem_interval_vector_append` | 115 | c/test_a_fm_index.c::region_profile_strategies (rp_smem: occurring, left-maximal, covering) and cli/test_a_candidate_generation.py (SMEM, substitution-only reads) |
| x | `region_profile_compute_smem` | 127 | c/test_a_fm_index.c::region_profile_strategies (rp_smem: occurring, left-maximal, covering) and cli/test_a_candidate_generation.py (SMEM, substitution-only reads) |
| x | `region_profile_generate_smem` | 203 | c/test_a_fm_index.c::region_profile_strategies (rp_smem: occurring, left-maximal, covering) and cli/test_a_candidate_generation.py (SMEM, substitution-only reads) |

### `src/filtering/region_profile/region_profile_optimum.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `region_profile_generate_ops_fetch_region` | 55 | c/test_a_fm_index.c::region_profile_strategies (rp_ops vs brute-force optimum) and cli/test_a_candidate_generation.py (OPS) |
| x | `region_profile_generate_optimum_fixed` | 76 | c/test_a_fm_index.c::region_profile_strategies (rp_ops) and cli/test_a_candidate_generation.py (OPS) |
|  | `region_profile_optimum_compute_occ` | 181 | partial: c/test_a_fm_index.c::region_profile_strategies (rp_opp: partition and optimum total vs brute force); CLI not tested: OPP misplaces reads (known-bugs/A-candidate-generation-strategies.md) |
|  | `region_profile_generate_optimum_variable` | 208 | partial: c/test_a_fm_index.c::region_profile_strategies (rp_opp); CLI not tested (known bug, as above) |

### `src/filtering/region_profile/region_profile_schedule.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `region_sort_cmp_value` | 36 | c/test_a_fm_index.c::region_profile_strategies (rp_cks selection) |
| x | `region_sort_sort` | 41 | c/test_a_fm_index.c::region_profile_strategies (rp_cks selection) |
| x | `region_profile_schedule_exact_matches` | 50 | c/test_a_fm_index.c::region_profile_strategies (whole-key regions), cli/test_a_candidate_generation.py |
| x | `region_profile_schedule_exact` | 70 | c/test_a_fm_index.c::region_profile_strategies (rp_check_schedule for MEM, OPS, adaptive) |
| x | `region_profile_schedule_exact_best` | 105 | c/test_a_fm_index.c::region_profile_strategies (rp_cks) |
| x | `region_profile_schedule_approximate` | 157 | c/test_a_fm_index.c::region_profile_strategies (rp_check_schedule for fixed, factors) |

### `src/fm_index/bwt/bwt_sbasic.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `bwt_sbasic_read_mem` | 308 | c/test_a_fm_index.c index_* (check_text/check_bwt/check_decode/check_bsearch/check_rank_mtable/check_2query/check_locator) |
| x | `bwt_sbasic_delete` | 331 | c/test_a_fm_index.c index_* (check_text/check_bwt/check_decode/check_bsearch/check_rank_mtable/check_2query/check_locator) |
|  | `bwt_sbasic_get_size` | 347 | partial: cli/test_a_index.py::test_verbose_dev_prints_the_loaded_index (printed, value not checked) |
| x | `bwt_sbasic_get_block_location` | 378 | c/test_a_fm_index.c check_bwt (every position and character A,C,G,T,N vs naive BWT ranks) |
| x | `bwt_sbasic_char_` | 386 | c/test_a_fm_index.c check_bwt (every position and character A,C,G,T,N vs naive BWT ranks) |
| x | `bwt_sbasic_sampling_is_sampled_` | 394 | c/test_a_fm_index.c check_bwt (LF of unsampled positions vs naive; variants agree) and check_decode (sampled path: decode = naive SA at 3 samplings) |
| x | `bwt_sbasic_sampling_erank_` | 397 | c/test_a_fm_index.c check_bwt (LF of unsampled positions vs naive; variants agree) and check_decode (sampled path: decode = naive SA at 3 samplings) |
| x | `bwt_sbasic_erank_` | 403 | c/test_a_fm_index.c check_bwt (every position and character A,C,G,T,N vs naive BWT ranks) |
| x | `bwt_sbasic_precompute_` | 434 | c/test_a_fm_index.c check_bwt (every position and character A,C,G,T,N vs naive BWT ranks) |
| x | `bwt_sbasic_erank` | 472 | c/test_a_fm_index.c check_bwt (every position and character A,C,G,T,N vs naive BWT ranks) |
| x | `bwt_sbasic_prefetch` | 495 | c/test_a_fm_index.c check_bwt (every position and character A,C,G,T,N vs naive BWT ranks) |
| x | `bwt_sbasic_precompute` | 526 | c/test_a_fm_index.c check_bwt (every position and character A,C,G,T,N vs naive BWT ranks) |
| x | `bwt_sbasic_precomputed_erank` | 551 | c/test_a_fm_index.c check_bwt (every position and character A,C,G,T,N vs naive BWT ranks) |
| x | `bwt_sbasic_LF` | 590 | c/test_a_fm_index.c check_bwt (LF of unsampled positions vs naive; variants agree) and check_decode (sampled path: decode = naive SA at 3 samplings) |
| x | `bwt_sbasic_prefetched_LF` | 595 | c/test_a_fm_index.c check_bwt (LF of unsampled positions vs naive; variants agree) and check_decode (sampled path: decode = naive SA at 3 samplings) |
| x | `bwt_sbasic_LF__enc` | 601 | c/test_a_fm_index.c check_bwt (LF of unsampled positions vs naive; variants agree) and check_decode (sampled path: decode = naive SA at 3 samplings) |
| x | `bwt_sbasic_prefetched_LF__enc` | 620 | c/test_a_fm_index.c check_bwt (LF of unsampled positions vs naive; variants agree) and check_decode (sampled path: decode = naive SA at 3 samplings) |
|  | `bwt_sbasic_print_` | 638 | partial: cli/test_a_index.py::test_verbose_dev_prints_the_loaded_index (section present, content not checked) |
|  | `bwt_sbasic_print` | 693 | partial: cli/test_a_index.py::test_verbose_dev_prints_the_loaded_index (section present, content not checked) |

### `src/fm_index/fm_index.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `fm_index_read_mem` | 82 | c/test_a_fm_index.c index_* (check_text/check_bwt/check_decode/check_bsearch/check_rank_mtable/check_2query/check_locator) |
| x | `fm_index_delete` | 99 | c/test_a_fm_index.c index_* (check_text/check_bwt/check_decode/check_bsearch/check_rank_mtable/check_2query/check_locator) |
| x | `fm_index_get_length` | 113 | c/test_a_fm_index.c check_text (= naive text length), index_run_length |
| x | `fm_index_get_proper_length` | 116 | indirectly: adaptive profiles in c/test_a_fm_index.c::region_profile_strategies and cli/test_a_candidate_generation.py |
|  | `fm_index_get_size` | 119 | partial: cli/test_a_index.py::test_verbose_dev_prints_the_loaded_index (printed, value not checked) |
| x | `fm_index_decode` | 128 | c/test_a_fm_index.c check_decode (every BWT position = naive SA; 3 samplings + RL); cli/test_a_index.py::test_index_sampling_does_not_change_output |
| x | `fm_index_print` | 215 | cli/test_a_index.py::test_verbose_dev_prints_the_loaded_index (sampling, packed length, rank table content checked) |

### `src/fm_index/fm_index_query.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `fm_index_2query_init` | 32 | c/test_a_fm_index.c check_2query (random forward/backward extension orders, direct and precomputed: backward = naive range of P, forward = of revcomp(P)); SMEM in cli/test_a_candidate_generation.py |
| x | `fm_index_precompute_2erank` | 44 | c/test_a_fm_index.c check_2query (random forward/backward extension orders, direct and precomputed: backward = naive range of P, forward = of revcomp(P)); SMEM in cli/test_a_candidate_generation.py |
| x | `fm_2erank_elms_upper_acc_erank` | 76 | c/test_a_fm_index.c check_2query (random forward/backward extension orders, direct and precomputed: backward = naive range of P, forward = of revcomp(P)); SMEM in cli/test_a_candidate_generation.py |
| x | `fm_2erank_elms_lower_acc_erank` | 107 | c/test_a_fm_index.c check_2query (random forward/backward extension orders, direct and precomputed: backward = naive range of P, forward = of revcomp(P)); SMEM in cli/test_a_candidate_generation.py |
| x | `fm_index_precomputed_2query_forward` | 138 | c/test_a_fm_index.c check_2query (random forward/backward extension orders, direct and precomputed: backward = naive range of P, forward = of revcomp(P)); SMEM in cli/test_a_candidate_generation.py |
| x | `fm_index_precomputed_2query_backward` | 154 | c/test_a_fm_index.c check_2query (random forward/backward extension orders, direct and precomputed: backward = naive range of P, forward = of revcomp(P)); SMEM in cli/test_a_candidate_generation.py |
| x | `fm_index_2query_forward_query` | 173 | c/test_a_fm_index.c check_2query (random forward/backward extension orders, direct and precomputed: backward = naive range of P, forward = of revcomp(P)); SMEM in cli/test_a_candidate_generation.py |
| x | `fm_index_2query_backward_query` | 186 | c/test_a_fm_index.c check_2query (random forward/backward extension orders, direct and precomputed: backward = naive range of P, forward = of revcomp(P)); SMEM in cli/test_a_candidate_generation.py |
| x | `fm_index_2query_forward_precompute` | 202 | c/test_a_fm_index.c check_2query (random forward/backward extension orders, direct and precomputed: backward = naive range of P, forward = of revcomp(P)); SMEM in cli/test_a_candidate_generation.py |
| x | `fm_index_2query_backward_precompute` | 211 | c/test_a_fm_index.c check_2query (random forward/backward extension orders, direct and precomputed: backward = naive range of P, forward = of revcomp(P)); SMEM in cli/test_a_candidate_generation.py |
| x | `fm_index_2query_precomputed_forward_query` | 220 | c/test_a_fm_index.c check_2query (random forward/backward extension orders, direct and precomputed: backward = naive range of P, forward = of revcomp(P)); SMEM in cli/test_a_candidate_generation.py |
| x | `fm_index_2query_precomputed_backward_query` | 231 | c/test_a_fm_index.c check_2query (random forward/backward extension orders, direct and precomputed: backward = naive range of P, forward = of revcomp(P)); SMEM in cli/test_a_candidate_generation.py |

### `src/fm_index/fm_index_search.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `fm_index_bsearch` | 34 | c/test_a_fm_index.c check_bsearch (6000 keys incl. N-containing, vs naive) and index_run_length |

### `src/fm_index/rank_mtable.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `rank_mtable_init_levels` | 33 | c/test_a_fm_index.c check_rank_mtable; cli/test_a_index.py::test_verbose_dev_prints_the_loaded_index |
| x | `rank_mtable_read_mem` | 62 | c/test_a_fm_index.c check_rank_mtable; cli/test_a_index.py::test_verbose_dev_prints_the_loaded_index |
| x | `rank_mtable_delete` | 80 | c/test_a_fm_index.c check_rank_mtable; cli/test_a_index.py::test_verbose_dev_prints_the_loaded_index |
|  | `rank_mtable_get_size` | 89 | partial: cli/test_a_index.py::test_verbose_dev_prints_the_loaded_index (printed, value not checked) |
| x | `rank_mquery_new` | 95 | c/test_a_fm_index.c check_rank_mtable (every k-mer up to 5, random up to 11, vs naive) |
| x | `rank_mquery_add_char` | 99 | c/test_a_fm_index.c check_rank_mtable (every k-mer up to 5, random up to 11, vs naive) |
| x | `rank_mquery_is_exhausted` | 110 | c/test_a_fm_index.c check_rank_mtable (every k-mer up to 5, random up to 11, vs naive) |
| x | `rank_mtable_prefetch` | 116 | c/test_a_fm_index.c check_rank_mtable (every k-mer up to 5, random up to 11, vs naive) |
| x | `rank_mtable_fetch` | 121 | c/test_a_fm_index.c check_rank_mtable (every k-mer up to 5, random up to 11, vs naive) |
| x | `rank_mtable_print_content_rec` | 132 | cli/test_a_index.py::test_verbose_dev_prints_the_loaded_index (all 84 printed k-mer intervals vs brute force) |
| x | `rank_mtable_print` | 168 | cli/test_a_index.py::test_verbose_dev_prints_the_loaded_index |

### `src/fm_index/sampled_sa.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `sampled_sa_read_mem` | 42 | c/test_a_fm_index.c index_* (check_text/check_bwt/check_decode/check_bsearch/check_rank_mtable/check_2query/check_locator) |
| x | `sampled_sa_delete` | 66 | c/test_a_fm_index.c index_* (check_text/check_bwt/check_decode/check_bsearch/check_rank_mtable/check_2query/check_locator) |
|  | `sampled_sa_get_size` | 164 | partial: cli/test_a_index.py::test_verbose_dev_prints_the_loaded_index (printed, value not checked) |
| x | `sampled_sa_get_sa_sampling_rate` | 167 | c/test_a_fm_index.c check_all (equals the build rate) |
| x | `sampled_sa_get_text_sampling_rate` | 170 | c/test_a_fm_index.c check_all |
|  | `sampled_sa_prefetch_sample` | 173 | partial: indirectly: candidate decoding in every mapping CLI test (prefetch has no observable effect) |
| x | `sampled_sa_get_sample` | 178 | c/test_a_fm_index.c check_decode |
| x | `sampled_sa_print_` | 192 | cli/test_a_index.py::test_verbose_dev_prints_the_loaded_index (text sampling rate and packed array length checked) |
| x | `sampled_sa_print` | 213 | cli/test_a_index.py::test_verbose_dev_prints_the_loaded_index |

### `src/io/buffered_input_file.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `buffered_input_file_new` | 46 | c/test_a_input.c::fastq_plain_all_buffer_sizes, fastq_compressed_streams_and_pipes, buffered_input_forced_lines, fastq_header_forms; cli/test_a_input.py |
| x | `buffered_input_file_close` | 72 | c/test_a_input.c::fastq_plain_all_buffer_sizes, fastq_compressed_streams_and_pipes, buffered_input_forced_lines, fastq_header_forms; cli/test_a_input.py |
| x | `buffered_input_file_get_file_name` | 79 | c/test_a_input.c parse_and_compare (file, 'stream', pipe command) |
|  | `buffered_input_file_get_current_line_num` | 88 | partial: c/test_a_input.c parse_and_compare (differences only: the counter starts uninitialized, known-bugs/A-input-line-counter-uninitialized.md) |
| x | `buffered_input_file_eob` | 91 | c/test_a_input.c::fastq_plain_all_buffer_sizes, fastq_compressed_streams_and_pipes, buffered_input_forced_lines, fastq_header_forms; cli/test_a_input.py |
| x | `buffered_input_file_read_buffer` | 98 | c/test_a_input.c::fastq_plain_all_buffer_sizes, fastq_compressed_streams_and_pipes, buffered_input_forced_lines, fastq_header_forms; cli/test_a_input.py |
| x | `buffered_input_file_load_next_chunk` | 136 | c/test_a_input.c::fastq_plain_all_buffer_sizes, fastq_compressed_streams_and_pipes, buffered_input_forced_lines, fastq_header_forms; cli/test_a_input.py |
| x | `buffered_input_file_get_line` | 165 | c/test_a_input.c::fastq_plain_all_buffer_sizes, fastq_compressed_streams_and_pipes, buffered_input_forced_lines, fastq_header_forms; cli/test_a_input.py |
| x | `buffered_input_file_reload` | 206 | c/test_a_input.c::fastq_plain_all_buffer_sizes, fastq_compressed_streams_and_pipes, buffered_input_forced_lines, fastq_header_forms; cli/test_a_input.py |
| x | `buffered_input_file_attach_buffered_output` | 232 | indirectly: every CLI mapping test (output buffers follow input blocks) |

### `src/io/input_buffer.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `input_buffer_new` | 39 | c/test_a_input.c::fastq_plain_all_buffer_sizes, fastq_compressed_streams_and_pipes, buffered_input_forced_lines, fastq_header_forms; cli/test_a_input.py |
| x | `input_buffer_delete` | 55 | c/test_a_input.c::fastq_plain_all_buffer_sizes, fastq_compressed_streams_and_pipes, buffered_input_forced_lines, fastq_header_forms; cli/test_a_input.py |
| x | `input_buffer_annotate_lines` | 63 | c/test_a_input.c::fastq_plain_all_buffer_sizes, fastq_compressed_streams_and_pipes, buffered_input_forced_lines, fastq_header_forms; cli/test_a_input.py |
| x | `input_buffer_get_num_lines` | 82 | c/test_a_input.c::fastq_plain_all_buffer_sizes, fastq_compressed_streams_and_pipes, buffered_input_forced_lines, fastq_header_forms; cli/test_a_input.py |

### `src/io/input_fasta.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `input_fasta_parse_sequence_tag_chomp_end_info` | 57 | c/test_a_input.c::fastq_header_forms and random headers in parse_and_compare; cli/test_a_input.py::test_read_names |
| x | `input_fasta_parse_sequence_tag_extra` | 76 | c/test_a_input.c::fastq_header_forms and random headers in parse_and_compare; cli/test_a_input.py::test_read_names |
| x | `input_fasta_parse_sequence_tag` | 141 | c/test_a_input.c::fastq_header_forms and random headers in parse_and_compare; cli/test_a_input.py::test_read_names |
| x | `input_fasta_parse_sequence_read` | 174 | c/test_a_input.c (IUPAC/odd characters -> N); cli/test_a_input.py::test_read_characters; last character unchecked is a known bug (known-bugs/A-fastq-last-character-unchecked.md) |
| x | `input_fasta_parse_sequence_qualities` | 207 | c/test_a_input.c; cli/test_a_input.py::test_malformed_fastq_is_rejected[low_quality] |
| x | `input_fasta_parse_sequence_components` | 231 | c/test_a_input.c::fastq_plain_all_buffer_sizes, fastq_compressed_streams_and_pipes, buffered_input_forced_lines, fastq_header_forms; cli/test_a_input.py; malformed cases in cli/test_a_input.py |
| x | `input_fasta_parse_sequence` | 261 | c/test_a_input.c::fastq_plain_all_buffer_sizes, fastq_compressed_streams_and_pipes, buffered_input_forced_lines, fastq_header_forms; cli/test_a_input.py |
| x | `input_fasta_parser_prompt_error` | 288 | cli/test_a_input.py::test_malformed_fastq_is_rejected (all four error kinds; line numbers not checked) |

### `src/io/input_file_sliced.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `input_file_sliced_initialize` | 32 | c/test_a_input.c::fastq_plain_all_buffer_sizes, fastq_compressed_streams_and_pipes, buffered_input_forced_lines, fastq_header_forms; cli/test_a_input.py |
| x | `input_stream_sliced_open` | 58 | c/test_a_input.c::fastq_compressed_streams_and_pipes; cli stdin tests |
| x | `input_gzip_stream_sliced_open` | 70 | c/test_a_input.c::fastq_compressed_streams_and_pipes; cli stdin_gzip (-z) |
| x | `input_bzip_stream_sliced_open` | 82 | c/test_a_input.c::fastq_compressed_streams_and_pipes; cli stdin_bzip2 (-j) |
| x | `input_file_sliced_open` | 94 | c/test_a_input.c::fastq_plain_all_buffer_sizes, fastq_compressed_streams_and_pipes, buffered_input_forced_lines, fastq_header_forms; cli/test_a_input.py |
| x | `input_file_sliced_popen` | 106 | c/test_a_input.c (pipe); cli pipe, pipe_gzip, failing pipe |
| x | `input_file_sliced_close` | 118 | c/test_a_input.c::fastq_plain_all_buffer_sizes, fastq_compressed_streams_and_pipes, buffered_input_forced_lines, fastq_header_forms; cli/test_a_input.py |
| x | `input_file_sliced_get_file_name` | 135 | c/test_a_input.c parse_and_compare |
| x | `input_file_sliced_get_next_id` | 140 | c/test_a_input.c::fastq_plain_all_buffer_sizes, fastq_compressed_streams_and_pipes, buffered_input_forced_lines, fastq_header_forms; cli/test_a_input.py (blocks of 64 B-1 MB; small blocks only with large buffers or forced lines, known-bugs/A-input-model-sma... |
| x | `input_file_sliced_request_buffer_empty` | 149 | c/test_a_input.c::fastq_plain_all_buffer_sizes, fastq_compressed_streams_and_pipes, buffered_input_forced_lines, fastq_header_forms; cli/test_a_input.py (blocks of 64 B-1 MB; small blocks only with large buffers or forced lines, known-bugs/A-input-model-sma... |
| x | `input_file_sliced_find_buffer` | 167 | c/test_a_input.c::fastq_plain_all_buffer_sizes, fastq_compressed_streams_and_pipes, buffered_input_forced_lines, fastq_header_forms; cli/test_a_input.py (blocks of 64 B-1 MB; small blocks only with large buffers or forced lines, known-bugs/A-input-model-sma... |
| x | `input_file_sliced_input_buffer_read` | 185 | c/test_a_input.c::fastq_plain_all_buffer_sizes, fastq_compressed_streams_and_pipes, buffered_input_forced_lines, fastq_header_forms; cli/test_a_input.py (blocks of 64 B-1 MB; small blocks only with large buffers or forced lines, known-bugs/A-input-model-sma... |
| x | `input_file_sliced_input_buffer_process` | 193 | c/test_a_input.c::fastq_plain_all_buffer_sizes, fastq_compressed_streams_and_pipes, buffered_input_forced_lines, fastq_header_forms; cli/test_a_input.py (blocks of 64 B-1 MB; small blocks only with large buffers or forced lines, known-bugs/A-input-model-sma... |
| x | `input_file_sliced_process` | 210 | c/test_a_input.c::fastq_plain_all_buffer_sizes, fastq_compressed_streams_and_pipes, buffered_input_forced_lines, fastq_header_forms; cli/test_a_input.py (blocks of 64 B-1 MB; small blocks only with large buffers or forced lines, known-bugs/A-input-model-sma... |
| x | `input_file_sliced_discard_exhausted_buffers` | 223 | c/test_a_input.c::fastq_plain_all_buffer_sizes, fastq_compressed_streams_and_pipes, buffered_input_forced_lines, fastq_header_forms; cli/test_a_input.py (blocks of 64 B-1 MB; small blocks only with large buffers or forced lines, known-bugs/A-input-model-sma... |
| x | `input_file_sliced_read_lines` | 249 | c/test_a_input.c::fastq_plain_all_buffer_sizes, fastq_compressed_streams_and_pipes, buffered_input_forced_lines, fastq_header_forms; cli/test_a_input.py (blocks of 64 B-1 MB; small blocks only with large buffers or forced lines, known-bugs/A-input-model-sma... |
| x | `input_file_sliced_input_buffer_get_current` | 299 | c/test_a_input.c::fastq_plain_all_buffer_sizes, fastq_compressed_streams_and_pipes, buffered_input_forced_lines, fastq_header_forms; cli/test_a_input.py (blocks of 64 B-1 MB; small blocks only with large buffers or forced lines, known-bugs/A-input-model-sma... |

### `src/io/input_text.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `input_text_parse_integer` | 41 | c/test_a_input.c::input_text_parsers (hand-worked values) |
| x | `input_text_parse_double` | 74 | c/test_a_input.c::input_text_parsers (hand-worked values) |
| x | `input_text_parse_size` | 145 | c/test_a_input.c::input_text_parsers (hand-worked values) |
| x | `input_text_parse_csv_arguments` | 187 | c/test_a_input.c::input_text_parsers (hand-worked values) |
| x | `input_text_parse_csv_array` | 206 | c/test_a_input.c::input_text_parsers (hand-worked values) |
| x | `input_text_parse_extended_uint64` | 221 | c/test_a_input.c::input_text_parsers (hand-worked values) |
| x | `input_text_parse_extended_double` | 249 | c/test_a_input.c::input_text_parsers (hand-worked values) |
| x | `input_text_parse_extended_bool` | 263 | c/test_a_input.c::input_text_parsers (hand-worked values) |
| x | `input_text_parse_count_colons_in_field` | 277 | c/test_a_input.c::input_text_parsers (hand-worked values) |

### `src/system/fm.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `fm_initialize` | 124 | c/test_a_input.c::file_and_memory_managers (plain/gzip/bzip2 byte-exact), parse tests; cli gzip/bzip2 inputs |
| x | `fm_check_file_type` | 163 | c/test_a_input.c::file_and_memory_managers (plain/gzip/bzip2 byte-exact), parse tests; cli gzip/bzip2 inputs |
| x | `fm_open_file` | 188 | c/test_a_input.c::file_and_memory_managers (plain/gzip/bzip2 byte-exact), parse tests; cli gzip/bzip2 inputs |
| x | `fm_open_popen` | 218 | c/test_a_input.c (pipe); cli pipe tests |
| x | `fm_open_FILE` | 271 | c/test_a_input.c (stream); cli stdin |
| x | `fm_open_gzFILE` | 293 | c/test_a_input.c (gzip stream); cli -z |
| x | `fm_open_bzFILE` | 317 | c/test_a_input.c (bzip2 stream); cli -j; multi-stream bzip2 is a known bug (known-bugs/A-bzip2-multistream-truncated.md) |
| x | `fm_close` | 342 | c/test_a_input.c::file_and_memory_managers (plain/gzip/bzip2 byte-exact), parse tests; cli gzip/bzip2 inputs |
| x | `fm_eof` | 387 | c/test_a_input.c::file_and_memory_managers (plain/gzip/bzip2 byte-exact), parse tests; cli gzip/bzip2 inputs |
| x | `fm_get_file_name` | 415 | c/test_a_input.c::file_and_memory_managers (plain/gzip/bzip2 byte-exact), parse tests; cli gzip/bzip2 inputs |
| x | `fm_read_mem` | 595 | c/test_a_input.c::file_and_memory_managers (plain/gzip/bzip2 byte-exact), parse tests; cli gzip/bzip2 inputs |
|  | `fm_prefetch_next` | 651 | partial: called on every block read in the parse tests; no observable effect |
| x | `fm_write_mem` | 680 | c/test_a_input.c::packed_integer_array_roundtrip (written with fm_write_uint64, read back); the C index fixtures write .gem files through it |
| x | `fm_bulk_read_fd` | 713 | c/test_a_input.c::file_and_memory_managers |
| x | `fm_bulk_read_file` | 729 | c/test_a_input.c::file_and_memory_managers (with and without offset) |
| x | `gem_access` | 768 | c/test_a_input.c::file_and_memory_managers |
| x | `vfmprintf` | 795 | c/test_a_input.c::file_and_memory_managers (fmprintf to a plain file; the gzip/bzip2 output branch is area D's --gzip-output/--bzip-output) |

### `src/system/mm.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `mm_get_tmp_folder` | 109 | c/test_a_input.c::file_and_memory_managers |
| x | `mm_set_tmp_folder` | 112 | c/test_a_input.c::file_and_memory_managers |
| x | `mm_malloc_nothrow` | 121 | c/test_a_input.c::file_and_memory_managers |
| x | `mm_malloc_` | 133 | indirectly: every test |
| x | `mm_realloc_nothrow` | 149 | c/test_a_input.c::file_and_memory_managers |
|  | `mm_realloc` | 154 | partial: reached when string buffers grow (src/utils/string_buffer.c); not checked directly |
| x | `mm_free` | 161 | indirectly: every test |
| x | `mm_bulk_malloc` | 169 | c/test_a_input.c::file_and_memory_managers (zero-filled with init_mem) |
| x | `mm_bulk_mmalloc` | 188 | c/test_a_input.c::file_and_memory_managers |
| x | `mm_bulk_mmalloc_temp` | 229 | c/test_a_input.c::file_and_memory_managers (in a set temporary folder) |
| x | `mm_bulk_free` | 271 | c/test_a_input.c::file_and_memory_managers |
| x | `mm_bulk_mload_file` | 351 | c/test_a_input.c::file_and_memory_managers (bytes equal the file); index loading everywhere |
| x | `mm_get_mem` | 367 | c/test_a_input.c::file_and_memory_managers (cursor) |
| x | `mm_skip_align_4KB` | 502 | c/test_a_input.c::file_and_memory_managers |
| x | `mm_read_uint64` | 520 | c/test_a_input.c::file_and_memory_managers |
| x | `mm_read_mem` | 544 | c/test_a_input.c::file_and_memory_managers |
|  | `mm_get_stat_meminfo` | 618 | partial: c/test_a_input.c::file_and_memory_managers (through the getters: positive, total <= MemTotal) |
|  | `mm_get_mem_available_cached` | 651 | partial: c/test_a_input.c::file_and_memory_managers (sign and bound only; values change while running) |
|  | `mm_get_mem_available_free` | 656 | partial: c/test_a_input.c::file_and_memory_managers (sign and bound only; values change while running) |
|  | `mm_get_mem_available_total` | 661 | partial: c/test_a_input.c::file_and_memory_managers (sign and bound only; values change while running) |

### `src/text/dna_text.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `dna_text_read_mem` | 226 | c/test_a_fm_index.c index_* (check_text/check_bwt/check_decode/check_bsearch/check_rank_mtable/check_2query/check_locator) |
| x | `dna_text_delete` | 244 | c/test_a_fm_index.c index_* (check_text/check_bwt/check_decode/check_bsearch/check_rank_mtable/check_2query/check_locator) |
| x | `dna_text_get_length` | 298 | c/test_a_fm_index.c check_text |
|  | `dna_text_get_size` | 306 | partial: cli/test_a_index.py::test_verbose_dev_prints_the_loaded_index (printed, value not checked) |
| x | `dna_text_retrieve_sequence` | 325 | c/test_a_fm_index.c check_text (stored forward text = FASTA-derived) |
| x | `dna_text_print` | 341 | cli/test_a_index.py::test_verbose_dev_prints_the_loaded_index (Text.Length checked) |

### `src/text/sequence.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `sequence_init` | 31 | c/test_a_input.c parse_and_compare, fastq_header_forms |
| x | `sequence_clear` | 49 | c/test_a_input.c parse_and_compare, fastq_header_forms |
| x | `sequence_get_length` | 67 | c/test_a_input.c parse_and_compare, fastq_header_forms |
| x | `sequence_get_tag` | 70 | c/test_a_input.c parse_and_compare, fastq_header_forms |
| x | `sequence_get_read` | 79 | c/test_a_input.c parse_and_compare, fastq_header_forms |
| x | `sequence_get_qualities` | 88 | c/test_a_input.c parse_and_compare, fastq_header_forms |
| x | `sequence_has_qualities` | 97 | c/test_a_input.c parse_and_compare, fastq_header_forms |
| x | `sequence_has_casava_tag` | 100 | c/test_a_input.c parse_and_compare, fastq_header_forms |
| x | `sequence_has_extra_tag` | 103 | c/test_a_input.c parse_and_compare, fastq_header_forms |
| x | `sequence_get_end_info` | 106 | c/test_a_input.c parse_and_compare, fastq_header_forms |
| x | `sequence_print` | 152 | cli/test_a_input.py::test_malformed_fastq_is_rejected[short_qualities] (record in the error report) |

### `src/text/sequence_qualities_model.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `sequence_qualities_model_process` | 31 | c/test_a_input.c::qualities_model (hand-worked); cli/test_a_clipping_quality.py::test_quality_offsets |

### `src/text/text_trace.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `text_trace_destroy` | 30 | c/test_a_input.c::text_trace_padding |
| x | `text_trace_compose_padded_text` | 44 | c/test_a_input.c::text_trace_padding (hand-worked) |

### `src/utils/packed_integer_array.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `packed_integer_array_delete` | 48 | c/test_a_input.c::packed_integer_array_roundtrip (1-64 bits); sampled SA decoding in c/test_a_fm_index.c |
| x | `packed_integer_array_read_mem` | 71 | c/test_a_input.c::packed_integer_array_roundtrip (1-64 bits); sampled SA decoding in c/test_a_fm_index.c |
| x | `packed_integer_array_get_size` | 270 | c/test_a_input.c::packed_integer_array_roundtrip (1-64 bits); sampled SA decoding in c/test_a_fm_index.c |
| x | `packed_integer_array_prefetch` | 281 | c/test_a_input.c::packed_integer_array_roundtrip (1-64 bits); sampled SA decoding in c/test_a_fm_index.c |
| x | `packed_integer_array_load` | 289 | c/test_a_input.c::packed_integer_array_roundtrip (1-64 bits); sampled SA decoding in c/test_a_fm_index.c |
| x | `packed_integer_array_print` | 317 | cli/test_a_index.py::test_verbose_dev_prints_the_loaded_index (PackedArray.Length checked) |

## Area B: Search and alignment (334 functions in 61 files)

### `include/utils/vector_sort.h`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `vector_sort_selection_align_distance` | 64 | C test_b_sort_queue (template instantiated on int64 vs qsort; gem's instances by order and permutation) |
| x | `vector_sort_selection_filtering_positions` | 64 | C test_b_sort_queue (template instantiated on int64 vs qsort; gem's instances by order and permutation) |
| x | `vector_sort_selection_match_alignment_region` | 64 | C test_b_sort_queue (template instantiated on int64 vs qsort; gem's instances by order and permutation) |
| x | `vector_sort_selection_region_search_by_position` | 64 | C test_b_sort_queue (template instantiated on int64 vs qsort; gem's instances by order and permutation) |
| x | `vector_sort_selection_scaffold_coverage` | 64 | C test_b_sort_queue (template instantiated on int64 vs qsort; gem's instances by order and permutation) |
| x | `vector_sort_heapify_align_distance` | 84 | C test_b_sort_queue (template instantiated on int64 vs qsort; gem's instances by order and permutation) |
| x | `vector_sort_heapify_filtering_positions` | 84 | C test_b_sort_queue (template instantiated on int64 vs qsort; gem's instances by order and permutation) |
| x | `vector_sort_heapify_match_alignment_region` | 84 | C test_b_sort_queue (template instantiated on int64 vs qsort; gem's instances by order and permutation) |
| x | `vector_sort_heapify_region_search_by_position` | 84 | C test_b_sort_queue (template instantiated on int64 vs qsort; gem's instances by order and permutation) |
| x | `vector_sort_heapify_scaffold_coverage` | 84 | C test_b_sort_queue (template instantiated on int64 vs qsort; gem's instances by order and permutation) |
| x | `vector_sort_heapsort_align_distance` | 102 | C test_b_sort_queue (template instantiated on int64 vs qsort; gem's instances by order and permutation) |
| x | `vector_sort_heapsort_filtering_positions` | 102 | C test_b_sort_queue (template instantiated on int64 vs qsort; gem's instances by order and permutation) |
| x | `vector_sort_heapsort_match_alignment_region` | 102 | C test_b_sort_queue (template instantiated on int64 vs qsort; gem's instances by order and permutation) |
| x | `vector_sort_heapsort_region_search_by_position` | 102 | C test_b_sort_queue (template instantiated on int64 vs qsort; gem's instances by order and permutation) |
| x | `vector_sort_heapsort_scaffold_coverage` | 102 | C test_b_sort_queue (template instantiated on int64 vs qsort; gem's instances by order and permutation) |
| x | `vector_sort_quicksort_partition_align_distance` | 118 | C test_b_sort_queue (template instantiated on int64 vs qsort; gem's instances by order and permutation) |
| x | `vector_sort_quicksort_partition_filtering_positions` | 118 | C test_b_sort_queue (template instantiated on int64 vs qsort; gem's instances by order and permutation) |
| x | `vector_sort_quicksort_partition_match_alignment_region` | 118 | C test_b_sort_queue (template instantiated on int64 vs qsort; gem's instances by order and permutation) |
| x | `vector_sort_quicksort_partition_region_search_by_position` | 118 | C test_b_sort_queue (template instantiated on int64 vs qsort; gem's instances by order and permutation) |
| x | `vector_sort_quicksort_partition_scaffold_coverage` | 118 | C test_b_sort_queue (template instantiated on int64 vs qsort; gem's instances by order and permutation) |
| x | `buffer_sort_align_distance` | 135 | C test_b_sort_queue (template instantiated on int64 vs qsort; gem's instances by order and permutation) |
| x | `buffer_sort_filtering_positions` | 135 | C test_b_sort_queue (template instantiated on int64 vs qsort; gem's instances by order and permutation) |
| x | `buffer_sort_match_alignment_region` | 135 | C test_b_sort_queue (template instantiated on int64 vs qsort; gem's instances by order and permutation) |
| x | `buffer_sort_region_search_by_position` | 135 | C test_b_sort_queue (template instantiated on int64 vs qsort; gem's instances by order and permutation) |
| x | `buffer_sort_scaffold_coverage` | 135 | C test_b_sort_queue (template instantiated on int64 vs qsort; gem's instances by order and permutation) |
| x | `vector_sort_align_distance` | 176 | C test_b_sort_queue (template instantiated on int64 vs qsort; gem's instances by order and permutation) |
| x | `vector_sort_filtering_positions` | 176 | C test_b_sort_queue (template instantiated on int64 vs qsort; gem's instances by order and permutation) |
| x | `vector_sort_scaffold_coverage` | 176 | C test_b_sort_queue (template instantiated on int64 vs qsort; gem's instances by order and permutation) |

### `src/align/align_bpm.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `align_bpm_compute_matrix` | 39 | C test_b_bpm_align (brute-force edit distance) |
| x | `align_bpm_backtrace_matrix` | 182 | C test_b_bpm_align (brute-force edit distance) |
| x | `align_bpm_match` | 292 | C test_b_bpm_align (brute-force edit distance) |

### `src/align/align_bpm_distance.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `bpm_reset_search_cutoff` | 59 | C test_bpm_distance (pilot) and C test_b_bpm_align |
| x | `bpm_compute_edit_distance` | 172 | C test_bpm_distance (pilot) and C test_b_bpm_align |
|  | `bpm_compute_edit_distance_all` | 283 | partial: paired-end only; executed by CLI test_check_alignments_correct_paired and the pilot PE tests (recovery extension runs), but no input was found where the extension changes the output (mates with 10 spaced substitutions and no exact 10-mer are not re... |

### `src/align/align_bpm_pattern.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `bpm_pattern_compile` | 47 | C test_bpm_distance (pilot), C test_b_bpm_align (tiles vs substring patterns) |
| x | `bpm_pattern_compile_tiles` | 122 | C test_b_bpm_align (brute-force edit distance) |

### `src/align/align_kmer_filter_nway.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `kmer_counting_compile_nway_allocate_tables` | 70 | C test_b_kmer_filter (hand values, bound <= brute-force distance); CLI test_qgram_filter_never_changes_output |
| x | `kmer_counting_compile_nway` | 95 | C test_b_kmer_filter (hand values, bound <= brute-force distance); CLI test_qgram_filter_never_changes_output |
| x | `kmer_counting_prepare_tiling` | 177 | C test_b_kmer_filter (hand values, bound <= brute-force distance); CLI test_qgram_filter_never_changes_output |
| x | `kmer_counting_increment_kmer_count` | 242 | C test_b_kmer_filter (hand values, bound <= brute-force distance); CLI test_qgram_filter_never_changes_output |
| x | `kmer_counting_min_bound_nway_ks1` | 288 | C test_b_kmer_filter (hand values, bound <= brute-force distance); CLI test_qgram_filter_never_changes_output |
| x | `kmer_counting_min_bound_nway_ks2` | 341 | C test_b_kmer_filter (hand values, bound <= brute-force distance); CLI test_qgram_filter_never_changes_output |
| x | `kmer_counting_min_bound_nway_ks3` | 390 | C test_b_kmer_filter (hand values, bound <= brute-force distance); CLI test_qgram_filter_never_changes_output |
| x | `kmer_counting_min_bound_nway_ks4` | 442 | C test_b_kmer_filter (hand values, bound <= brute-force distance); CLI test_qgram_filter_never_changes_output |
| x | `kmer_counting_min_bound_nway_ksn` | 500 | C test_b_kmer_filter (hand values, bound <= brute-force distance); CLI test_qgram_filter_never_changes_output |
| x | `kmer_counting_min_bound_nway` | 573 | C test_b_kmer_filter (hand values, bound <= brute-force distance); CLI test_qgram_filter_never_changes_output |

### `src/align/align_swg.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `align_swg_traceback` | 33 | C test_b_swg (brute-force Gotoh, hand values, --help defaults) |
| x | `align_swg_allocate_table` | 149 | C test_b_swg (brute-force Gotoh, hand values, --help defaults) |
| x | `align_swg_init_table` | 171 | C test_b_swg (brute-force Gotoh, hand values, --help defaults) |
| x | `align_swg_init_table_banded` | 190 | C test_b_swg (brute-force Gotoh, hand values, --help defaults) |
|  | `align_swg_allocate__init_column_banded` | 233 | indirect: CLI test_b_* (every mapping) |
| x | `align_swg_base` | 287 | C test_b_swg (brute-force Gotoh, hand values, --help defaults) |
| x | `align_swg` | 429 | C test_b_swg (brute-force Gotoh, hand values, --help defaults) |

### `src/align/align_swg_banded.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `align_swg_banded_compute_column` | 33 | C test_b_swg (brute-force Gotoh, hand values, --help defaults) |
|  | `align_swg_banded_compute_column__max_score` | 79 | indirect: CLI test_b_* (every mapping) |
| x | `align_swg_banded` | 136 | C test_b_swg (brute-force Gotoh, hand values, --help defaults) |
|  | `align_swg_banded_extend` | 216 | indirect: CLI test_b_* (every mapping) |

### `src/align/align_swg_score.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `align_swg_score_deletion` | 32 | C test_b_swg (brute-force Gotoh, hand values, --help defaults) |
| x | `align_swg_score_insertion` | 37 | C test_b_swg (brute-force Gotoh, hand values, --help defaults) |
| x | `align_swg_score_mismatch` | 42 | C test_b_swg (brute-force Gotoh, hand values, --help defaults) |
| x | `align_swg_score_match` | 45 | C test_b_swg (brute-force Gotoh, hand values, --help defaults) |
| x | `align_swg_score_cigar_element` | 48 | C test_b_swg (brute-force Gotoh, hand values, --help defaults) |
| x | `align_swg_score_cigar_excluding_clipping` | 106 | C test_b_swg (brute-force Gotoh, hand values, --help defaults) |
| x | `align_swg_score_compute_max_score_bound` | 138 | C test_b_swg (brute-force Gotoh, hand values, --help defaults) |
| x | `align_swg_score_compute_max_edit_bound` | 150 | C test_b_swg (brute-force Gotoh, hand values, --help defaults) |

### `src/align/alignment.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `alignment_init` | 37 | C test_b_bpm_align (brute-force edit distance) |
| x | `alignment_verify_edit_kmer` | 129 | C test_b_kmer_filter (hand values, bound <= brute-force distance); CLI test_qgram_filter_never_changes_output |
| x | `alignment_verify_edit_bpm_tile` | 161 | C test_b_bpm_align (brute-force edit distance) |
| x | `alignment_verify_edit_bpm` | 177 | C test_b_bpm_align (brute-force edit distance) |
|  | `alignment_check` | 242 | partial: CLI test_check_alignments_correct(_paired) (-c correct finds no fault on correct output; the failure branch is not reached) |
|  | `alignment_print_pretty` | 335 | not covered: debug/-c failure printing only |

### `src/align/pattern/pattern.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `pattern_init_encode` | 36 | indirect: CLI test_b_* (every mapping) |
|  | `pattern_init_encode_rl` | 112 | indirect: CLI test_b_* (every mapping) |
|  | `pattern_init` | 138 | indirect: CLI test_b_* (every mapping) |
| x | `pattern_tiling_init` | 230 | C test_b_bpm_align (brute-force edit distance) |
| x | `pattern_tiling_next` | 257 | C test_b_bpm_align (brute-force edit distance) |

### `src/align/pattern/pattern_tiled.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `pattern_tiled_compile` | 31 | C test_b_bpm_align (brute-force edit distance) |

### `src/approximate_search/approximate_search.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `approximate_search_init` | 67 | indirect: CLI test_b_mapping_modes (every mode vs truth, planted genome vs brute-force loci within -E) |
|  | `approximate_search_prepare` | 104 | indirect: CLI test_b_mapping_modes (every mode vs truth, planted genome vs brute-force loci within -E) |
|  | `approximate_search` | 154 | indirect: CLI test_b_mapping_modes (every mode vs truth, planted genome vs brute-force loci within -E) |

### `src/approximate_search/approximate_search_filtering_adaptive.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `as_filtering_control_begin` | 34 | indirect: gcov 0%%, inlined on the main path): CLI test_b_mapping_modes (every mode vs truth, planted genome vs brute-force loci within -E |
|  | `as_filtering_control_filtering_adaptive_next_state` | 43 | indirect: CLI test_b_mapping_modes (every mode vs truth, planted genome vs brute-force loci within -E) |
|  | `approximate_search_filtering_adaptive` | 68 | indirect: CLI test_b_mapping_modes (every mode vs truth, planted genome vs brute-force loci within -E) |

### `src/approximate_search/approximate_search_generate_candidates.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `approximate_search_generate_candidates_limit_exact_matches` | 59 | not covered: needs more exact matches than the limit (gcov: not executed) |
|  | `approximate_search_generate_candidates` | 73 | indirect: CLI test_b_mapping_modes (every mode vs truth, planted genome vs brute-force loci within -E) |

### `src/approximate_search/approximate_search_hybrid.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `as_hybrid_control_begin` | 40 | indirect: gcov 0%%, inlined on the main path): CLI test_b_mapping_modes (every mode vs truth, planted genome vs brute-force loci within -E |
|  | `as_hybrid_control_filtering_adaptive_next_state` | 49 | indirect: CLI test_b_mapping_modes (every mode vs truth, planted genome vs brute-force loci within -E) |
|  | `as_hybrid_control_neighborhood_next_state` | 86 | indirect: gcov 0%%, inlined on the main path): CLI test_b_mapping_modes (every mode vs truth, planted genome vs brute-force loci within -E |
|  | `approximate_search_hybrid` | 107 | indirect: CLI test_b_mapping_modes (every mode vs truth, planted genome vs brute-force loci within -E) |

### `src/approximate_search/approximate_search_neighborhood.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `approximate_search_neighborhood_search_brute_force` | 78 | indirect: CLI test_b_mapping_modes (every mode vs truth, planted genome vs brute-force loci within -E) |
|  | `approximate_search_neighborhood_search_partition` | 109 | indirect: CLI test_b_mapping_modes (every mode vs truth, planted genome vs brute-force loci within -E) |
|  | `approximate_search_neighborhood_search_partition_preconditioned` | 127 | indirect: CLI test_b_mapping_modes (every mode vs truth, planted genome vs brute-force loci within -E) |

### `src/approximate_search/approximate_search_stages.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `approximate_search_exact_filtering_adaptive` | 71 | indirect: CLI test_b_mapping_modes (every mode vs truth, planted genome vs brute-force loci within -E) |
|  | `approximate_search_neighborhood` | 182 | indirect: CLI test_b_mapping_modes (every mode vs truth, planted genome vs brute-force loci within -E) |
|  | `approximate_search_align_local` | 200 | indirect: CLI test_b_mapping_modes (every mode vs truth, planted genome vs brute-force loci within -E) |
|  | `approximate_search_end` | 218 | indirect: CLI test_b_mapping_modes (every mode vs truth, planted genome vs brute-force loci within -E) |

### `src/approximate_search/approximate_search_verify_candidates.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `approximate_search_verify_candidates` | 49 | indirect: CLI test_b_mapping_modes (every mode vs truth, planted genome vs brute-force loci within -E) |
|  | `approximate_search_verify_extend_candidate` | 68 | partial: paired-end only; executed by CLI test_check_alignments_correct_paired and the pilot PE tests (recovery extension runs), but no input was found where the extension changes the output (mates with 10 spaced substitutions and no exact 10-mer are not re... |

### `src/archive/search/archive_check.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `archive_check_se_match_retrieve_text` | 52 | partial: CLI test_check_alignments_correct (passing case only) |
|  | `archive_check_se_match_check_optimum` | 88 | not covered: reached only by -c first-best\|best\|all-best\|complete, known-bugs/B-check-alignments-crash.md |
|  | `archive_check_se_match_print` | 121 | not covered: -c failure report; the only way to reach it is the -c crash (B-check-alignments-crash.md) |
|  | `archive_check_se_match_print_incorrect` | 160 | not covered: needs an incorrect alignment under -c correct; none found |
|  | `archive_check_se_match_print_suboptimum` | 170 | not covered: known-bugs/B-check-alignments-crash.md |
|  | `archive_check_se_matches` | 207 | partial: CLI test_check_alignments_correct (-c correct, none) |
|  | `archive_check_pe_matches` | 281 | partial: CLI test_check_alignments_correct_paired |

### `src/archive/search/archive_search.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `archive_search_configure` | 50 | indirect: CLI test_b_* (every SE mapping) |
|  | `archive_search_se_new` | 60 | indirect: CLI test_b_* (every SE mapping) |
|  | `archive_search_pe_new` | 69 | partial: paired-end only; executed by CLI test_check_alignments_correct_paired and the pilot PE tests (recovery extension runs), but no input was found where the extension changes the output (mates with 10 spaced substitutions and no exact 10-mer are not re... |
|  | `archive_search_delete` | 92 | indirect: CLI test_b_* (every SE mapping) |
|  | `archive_search_prepare_sequence` | 100 | indirect: CLI test_b_* (every SE mapping) |
|  | `archive_search_inject_handlers` | 115 | indirect: CLI test_b_* (every SE mapping) |

### `src/archive/search/archive_search_handlers.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `archive_search_handlers_new` | 34 | indirect: CLI test_b_* (every SE mapping) |
|  | `archive_search_handlers_clear` | 50 | indirect: CLI test_b_* (every SE mapping) |
|  | `archive_search_handlers_delete` | 55 | indirect: CLI test_b_* (every SE mapping) |
|  | `archive_search_handlers_prepare_se` | 66 | indirect: CLI test_b_* (every SE mapping) |
|  | `archive_search_handlers_prepare_pe` | 79 | partial: paired-end only; executed by CLI test_check_alignments_correct_paired and the pilot PE tests (recovery extension runs), but no input was found where the extension changes the output (mates with 10 spaced substitutions and no exact 10-mer are not re... |

### `src/archive/search/archive_search_se.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `archive_search_se` | 57 | indirect: CLI test_b_* (every SE mapping) |

### `src/archive/search/archive_search_se_parameters.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `search_parameters_init_alignment_swg_penalties` | 38 | C test_b_swg (brute-force Gotoh, hand values, --help defaults) |
| x | `search_parameters_init_alignment` | 69 | C test_b_swg (brute-force Gotoh, hand values, --help defaults) |
| x | `select_parameters_init_mapq` | 97 | C test_b_swg (brute-force Gotoh, hand values, --help defaults) |
| x | `select_parameters_init_gpu` | 101 | C test_b_swg (brute-force Gotoh, hand values, --help defaults) |
| x | `search_parameters_init` | 108 | C test_b_swg (brute-force Gotoh, hand values, --help defaults) |
| x | `search_configure_alignment_match_scores` | 157 | C test_b_swg (brute-force Gotoh, hand values, --help defaults) |
|  | `search_configure_alignment_mismatch_scores` | 167 | partial: C test_b_swg checks the matrix; its generic_mismatch_score is known-bugs/B-mismatch-penalty-sign.md (-B, --gap-affine-penalties untested) |
| x | `search_instantiate_values` | 203 | C test_b_swg (brute-force Gotoh, hand values, --help defaults) |

### `src/archive/search/archive_search_se_stepwise.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `archive_search_se_stepwise_debug_prologue` | 48 | not covered: GEM_DEEP_DEBUG only |
|  | `archive_search_se_stepwise_finish_search` | 235 | partial: paired-end only; executed by CLI test_check_alignments_correct_paired and the pilot PE tests (recovery extension runs), but no input was found where the extension changes the output (mates with 10 spaced substitutions and no exact 10-mer are not re... |

### `src/filtering/candidates/filtering_candidates.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `filtering_candidates_init` | 54 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_init_alignment` | 66 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_clear` | 83 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_destroy` | 108 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_inject_handlers` | 121 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_allocate_position` | 136 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_free_position` | 143 | not covered: memory release, no observable output (gcov: not executed) |
|  | `filtering_candidates_allocate_region` | 148 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_free_region` | 155 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_allocate_alignment_tiles` | 170 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_free_alignment_tiles` | 177 | not covered: memory release, no observable output (gcov: not executed) |
|  | `filtering_candidates_get_num_positions` | 187 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_get_positions` | 196 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_clear_positions` | 200 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_get_num_regions` | 213 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_set_num_regions` | 217 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_get_regions` | 222 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_get_num_discarded_regions` | 239 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_add_num_discarded_regions` | 248 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_get_discarded_regions` | 253 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_reserve_discarded_regions` | 257 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_clear_discarded_regions` | 263 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
| x | `filtering_position_cmp_position` | 276 | C test_b_sort_queue (sort_instances_order_and_permutation) |
| x | `filtering_region_cmp_align_distance` | 279 | C test_b_sort_queue (sort_instances_order_and_permutation) |
| x | `filtering_region_cmp_scaffold_coverage` | 285 | C test_b_sort_queue (sort_instances_order_and_permutation) |
|  | `filtering_candidates_sort_positions` | 295 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_sort_regions_by_align_distance` | 303 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_sort_regions_by_scaffold_coverage` | 311 | indirect: (gcov 0%%, inlined on the main path): CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_sort_discarded_by_scaffold_coverage` | 315 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |

### `src/filtering/candidates/filtering_candidates_accessors.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `filtering_candidates_compose_filtering_region_from_positions_exact` | 59 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_compose_filtering_region_from_positions` | 110 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_add_positions_from_interval` | 171 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_add_region_from_group_positions` | 212 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_add_region_verified` | 258 | partial: paired-end only; executed by CLI test_check_alignments_correct_paired and the pilot PE tests (recovery extension runs), but no input was found where the extension changes the output (mates with 10 spaced substitutions and no exact 10-mer are not re... |

### `src/filtering/candidates/filtering_candidates_align.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `filtering_candidates_align_search_cache` | 45 | not covered: gcov shows it is not executed by the B tests (cache-hit path) |
|  | `filtering_candidates_align_region` | 62 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_align_local_region` | 99 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_align_extended_region` | 127 | partial: paired-end only; executed by CLI test_check_alignments_correct_paired and the pilot PE tests (recovery extension runs), but no input was found where the extension changes the output (mates with 10 spaced substitutions and no exact 10-mer are not re... |
|  | `filtering_candidates_align_candidates` | 156 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_align_extended_candidates` | 227 | partial: paired-end only; executed by CLI test_check_alignments_correct_paired and the pilot PE tests (recovery extension runs), but no input was found where the extension changes the output (mates with 10 spaced substitutions and no exact 10-mer are not re... |

### `src/filtering/candidates/filtering_candidates_align_local.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `filtering_candidates_align_local_discarded` | 46 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_align_local` | 73 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |

### `src/filtering/candidates/filtering_candidates_buffered.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `filtering_candidates_buffered_inject_handlers` | 34 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |

### `src/filtering/candidates/filtering_candidates_classify.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `filtering_candidates_classify_subdominant_match_edit` | 30 | not covered: gcov: not executed by SE tests |
|  | `filtering_candidates_classify_subdominant_match_swg` | 46 | not covered: gcov: not executed by SE tests |
|  | `filtering_candidates_classify_subdominant_match` | 65 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |

### `src/filtering/candidates/filtering_candidates_extend.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `filtering_candidates_compute_extension_region` | 41 | partial: paired-end only; executed by CLI test_check_alignments_correct_paired and the pilot PE tests (recovery extension runs), but no input was found where the extension changes the output (mates with 10 spaced substitutions and no exact 10-mer are not re... |
|  | `filtering_candidates_extend_match` | 69 | partial: paired-end only; executed by CLI test_check_alignments_correct_paired and the pilot PE tests (recovery extension runs), but no input was found where the extension changes the output (mates with 10 spaced substitutions and no exact 10-mer are not re... |

### `src/filtering/candidates/filtering_candidates_process.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `filtering_candidates_compose_filtering_regions` | 48 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_compute_text_boundaries` | 131 | indirect: (gcov 0%%, inlined on the main path): CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_compute_text_coordinates` | 170 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_decode_filtering_positions` | 215 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_decode_filtering_position_batch_prefetched` | 236 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_candidates_process_candidates` | 340 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |

### `src/filtering/candidates/filtering_candidates_verify.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `filtering_candidates_extend_candidates` | 45 | partial: paired-end only; executed by CLI test_check_alignments_correct_paired and the pilot PE tests (recovery extension runs), but no input was found where the extension changes the output (mates with 10 spaced substitutions and no exact 10-mer are not re... |
|  | `filtering_candidates_extend_discarded_candidates` | 71 | indirect: CLI test_b_alignment_options (local chimeras: discarded candidates extended for local alignment), test_b_pins |
|  | `filtering_candidates_verify_candidates` | 100 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |

### `src/filtering/region/filtering_region.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `filtering_region_retrieve_text` | 53 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_region_compute_key_trims` | 108 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |

### `src/filtering/region/filtering_region_align.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `filtering_region_align_debug_prologue` | 40 | not covered: GEM_DEEP_DEBUG only |
|  | `filtering_region_align_debug_epilogue` | 49 | not covered: GEM_DEEP_DEBUG only |
|  | `filtering_region_align_clone` | 87 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_region_align` | 130 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_region_align_local` | 205 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |

### `src/filtering/region/filtering_region_cache.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `filtering_region_cache_init` | 43 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_region_cache_clear` | 47 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_region_transient_cache_add` | 58 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_region_cache_cmp_regions` | 68 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_region_transient_cache_search` | 80 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |

### `src/filtering/region/filtering_region_verify.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `filtering_region_verify_hamming_text` | 46 | indirect: (gcov 0%%, inlined on the main path): CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_region_verify_hamming` | 67 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_region_verify_levenshtein` | 104 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_region_verify` | 146 | indirect: CLI test_b_mapping_modes, test_b_alignment_options (verification strategies, -e, q-gram), test_b_pins |
|  | `filtering_region_verify_extension` | 195 | partial: paired-end only; executed by CLI test_check_alignments_correct_paired and the pilot PE tests (recovery extension runs), but no input was found where the extension changes the output (mates with 10 spaced substitutions and no exact 10-mer are not re... |

### `src/matches/align/match_align.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `match_aling_add_clipping` | 42 | indirect: CLI test_b_alignment_options (AS vs brute-force Gotoh under -A/-O/-X, local chimeras, curation, bandwidth), test_b_mapping_modes (models), test_b_pins |
|  | `match_align_exact` | 99 | indirect: CLI test_b_alignment_options (AS vs brute-force Gotoh under -A/-O/-X, local chimeras, curation, bandwidth), test_b_mapping_modes (models), test_b_pins |
|  | `match_align_pseudoalignment` | 136 | indirect: CLI test_b_alignment_options (AS vs brute-force Gotoh under -A/-O/-X, local chimeras, curation, bandwidth), test_b_mapping_modes (models), test_b_pins |
|  | `match_align_hamming` | 181 | indirect: CLI test_b_alignment_options (AS vs brute-force Gotoh under -A/-O/-X, local chimeras, curation, bandwidth), test_b_mapping_modes (models), test_b_pins |
|  | `match_align_levenshtein` | 241 | indirect: CLI test_b_alignment_options (AS vs brute-force Gotoh under -A/-O/-X, local chimeras, curation, bandwidth), test_b_mapping_modes (models), test_b_pins |
|  | `match_align_smith_waterman_gotoh` | 302 | indirect: CLI test_b_alignment_options (AS vs brute-force Gotoh under -A/-O/-X, local chimeras, curation, bandwidth), test_b_mapping_modes (models), test_b_pins |
|  | `match_align_smith_waterman_gotoh_local` | 355 | indirect: CLI test_b_alignment_options (AS vs brute-force Gotoh under -A/-O/-X, local chimeras, curation, bandwidth), test_b_mapping_modes (models), test_b_pins |

### `src/matches/align/match_align_normalize.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `match_align_normalize_cigar_trim` | 31 | indirect: CLI test_b_alignment_options (AS vs brute-force Gotoh under -A/-O/-X, local chimeras, curation, bandwidth), test_b_mapping_modes (models), test_b_pins |
|  | `match_align_normalize_cigar` | 64 | indirect: CLI test_b_alignment_options (AS vs brute-force Gotoh under -A/-O/-X, local chimeras, curation, bandwidth), test_b_mapping_modes (models), test_b_pins |
|  | `match_align_normalize` | 172 | indirect: CLI test_b_alignment_options (AS vs brute-force Gotoh under -A/-O/-X, local chimeras, curation, bandwidth), test_b_mapping_modes (models), test_b_pins |

### `src/matches/align/match_align_swg.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `match_align_swg_compute_alignment_type` | 70 | indirect: CLI test_b_alignment_options (AS vs brute-force Gotoh under -A/-O/-X, local chimeras, curation, bandwidth), test_b_mapping_modes (models), test_b_pins |
|  | `match_align_swg` | 111 | indirect: CLI test_b_alignment_options (AS vs brute-force Gotoh under -A/-O/-X, local chimeras, curation, bandwidth), test_b_mapping_modes (models), test_b_pins |

### `src/matches/align/match_align_swg_chained.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `match_align_swg_chained_add_region` | 36 | indirect: CLI test_b_alignment_options (AS vs brute-force Gotoh under -A/-O/-X, local chimeras, curation, bandwidth), test_b_mapping_modes (models), test_b_pins |
|  | `match_align_swg_chained_gap_full` | 95 | indirect: CLI test_b_alignment_options (AS vs brute-force Gotoh under -A/-O/-X, local chimeras, curation, bandwidth), test_b_mapping_modes (models), test_b_pins |
|  | `match_align_swg_chained_gap_extending` | 142 | indirect: CLI test_b_alignment_options (AS vs brute-force Gotoh under -A/-O/-X, local chimeras, curation, bandwidth), test_b_mapping_modes (models), test_b_pins |
|  | `match_align_swg_end_extending` | 238 | indirect: CLI test_b_alignment_options (AS vs brute-force Gotoh under -A/-O/-X, local chimeras, curation, bandwidth), test_b_mapping_modes (models), test_b_pins |
|  | `match_align_swg_chained_to_first` | 289 | indirect: CLI test_b_alignment_options (AS vs brute-force Gotoh under -A/-O/-X, local chimeras, curation, bandwidth), test_b_mapping_modes (models), test_b_pins |
|  | `match_align_swg_chained_from_last` | 326 | indirect: CLI test_b_alignment_options (AS vs brute-force Gotoh under -A/-O/-X, local chimeras, curation, bandwidth), test_b_mapping_modes (models), test_b_pins |
|  | `match_align_swg_chained_gap` | 363 | indirect: CLI test_b_alignment_options (AS vs brute-force Gotoh under -A/-O/-X, local chimeras, curation, bandwidth), test_b_mapping_modes (models), test_b_pins |
|  | `match_align_swg_chained` | 411 | indirect: CLI test_b_alignment_options (AS vs brute-force Gotoh under -A/-O/-X, local chimeras, curation, bandwidth), test_b_mapping_modes (models), test_b_pins |

### `src/matches/align/match_align_swg_local.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `match_align_swg_local_alignment_add_match_compose_cigar` | 43 | indirect: CLI test_b_alignment_options (AS vs brute-force Gotoh under -A/-O/-X, local chimeras, curation, bandwidth), test_b_mapping_modes (models), test_b_pins |
|  | `match_align_swg_local_alignment_add_match` | 86 | indirect: CLI test_b_alignment_options (AS vs brute-force Gotoh under -A/-O/-X, local chimeras, curation, bandwidth), test_b_mapping_modes (models), test_b_pins |
|  | `match_align_swg_local_alignment` | 132 | indirect: CLI test_b_alignment_options (AS vs brute-force Gotoh under -A/-O/-X, local chimeras of 40, 50 and 60 genomic bases on the --alignment-local-min-identity/-min-score thresholds, curation, bandwidth), test_b_mapping_modes, test_b_pins |

### `src/matches/align/match_alignment.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `match_alignment_print_pretty` | 52 | not covered: debug/-c failure printing only |

### `src/matches/align/match_alignment_region.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `match_alignment_region_init` | 31 | C test_b_sort_queue (match_alignment_region_accessors) |
| x | `match_alignment_region_get_type` | 50 | C test_b_sort_queue (match_alignment_region_accessors) |
| x | `match_alignment_region_set_type` | 54 | C test_b_sort_queue (match_alignment_region_accessors) |
| x | `match_alignment_region_get_error` | 63 | C test_b_sort_queue (match_alignment_region_accessors) |
| x | `match_alignment_region_get_cigar_buffer_offset` | 72 | C test_b_sort_queue (match_alignment_region_accessors) |
| x | `match_alignment_region_get_cigar_length` | 76 | C test_b_sort_queue (match_alignment_region_accessors) |
| x | `match_alignment_region_get_key_begin` | 80 | C test_b_sort_queue (match_alignment_region_accessors) |
| x | `match_alignment_region_get_key_end` | 84 | C test_b_sort_queue (match_alignment_region_accessors) |
| x | `match_alignment_region_get_text_begin` | 88 | C test_b_sort_queue (match_alignment_region_accessors) |
| x | `match_alignment_region_get_text_end` | 92 | C test_b_sort_queue (match_alignment_region_accessors) |
| x | `match_alignment_region_set_cigar_buffer_offset` | 96 | C test_b_sort_queue (match_alignment_region_accessors) |
| x | `match_alignment_region_set_cigar_length` | 101 | C test_b_sort_queue (match_alignment_region_accessors) |
| x | `match_alignment_region_set_key_begin` | 106 | C test_b_sort_queue (match_alignment_region_accessors) |
| x | `match_alignment_region_set_key_end` | 111 | C test_b_sort_queue (match_alignment_region_accessors) |
| x | `match_alignment_region_set_text_begin` | 116 | C test_b_sort_queue (match_alignment_region_accessors) |
| x | `match_alignment_region_set_text_end` | 121 | C test_b_sort_queue (match_alignment_region_accessors) |
| x | `match_alignment_region_text_coverage` | 129 | C test_b_sort_queue (match_alignment_region_accessors) |
| x | `match_alignment_region_text_distance` | 133 | C test_b_sort_queue (match_alignment_region_accessors) |
| x | `match_alignment_region_text_overlap` | 138 | C test_b_sort_queue (match_alignment_region_accessors) |
| x | `match_alignment_region_cmp_key_position` | 147 | C test_b_sort_queue (match_alignment_region_accessors) |
| x | `match_alignment_region_cmp_text_position` | 152 | C test_b_sort_queue (match_alignment_region_accessors) |
|  | `match_alignment_region_print` | 196 | not covered: debug printing only |

### `src/matches/align/match_alignment_region_rl.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `match_alignment_region_rl_translate_cigar` | 31 | not covered: run-length (RL) indexes only; the suite's indexes are not RL |
|  | `match_alignment_region_rl_translate` | 77 | not covered: RL indexes only |

### `src/matches/scaffold/match_scaffold.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `match_scaffold_init` | 45 | indirect: CLI test_b_alignment_options (local chimeras, indel reads), test_b_pins (scaffolding-3/-30) |
|  | `match_scaffold_destroy` | 57 | indirect: CLI test_b_alignment_options (local chimeras, indel reads), test_b_pins (scaffolding-3/-30) |
|  | `match_scaffold_allocate_alignment_region` | 71 | indirect: CLI test_b_alignment_options (local chimeras, indel reads), test_b_pins (scaffolding-3/-30) |
|  | `match_scaffold_free_alignment_region` | 77 | indirect: CLI test_b_alignment_options (local chimeras, indel reads), test_b_pins (scaffolding-3/-30) |
|  | `match_scaffold_rl_translate_regions` | 86 | not covered: RL indexes only |
|  | `match_scaffold_region_chain` | 125 | indirect: CLI test_b_alignment_options (local chimeras, indel reads), test_b_pins (scaffolding-3/-30) |
|  | `match_scaffold_adaptive` | 151 | indirect: CLI test_b_alignment_options (local chimeras, indel reads), test_b_pins (scaffolding-3/-30) |
|  | `match_scaffold_sort_alignment_regions` | 224 | indirect: CLI test_b_alignment_options (local chimeras, indel reads), test_b_pins (scaffolding-3/-30) |
|  | `match_scaffold_print` | 231 | not covered: debug printing only |

### `src/matches/scaffold/match_scaffold_chain.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `match_scaffold_lis_cmp` | 39 | not covered: gcov: not executed (qsort comparator of the chaining LIS) |
|  | `match_scaffold_chain_is_compatible` | 49 | indirect: CLI test_b_alignment_options (local chimeras, indel reads), test_b_pins (scaffolding-3/-30) |
|  | `match_scaffold_chain_compute_lis` | 59 | indirect: CLI test_b_alignment_options (local chimeras, indel reads), test_b_pins (scaffolding-3/-30) |
|  | `match_scaffold_chain_backtrace_lis` | 109 | indirect: CLI test_b_alignment_options (local chimeras, indel reads), test_b_pins (scaffolding-3/-30) |
|  | `match_scaffold_chain_alignment_regions` | 126 | indirect: CLI test_b_alignment_options (local chimeras, indel reads), test_b_pins (scaffolding-3/-30) |
|  | `match_scaffold_chain_exact_extend_right` | 149 | indirect: CLI test_b_alignment_options (local chimeras, indel reads), test_b_pins (scaffolding-3/-30) |
|  | `match_scaffold_chain_exact_extend_left` | 182 | indirect: CLI test_b_alignment_options (local chimeras, indel reads), test_b_pins (scaffolding-3/-30) |
|  | `match_scaffold_chain_exact_extend` | 218 | indirect: CLI test_b_alignment_options (local chimeras, indel reads), test_b_pins (scaffolding-3/-30) |
|  | `match_scaffold_chain` | 256 | indirect: CLI test_b_alignment_options (local chimeras, indel reads), test_b_pins (scaffolding-3/-30) |

### `src/matches/scaffold/match_scaffold_compose.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `match_scaffold_compose_add_exact_match` | 47 | indirect: CLI test_b_alignment_options (local chimeras, indel reads), test_b_pins (scaffolding-3/-30) |
|  | `match_scaffold_compose_add_mismatch` | 69 | indirect: CLI test_b_alignment_options (local chimeras, indel reads), test_b_pins (scaffolding-3/-30) |

### `src/matches/scaffold/match_scaffold_levenshtein.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `match_scaffold_levenshtein_allocate` | 44 | indirect: CLI test_b_alignment_options (local chimeras, indel reads), test_b_pins (scaffolding-3/-30) |
|  | `match_scaffold_levenshtein_compose_alignment` | 63 | indirect: CLI test_b_alignment_options (local chimeras, indel reads), test_b_pins (scaffolding-3/-30) |
|  | `match_scaffold_levenshtein_tile_align_bpm` | 127 | indirect: CLI test_b_alignment_options (local chimeras, indel reads), test_b_pins (scaffolding-3/-30) |
|  | `match_scaffold_levenshtein_tile` | 178 | indirect: CLI test_b_alignment_options (local chimeras, indel reads), test_b_pins (scaffolding-3/-30) |
|  | `match_scaffold_levenshtein` | 218 | indirect: CLI test_b_alignment_options (local chimeras, indel reads), test_b_pins (scaffolding-3/-30) |

### `src/neighborhood_search/nsearch_filtering.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `nsearch_filtering` | 34 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |

### `src/neighborhood_search/nsearch_hamming.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `nsearch_hamming_query` | 35 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_hamming_directional_query` | 52 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_hamming_terminate` | 74 | not covered: gcov: not executed |
|  | `nsearch_hamming_brute_force_step` | 98 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_hamming_brute_force` | 124 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_hamming_scheduled_search_operation_next` | 154 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_hamming_scheduled_search` | 217 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_hamming` | 238 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_hamming_preconditioned` | 262 | not covered: hybrid complete + Hamming reaches it only when filtering leaves strata unsearched; no such input found (gcov: not executed) |

### `src/neighborhood_search/nsearch_levenshtein.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `nsearch_levenshtein_base_step` | 34 | not covered: known bug: reached only by --candidate-generation fixed or factors with region_error > 0 (A-candidate-generation-strategies.md) |
|  | `nsearch_levenshtein_base` | 75 | not covered: known bug: reached only by --candidate-generation fixed or factors with region_error > 0 (A-candidate-generation-strategies.md) |
|  | `nsearch_levenshtein_brute_force_step` | 100 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_levenshtein_brute_force` | 156 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_levenshtein` | 194 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_levenshtein_preconditioned` | 218 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |

### `src/neighborhood_search/nsearch_levenshtein_control.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `nsearch_levenshtein_candidates_cutoff` | 35 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_levenshtein_control_limit_interval` | 66 | not covered: gcov: not executed |
|  | `nsearch_levenshtein_terminate` | 81 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_levenshtein_scheduled_terminate` | 112 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |

### `src/neighborhood_search/nsearch_levenshtein_query.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `nsearch_query_init` | 29 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_levenshtein_query` | 47 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_levenshtein_scheduled_precompute_query` | 68 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_levenshtein_scheduled_query_exact` | 88 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_levenshtein_scheduled_query` | 118 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |

### `src/neighborhood_search/nsearch_levenshtein_scheduled.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `nsearch_levenshtein_scheduled_align_distance_pass` | 38 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_levenshtein_scheduled_operation_step_next` | 69 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_levenshtein_scheduled_operation_step_query` | 103 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_levenshtein_scheduled_operation_step_query_exact` | 167 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_levenshtein_scheduled_search_next` | 208 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_levenshtein_scheduled_search` | 228 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |

### `src/neighborhood_search/nsearch_levenshtein_state.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `nsearch_levenshtein_state_init` | 29 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_levenshtein_state_prepare` | 52 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_levenshtein_state_prepare_column` | 68 | not covered: gcov: not executed |
|  | `nsearch_levenshtein_state_get_global_align_distance` | 82 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_levenshtein_state_get_local_align_distance` | 95 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_levenshtein_state_compute_chararacter` | 137 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_levenshtein_state_compute_text_banded` | 182 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_levenshtein_state_compute_chararacter_banded` | 206 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |

### `src/neighborhood_search/nsearch_operation.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `nsearch_operation_init` | 30 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_operation_chained_prepare_sequence_forward` | 47 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_operation_chained_prepare_sequence_reverse` | 57 | not covered: gcov: not executed |
|  | `nsearch_operation_chained_prepare` | 69 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |

### `src/neighborhood_search/nsearch_parameters.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `nsearch_parameters_init` | 26 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |

### `src/neighborhood_search/nsearch_partition.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `nsearch_partition_compute` | 29 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_partition_preconditioned_compute` | 42 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_partition_compute_error` | 66 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |

### `src/neighborhood_search/nsearch_schedule.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `nsearch_schedule_init` | 35 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_schedule_inject_mm` | 76 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_schedule_add_pending_search` | 84 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_schedule_search_step` | 121 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_schedule_search` | 190 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_schedule_search_preconditioned_step` | 211 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_schedule_search_preconditioned` | 278 | indirect: CLI test_b_mapping_modes::test_neighborhood_search_reports_exactly_loci_within_E (exact set vs brute-force Myers/Hamming scan) |
|  | `nsearch_schedule_compute_min_error` | 305 | not covered: gcov: not executed |

### `src/utils/priority_queue.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `pqueue_new` | 37 | C test_b_sort_queue (pqueue_orders_by_priority) |
| x | `pqueue_delete` | 51 | C test_b_sort_queue (pqueue_orders_by_priority) |
| x | `pqueue_top_priority` | 60 | C test_b_sort_queue (pqueue_orders_by_priority) |
| x | `pqueue_push_` | 72 | C test_b_sort_queue (pqueue_orders_by_priority) |
| x | `pqueue_pop_` | 87 | C test_b_sort_queue (pqueue_orders_by_priority) |

## Area C: Paired-end, bisulfite, matches and MAPQ (149 functions in 20 files)

### `src/archive/score/archive_score_logit.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `archive_score_logit` | 30 | test_c_matches:score_probability_and_logit, random coefficients vs a recomputed logistic |
| x | `archive_score_logit_unique` | 106 | test_c_matches:score_probability_and_logit, score_se_models |
| x | `archive_score_logit_mmap` | 112 | test_c_matches:score_probability_and_logit, score_se_models |

### `src/archive/score/archive_score_pe.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `archive_score_matches_pe_default` | 80 | test_c_paired:score_pe_models, all branches; test_c_paired.py orientation/discordant tests, pins |
| x | `archive_score_matches_pe_classify` | 130 | test_c_paired:score_pe_models; test_c_paired.py::test_pe_mapq_model_none_and_classify |
| x | `archive_score_matches_pe` | 182 | test_c_paired:score_pe_models model switch incl. none and no alignment model; test_c_paired.py MAPQ tests, test_pe_dump_predictors |

### `src/archive/score/archive_score_se.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `archive_score_probability_to_mapq` | 49 | test_c_matches:score_probability_and_logit, hand values |
| x | `archive_score_matches_se_get_model` | 55 | test_c_matches:score_se_models |
| x | `archive_score_matches_se_logit_unique` | 63 | test_c_matches:score_se_models, cap at 59 |
| x | `archive_score_matches_se_logit_mmap` | 76 | test_c_matches:score_se_models via the default model's mmap branch |
| x | `archive_score_matches_se_default` | 92 | test_c_matches:score_se_models; test_c_select.py::test_mapq_ties_and_unique_reads |
| x | `archive_score_matches_se_classify` | 129 | test_c_matches:score_se_models; test_c_select.py::test_mapq_model_classify_classes |
| x | `archive_score_matches_se` | 190 | test_c_matches:score_se_models; test_c_select.py MAPQ model tests incl. dump-predictors |

### `src/archive/search/archive_search_pe.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `archive_search_pe_continue` | 50 | every test in test_c_paired.py; state machine reached end to end |
| x | `archive_search_pe` | 83 | every test in test_c_paired.py, test_c_bisulfite.py::test_pe_directional_pairs |

### `src/archive/search/archive_search_pe_extend.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `archive_search_pe_extend_matches_ends` | 42 | test_c_paired.py::test_mate_rescue_in_high_copy_repeat; the no-recovery mutant fails it and 3 pins |
| x | `archive_search_pe_extend_matches` | 111 | same |

### `src/archive/search/archive_search_pe_parameters.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `search_paired_parameters_init` | 31 | test_c_paired:search_paired_parameters_defaults, checked against --help defaults |

### `src/archive/search/archive_search_pe_stages.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `archive_search_pe_begin` | 54 | every paired CLI test |
| x | `archive_search_pe_search_end1` | 67 | every paired CLI test |
| x | `archive_search_pe_search_end2` | 100 | every paired CLI test |
| x | `archive_search_pe_use_shortcut_extension` | 113 | test_c_paired:pe_extension_switches; always false in this version |
| x | `archive_search_pe_use_recovery_extension` | 129 | test_c_paired:pe_extension_switches, key length 0/300/301 and the switch; test_mate_rescue_in_high_copy_repeat |
| x | `archive_search_pe_recovery` | 153 | test_c_paired.py::test_mate_rescue_in_high_copy_repeat |
| x | `archive_search_pe_cross_pair_ends` | 185 | every paired CLI test; test_unmapped_and_cross_contig_mates for an end without matches |
| x | `archive_search_pe_find_pairs` | 200 | every paired CLI test; discordant modes in test_orientation_and_discordant_search, test_discordant_search_only_without_concordant |
| x | `archive_search_pe_end` | 230 | test_c_paired.py::test_pe_dump_predictors checks the read length it sets through the edit and swg predictors |

### `src/archive/search/archive_select.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `archive_select_configure_se` | 37 | test_c_matches:select_parameters |
| x | `archive_select_configure_pe` | 43 | test_c_matches:select_parameters |
| x | `archive_select_se_matches_discard` | 52 | test_c_matches:select_se_matches; test_c_select.py -m tests |
| x | `archive_select_se_matches` | 72 | test_c_matches:select_se_matches; test_c_select.py -M/-m tests, pins |
| x | `archive_select_pe_matches` | 100 | test_c_paired:select_pe_matches incl. template-length sampling; test_c_paired.py::test_pair_tie_between_two_copies with -M 1, pins |

### `src/archive/search/archive_select_parameters.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `select_parameters_init` | 32 | test_c_matches:select_parameters; test_c_select.py::test_max_reported_matches_default_is_five |
| x | `select_parameters_instantiate_values` | 47 | test_c_matches:select_parameters, fraction and count; test_c_select.py::test_min_reported_strata_fraction_of_read_length; pin m0.03 |

### `src/matches/classify/matches_classify.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `matches_classify` | 29 | test_c_matches:classify_hand |
| x | `matches_classify_pattern_viable` | 72 | test_c_matches:classify_fallbacks |
| x | `matches_classify_max_matches_searched` | 79 | test_c_matches:classify_fallbacks |
| x | `matches_classify_min_depth_searched` | 113 | test_c_matches:classify_min_depth_searched, truth table |
| x | `matches_classify_neighbourhood_fallback` | 156 | test_c_matches:classify_fallbacks |
| x | `matches_classify_local_alignment_fallback` | 169 | test_c_matches:classify_fallbacks |
| x | `matches_classify_adjust_max_error_by_strata_after_best` | 178 | test_c_matches:classify_fallbacks |
| x | `matches_classify_compute_max_search_error` | 196 | test_c_matches:classify_max_search_error, every class |

### `src/matches/classify/matches_metrics.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `matches_metrics_init` | 29 | test_c_matches:metrics_update |
| x | `matches_metrics_set_proper_length` | 43 | test_c_matches:metrics_update, predictors |
| x | `matches_metrics_set_read_length` | 48 | test_c_matches:metrics_update; test_pe_dump_predictors |
| x | `matches_metrics_set_swg_match_score` | 53 | test_c_matches:metrics_update, predictors |
| x | `matches_metrics_set_region_profile_metrics` | 58 | test_c_matches:metrics_update |
| x | `matches_metrics_add_accepted_candidates` | 67 | test_c_matches:metrics_update |
| x | `matches_metrics_update` | 80 | test_c_matches:metrics_update, random multisets vs the two smallest/largest |

### `src/matches/classify/matches_predictors.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `matches_predictors_compute_unique` | 38 | test_c_matches:predictors |
| x | `matches_predictors_compute_mmaps` | 62 | test_c_matches:predictors |
| x | `matches_predictors_compute_se` | 89 | test_c_matches:predictors; test_c_select.py::test_mapq_model_dump_predictors |
| x | `matches_predictors_compute_pe` | 129 | test_c_paired:predictors_pe; test_c_paired.py::test_pe_dump_predictors |
| x | `matches_predictors_print` | 172 | test_c_select.py::test_mapq_model_dump_predictors, test_c_paired.py::test_pe_dump_predictors; fields recomputed from the SAM records |
| x | `matches_predictors_se_print` | 238 | test_c_select.py::test_mapq_model_dump_predictors, one line per mapped read |
| x | `matches_predictors_pe_print` | 248 | test_c_paired.py::test_pe_dump_predictors, one line per paired read |

### `src/matches/classify/paired_matches_classify.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `paired_matches_classify` | 30 | test_c_paired:paired_classify |
| x | `paired_matches_classify_search_accomplished` | 68 | test_c_paired:paired_classify |
| x | `paired_matches_classify_subdominant_end` | 82 | test_c_paired:paired_classify |

### `src/matches/match_trace.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `match_trace_get_cigar_buffer` | 29 | test_c_matches:matches_recompute_and_accessors |
| x | `match_trace_get_cigar_length` | 34 | test_c_matches:matches_recompute_and_accessors |
| x | `match_trace_cmp_swg_score` | 43 | test_c_matches:matches_insert_sorted, vs the documented order |
| x | `match_trace_cmp_genomic_position` | 61 | test_c_matches:matches_sort_by_genomic_position |

### `src/matches/matches.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `matches_new` | 38 | test_c_matches, all matches cases |
| x | `matches_clear` | 63 | test_c_matches:matches_index_lookup, matches_strata_bookkeeping |
| x | `matches_delete` | 78 | test_c_matches |
| x | `matches_is_mapped` | 94 | test_c_matches:matches_insert_sorted |
| x | `matches_recompute_metrics` | 97 | test_c_matches:matches_recompute_and_accessors |
| x | `matches_update_mcs` | 125 | test_c_matches:matches_strata_bookkeeping |
| x | `matches_update_limited_exact_matches` | 130 | test_c_matches:matches_strata_bookkeeping |
| x | `matches_index_match` | 142 | test_c_matches:matches_index_lookup |
| x | `matches_lookup_match` | 152 | test_c_matches:matches_index_lookup, by begin and by end |
| x | `matches_get_num_match_traces` | 166 | test_c_matches |
| x | `matches_get_num_match_traces_extended` | 169 | test_c_matches:matches_local_pending |
| x | `matches_get_primary_match` | 172 | test_c_matches:matches_insert_sorted |
| x | `matches_get_match_traces` | 175 | test_c_matches |
| x | `matches_traces_sort_by_genomic_position` | 181 | test_c_matches:matches_sort_by_genomic_position; test_c_paired:find_pairs_all_combinations |
| x | `matches_add_match_trace_insert_sorted` | 188 | test_c_matches:matches_insert_sorted, stable for equal keys |
| x | `matches_add_match_trace_preserve_sorted` | 210 | test_c_matches:matches_preserve_sorted |
|  | `match_trace_locate` | 232 | indirect: test_c_select.py::test_max_reported_all_finds_every_designed_copy: every designed copy's position and strand, both read strands; needs an index's locator, so no C case |
| x | `match_trace_copy` | 250 | test_c_matches:matches_recompute_and_accessors |
|  | `matches_add_match_trace` | 273 | indirect: test_c_select.py -M all: each designed placement exactly once, no duplicates; every CLI test; no C case, needs a locator |
|  | `matches_add_match_trace_extended` | 327 | indirect: test_c_paired.py::test_mate_rescue_in_high_copy_repeat |
| x | `matches_add_match_trace_local_pending` | 343 | test_c_matches:matches_local_pending; local alignments in test_c_select.py mixed_se reads |
|  | `matches_local_pending_add_to_regular_matches` | 350 | indirect: test_c_select.py: reads with 14 substitutions map as local alignments; classify MAPQ 1, dump-predictors, pins |
| x | `matches_local_pending_add_to_extended_matches` | 363 | test_c_matches:matches_local_pending |

### `src/matches/matches_cigar.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `matches_cigar_buffer_add_cigar_element__attr` | 35 | test_c_matches:cigar_buffer_sentinel |
| x | `matches_cigar_buffer_add_cigar_element` | 50 | test_c_matches:cigar_buffer_sentinel |
| x | `matches_cigar_buffer_add_mismatch` | 57 | test_c_matches:cigar_buffer_sentinel |
| x | `matches_cigar_vector_append_insertion` | 68 | test_c_matches:cigar_append_merges_runs, cigar_append_starts_new_cigar |
| x | `matches_cigar_vector_append_deletion` | 90 | same |
| x | `matches_cigar_vector_append_match` | 112 | same |
| x | `matches_cigar_vector_append_mismatch` | 135 | same |
| x | `matches_cigar_vector_append_cigar_element` | 148 | test_c_matches:cigar_append_merges_runs |
| x | `matches_cigar_vector_insert_cigar_element` | 170 | test_c_matches:cigar_insert_element |
| x | `matches_cigar_compute_event_distance` | 224 | test_c_matches:cigar_distances, column expansion reference; CLI dump-predictors event field |
| x | `matches_cigar_compute_edit_distance` | 247 | test_c_matches:cigar_distances; CLI NM checks |
| x | `matches_cigar_compute_edit_distance__excluding_clipping` | 271 | test_c_matches:cigar_distances |
| x | `matches_cigar_compute_matching_bases` | 297 | test_c_matches:cigar_distances |
| x | `matches_cigar_element_effective_length` | 320 | test_c_matches:cigar_distances |
| x | `matches_cigar_compute_effective_length` | 340 | test_c_matches:cigar_distances; CLI TLEN checks |
| x | `matches_cigar_compute_error_quality` | 354 | test_c_matches:cigar_error_quality, hand-worked and random vs per-column reference |
| x | `match_cigar_print` | 434 | test_c_matches:cigar_print, hand-worked |

### `src/matches/matches_counters.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `matches_counters_new` | 34 | test_c_matches:counters_random_ops |
| x | `matches_counters_clear` | 43 | same |
| x | `matches_counters_delete` | 47 | same |
| x | `matches_counters_get_num_counters` | 54 | same |
| x | `matches_counters_get_counts` | 58 | same |
| x | `matches_counters_add` | 71 | same |
| x | `matches_counters_sub` | 85 | same |
| x | `matches_counters_compact` | 96 | same |
| x | `matches_counters_compute_matches_to_report` | 103 | test_c_matches:counters_matches_to_report, ten hand-worked cases; test_c_select.py -m tests |
| x | `matches_counters_count_first_subdominant` | 136 | test_c_matches:counters_random_ops, test_c_matches:predictors |

### `src/matches/paired_map.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `paired_map_compute_event_distance` | 29 | test_c_paired:paired_map_scores |
| x | `paired_map_compute_edit_distance` | 34 | test_c_paired:paired_map_scores, find_pairs_all_combinations |
| x | `paired_map_compute_swg_score` | 39 | same |
| x | `paired_map_compute_error_quality` | 44 | test_c_paired:paired_map_scores |

### `src/matches/paired_matches.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `paired_matches_new` | 49 | test_c_paired |
| x | `paired_matches_clear` | 68 | test_c_paired:add_and_cross_pair, with and without clearing the ends |
| x | `paired_matches_delete` | 84 | test_c_paired |
| x | `paired_matches_is_mapped` | 102 | test_c_paired:find_pairs_all_combinations, paired_map_order |
| x | `paired_matches_get_max_complete_stratum` | 115 | test_c_paired:paired_classify, limited and extended ends |
| x | `paired_matches_get_primary_map` | 125 | test_c_paired:paired_map_order |
| x | `paired_matches_get_num_maps` | 131 | test_c_paired |
| x | `paired_matches_get_maps` | 134 | test_c_paired |
| x | `paired_matches_get_num_discordant_maps` | 137 | test_c_paired:find_pairs_all_combinations, add_and_cross_pair |
| x | `paired_matches_get_discordant_maps` | 140 | same |
| x | `paired_matches_limit_maps` | 143 | test_c_paired:paired_map_order, select_pe_matches |
| x | `paired_matches_recompute_metrics` | 149 | test_c_paired:paired_map_order |
| x | `paired_matches_cmp_swg_score` | 167 | test_c_paired:paired_map_order (hand-worked order); the order of tied pairs (layout tie-break returns 0, known-bugs/r2-pe-tied-pair-order.md) is pinned by test_c_paired.py::test_tandem_array_pairs_are_pinned and test_d_output.py tie pins |
| x | `paired_matches_add_map_insert_sorted` | 190 | test_c_paired:paired_map_order |
| x | `paired_matches_add_map` | 212 | test_c_paired:paired_map_order, counters and metrics |
| x | `paired_matches_add` | 224 | test_c_paired:add_and_cross_pair |
| x | `paired_matches_compute_orientation` | 274 | test_c_paired:orientation_layout_template_length, 20,000 random pairs |
|  | `paired_matches_compute_layout` | 299 | partial: test_c_paired:orientation_layout_template_length; read 2 containing read 1 is a known bug, C-pair-layout-end2-contains-end1.md |
|  | `paired_matches_compute_template_length` | 326 | partial: test_c_paired:orientation_layout_template_length vs the SAM span; RF overlap is a known bug, C-tlen-sam-spec.md; CLI TLEN checks |
| x | `paired_matches_compute_relation` | 382 | test_c_paired:relation_tables_and_limits, random tables, -l/-L, distribution; CLI -l/-L and long-template tests |
| x | `paired_matches_cross_pair` | 423 | test_c_paired:add_and_cross_pair, find_pairs_all_combinations |
| x | `paired_matches_find_pairs_locate` | 461 | test_c_paired:find_pairs_all_combinations |
| x | `paired_matches_find_pairs_sorted` | 471 | test_c_paired:find_pairs_all_combinations, every same-sequence combination |
| x | `paired_matches_find_pairs` | 505 | same |
| x | `paired_matches_find_discordant_pairs` | 536 | test_c_paired:find_pairs_all_combinations, all three modes; test_c_paired.py discordant tests |

### `src/text/sequence_bisulfite.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `sequence_bisulfite_convert` | 32 | test_c_bisulfite:convert_tables |
| x | `sequence_bisulfite_check_cg_depletion_se` | 46 | test_c_bisulfite:cg_depletion; test_c_bisulfite.py non-stranded tests |
|  | `sequence_bisulfite_check_cg_depletion_pe` | 57 | partial: test_c_bisulfite:cg_depletion); not reachable from the CLI (known bug pe-bisulfite-non-stranded-ignored.md |
|  | `sequence_bisulfite_process_se` | 75 | partial: test_c_bisulfite:process_se_modes for inferred, 1, 2, non-stranded; CLI for inferred, 1, non-stranded; '2' and 'interleaved' from the CLI are a known bug, bisulfite-read-2-ignored.md |
| x | `sequence_bisulfite_process_pe` | 101 | test_c_bisulfite:process_pe_modes incl. the non-stranded flip; CLI directional pairs |
| x | `sequence_bisulfite_restore_se` | 123 | test_c_bisulfite:process_se_modes; CLI SEQ as sequenced |
| x | `sequence_bisulfite_restore_pe` | 128 | test_c_bisulfite:process_pe_modes; test_c_bisulfite.py::test_pe_directional_pairs |

## Area D: Output, CLI, threads and runtime (276 functions in 27 files)

### `src/io/buffered_output_file.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `buffered_output_file_new` | 50 | Cfmt buffered_output_file_printers_reach_the_stream, Cout all cases; every CLI run |
| x | `buffered_output_file_close` | 57 | Cfmt, Cout |
| x | `buffered_output_file_request_buffer` | 64 | Cout: block ids from 1-16 threads |
| x | `buffered_output_file_dump_buffer` | 69 | Cout |
| x | `buffered_output_file_safety_dump_buffer` | 78 | Cout 512-2048 B buffers force extensions; out test_output_model_does_not_change_output with 512 B / 1 KB buffers |
| x | `buffered_output_file_reserve` | 82 | Cout, Cfmt, out output-model tests |
| x | `bofprintf_uint64` | 136 | Cfmt against snprintf; Cout; all SAM/MAP checks |
| x | `bofprintf_int64` | 139 | Cfmt against snprintf; Cout; all SAM/MAP checks |
| x | `bofprintf_char` | 142 | Cfmt against snprintf; Cout; all SAM/MAP checks |
| x | `bofprintf_string` | 145 | Cfmt against snprintf; Cout; all SAM/MAP checks |

### `src/io/output_buffer.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `output_buffer_new` | 31 | Cfmt output_buffer_printers_match_printf |
| x | `output_buffer_clear` | 38 | Cfmt output_buffer_printers_match_printf |
| x | `output_buffer_delete` | 47 | Cfmt output_buffer_printers_match_printf |
| x | `output_buffer_set_state` | 54 | Cfmt |
| x | `output_buffer_get_state` | 57 | Cfmt |
|  | `output_buffer_set_incomplete` | 60 | indirect: Cout: extension buffers must be written in minor-id order; mutants m19/m20 |
| x | `output_buffer_get_used` | 63 | Cfmt |
| x | `output_buffer_get_buffer` | 66 | Cfmt |
| x | `bprintf_uint64` | 72 | Cfmt against snprintf, incl. negatives and INT64_MIN+1 |
| x | `bprintf_int64` | 77 | Cfmt against snprintf, incl. negatives and INT64_MIN+1 |
| x | `bprintf_char` | 89 | Cfmt against snprintf, incl. negatives and INT64_MIN+1 |
| x | `bprintf_buffer` | 94 | Cfmt against snprintf, incl. negatives and INT64_MIN+1 |

### `src/io/output_file.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `output_file_init_buffers` | 34 | Cout, all buffer counts 1-16 |
| x | `output_stream_new` | 67 | every CLI run to stdout or -o; out test_output_file_matches_stdout |
| x | `output_gzip_stream_new` | 78 | Cout blocks_in_order_gzip; out test_compressed_output[gzip] |
| x | `output_bzip_stream_new` | 89 | Cout blocks_in_order_bzip2; out test_compressed_output[bzip] |
| x | `output_file_close` | 100 | Cout, Cfmt |
| x | `output_file_serve_buffer_cond` | 120 | Cout: 1-16 threads, fewer buffers than threads; out output-model x threads |
| x | `output_file_eligible_request_cond` | 136 | Cout: 1-16 threads, fewer buffers than threads; out output-model x threads |
| x | `output_file_write_victim_cond` | 139 | Cout: 1-16 threads, fewer buffers than threads; out output-model x threads |
| x | `output_file_get_free_buffer` | 148 | Cout |
| x | `output_file_get_next_buffer_to_write` | 166 | Cout |
| x | `output_file_release_buffer` | 197 | Cout |
| x | `output_file_request_buffer` | 205 | Cout |
| x | `output_file_request_buffer_extension` | 243 | Cout |
| x | `output_file_next_block_id` | 259 | Cout |
| x | `output_file_return_buffer` | 269 | Cout |
| x | `vofprintf` | 323 | Cout header line; sam test_header_lines |
| x | `ofprintf` | 332 | Cout header line; sam test_header_lines |

### `src/io/output_map.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `output_map_parameters_set_defaults` | 44 | map test_map_format_3_differs_only_in_homopolymer_marks: default == --map-format 2 |
| x | `output_map_print_cigar_mapv2` | 50 | map test_se_map_matches_sam, test_pe_map_matches_sam: CIGAR rebuilt from SAM CIGAR + read + reference |
|  | `output_map_print_cigar_mapv3` | 89 | partial: map: v3 == v2 bytes; the '&' homopolymer branches are set only by the run-length text path and were not reached, also not with a run-length index |
| x | `output_map_print_match` | 138 | map |
| x | `output_map_print_paired_match` | 185 | map |
| x | `output_map_print_separator` | 228 | map |
|  | `output_map_single_end_print_tag` | 237 | partial: map: name, /1 and /2; CASAVA and extra-tag branches not tested, known bug D-map-tag-comments.md |
| x | `output_map_paired_end_print_tag` | 272 | map test_pe_map_matches_sam; comment branches as above |
| x | `output_map_single_end_print_read__qualities` | 296 | test_d_map.py::test_se_map_matches_sam (with qualities), test_fasta_input_map_lines_have_no_qualities (FASTA, no qualities column) |
|  | `output_map_paired_end_print_read__qualities` | 311 | partial: with qualities only |
|  | `output_map_print_counters_account_mcs` | 337 | partial: map: counters syntax, at most one '+', sum >= reported matches; values pinned by map pins and X4 == MAP column; no independent stratum count |
|  | `output_map_print_counters` | 359 | partial: map: counters syntax, at most one '+', sum >= reported matches; values pinned by map pins and X4 == MAP column; no independent stratum count |
| x | `output_map_single_end_matches` | 416 | map |
| x | `output_map_paired_end_matches` | 453 | map |

### `src/io/output_sam.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `output_sam_parameters_set_defaults` | 55 | sam test_sam_compact_explicit_true_is_default; cli test_documented_defaults_are_the_defaults |
| x | `output_sam_print_header` | 84 | sam test_header_lines, test_bisulfite_index_header_lists_plain_contigs, test_command_line_escaped_in_pg |
|  | `output_sam_print_cigar_element` | 188 | partial: I/D/S through the reduced printer (sam CIGAR checks); its print_mismatches branch is unreachable (no option sets print_mismatches) |
| x | `output_sam_print_reduced_cigar_element` | 217 | sam CIGAR validity, NM/MD/AS recomputed |
| x | `output_sam_print_md_cigar_element` | 234 | sam: MD recomputed from the reference, both strands |
| x | `output_sam_parse_read_group_header` | 285 | sam test_read_group; cli errors rg-prefix, rg-id, rg-long-id; test_rg_id_of_255_characters_is_accepted |
|  | `output_sam_print_opt_field_tag_AS` | 336 | partial: gap-affine: AS recomputed; hamming/edit models are B's options |
| x | `output_sam_print_opt_field_tag_MD` | 380 | sam: MD recomputed from the reference, both strands |
| x | `output_sam_print_opt_field_tag_NM` | 417 | sam |
| x | `output_sam_print_opt_field_tag_RG` | 441 | sam test_read_group, SE and PE |
| x | `output_sam_print_opt_field_tag_XB` | 492 | sam test_bisulfite_records_carry_xb_not_md; rep test_report_bisulfite_pairs counts XB per end; pilot checks XB against truth |
|  | `output_sam_print_opt_field_tag_XS` | 513 | partial: values pinned by test_c_paired.py pin rg-Mall-full, where each end has as many alignments as pairs; otherwise XS reads stale memory or crashes (known bug D-pe-sam-compact-false-crash.md) |
| x | `output_sam_calculate_flag_se` | 568 | sam test_se_records_follow_spec, test_pe_records_follow_spec: every bit from SAMv1 |
| x | `output_sam_calculate_flag_pe` | 590 | sam test_se_records_follow_spec, test_pe_records_follow_spec: every bit from SAMv1 |
| x | `output_sam_calculate_flag_se_match` | 633 | sam test_se_records_follow_spec, test_pe_records_follow_spec: every bit from SAMv1 |
| x | `output_sam_calculate_flag_pe_map` | 645 | sam test_se_records_follow_spec, test_pe_records_follow_spec: every bit from SAMv1 |
| x | `output_sam_print_match_cigar` | 668 | sam |
| x | `output_sam_print_opt_field_tag_X4` | 715 | sam test_gem_compatible_adds_x4_x5; map test_gem_compatible_x4_x5_equal_map_columns |
| x | `output_sam_print_opt_field_tag_X4_pe` | 724 | sam test_gem_compatible_adds_x4_x5; map test_gem_compatible_x4_x5_equal_map_columns |
| x | `output_sam_print_opt_field_tag_X5` | 735 | sam test_gem_compatible_adds_x4_x5; map test_gem_compatible_x4_x5_equal_map_columns |
| x | `output_sam_print_opt_field_tag_X5_pe` | 758 | sam test_gem_compatible_adds_x4_x5; map test_gem_compatible_x4_x5_equal_map_columns |
| x | `output_sam_print_opt_field_tag_XA_match` | 804 | sam test_se_compact_xa_lists_the_secondary_records, test_pe_compact_xa_per_end; map); XA distance of clipped entries is a known bug (D-xa-distance-counts-clipping.md |
| x | `output_sam_print_opt_field_tag_XA_se` | 832 | sam test_se_compact_xa_lists_the_secondary_records, test_pe_compact_xa_per_end; map); XA distance of clipped entries is a known bug (D-xa-distance-counts-clipping.md |
| x | `output_sam_print_opt_field_tag_XA_pe_end1` | 850 | sam test_se_compact_xa_lists_the_secondary_records, test_pe_compact_xa_per_end; map); XA distance of clipped entries is a known bug (D-xa-distance-counts-clipping.md |
| x | `output_sam_print_opt_field_tag_XA_pe_end2` | 868 | sam test_se_compact_xa_lists_the_secondary_records, test_pe_compact_xa_per_end; map); XA distance of clipped entries is a known bug (D-xa-distance-counts-clipping.md |
| x | `output_sam_print_qname` | 886 | sam test_qname_drops_comment_and_mate_suffix |
| x | `output_sam_print_seq__qualities` | 900 | sam: forward/reverse, secondary '*', FASTA QUAL '*' in test_fasta_input_prints_star_qual |
| x | `output_sam_print_core_fields_pe` | 952 | sam |
| x | `output_sam_print_core_fields_se` | 1035 | sam |
| x | `output_sam_print_optional_fields_se` | 1086 | sam |
| x | `output_sam_print_optional_fields_pe` | 1125 | sam |
| x | `output_sam_single_end_matches` | 1169 | sam; pins); RNEXT of an unmapped mate accepted as '=' or name (known bug D-rnext-unmapped-mate.md |
| x | `output_sam_paired_end_matches_unpaired` | 1231 | sam; pins); RNEXT of an unmapped mate accepted as '=' or name (known bug D-rnext-unmapped-mate.md |
| x | `output_sam_paired_end_matches_paired` | 1292 | sam; pins); RNEXT of an unmapped mate accepted as '=' or name (known bug D-rnext-unmapped-mate.md |
| x | `output_sam_paired_end_matches` | 1368 | sam; pins); RNEXT of an unmapped mate accepted as '=' or name (known bug D-rnext-unmapped-mate.md |

### `src/mapper/mapper.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `mapper_se_thread` | 144 | out test_thread_count_does_not_change_output at -t 1,2,3,4,8 over many input blocks, SE/PE/bisulfite; input order kept |
| x | `mapper_pe_thread` | 221 | out test_thread_count_does_not_change_output at -t 1,2,3,4,8 over many input blocks, SE/PE/bisulfite; input order kept |
| x | `mapper_run` | 312 | out test_thread_count_does_not_change_output at -t 1,2,3,4,8 over many input blocks, SE/PE/bisulfite; input order kept |
| x | `mapper_se_run` | 399 | out test_thread_count_does_not_change_output at -t 1,2,3,4,8 over many input blocks, SE/PE/bisulfite; input order kept |
| x | `mapper_pe_run` | 402 | out test_thread_count_does_not_change_output at -t 1,2,3,4,8 over many input blocks, SE/PE/bisulfite; input order kept |

### `src/mapper/mapper_io.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `mapper_load_index` | 89 | every run; cli test_verbose_user_log logs it |
| x | `mapper_io_handler_new_se` | 102 | all runs |
| x | `mapper_io_handler_new_pe` | 126 | all runs |
| x | `mapper_io_handler_delete` | 166 | all runs |
| x | `mapper_se_reload_buffer` | 180 | out: 4 KB input buffers, many blocks |
| x | `mapper_se_parse_sequence` | 184 | out: 4 KB input buffers, many blocks |
| x | `mapper_pe_reload_buffer` | 189 | out test_paired_output_follows_input_buffers_not_threads: tie-prone pairs pinned with default and 4 KB buffers, the per-buffer template-length reset; thread-count invariance); its unsynchronized-input errors not tested (known bug D-unsynced-pairs-hang.md |
| x | `mapper_pe_parse_sequence` | 235 | split and bisulfite pairs |
| x | `mapper_read_sequence` | 261 | out: 4 KB input buffers, many blocks |
| x | `mapper_read_paired_sequence` | 295 | split and bisulfite pairs |
| x | `mapper_io_handler_output_matches` | 328 | SAM and MAP, with and without --report-file |
| x | `mapper_io_handler_output_paired_matches` | 361 | SAM and MAP, with and without --report-file |

### `src/mapper/mapper_parameters.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `mapper_parameters_set_defaults_io` | 32 | cli test_documented_defaults_are_the_defaults: explicit --help defaults change no byte |
| x | `mapper_parameters_set_defaults_system` | 58 | cli test_documented_defaults_are_the_defaults: explicit --help defaults change no byte |
|  | `mapper_parameters_set_defaults_cuda` | 65 | indirect: runs; no effect in a --disable-cuda build |
| x | `mapper_parameters_set_defaults_misc` | 89 | cli test_documented_defaults_are_the_defaults: explicit --help defaults change no byte |
| x | `mapper_parameters_set_defaults` | 98 | cli test_documented_defaults_are_the_defaults: explicit --help defaults change no byte |

### `src/mapper/mapper_profile.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `mapper_profile_print_mapper_se` | 166 | partial: cli test_profile calls them; empty in release builds, their GEM_PROFILE output is not checked |
|  | `mapper_profile_print_mapper_pe` | 170 | partial: cli test_profile calls them; empty in release builds, their GEM_PROFILE output is not checked |

### `src/mapper/mapper_stats.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `mapper_stats_new` | 45 | Cst template_length_estimate: pooled mean and sum of squares, getter formulas |
| x | `mapper_stats_clear` | 53 | Cst template_length_estimate: pooled mean and sum of squares, getter formulas |
| x | `mapper_stats_delete` | 64 | Cst template_length_estimate: pooled mean and sum of squares, getter formulas |
| x | `mapper_stats_template_init` | 70 | Cst template_length_estimate: pooled mean and sum of squares, getter formulas |
| x | `mapper_stats_template_length_sample` | 91 | Cst template_length_estimate: pooled mean and sum of squares, getter formulas |
| x | `mapper_stats_template_length_get_num_samples` | 103 | Cst template_length_estimate: pooled mean and sum of squares, getter formulas |
| x | `mapper_stats_template_length_get_ci_min_samples` | 106 | Cst template_length_estimate: pooled mean and sum of squares, getter formulas |
| x | `mapper_stats_template_length_is_reliable` | 116 | Cst template_length_estimate: pooled mean and sum of squares, getter formulas |
| x | `mapper_stats_template_length_get_mean` | 122 | Cst template_length_estimate: pooled mean and sum of squares, getter formulas |
| x | `mapper_stats_template_length_get_stddev` | 125 | Cst template_length_estimate: pooled mean and sum of squares, getter formulas |
| x | `mapper_stats_template_length_get_expected_max` | 128 | Cst template_length_estimate: pooled mean and sum of squares, getter formulas |
| x | `mapper_stats_template_length_get_expected_min` | 134 | Cst template_length_estimate: pooled mean and sum of squares, getter formulas |
| x | `mapper_stats_template_length_get_sigma_dev` | 140 | Cst template_length_estimate: pooled mean and sum of squares, getter formulas |

### `src/profiler/profiler_timer.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `TIMER_START` | 29 | Cst timers_measure_elapsed_time against nanosleep |
| x | `TIMER_STOP` | 33 | Cst timers_measure_elapsed_time against nanosleep |
| x | `TIMER_PAUSE` | 37 | Cst timers_measure_elapsed_time against nanosleep |
| x | `TIMER_CONTINUE` | 41 | Cst timers_measure_elapsed_time against nanosleep |
| x | `TIMER_RESET` | 44 | Cst timers_measure_elapsed_time against nanosleep |
| x | `TIMER_RESTART` | 47 | Cst timers_measure_elapsed_time against nanosleep |
| x | `TIMER_GET_TOTAL_NS` | 51 | Cst timers_measure_elapsed_time against nanosleep |

### `src/profiler/profiler_vtune.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `PROFILE_VTUNE_START` | 42 | not covered: empty without VTune, no behaviour to test |
|  | `PROFILE_VTUNE_STOP` | 43 | not covered: empty without VTune, no behaviour to test |

### `src/stats/report_stats.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `init_mapping_stats` | 30 | Cst merge_mapping_stats_sums_threads; rep test_report_same_at_any_thread_count |
| x | `merge_ihash` | 46 | Cst merge_mapping_stats_sums_threads; rep test_report_same_at_any_thread_count |
| x | `merge_mapping_stats` | 63 | Cst merge_mapping_stats_sums_threads; rep test_report_same_at_any_thread_count |
| x | `update_counts` | 88 | Cst; rep BaseCounts A/C/G/T and HistReadLen); N counts excluded (known bug D-report-basecounts-extra-n.md |
| x | `update_distance_counts` | 107 | rep HistMismatch recomputed from reads and reference |
| x | `update_conversion_counts` | 130 | rep test_report_bisulfite_pairs: General and conversion-control C2T/G2A counts |
| x | `get_read_type` | 144 | Cst read_type_by_contig_name; rep: genome with the three control contigs |
| x | `collect_se_mapping_stats` | 152 | rep, SE/PE/bisulfite |
| x | `collect_pe_mapping_stats` | 177 | rep, SE/PE/bisulfite |
| x | `output_json_uint_element` | 257 | Cst json_printers; rep test_report_json_layout |
| x | `output_json_uint_array` | 265 | Cst json_printers; rep test_report_json_layout |
| x | `output_mapping_stats` | 295 | rep, SE/PE/bisulfite |

### `src/stats/stats_counter.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `COUNTER_RESET` | 30 | Cst counter_mean_variance_min_max against a two-pass reference |
| x | `COUNTER_ADD` | 33 | Cst counter_mean_variance_min_max against a two-pass reference |
| x | `COUNTER_GET_TOTAL` | 53 | Cst counter_mean_variance_min_max against a two-pass reference |
| x | `COUNTER_GET_NUM_SAMPLES` | 56 | Cst counter_mean_variance_min_max against a two-pass reference |
| x | `COUNTER_GET_MEAN` | 65 | Cst counter_mean_variance_min_max against a two-pass reference |
| x | `COUNTER_GET_VARIANCE` | 68 | Cst counter_mean_variance_min_max against a two-pass reference |
| x | `COUNTER_GET_STDDEV` | 71 | Cst counter_mean_variance_min_max against a two-pass reference |

### `src/system/commons.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `integer_num_ciphers` | 97 | Cfmt integer_to_ascii_matches_printf: 0..199999, powers of ten +-1, 200k random 64-bit values |
| x | `integer_proportion` | 118 | Cfmt integer_proportion_cases |
| x | `integer_to_ascii` | 134 | Cfmt integer_to_ascii_matches_printf: 0..199999, powers of ten +-1, 200k random 64-bit values |
| x | `gem_log2` | 181 | Cfmt gem_log2_matches_log2 |
| x | `gem_loge` | 193 | Cfmt gem_loge_approximation: the quadratic written out independently, error against log |
| x | `gem_rand_IID` | 206 | Cfmt gem_rand_iid_in_range; range only |
| x | `system_get_num_processors` | 273 | Csys system_information_matches_libc; sam test_header_lines @CO; cli test_profile |
| x | `system_get_cwd` | 277 | Csys system_information_matches_libc; sam test_header_lines @CO; cli test_profile |
| x | `system_get_hostname` | 282 | Csys system_information_matches_libc; sam test_header_lines @CO; cli test_profile |
| x | `system_get_user_name` | 287 | Csys system_information_matches_libc; sam test_header_lines @CO; cli test_profile |
| x | `system_print_info` | 293 | Csys system_information_matches_libc; sam test_header_lines @CO; cli test_profile |
| x | `system_get_time` | 320 | Csys system_information_matches_libc; sam test_header_lines @CO; cli test_profile |

### `src/system/errors.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `gem_print_stack_trace` | 50 | partial: runs in GEM_DEBUG builds inside the handler test; its text is not checked; empty in release |
| x | `gem_handle_error_signals` | 76 | Csys error_signal_handler_reports_and_exits_1: SIGSEGV in a child -> report function, message, exit 1 |

### `src/system/gruntime.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `gruntime_init` | 33 | indirect: every C driver and run; the handler they install is tested in Csys |
|  | `gruntime_destroy` | 45 | indirect: every C driver and run; the handler they install is tested in Csys |

### `src/system/gthread.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `gem_thread_register_id` | 33 | Csys thread_ids_are_per_thread, 8 threads |
|  | `gem_thread_cleanup` | 51 | partial: runs at every driver exit via gruntime_destroy; nothing observable to check |

### `src/system/mm_allocator.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `mm_allocator_segment_allocate` | 59 | Cmem: disjoint, kept, zeroed, 8-aligned blocks over several segments |
| x | `mm_allocator_segment_clear` | 75 | Cmem: disjoint, kept, zeroed, 8-aligned blocks over several segments |
| x | `mm_allocator_segment_free` | 82 | Cmem: disjoint, kept, zeroed, 8-aligned blocks over several segments |
| x | `mm_allocator_new_segment` | 91 | Cmem: disjoint, kept, zeroed, 8-aligned blocks over several segments |
| x | `mm_allocator_new` | 102 | Cmem: disjoint, kept, zeroed, 8-aligned blocks over several segments |
| x | `mm_allocator_clear` | 125 | Cmem: disjoint, kept, zeroed, 8-aligned blocks over several segments |
| x | `mm_allocator_delete` | 139 | Cmem: disjoint, kept, zeroed, 8-aligned blocks over several segments |
| x | `mm_allocator_get_segment` | 160 | Cmem: disjoint, kept, zeroed, 8-aligned blocks over several segments |
| x | `mm_allocator_get_num_segments` | 165 | Cmem: disjoint, kept, zeroed, 8-aligned blocks over several segments |
| x | `mm_allocator_allocate_reserve` | 203 | Cmem: disjoint, kept, zeroed, 8-aligned blocks over several segments |
| x | `mm_allocator_allocate` | 231 | Cmem: disjoint, kept, zeroed, 8-aligned blocks over several segments |
| x | `mm_allocator_free_recycle_segment` | 341 | Cmem: random-order frees, survivors intact |
| x | `mm_allocator_free_contiguous_segment` | 359 | Cmem: random-order frees, survivors intact |
| x | `mm_allocator_free_contiguous_memory` | 372 | Cmem: random-order frees, survivors intact |
|  | `mm_allocator_free_malloc` | 381 | not covered: reached only with MM_ALLOCATOR_MALLOC defined (malloc-backed blocks), not in any build here |
| x | `mm_allocator_free` | 415 | Cmem: random-order frees, survivors intact |
| x | `mm_allocator_push_memory_state` | 443 | Cmem allocator_push_pop_releases_later_blocks, nested states |
| x | `mm_allocator_pop_memory_state` | 464 | Cmem allocator_push_pop_releases_later_blocks, nested states |

### `src/system/mm_slab.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `mm_slab_add_new_segment` | 46 | Cmem slab_units_are_disjoint_and_reusable |
| x | `mm_slab_new_` | 85 | Cmem slab_units_are_disjoint_and_reusable |
| x | `mm_slab_delete` | 113 | Cmem slab_units_are_disjoint_and_reusable |
| x | `mm_slab_request` | 128 | Cmem slab_units_are_disjoint_and_reusable |
| x | `mm_slab_return` | 140 | Cmem slab_units_are_disjoint_and_reusable |
| x | `mm_slab_get_slab_size` | 147 | Cmem slab_units_are_disjoint_and_reusable |

### `src/system/report.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `gem_error_get_report_function` | 41 | Csys |
| x | `gem_error_set_report_function` | 44 | Csys |
| x | `gem_error_get_stream` | 47 | Csys |
| x | `gem_log_get_stream` | 53 | Csys |
| x | `gem_is_mute_error_stream` | 64 | Csys |
| x | `gem_is_mute_log_stream` | 67 | Csys |
| x | `vtfprintf` | 77 | Csys log_streams_and_timestamped_print; cli test_verbose_user_log |
| x | `tfprintf` | 104 | Csys log_streams_and_timestamped_print; cli test_verbose_user_log |
| x | `tab_fprintf` | 127 | Csys tabulated_printing |
| x | `fprintf_tabs` | 172 | Csys tabulated_printing |
| x | `tab_global_print` | 181 | Csys tabulated_printing |
| x | `tab_global_inc` | 184 | Csys tabulated_printing |
| x | `tab_global_dec` | 190 | Csys tabulated_printing |
| x | `ticker_count_reset` | 227 | Csys count_ticker; cli test_progress_step, test_verbose_user_log |
| x | `ticker_percentage_finish` | 252 | Csys percentage_ticker_finishes_at_max |
| x | `ticker_count_finish` | 266 | Csys count_ticker; cli test_progress_step, test_verbose_user_log |
| x | `ticker_update` | 280 | Csys count_ticker; cli test_progress_step, test_verbose_user_log |
| x | `ticker_update_mutex` | 301 | Csys count_ticker; cli test_progress_step, test_verbose_user_log |
| x | `ticker_finish` | 306 | Csys count_ticker; cli test_progress_step, test_verbose_user_log |
| x | `ticker_add_process_label` | 322 | Csys count_ticker; cli test_progress_step, test_verbose_user_log |
| x | `ticker_add_finish_label` | 326 | Csys count_ticker; cli test_progress_step, test_verbose_user_log |
| x | `ticker_mutex_enable` | 339 | Csys count_ticker; cli test_progress_step, test_verbose_user_log |
| x | `ticker_mutex_cleanup` | 342 | Csys count_ticker; cli test_progress_step, test_verbose_user_log |

### `src/tools/gem-mapper.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `gem_mapper_open_input_file` | 45 | file, split pairs, stdin in out test_standard_input; '\|' pipes and -z/-j are A's rows |
| x | `gem_mapper_open_input` | 79 | file, split pairs, stdin in out test_standard_input; '\|' pipes and -z/-j are A's rows |
| x | `gem_mapper_close_input` | 96 | file, split pairs, stdin in out test_standard_input; '\|' pipes and -z/-j are A's rows |
| x | `gem_mapper_open_output` | 104 | out: stdout, -o, gzip, bzip2; cli error output-dir |
| x | `gem_mapper_close_output` | 131 | out: stdout, -o, gzip, bzip2; cli error output-dir |
| x | `gem_mapper_print_profile` | 134 | cli test_profile, all five reduce types; system information checked |
| x | `main` | 182 | all CLI tests |

### `src/tools/interface/mapper_arguments.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `parse_integer_system` | 136 | Cfmt parse_integer_system_forms: n, c, xc); 0 and non-numbers are a known bug (D-zero-or-non-numeric-counts.md |
| x | `gem_mapper_print_usage` | 158 | cli test_help_lists_options_by_visibility, test_help_spellings_agree |
|  | `gem_mapper_parameters_presets` | 165 | indirect: runs on every call; the 'exact' preset is B's --mapping-mode row |
| x | `gem_mapper_parameters_check` | 181 | cli test_argument_errors: no index, same files, missing ends, -m with -M, -M 0 |
| x | `gem_mapper_parse_arguments_io` | 239 | cli errors and defaults; out/sam/map option tests; CUDA options rejected |
|  | `gem_mapper_parse_arguments_single_end` | 370 | partial: cli test_documented_defaults_are_the_defaults parses each with its default; behaviour belongs to areas B and C |
|  | `gem_mapper_parse_arguments_paired_end` | 590 | partial: cli test_documented_defaults_are_the_defaults parses each with its default; behaviour belongs to areas B and C |
|  | `gem_mapper_parse_arguments_bisulfite` | 738 | partial: cli test_documented_defaults_are_the_defaults parses each with its default; behaviour belongs to areas B and C |
|  | `gem_mapper_parse_arguments_alignment_score` | 784 | partial: cli test_documented_defaults_are_the_defaults parses each with its default; behaviour belongs to areas B and C |
|  | `gem_mapper_parse_arguments_mapq_score` | 855 | partial: cli test_documented_defaults_are_the_defaults parses each with its default; behaviour belongs to areas B and C |
|  | `gem_mapper_parse_arguments_reporting` | 882 | partial: cli test_documented_defaults_are_the_defaults parses each with its default; behaviour belongs to areas B and C |
| x | `gem_mapper_parse_arguments_output_format` | 904 | cli errors and defaults; out/sam/map option tests; CUDA options rejected |
| x | `gem_mapper_parse_arguments_system` | 947 | cli errors and defaults; out/sam/map option tests; CUDA options rejected |
| x | `gem_mapper_parse_arguments_gpu` | 971 | cli errors and defaults; out/sam/map option tests; CUDA options rejected |
| x | `gem_mapper_parse_arguments_misc` | 1053 | cli errors and defaults; out/sam/map option tests; CUDA options rejected |
| x | `gem_mapper_parse_arguments` | 1139 | cli errors and defaults; out/sam/map option tests; CUDA options rejected |

### `src/utils/ihash.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `ihash_new` | 44 | Cfmt ihash_insert_get_sort; rep histograms sorted by key |
| x | `ihash_clear` | 50 | Cfmt ihash_insert_get_sort; rep histograms sorted by key |
| x | `ihash_delete` | 61 | Cfmt ihash_insert_get_sort; rep histograms sorted by key |
| x | `ihash_get_ihash_element` | 68 | Cfmt ihash_insert_get_sort; rep histograms sorted by key |
| x | `ihash_insert_element` | 73 | Cfmt ihash_insert_get_sort; rep histograms sorted by key |
| x | `ihash_get_element` | 91 | Cfmt ihash_insert_get_sort; rep histograms sorted by key |
| x | `ihash_cmp_keys` | 115 | Cfmt ihash_insert_get_sort; rep histograms sorted by key |
| x | `ihash_sort_by_key` | 124 | Cfmt ihash_insert_get_sort; rep histograms sorted by key |

### `src/utils/options_menu.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `options_get_num_options` | 26 | Cfmt options_menu_adaptors with getopt_long; cli help tests |
| x | `options_adaptor_getopt` | 31 | Cfmt options_menu_adaptors with getopt_long; cli help tests |
| x | `options_adaptor_getopt_short` | 48 | Cfmt options_menu_adaptors with getopt_long; cli help tests |
| x | `options_fprint_menu` | 69 | Cfmt options_menu_adaptors with getopt_long; cli help tests |

### `src/utils/string_buffer.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `string_init` | 35 | Cfmt string_buffer_operations, with and without an mm_allocator |
| x | `string_resize` | 50 | Cfmt string_buffer_operations, with and without an mm_allocator |
| x | `string_clear` | 71 | Cfmt string_buffer_operations, with and without an mm_allocator |
| x | `string_destroy` | 75 | Cfmt string_buffer_operations, with and without an mm_allocator |
| x | `string_set_buffer` | 85 | Cfmt string_buffer_operations, with and without an mm_allocator |
| x | `string_get_buffer` | 93 | Cfmt string_buffer_operations, with and without an mm_allocator |
| x | `string_set_length` | 96 | Cfmt string_buffer_operations, with and without an mm_allocator |
| x | `string_get_length` | 101 | Cfmt string_buffer_operations, with and without an mm_allocator |
| x | `string_append_char` | 113 | Cfmt string_buffer_operations, with and without an mm_allocator |
| x | `string_append_eos` | 120 | Cfmt string_buffer_operations, with and without an mm_allocator |
| x | `string_right_append_buffer` | 166 | Cfmt string_buffer_operations, with and without an mm_allocator |
| x | `string_right_append_string` | 176 | Cfmt string_buffer_operations, with and without an mm_allocator |
| x | `string_is_null` | 230 | Cfmt string_buffer_operations, with and without an mm_allocator |
| x | `string_copy` | 284 | Cfmt string_buffer_operations, with and without an mm_allocator |
| x | `gem_strncpy` | 340 | Cfmt gem_string_helpers_match_libc |
| x | `gem_strcmp` | 352 | Cfmt gem_string_helpers_match_libc |
| x | `gem_streq` | 358 | Cfmt gem_string_helpers_match_libc |
| x | `gem_strcaseeq` | 361 | Cfmt gem_string_helpers_match_libc |
| x | `gem_strcat` | 373 | Cfmt gem_string_helpers_match_libc |
| x | `gem_strlen` | 381 | Cfmt gem_string_helpers_match_libc |

### `src/utils/vector.c`

| Done | Function | Line | Test |
|---|---|---|---|
| x | `vector_new_` | 36 | Cfmt vector_grows_and_keeps_elements |
| x | `vector_reserve` | 46 | Cfmt vector_grows_and_keeps_elements |
| x | `vector_delete` | 78 | Cfmt vector_grows_and_keeps_elements |

## CUDA-only (not reachable in `--disable-cuda` builds; no tests needed)

- `gpu_buffer_collection_new` (`src/gpu/gpu_buffer_collection.c`:137)
- `gpu_buffer_collection_delete` (`src/gpu/gpu_buffer_collection.c`:146)
- `gpu_supported` (`src/gpu/gpu_config.c`:44)
- `mapper_load_gpu_index` (`src/mapper/mapper_cuda.c`:49)
- `mapper_cuda_prepare_ticker` (`src/mapper/mapper_cuda.c`:78)
- `mapper_cuda_setup_thread` (`src/mapper/mapper_cuda.c`:89)
- `mapper_cuda_run` (`src/mapper/mapper_cuda.c`:110)
- `mapper_cuda_se_run` (`src/mapper/mapper_cuda.c`:193)
- `mapper_cuda_pe_run` (`src/mapper/mapper_cuda.c`:199)
- `mapper_profile_print_mapper_se_cuda` (`src/mapper/mapper_profile_cuda.c`:439)
- `mapper_profile_print_mapper_pe_cuda` (`src/mapper/mapper_profile_cuda.c`:444)
