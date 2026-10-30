# gem3-mapper test targets

Everything the new suite has to cover, listed before the tests were written.
Area writers tick the **Done** column (`x`) and name the covering test in
**Test** when a test exercises the row through its inputs and outputs. A row
may stay unticked with a reason in **Test** (for example `known bug`,
`unreachable`, `covered indirectly by ...`).

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
|  | `--index` | `-I` | required | user | A | memory-mapped vs loaded index |  |
|  | `--input` | `-i` | required | user | A | FASTA/FASTQ, stdin when absent; a trailing '\|' opens a pipe |  |
|  | `--i1` | `-1` | required | user | A | split paired input, end 1 (sets paired mode) |  |
|  | `--i2` | `-2` | required | user | A | split paired input, end 2 |  |
|  | `--gzip-input` | `-z` | none | user | A | only consulted for stdin |  |
|  | `--bzip-input` | `-j` | none | user | A | only consulted for stdin |  |
|  | `--input-model` |  | required | developer | A | block_size,num_blocks,buffer_size; the buffer size changes PE output (README); known bug: PE hangs when num_blocks < threads |  |
|  | `--clipping` |  | optional (=) | advanced | A | none, uncalled, masked, fixed,<l>,<r> |  |
|  | `--quality-format` | `-q` | required | advanced | A | ignore, offset-33, offset-64 |  |
|  | `--candidate-generation` |  | required | developer | A | region-profile strategy: fixed, CKS, OPS, factors, MEM, SMEM, OPP, test; known bug: 'adaptive' and 'adaptive-limited' rejected |  |
|  | `--candidate-generation-adaptive` |  | required | developer | A | <app_threshold>,<app_steps>,<app_dec> |  |
|  | `--mapping-mode` |  | required | user | B | fast, sensitive\|hybrid, customed\|complete, exact, complete-pure, complete-brute-force, complete-partition, complete-partition-filter |  |
|  | `--complete-search-error` | `-E` | required | advanced | B |  |  |
|  | `--complete-strata-after-best` | `-s` | required | advanced | B |  |  |
|  | `--alignment-max-error` | `-e` | required | user | B | number or fraction |  |
|  | `--alignment-max-bandwidth` |  | required | advanced | B |  |  |
|  | `--alignment-global-min-identity` |  | required | user | B |  |  |
|  | `--alignment-global-min-score` |  | required | user | B |  |  |
|  | `--alignment-local` |  | required | user | B | never, if-unmapped |  |
|  | `--alignment-local-min-identity` |  | required | user | B |  |  |
|  | `--alignment-local-min-score` |  | required | user | B |  |  |
|  | `--alignment-force-full-swg` |  | optional (=) | advanced | B |  |  |
|  | `--alignment-scaffolding-min-matching_length` |  | required | developer | B |  |  |
|  | `--alignment-curation` |  | optional (=) | advanced | B |  |  |
|  | `--alignment-curation-min-end-context` |  | required | developer | B |  |  |
|  | `--candidate-verification` |  | required | developer | B | BPM, chained |  |
|  | `--qgram-filter` |  | required | developer | B | <kmer_tiles>,<qgram_length> |  |
|  | `--alignment-model` |  | required | advanced | B | none\|pseudoalignment, hamming, edit\|levenshtein, gap-affine |  |
|  | `--gap-affine-penalties` |  | required | user | B | A,B,O,X |  |
|  | `--matching-score` | `-A` | required | advanced | B |  |  |
|  | `--mismatch-penalty` | `-B` | required | advanced | B |  |  |
|  | `--gap-open-penalty` | `-O` | required | advanced | B |  |  |
|  | `--gap-extension-penalty` | `-X` | required | advanced | B |  |  |
|  | `--paired-end-alignment` | `-p` | none | user | C | interleaved paired input |  |
|  | `--min-template-length` | `-l` | required | user | C |  |  |
|  | `--max-template-length` | `-L` | required | user | C |  |  |
|  | `--discordant-pair-search` |  | required | user | C | always, if-no-concordant, never |  |
|  | `--pair-orientation` |  | required | advanced | C | known bug: hangs on every valid value |  |
|  | `--discordant-pair-orientation` |  | required | advanced | C | known bug: hangs |  |
|  | `--pair-layout` |  | required | advanced | C | known bug: hangs |  |
|  | `--discordant-pair-layout` |  | required | advanced | C | known bug: hangs |  |
|  | `--pe-template-length` |  | required | advanced | C | known bug: 3 values rejected, 2 values crash |  |
|  | `--bisulfite-read` |  | required | user | C | inferred, 1, 2, interleaved, non-stranded; known bugs: SE '2' and 'interleaved' act as '1', PE 'non-stranded' ignored |  |
|  | `--underconversion_sequence` |  | required | user | C |  |  |
|  | `--overconversion_sequence` |  | required | user | C |  |  |
|  | `--control_sequence` |  | required | user | C |  |  |
|  | `--mapq-model` |  | required | advanced | C | none, gem, classify, dump-predictors |  |
|  | `--min-reported-strata` | `-m` | required | advanced | C | cannot combine with -M |  |
|  | `--max-reported-matches` | `-M` | required | user | C | number or 'all'; 0 rejected |  |
|  | `--output` | `-o` | required | user | D |  |  |
|  | `--gzip-output` |  | none | user | D |  |  |
|  | `--bzip-output` |  | none | user | D |  |  |
|  | `--output-model` |  | required | developer | D | buffer_size,num_buffers |  |
|  | `--report-file` |  | required | user | D | JSON mapping statistics |  |
|  | `--output-format` | `-F` | required | user | D | MAP, SAM |  |
|  | `--sam-compact` |  | optional (=) | user | D | known bug (upstream #37): a space-separated value is ignored |  |
|  | `--sam-read-group-header` | `-r` | required | user | D |  |  |
|  | `--sam-gem-compatible` |  | optional (=) | developer | D |  |  |
|  | `--map-format` |  | required | developer | D | 1, 2, 3 |  |
|  | `--threads` | `-t` | required | user | D | number or <x>c (cores) |  |
|  | `--max-memory` |  | required | developer | D |  |  |
|  | `--tmp-folder` |  | required | developer | D |  |  |
|  | `--check-alignments` | `-c` | required | developer | B | correct, first-best\|best, all-best, complete, none (archive_check.c) |  |
|  | `--profile` |  | optional (=) | developer | D | sum, min, max, mean, sample |  |
|  | `--verbose` | `-v` | optional (=) | user | D | quiet, user, dev\|debug |  |
|  | `--progress-step` |  | required | developer | D |  |  |
|  | `--version` |  | none | user | D |  |  |
|  | `--help` | `-h` | optional (=) | user | D | user, advanced, developer |  |
|  | `(no option)` |  |  |  | D | argument errors: no -I, same file twice, -m with -M, unknown option, positional arguments |  |

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
|  | `approximate_search_region_profile_stats` | 59 |  |
|  | `approximate_search_region_profile` | 109 |  |

### `src/archive/archive.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `archive_read_mem` | 43 |  |
|  | `archive_read` | 65 |  |
|  | `archive_delete` | 71 |  |
|  | `archive_print` | 92 |  |

### `src/archive/archive_text.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `archive_text_read_mem` | 64 |  |
|  | `archive_text_delete` | 81 |  |
|  | `archive_text_get_size` | 91 |  |
|  | `archive_text_get_position_strand` | 97 |  |
|  | `archive_text_get_projection` | 107 |  |
|  | `archive_text_retrieve` | 116 |  |
|  | `archive_text_print` | 186 |  |

### `src/archive/archive_text_rl.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `archive_text_rl_encode` | 33 |  |
|  | `archive_text_rl_position_translate` | 67 |  |
|  | `archive_text_rl_get_run_length` | 95 |  |
|  | `archive_text_rl_get_decoded_offset_exl` | 105 |  |

### `src/archive/locator.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `locator_intervals_init` | 45 |  |
|  | `locator_read_mem` | 61 |  |
|  | `locator_delete` | 85 |  |
|  | `locator_get_size` | 96 |  |
|  | `locator_interval_get_tag` | 112 |  |
|  | `locator_lookup_interval_index_4way` | 150 |  |
|  | `locator_lookup_interval` | 187 |  |
|  | `locator_lookup_rl_interval_index` | 195 |  |
|  | `locator_lookup_rl_interval` | 216 |  |
|  | `locator_map` | 225 |  |
|  | `locator_inverse_map` | 245 |  |
|  | `locator_inverse_map_position` | 271 |  |
|  | `locator_interval_print` | 283 |  |
|  | `locator_print_summary` | 315 |  |
|  | `locator_print` | 327 |  |

### `src/archive/sampled_rl.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `sampled_rl_read_mem` | 56 |  |
|  | `sampled_rl_delete` | 78 |  |
|  | `sampled_rl_get_size` | 87 |  |
|  | `sampled_rl_get_sample` | 97 |  |
|  | `sampled_rl_print` | 104 |  |

### `src/filtering/region_profile/region_profile.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `region_profile_init` | 39 |  |
|  | `region_profile_model_init` | 48 |  |
|  | `region_profile_clear` | 63 |  |
|  | `region_profile_inject_mm` | 72 |  |
|  | `region_profile_allocate_regions` | 86 |  |
|  | `region_profile_has_exact_matches` | 103 |  |
|  | `region_profile_query_character` | 114 |  |
|  | `region_profile_query_regions` | 179 |  |
|  | `region_profile_extend_last_region` | 190 |  |
|  | `region_profile_fill_gaps_add` | 215 |  |
|  | `region_profile_fill_gaps` | 226 |  |
|  | `region_profile_merge_small_regions` | 276 |  |
|  | `region_profile_compute_error_limits` | 346 |  |
|  | `region_profile_compute_kmer_frequency` | 363 |  |
|  | `region_search_cmp_position` | 446 |  |
|  | `region_profile_sort_by_position` | 453 |  |
|  | `region_profile_print_region_pretty` | 483 |  |
|  | `region_profile_print_pretty` | 524 |  |

### `src/filtering/region_profile/region_profile_adaptive.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `region_profile_generator_save_cut_point` | 54 |  |
|  | `region_profile_generator_restart` | 59 |  |
|  | `region_profile_generator_init` | 73 |  |
|  | `region_profile_generator_close_region` | 96 |  |
|  | `region_profile_generator_close_profile` | 110 |  |
|  | `region_profile_generator_add_character` | 140 |  |
|  | `region_profile_generator_disallow_character` | 192 |  |
|  | `region_profile_generate_adaptive` | 248 |  |
|  | `region_profile_generate_adaptive_limited` | 304 |  |

### `src/filtering/region_profile/region_profile_fixed.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `region_profile_generate_fixed` | 97 |  |
|  | `region_profile_generate_cks` | 134 |  |
|  | `region_profile_generate_factors` | 174 |  |

### `src/filtering/region_profile/region_profile_mem.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `region_profile_compute_mem` | 39 |  |
|  | `region_profile_generate_mem` | 69 |  |
|  | `smem_interval_vector_append` | 115 |  |
|  | `region_profile_compute_smem` | 127 |  |
|  | `region_profile_generate_smem` | 203 |  |

### `src/filtering/region_profile/region_profile_optimum.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `region_profile_generate_ops_fetch_region` | 55 |  |
|  | `region_profile_generate_optimum_fixed` | 76 |  |
|  | `region_profile_optimum_compute_occ` | 181 |  |
|  | `region_profile_generate_optimum_variable` | 208 |  |

### `src/filtering/region_profile/region_profile_schedule.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `region_sort_cmp_value` | 36 |  |
|  | `region_sort_sort` | 41 |  |
|  | `region_profile_schedule_exact_matches` | 50 |  |
|  | `region_profile_schedule_exact` | 70 |  |
|  | `region_profile_schedule_exact_best` | 105 |  |
|  | `region_profile_schedule_approximate` | 157 |  |

### `src/fm_index/bwt/bwt_sbasic.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `bwt_sbasic_read_mem` | 308 |  |
|  | `bwt_sbasic_delete` | 331 |  |
|  | `bwt_sbasic_get_size` | 347 |  |
|  | `bwt_sbasic_get_block_location` | 378 |  |
|  | `bwt_sbasic_char_` | 386 |  |
|  | `bwt_sbasic_sampling_is_sampled_` | 394 |  |
|  | `bwt_sbasic_sampling_erank_` | 397 |  |
|  | `bwt_sbasic_erank_` | 403 |  |
|  | `bwt_sbasic_precompute_` | 434 |  |
|  | `bwt_sbasic_erank` | 472 |  |
|  | `bwt_sbasic_prefetch` | 495 |  |
|  | `bwt_sbasic_precompute` | 526 |  |
|  | `bwt_sbasic_precomputed_erank` | 551 |  |
|  | `bwt_sbasic_LF` | 590 |  |
|  | `bwt_sbasic_prefetched_LF` | 595 |  |
|  | `bwt_sbasic_LF__enc` | 601 |  |
|  | `bwt_sbasic_prefetched_LF__enc` | 620 |  |
|  | `bwt_sbasic_print_` | 638 |  |
|  | `bwt_sbasic_print` | 693 |  |

### `src/fm_index/fm_index.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `fm_index_read_mem` | 82 |  |
|  | `fm_index_delete` | 99 |  |
|  | `fm_index_get_length` | 113 |  |
|  | `fm_index_get_proper_length` | 116 |  |
|  | `fm_index_get_size` | 119 |  |
|  | `fm_index_decode` | 128 |  |
|  | `fm_index_print` | 215 |  |

### `src/fm_index/fm_index_query.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `fm_index_2query_init` | 32 |  |
|  | `fm_index_precompute_2erank` | 44 |  |
|  | `fm_2erank_elms_upper_acc_erank` | 76 |  |
|  | `fm_2erank_elms_lower_acc_erank` | 107 |  |
|  | `fm_index_precomputed_2query_forward` | 138 |  |
|  | `fm_index_precomputed_2query_backward` | 154 |  |
|  | `fm_index_2query_forward_query` | 173 |  |
|  | `fm_index_2query_backward_query` | 186 |  |
|  | `fm_index_2query_forward_precompute` | 202 |  |
|  | `fm_index_2query_backward_precompute` | 211 |  |
|  | `fm_index_2query_precomputed_forward_query` | 220 |  |
|  | `fm_index_2query_precomputed_backward_query` | 231 |  |

### `src/fm_index/fm_index_search.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `fm_index_bsearch` | 34 |  |

### `src/fm_index/rank_mtable.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `rank_mtable_init_levels` | 33 |  |
|  | `rank_mtable_read_mem` | 62 |  |
|  | `rank_mtable_delete` | 80 |  |
|  | `rank_mtable_get_size` | 89 |  |
|  | `rank_mquery_new` | 95 |  |
|  | `rank_mquery_add_char` | 99 |  |
|  | `rank_mquery_is_exhausted` | 110 |  |
|  | `rank_mtable_prefetch` | 116 |  |
|  | `rank_mtable_fetch` | 121 |  |
|  | `rank_mtable_print_content_rec` | 132 |  |
|  | `rank_mtable_print` | 168 |  |

### `src/fm_index/sampled_sa.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `sampled_sa_read_mem` | 42 |  |
|  | `sampled_sa_delete` | 66 |  |
|  | `sampled_sa_get_size` | 164 |  |
|  | `sampled_sa_get_sa_sampling_rate` | 167 |  |
|  | `sampled_sa_get_text_sampling_rate` | 170 |  |
|  | `sampled_sa_prefetch_sample` | 173 |  |
|  | `sampled_sa_get_sample` | 178 |  |
|  | `sampled_sa_print_` | 192 |  |
|  | `sampled_sa_print` | 213 |  |

### `src/io/buffered_input_file.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `buffered_input_file_new` | 46 |  |
|  | `buffered_input_file_close` | 72 |  |
|  | `buffered_input_file_get_file_name` | 79 |  |
|  | `buffered_input_file_get_current_line_num` | 88 |  |
|  | `buffered_input_file_eob` | 91 |  |
|  | `buffered_input_file_read_buffer` | 98 |  |
|  | `buffered_input_file_load_next_chunk` | 136 |  |
|  | `buffered_input_file_get_line` | 165 |  |
|  | `buffered_input_file_reload` | 206 |  |
|  | `buffered_input_file_attach_buffered_output` | 232 |  |

### `src/io/input_buffer.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `input_buffer_new` | 39 |  |
|  | `input_buffer_delete` | 55 |  |
|  | `input_buffer_annotate_lines` | 63 |  |
|  | `input_buffer_get_num_lines` | 82 |  |

### `src/io/input_fasta.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `input_fasta_parse_sequence_tag_chomp_end_info` | 57 |  |
|  | `input_fasta_parse_sequence_tag_extra` | 76 |  |
|  | `input_fasta_parse_sequence_tag` | 141 |  |
|  | `input_fasta_parse_sequence_read` | 174 |  |
|  | `input_fasta_parse_sequence_qualities` | 207 |  |
|  | `input_fasta_parse_sequence_components` | 231 |  |
|  | `input_fasta_parse_sequence` | 261 |  |
|  | `input_fasta_parser_prompt_error` | 288 |  |

### `src/io/input_file_sliced.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `input_file_sliced_initialize` | 32 |  |
|  | `input_stream_sliced_open` | 58 |  |
|  | `input_gzip_stream_sliced_open` | 70 |  |
|  | `input_bzip_stream_sliced_open` | 82 |  |
|  | `input_file_sliced_open` | 94 |  |
|  | `input_file_sliced_popen` | 106 |  |
|  | `input_file_sliced_close` | 118 |  |
|  | `input_file_sliced_get_file_name` | 135 |  |
|  | `input_file_sliced_get_next_id` | 140 |  |
|  | `input_file_sliced_request_buffer_empty` | 149 |  |
|  | `input_file_sliced_find_buffer` | 167 |  |
|  | `input_file_sliced_input_buffer_read` | 185 |  |
|  | `input_file_sliced_input_buffer_process` | 193 |  |
|  | `input_file_sliced_process` | 210 |  |
|  | `input_file_sliced_discard_exhausted_buffers` | 223 |  |
|  | `input_file_sliced_read_lines` | 249 |  |
|  | `input_file_sliced_input_buffer_get_current` | 299 |  |

### `src/io/input_text.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `input_text_parse_integer` | 41 |  |
|  | `input_text_parse_double` | 74 |  |
|  | `input_text_parse_size` | 145 |  |
|  | `input_text_parse_csv_arguments` | 187 |  |
|  | `input_text_parse_csv_array` | 206 |  |
|  | `input_text_parse_extended_uint64` | 221 |  |
|  | `input_text_parse_extended_double` | 249 |  |
|  | `input_text_parse_extended_bool` | 263 |  |
|  | `input_text_parse_count_colons_in_field` | 277 |  |

### `src/system/fm.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `fm_initialize` | 124 |  |
|  | `fm_check_file_type` | 163 |  |
|  | `fm_open_file` | 188 |  |
|  | `fm_open_popen` | 218 |  |
|  | `fm_open_FILE` | 271 |  |
|  | `fm_open_gzFILE` | 293 |  |
|  | `fm_open_bzFILE` | 317 |  |
|  | `fm_close` | 342 |  |
|  | `fm_eof` | 387 |  |
|  | `fm_get_file_name` | 415 |  |
|  | `fm_read_mem` | 595 |  |
|  | `fm_prefetch_next` | 651 |  |
|  | `fm_write_mem` | 680 |  |
|  | `fm_bulk_read_fd` | 713 |  |
|  | `fm_bulk_read_file` | 729 |  |
|  | `gem_access` | 768 |  |
|  | `vfmprintf` | 795 |  |

### `src/system/mm.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `mm_get_tmp_folder` | 109 |  |
|  | `mm_set_tmp_folder` | 112 |  |
|  | `mm_malloc_nothrow` | 121 |  |
|  | `mm_malloc_` | 133 |  |
|  | `mm_realloc_nothrow` | 149 |  |
|  | `mm_realloc` | 154 |  |
|  | `mm_free` | 161 |  |
|  | `mm_bulk_malloc` | 169 |  |
|  | `mm_bulk_mmalloc` | 188 |  |
|  | `mm_bulk_mmalloc_temp` | 229 |  |
|  | `mm_bulk_free` | 271 |  |
|  | `mm_bulk_mload_file` | 351 |  |
|  | `mm_get_mem` | 367 |  |
|  | `mm_skip_align_4KB` | 502 |  |
|  | `mm_read_uint64` | 520 |  |
|  | `mm_read_mem` | 544 |  |
|  | `mm_get_stat_meminfo` | 618 |  |
|  | `mm_get_mem_available_cached` | 651 |  |
|  | `mm_get_mem_available_free` | 656 |  |
|  | `mm_get_mem_available_total` | 661 |  |

### `src/text/dna_text.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `dna_text_read_mem` | 226 |  |
|  | `dna_text_delete` | 244 |  |
|  | `dna_text_get_length` | 298 |  |
|  | `dna_text_get_size` | 306 |  |
|  | `dna_text_retrieve_sequence` | 325 |  |
|  | `dna_text_print` | 341 |  |

### `src/text/sequence.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `sequence_init` | 31 |  |
|  | `sequence_clear` | 49 |  |
|  | `sequence_get_length` | 67 |  |
|  | `sequence_get_tag` | 70 |  |
|  | `sequence_get_read` | 79 |  |
|  | `sequence_get_qualities` | 88 |  |
|  | `sequence_has_qualities` | 97 |  |
|  | `sequence_has_casava_tag` | 100 |  |
|  | `sequence_has_extra_tag` | 103 |  |
|  | `sequence_get_end_info` | 106 |  |
|  | `sequence_print` | 152 |  |

### `src/text/sequence_qualities_model.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `sequence_qualities_model_process` | 31 |  |

### `src/text/text_trace.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `text_trace_destroy` | 30 |  |
|  | `text_trace_compose_padded_text` | 44 |  |

### `src/utils/packed_integer_array.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `packed_integer_array_delete` | 48 |  |
|  | `packed_integer_array_read_mem` | 71 |  |
|  | `packed_integer_array_get_size` | 270 |  |
|  | `packed_integer_array_prefetch` | 281 |  |
|  | `packed_integer_array_load` | 289 |  |
|  | `packed_integer_array_print` | 317 |  |

## Area B: Search and alignment (334 functions in 61 files)

### `include/utils/vector_sort.h`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `vector_sort_selection_align_distance` | 64 |  |
|  | `vector_sort_selection_filtering_positions` | 64 |  |
|  | `vector_sort_selection_match_alignment_region` | 64 |  |
|  | `vector_sort_selection_region_search_by_position` | 64 |  |
|  | `vector_sort_selection_scaffold_coverage` | 64 |  |
|  | `vector_sort_heapify_align_distance` | 84 |  |
|  | `vector_sort_heapify_filtering_positions` | 84 |  |
|  | `vector_sort_heapify_match_alignment_region` | 84 |  |
|  | `vector_sort_heapify_region_search_by_position` | 84 |  |
|  | `vector_sort_heapify_scaffold_coverage` | 84 |  |
|  | `vector_sort_heapsort_align_distance` | 102 |  |
|  | `vector_sort_heapsort_filtering_positions` | 102 |  |
|  | `vector_sort_heapsort_match_alignment_region` | 102 |  |
|  | `vector_sort_heapsort_region_search_by_position` | 102 |  |
|  | `vector_sort_heapsort_scaffold_coverage` | 102 |  |
|  | `vector_sort_quicksort_partition_align_distance` | 118 |  |
|  | `vector_sort_quicksort_partition_filtering_positions` | 118 |  |
|  | `vector_sort_quicksort_partition_match_alignment_region` | 118 |  |
|  | `vector_sort_quicksort_partition_region_search_by_position` | 118 |  |
|  | `vector_sort_quicksort_partition_scaffold_coverage` | 118 |  |
|  | `buffer_sort_align_distance` | 135 |  |
|  | `buffer_sort_filtering_positions` | 135 |  |
|  | `buffer_sort_match_alignment_region` | 135 |  |
|  | `buffer_sort_region_search_by_position` | 135 |  |
|  | `buffer_sort_scaffold_coverage` | 135 |  |
|  | `vector_sort_align_distance` | 176 |  |
|  | `vector_sort_filtering_positions` | 176 |  |
|  | `vector_sort_scaffold_coverage` | 176 |  |

### `src/align/align_bpm.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `align_bpm_compute_matrix` | 39 |  |
|  | `align_bpm_backtrace_matrix` | 182 |  |
|  | `align_bpm_match` | 292 |  |

### `src/align/align_bpm_distance.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `bpm_reset_search_cutoff` | 59 |  |
|  | `bpm_compute_edit_distance` | 172 |  |
|  | `bpm_compute_edit_distance_all` | 283 |  |

### `src/align/align_bpm_pattern.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `bpm_pattern_compile` | 47 |  |
|  | `bpm_pattern_compile_tiles` | 122 |  |

### `src/align/align_kmer_filter_nway.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `kmer_counting_compile_nway_allocate_tables` | 70 |  |
|  | `kmer_counting_compile_nway` | 95 |  |
|  | `kmer_counting_prepare_tiling` | 177 |  |
|  | `kmer_counting_increment_kmer_count` | 242 |  |
|  | `kmer_counting_min_bound_nway_ks1` | 288 |  |
|  | `kmer_counting_min_bound_nway_ks2` | 341 |  |
|  | `kmer_counting_min_bound_nway_ks3` | 390 |  |
|  | `kmer_counting_min_bound_nway_ks4` | 442 |  |
|  | `kmer_counting_min_bound_nway_ksn` | 500 |  |
|  | `kmer_counting_min_bound_nway` | 573 |  |

### `src/align/align_swg.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `align_swg_traceback` | 33 |  |
|  | `align_swg_allocate_table` | 149 |  |
|  | `align_swg_init_table` | 171 |  |
|  | `align_swg_init_table_banded` | 190 |  |
|  | `align_swg_allocate__init_column_banded` | 233 |  |
|  | `align_swg_base` | 287 |  |
|  | `align_swg` | 429 |  |

### `src/align/align_swg_banded.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `align_swg_banded_compute_column` | 33 |  |
|  | `align_swg_banded_compute_column__max_score` | 79 |  |
|  | `align_swg_banded` | 136 |  |
|  | `align_swg_banded_extend` | 216 |  |

### `src/align/align_swg_score.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `align_swg_score_deletion` | 32 |  |
|  | `align_swg_score_insertion` | 37 |  |
|  | `align_swg_score_mismatch` | 42 |  |
|  | `align_swg_score_match` | 45 |  |
|  | `align_swg_score_cigar_element` | 48 |  |
|  | `align_swg_score_cigar_excluding_clipping` | 106 |  |
|  | `align_swg_score_compute_max_score_bound` | 138 |  |
|  | `align_swg_score_compute_max_edit_bound` | 150 |  |

### `src/align/alignment.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `alignment_init` | 37 |  |
|  | `alignment_verify_edit_kmer` | 129 |  |
|  | `alignment_verify_edit_bpm_tile` | 161 |  |
|  | `alignment_verify_edit_bpm` | 177 |  |
|  | `alignment_check` | 242 |  |
|  | `alignment_print_pretty` | 335 |  |

### `src/align/pattern/pattern.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `pattern_init_encode` | 36 |  |
|  | `pattern_init_encode_rl` | 112 |  |
|  | `pattern_init` | 138 |  |
|  | `pattern_tiling_init` | 230 |  |
|  | `pattern_tiling_next` | 257 |  |

### `src/align/pattern/pattern_tiled.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `pattern_tiled_compile` | 31 |  |

### `src/approximate_search/approximate_search.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `approximate_search_init` | 67 |  |
|  | `approximate_search_prepare` | 104 |  |
|  | `approximate_search` | 154 |  |

### `src/approximate_search/approximate_search_filtering_adaptive.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `as_filtering_control_begin` | 34 |  |
|  | `as_filtering_control_filtering_adaptive_next_state` | 43 |  |
|  | `approximate_search_filtering_adaptive` | 68 |  |

### `src/approximate_search/approximate_search_generate_candidates.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `approximate_search_generate_candidates_limit_exact_matches` | 59 |  |
|  | `approximate_search_generate_candidates` | 73 |  |

### `src/approximate_search/approximate_search_hybrid.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `as_hybrid_control_begin` | 40 |  |
|  | `as_hybrid_control_filtering_adaptive_next_state` | 49 |  |
|  | `as_hybrid_control_neighborhood_next_state` | 86 |  |
|  | `approximate_search_hybrid` | 107 |  |

### `src/approximate_search/approximate_search_neighborhood.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `approximate_search_neighborhood_search_brute_force` | 78 |  |
|  | `approximate_search_neighborhood_search_partition` | 109 |  |
|  | `approximate_search_neighborhood_search_partition_preconditioned` | 127 |  |

### `src/approximate_search/approximate_search_stages.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `approximate_search_exact_filtering_adaptive` | 71 |  |
|  | `approximate_search_neighborhood` | 182 |  |
|  | `approximate_search_align_local` | 200 |  |
|  | `approximate_search_end` | 218 |  |

### `src/approximate_search/approximate_search_verify_candidates.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `approximate_search_verify_candidates` | 49 |  |
|  | `approximate_search_verify_extend_candidate` | 68 |  |

### `src/archive/search/archive_check.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `archive_check_se_match_retrieve_text` | 52 |  |
|  | `archive_check_se_match_check_optimum` | 88 |  |
|  | `archive_check_se_match_print` | 121 |  |
|  | `archive_check_se_match_print_incorrect` | 160 |  |
|  | `archive_check_se_match_print_suboptimum` | 170 |  |
|  | `archive_check_se_matches` | 207 |  |
|  | `archive_check_pe_matches` | 281 |  |

### `src/archive/search/archive_search.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `archive_search_configure` | 50 |  |
|  | `archive_search_se_new` | 60 |  |
|  | `archive_search_pe_new` | 69 |  |
|  | `archive_search_delete` | 92 |  |
|  | `archive_search_prepare_sequence` | 100 |  |
|  | `archive_search_inject_handlers` | 115 |  |

### `src/archive/search/archive_search_handlers.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `archive_search_handlers_new` | 34 |  |
|  | `archive_search_handlers_clear` | 50 |  |
|  | `archive_search_handlers_delete` | 55 |  |
|  | `archive_search_handlers_prepare_se` | 66 |  |
|  | `archive_search_handlers_prepare_pe` | 79 |  |

### `src/archive/search/archive_search_se.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `archive_search_se` | 57 |  |

### `src/archive/search/archive_search_se_parameters.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `search_parameters_init_alignment_swg_penalties` | 38 |  |
|  | `search_parameters_init_alignment` | 69 |  |
|  | `select_parameters_init_mapq` | 97 |  |
|  | `select_parameters_init_gpu` | 101 |  |
|  | `search_parameters_init` | 108 |  |
|  | `search_configure_alignment_match_scores` | 157 |  |
|  | `search_configure_alignment_mismatch_scores` | 167 |  |
|  | `search_instantiate_values` | 203 |  |

### `src/archive/search/archive_search_se_stepwise.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `archive_search_se_stepwise_debug_prologue` | 48 |  |
|  | `archive_search_se_stepwise_finish_search` | 235 |  |

### `src/filtering/candidates/filtering_candidates.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `filtering_candidates_init` | 54 |  |
|  | `filtering_candidates_init_alignment` | 66 |  |
|  | `filtering_candidates_clear` | 83 |  |
|  | `filtering_candidates_destroy` | 108 |  |
|  | `filtering_candidates_inject_handlers` | 121 |  |
|  | `filtering_candidates_allocate_position` | 136 |  |
|  | `filtering_candidates_free_position` | 143 |  |
|  | `filtering_candidates_allocate_region` | 148 |  |
|  | `filtering_candidates_free_region` | 155 |  |
|  | `filtering_candidates_allocate_alignment_tiles` | 170 |  |
|  | `filtering_candidates_free_alignment_tiles` | 177 |  |
|  | `filtering_candidates_get_num_positions` | 187 |  |
|  | `filtering_candidates_get_positions` | 196 |  |
|  | `filtering_candidates_clear_positions` | 200 |  |
|  | `filtering_candidates_get_num_regions` | 213 |  |
|  | `filtering_candidates_set_num_regions` | 217 |  |
|  | `filtering_candidates_get_regions` | 222 |  |
|  | `filtering_candidates_get_num_discarded_regions` | 239 |  |
|  | `filtering_candidates_add_num_discarded_regions` | 248 |  |
|  | `filtering_candidates_get_discarded_regions` | 253 |  |
|  | `filtering_candidates_reserve_discarded_regions` | 257 |  |
|  | `filtering_candidates_clear_discarded_regions` | 263 |  |
|  | `filtering_position_cmp_position` | 276 |  |
|  | `filtering_region_cmp_align_distance` | 279 |  |
|  | `filtering_region_cmp_scaffold_coverage` | 285 |  |
|  | `filtering_candidates_sort_positions` | 295 |  |
|  | `filtering_candidates_sort_regions_by_align_distance` | 303 |  |
|  | `filtering_candidates_sort_regions_by_scaffold_coverage` | 311 |  |
|  | `filtering_candidates_sort_discarded_by_scaffold_coverage` | 315 |  |

### `src/filtering/candidates/filtering_candidates_accessors.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `filtering_candidates_compose_filtering_region_from_positions_exact` | 59 |  |
|  | `filtering_candidates_compose_filtering_region_from_positions` | 110 |  |
|  | `filtering_candidates_add_positions_from_interval` | 171 |  |
|  | `filtering_candidates_add_region_from_group_positions` | 212 |  |
|  | `filtering_candidates_add_region_verified` | 258 |  |

### `src/filtering/candidates/filtering_candidates_align.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `filtering_candidates_align_search_cache` | 45 |  |
|  | `filtering_candidates_align_region` | 62 |  |
|  | `filtering_candidates_align_local_region` | 99 |  |
|  | `filtering_candidates_align_extended_region` | 127 |  |
|  | `filtering_candidates_align_candidates` | 156 |  |
|  | `filtering_candidates_align_extended_candidates` | 227 |  |

### `src/filtering/candidates/filtering_candidates_align_local.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `filtering_candidates_align_local_discarded` | 46 |  |
|  | `filtering_candidates_align_local` | 73 |  |

### `src/filtering/candidates/filtering_candidates_buffered.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `filtering_candidates_buffered_inject_handlers` | 34 |  |

### `src/filtering/candidates/filtering_candidates_classify.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `filtering_candidates_classify_subdominant_match_edit` | 30 |  |
|  | `filtering_candidates_classify_subdominant_match_swg` | 46 |  |
|  | `filtering_candidates_classify_subdominant_match` | 65 |  |

### `src/filtering/candidates/filtering_candidates_extend.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `filtering_candidates_compute_extension_region` | 41 |  |
|  | `filtering_candidates_extend_match` | 69 |  |

### `src/filtering/candidates/filtering_candidates_process.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `filtering_candidates_compose_filtering_regions` | 48 |  |
|  | `filtering_candidates_compute_text_boundaries` | 131 |  |
|  | `filtering_candidates_compute_text_coordinates` | 170 |  |
|  | `filtering_candidates_decode_filtering_positions` | 215 |  |
|  | `filtering_candidates_decode_filtering_position_batch_prefetched` | 236 |  |
|  | `filtering_candidates_process_candidates` | 340 |  |

### `src/filtering/candidates/filtering_candidates_verify.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `filtering_candidates_extend_candidates` | 45 |  |
|  | `filtering_candidates_extend_discarded_candidates` | 71 |  |
|  | `filtering_candidates_verify_candidates` | 100 |  |

### `src/filtering/region/filtering_region.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `filtering_region_retrieve_text` | 53 |  |
|  | `filtering_region_compute_key_trims` | 108 |  |

### `src/filtering/region/filtering_region_align.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `filtering_region_align_debug_prologue` | 40 |  |
|  | `filtering_region_align_debug_epilogue` | 49 |  |
|  | `filtering_region_align_clone` | 87 |  |
|  | `filtering_region_align` | 130 |  |
|  | `filtering_region_align_local` | 205 |  |

### `src/filtering/region/filtering_region_cache.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `filtering_region_cache_init` | 43 |  |
|  | `filtering_region_cache_clear` | 47 |  |
|  | `filtering_region_transient_cache_add` | 58 |  |
|  | `filtering_region_cache_cmp_regions` | 68 |  |
|  | `filtering_region_transient_cache_search` | 80 |  |

### `src/filtering/region/filtering_region_verify.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `filtering_region_verify_hamming_text` | 46 |  |
|  | `filtering_region_verify_hamming` | 67 |  |
|  | `filtering_region_verify_levenshtein` | 104 |  |
|  | `filtering_region_verify` | 146 |  |
|  | `filtering_region_verify_extension` | 195 |  |

### `src/matches/align/match_align.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `match_aling_add_clipping` | 42 |  |
|  | `match_align_exact` | 99 |  |
|  | `match_align_pseudoalignment` | 136 |  |
|  | `match_align_hamming` | 181 |  |
|  | `match_align_levenshtein` | 241 |  |
|  | `match_align_smith_waterman_gotoh` | 302 |  |
|  | `match_align_smith_waterman_gotoh_local` | 355 |  |

### `src/matches/align/match_align_normalize.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `match_align_normalize_cigar_trim` | 31 |  |
|  | `match_align_normalize_cigar` | 64 |  |
|  | `match_align_normalize` | 172 |  |

### `src/matches/align/match_align_swg.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `match_align_swg_compute_alignment_type` | 70 |  |
|  | `match_align_swg` | 111 |  |

### `src/matches/align/match_align_swg_chained.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `match_align_swg_chained_add_region` | 36 |  |
|  | `match_align_swg_chained_gap_full` | 95 |  |
|  | `match_align_swg_chained_gap_extending` | 142 |  |
|  | `match_align_swg_end_extending` | 238 |  |
|  | `match_align_swg_chained_to_first` | 289 |  |
|  | `match_align_swg_chained_from_last` | 326 |  |
|  | `match_align_swg_chained_gap` | 363 |  |
|  | `match_align_swg_chained` | 411 |  |

### `src/matches/align/match_align_swg_local.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `match_align_swg_local_alignment_add_match_compose_cigar` | 43 |  |
|  | `match_align_swg_local_alignment_add_match` | 86 |  |
|  | `match_align_swg_local_alignment` | 132 |  |

### `src/matches/align/match_alignment.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `match_alignment_print_pretty` | 52 |  |

### `src/matches/align/match_alignment_region.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `match_alignment_region_init` | 31 |  |
|  | `match_alignment_region_get_type` | 50 |  |
|  | `match_alignment_region_set_type` | 54 |  |
|  | `match_alignment_region_get_error` | 63 |  |
|  | `match_alignment_region_get_cigar_buffer_offset` | 72 |  |
|  | `match_alignment_region_get_cigar_length` | 76 |  |
|  | `match_alignment_region_get_key_begin` | 80 |  |
|  | `match_alignment_region_get_key_end` | 84 |  |
|  | `match_alignment_region_get_text_begin` | 88 |  |
|  | `match_alignment_region_get_text_end` | 92 |  |
|  | `match_alignment_region_set_cigar_buffer_offset` | 96 |  |
|  | `match_alignment_region_set_cigar_length` | 101 |  |
|  | `match_alignment_region_set_key_begin` | 106 |  |
|  | `match_alignment_region_set_key_end` | 111 |  |
|  | `match_alignment_region_set_text_begin` | 116 |  |
|  | `match_alignment_region_set_text_end` | 121 |  |
|  | `match_alignment_region_text_coverage` | 129 |  |
|  | `match_alignment_region_text_distance` | 133 |  |
|  | `match_alignment_region_text_overlap` | 138 |  |
|  | `match_alignment_region_cmp_key_position` | 147 |  |
|  | `match_alignment_region_cmp_text_position` | 152 |  |
|  | `match_alignment_region_print` | 196 |  |

### `src/matches/align/match_alignment_region_rl.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `match_alignment_region_rl_translate_cigar` | 31 |  |
|  | `match_alignment_region_rl_translate` | 77 |  |

### `src/matches/scaffold/match_scaffold.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `match_scaffold_init` | 45 |  |
|  | `match_scaffold_destroy` | 57 |  |
|  | `match_scaffold_allocate_alignment_region` | 71 |  |
|  | `match_scaffold_free_alignment_region` | 77 |  |
|  | `match_scaffold_rl_translate_regions` | 86 |  |
|  | `match_scaffold_region_chain` | 125 |  |
|  | `match_scaffold_adaptive` | 151 |  |
|  | `match_scaffold_sort_alignment_regions` | 224 |  |
|  | `match_scaffold_print` | 231 |  |

### `src/matches/scaffold/match_scaffold_chain.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `match_scaffold_lis_cmp` | 39 |  |
|  | `match_scaffold_chain_is_compatible` | 49 |  |
|  | `match_scaffold_chain_compute_lis` | 59 |  |
|  | `match_scaffold_chain_backtrace_lis` | 109 |  |
|  | `match_scaffold_chain_alignment_regions` | 126 |  |
|  | `match_scaffold_chain_exact_extend_right` | 149 |  |
|  | `match_scaffold_chain_exact_extend_left` | 182 |  |
|  | `match_scaffold_chain_exact_extend` | 218 |  |
|  | `match_scaffold_chain` | 256 |  |

### `src/matches/scaffold/match_scaffold_compose.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `match_scaffold_compose_add_exact_match` | 47 |  |
|  | `match_scaffold_compose_add_mismatch` | 69 |  |

### `src/matches/scaffold/match_scaffold_levenshtein.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `match_scaffold_levenshtein_allocate` | 44 |  |
|  | `match_scaffold_levenshtein_compose_alignment` | 63 |  |
|  | `match_scaffold_levenshtein_tile_align_bpm` | 127 |  |
|  | `match_scaffold_levenshtein_tile` | 178 |  |
|  | `match_scaffold_levenshtein` | 218 |  |

### `src/neighborhood_search/nsearch_filtering.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `nsearch_filtering` | 34 |  |

### `src/neighborhood_search/nsearch_hamming.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `nsearch_hamming_query` | 35 |  |
|  | `nsearch_hamming_directional_query` | 52 |  |
|  | `nsearch_hamming_terminate` | 74 |  |
|  | `nsearch_hamming_brute_force_step` | 98 |  |
|  | `nsearch_hamming_brute_force` | 124 |  |
|  | `nsearch_hamming_scheduled_search_operation_next` | 154 |  |
|  | `nsearch_hamming_scheduled_search` | 217 |  |
|  | `nsearch_hamming` | 238 |  |
|  | `nsearch_hamming_preconditioned` | 262 |  |

### `src/neighborhood_search/nsearch_levenshtein.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `nsearch_levenshtein_base_step` | 34 |  |
|  | `nsearch_levenshtein_base` | 75 |  |
|  | `nsearch_levenshtein_brute_force_step` | 100 |  |
|  | `nsearch_levenshtein_brute_force` | 156 |  |
|  | `nsearch_levenshtein` | 194 |  |
|  | `nsearch_levenshtein_preconditioned` | 218 |  |

### `src/neighborhood_search/nsearch_levenshtein_control.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `nsearch_levenshtein_candidates_cutoff` | 35 |  |
|  | `nsearch_levenshtein_control_limit_interval` | 66 |  |
|  | `nsearch_levenshtein_terminate` | 81 |  |
|  | `nsearch_levenshtein_scheduled_terminate` | 112 |  |

### `src/neighborhood_search/nsearch_levenshtein_query.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `nsearch_query_init` | 29 |  |
|  | `nsearch_levenshtein_query` | 47 |  |
|  | `nsearch_levenshtein_scheduled_precompute_query` | 68 |  |
|  | `nsearch_levenshtein_scheduled_query_exact` | 88 |  |
|  | `nsearch_levenshtein_scheduled_query` | 118 |  |

### `src/neighborhood_search/nsearch_levenshtein_scheduled.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `nsearch_levenshtein_scheduled_align_distance_pass` | 38 |  |
|  | `nsearch_levenshtein_scheduled_operation_step_next` | 69 |  |
|  | `nsearch_levenshtein_scheduled_operation_step_query` | 103 |  |
|  | `nsearch_levenshtein_scheduled_operation_step_query_exact` | 167 |  |
|  | `nsearch_levenshtein_scheduled_search_next` | 208 |  |
|  | `nsearch_levenshtein_scheduled_search` | 228 |  |

### `src/neighborhood_search/nsearch_levenshtein_state.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `nsearch_levenshtein_state_init` | 29 |  |
|  | `nsearch_levenshtein_state_prepare` | 52 |  |
|  | `nsearch_levenshtein_state_prepare_column` | 68 |  |
|  | `nsearch_levenshtein_state_get_global_align_distance` | 82 |  |
|  | `nsearch_levenshtein_state_get_local_align_distance` | 95 |  |
|  | `nsearch_levenshtein_state_compute_chararacter` | 137 |  |
|  | `nsearch_levenshtein_state_compute_text_banded` | 182 |  |
|  | `nsearch_levenshtein_state_compute_chararacter_banded` | 206 |  |

### `src/neighborhood_search/nsearch_operation.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `nsearch_operation_init` | 30 |  |
|  | `nsearch_operation_chained_prepare_sequence_forward` | 47 |  |
|  | `nsearch_operation_chained_prepare_sequence_reverse` | 57 |  |
|  | `nsearch_operation_chained_prepare` | 69 |  |

### `src/neighborhood_search/nsearch_parameters.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `nsearch_parameters_init` | 26 |  |

### `src/neighborhood_search/nsearch_partition.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `nsearch_partition_compute` | 29 |  |
|  | `nsearch_partition_preconditioned_compute` | 42 |  |
|  | `nsearch_partition_compute_error` | 66 |  |

### `src/neighborhood_search/nsearch_schedule.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `nsearch_schedule_init` | 35 |  |
|  | `nsearch_schedule_inject_mm` | 76 |  |
|  | `nsearch_schedule_add_pending_search` | 84 |  |
|  | `nsearch_schedule_search_step` | 121 |  |
|  | `nsearch_schedule_search` | 190 |  |
|  | `nsearch_schedule_search_preconditioned_step` | 211 |  |
|  | `nsearch_schedule_search_preconditioned` | 278 |  |
|  | `nsearch_schedule_compute_min_error` | 305 |  |

### `src/utils/priority_queue.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `pqueue_new` | 37 |  |
|  | `pqueue_delete` | 51 |  |
|  | `pqueue_top_priority` | 60 |  |
|  | `pqueue_push_` | 72 |  |
|  | `pqueue_pop_` | 87 |  |

## Area C: Paired-end, bisulfite, matches and MAPQ (149 functions in 20 files)

### `src/archive/score/archive_score_logit.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `archive_score_logit` | 30 |  |
|  | `archive_score_logit_unique` | 106 |  |
|  | `archive_score_logit_mmap` | 112 |  |

### `src/archive/score/archive_score_pe.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `archive_score_matches_pe_default` | 80 |  |
|  | `archive_score_matches_pe_classify` | 130 |  |
|  | `archive_score_matches_pe` | 182 |  |

### `src/archive/score/archive_score_se.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `archive_score_probability_to_mapq` | 49 |  |
|  | `archive_score_matches_se_get_model` | 55 |  |
|  | `archive_score_matches_se_logit_unique` | 63 |  |
|  | `archive_score_matches_se_logit_mmap` | 76 |  |
|  | `archive_score_matches_se_default` | 92 |  |
|  | `archive_score_matches_se_classify` | 129 |  |
|  | `archive_score_matches_se` | 190 |  |

### `src/archive/search/archive_search_pe.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `archive_search_pe_continue` | 50 |  |
|  | `archive_search_pe` | 83 |  |

### `src/archive/search/archive_search_pe_extend.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `archive_search_pe_extend_matches_ends` | 42 |  |
|  | `archive_search_pe_extend_matches` | 111 |  |

### `src/archive/search/archive_search_pe_parameters.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `search_paired_parameters_init` | 31 |  |

### `src/archive/search/archive_search_pe_stages.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `archive_search_pe_begin` | 54 |  |
|  | `archive_search_pe_search_end1` | 67 |  |
|  | `archive_search_pe_search_end2` | 100 |  |
|  | `archive_search_pe_use_shortcut_extension` | 113 |  |
|  | `archive_search_pe_use_recovery_extension` | 129 |  |
|  | `archive_search_pe_recovery` | 153 |  |
|  | `archive_search_pe_cross_pair_ends` | 185 |  |
|  | `archive_search_pe_find_pairs` | 200 |  |
|  | `archive_search_pe_end` | 230 |  |

### `src/archive/search/archive_select.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `archive_select_configure_se` | 37 |  |
|  | `archive_select_configure_pe` | 43 |  |
|  | `archive_select_se_matches_discard` | 52 |  |
|  | `archive_select_se_matches` | 72 |  |
|  | `archive_select_pe_matches` | 100 |  |

### `src/archive/search/archive_select_parameters.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `select_parameters_init` | 32 |  |
|  | `select_parameters_instantiate_values` | 47 |  |

### `src/matches/classify/matches_classify.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `matches_classify` | 29 |  |
|  | `matches_classify_pattern_viable` | 72 |  |
|  | `matches_classify_max_matches_searched` | 79 |  |
|  | `matches_classify_min_depth_searched` | 113 |  |
|  | `matches_classify_neighbourhood_fallback` | 156 |  |
|  | `matches_classify_local_alignment_fallback` | 169 |  |
|  | `matches_classify_adjust_max_error_by_strata_after_best` | 178 |  |
|  | `matches_classify_compute_max_search_error` | 196 |  |

### `src/matches/classify/matches_metrics.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `matches_metrics_init` | 29 |  |
|  | `matches_metrics_set_proper_length` | 43 |  |
|  | `matches_metrics_set_read_length` | 48 |  |
|  | `matches_metrics_set_swg_match_score` | 53 |  |
|  | `matches_metrics_set_region_profile_metrics` | 58 |  |
|  | `matches_metrics_add_accepted_candidates` | 67 |  |
|  | `matches_metrics_update` | 80 |  |

### `src/matches/classify/matches_predictors.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `matches_predictors_compute_unique` | 38 |  |
|  | `matches_predictors_compute_mmaps` | 62 |  |
|  | `matches_predictors_compute_se` | 89 |  |
|  | `matches_predictors_compute_pe` | 129 |  |
|  | `matches_predictors_print` | 172 |  |
|  | `matches_predictors_se_print` | 238 |  |
|  | `matches_predictors_pe_print` | 248 |  |

### `src/matches/classify/paired_matches_classify.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `paired_matches_classify` | 30 |  |
|  | `paired_matches_classify_search_accomplished` | 68 |  |
|  | `paired_matches_classify_subdominant_end` | 82 |  |

### `src/matches/match_trace.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `match_trace_get_cigar_buffer` | 29 |  |
|  | `match_trace_get_cigar_length` | 34 |  |
|  | `match_trace_cmp_swg_score` | 43 |  |
|  | `match_trace_cmp_genomic_position` | 61 |  |

### `src/matches/matches.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `matches_new` | 38 |  |
|  | `matches_clear` | 63 |  |
|  | `matches_delete` | 78 |  |
|  | `matches_is_mapped` | 94 |  |
|  | `matches_recompute_metrics` | 97 |  |
|  | `matches_update_mcs` | 125 |  |
|  | `matches_update_limited_exact_matches` | 130 |  |
|  | `matches_index_match` | 142 |  |
|  | `matches_lookup_match` | 152 |  |
|  | `matches_get_num_match_traces` | 166 |  |
|  | `matches_get_num_match_traces_extended` | 169 |  |
|  | `matches_get_primary_match` | 172 |  |
|  | `matches_get_match_traces` | 175 |  |
|  | `matches_traces_sort_by_genomic_position` | 181 |  |
|  | `matches_add_match_trace_insert_sorted` | 188 |  |
|  | `matches_add_match_trace_preserve_sorted` | 210 |  |
|  | `match_trace_locate` | 232 |  |
|  | `match_trace_copy` | 250 |  |
|  | `matches_add_match_trace` | 273 |  |
|  | `matches_add_match_trace_extended` | 327 |  |
|  | `matches_add_match_trace_local_pending` | 343 |  |
|  | `matches_local_pending_add_to_regular_matches` | 350 |  |
|  | `matches_local_pending_add_to_extended_matches` | 363 |  |

### `src/matches/matches_cigar.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `matches_cigar_buffer_add_cigar_element__attr` | 35 |  |
|  | `matches_cigar_buffer_add_cigar_element` | 50 |  |
|  | `matches_cigar_buffer_add_mismatch` | 57 |  |
|  | `matches_cigar_vector_append_insertion` | 68 |  |
|  | `matches_cigar_vector_append_deletion` | 90 |  |
|  | `matches_cigar_vector_append_match` | 112 |  |
|  | `matches_cigar_vector_append_mismatch` | 135 |  |
|  | `matches_cigar_vector_append_cigar_element` | 148 |  |
|  | `matches_cigar_vector_insert_cigar_element` | 170 |  |
|  | `matches_cigar_compute_event_distance` | 224 |  |
|  | `matches_cigar_compute_edit_distance` | 247 |  |
|  | `matches_cigar_compute_edit_distance__excluding_clipping` | 271 |  |
|  | `matches_cigar_compute_matching_bases` | 297 |  |
|  | `matches_cigar_element_effective_length` | 320 |  |
|  | `matches_cigar_compute_effective_length` | 340 |  |
|  | `matches_cigar_compute_error_quality` | 354 |  |
|  | `match_cigar_print` | 434 |  |

### `src/matches/matches_counters.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `matches_counters_new` | 34 |  |
|  | `matches_counters_clear` | 43 |  |
|  | `matches_counters_delete` | 47 |  |
|  | `matches_counters_get_num_counters` | 54 |  |
|  | `matches_counters_get_counts` | 58 |  |
|  | `matches_counters_add` | 71 |  |
|  | `matches_counters_sub` | 85 |  |
|  | `matches_counters_compact` | 96 |  |
|  | `matches_counters_compute_matches_to_report` | 103 |  |
|  | `matches_counters_count_first_subdominant` | 136 |  |

### `src/matches/paired_map.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `paired_map_compute_event_distance` | 29 |  |
|  | `paired_map_compute_edit_distance` | 34 |  |
|  | `paired_map_compute_swg_score` | 39 |  |
|  | `paired_map_compute_error_quality` | 44 |  |

### `src/matches/paired_matches.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `paired_matches_new` | 49 |  |
|  | `paired_matches_clear` | 68 |  |
|  | `paired_matches_delete` | 84 |  |
|  | `paired_matches_is_mapped` | 102 |  |
|  | `paired_matches_get_max_complete_stratum` | 115 |  |
|  | `paired_matches_get_primary_map` | 125 |  |
|  | `paired_matches_get_num_maps` | 131 |  |
|  | `paired_matches_get_maps` | 134 |  |
|  | `paired_matches_get_num_discordant_maps` | 137 |  |
|  | `paired_matches_get_discordant_maps` | 140 |  |
|  | `paired_matches_limit_maps` | 143 |  |
|  | `paired_matches_recompute_metrics` | 149 |  |
|  | `paired_matches_cmp_swg_score` | 167 |  |
|  | `paired_matches_add_map_insert_sorted` | 190 |  |
|  | `paired_matches_add_map` | 212 |  |
|  | `paired_matches_add` | 224 |  |
|  | `paired_matches_compute_orientation` | 274 |  |
|  | `paired_matches_compute_layout` | 299 |  |
|  | `paired_matches_compute_template_length` | 326 |  |
|  | `paired_matches_compute_relation` | 382 |  |
|  | `paired_matches_cross_pair` | 423 |  |
|  | `paired_matches_find_pairs_locate` | 461 |  |
|  | `paired_matches_find_pairs_sorted` | 471 |  |
|  | `paired_matches_find_pairs` | 505 |  |
|  | `paired_matches_find_discordant_pairs` | 536 |  |

### `src/text/sequence_bisulfite.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `sequence_bisulfite_convert` | 32 |  |
|  | `sequence_bisulfite_check_cg_depletion_se` | 46 |  |
|  | `sequence_bisulfite_check_cg_depletion_pe` | 57 |  |
|  | `sequence_bisulfite_process_se` | 75 |  |
|  | `sequence_bisulfite_process_pe` | 101 |  |
|  | `sequence_bisulfite_restore_se` | 123 |  |
|  | `sequence_bisulfite_restore_pe` | 128 |  |

## Area D: Output, CLI, threads and runtime (276 functions in 27 files)

### `src/io/buffered_output_file.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `buffered_output_file_new` | 50 |  |
|  | `buffered_output_file_close` | 57 |  |
|  | `buffered_output_file_request_buffer` | 64 |  |
|  | `buffered_output_file_dump_buffer` | 69 |  |
|  | `buffered_output_file_safety_dump_buffer` | 78 |  |
|  | `buffered_output_file_reserve` | 82 |  |
|  | `bofprintf_uint64` | 136 |  |
|  | `bofprintf_int64` | 139 |  |
|  | `bofprintf_char` | 142 |  |
|  | `bofprintf_string` | 145 |  |

### `src/io/output_buffer.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `output_buffer_new` | 31 |  |
|  | `output_buffer_clear` | 38 |  |
|  | `output_buffer_delete` | 47 |  |
|  | `output_buffer_set_state` | 54 |  |
|  | `output_buffer_get_state` | 57 |  |
|  | `output_buffer_set_incomplete` | 60 |  |
|  | `output_buffer_get_used` | 63 |  |
|  | `output_buffer_get_buffer` | 66 |  |
|  | `bprintf_uint64` | 72 |  |
|  | `bprintf_int64` | 77 |  |
|  | `bprintf_char` | 89 |  |
|  | `bprintf_buffer` | 94 |  |

### `src/io/output_file.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `output_file_init_buffers` | 34 |  |
|  | `output_stream_new` | 67 |  |
|  | `output_gzip_stream_new` | 78 |  |
|  | `output_bzip_stream_new` | 89 |  |
|  | `output_file_close` | 100 |  |
|  | `output_file_serve_buffer_cond` | 120 |  |
|  | `output_file_eligible_request_cond` | 136 |  |
|  | `output_file_write_victim_cond` | 139 |  |
|  | `output_file_get_free_buffer` | 148 |  |
|  | `output_file_get_next_buffer_to_write` | 166 |  |
|  | `output_file_release_buffer` | 197 |  |
|  | `output_file_request_buffer` | 205 |  |
|  | `output_file_request_buffer_extension` | 243 |  |
|  | `output_file_next_block_id` | 259 |  |
|  | `output_file_return_buffer` | 269 |  |
|  | `vofprintf` | 323 |  |
|  | `ofprintf` | 332 |  |

### `src/io/output_map.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `output_map_parameters_set_defaults` | 44 |  |
|  | `output_map_print_cigar_mapv2` | 50 |  |
|  | `output_map_print_cigar_mapv3` | 89 |  |
|  | `output_map_print_match` | 138 |  |
|  | `output_map_print_paired_match` | 185 |  |
|  | `output_map_print_separator` | 228 |  |
|  | `output_map_single_end_print_tag` | 237 |  |
|  | `output_map_paired_end_print_tag` | 272 |  |
|  | `output_map_single_end_print_read__qualities` | 296 |  |
|  | `output_map_paired_end_print_read__qualities` | 311 |  |
|  | `output_map_print_counters_account_mcs` | 337 |  |
|  | `output_map_print_counters` | 359 |  |
|  | `output_map_single_end_matches` | 416 |  |
|  | `output_map_paired_end_matches` | 453 |  |

### `src/io/output_sam.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `output_sam_parameters_set_defaults` | 55 |  |
|  | `output_sam_print_header` | 84 |  |
|  | `output_sam_print_cigar_element` | 188 |  |
|  | `output_sam_print_reduced_cigar_element` | 217 |  |
|  | `output_sam_print_md_cigar_element` | 234 |  |
|  | `output_sam_parse_read_group_header` | 285 |  |
|  | `output_sam_print_opt_field_tag_AS` | 336 |  |
|  | `output_sam_print_opt_field_tag_MD` | 380 |  |
|  | `output_sam_print_opt_field_tag_NM` | 417 |  |
|  | `output_sam_print_opt_field_tag_RG` | 441 |  |
|  | `output_sam_print_opt_field_tag_XB` | 492 |  |
|  | `output_sam_print_opt_field_tag_XS` | 513 |  |
|  | `output_sam_calculate_flag_se` | 568 |  |
|  | `output_sam_calculate_flag_pe` | 590 |  |
|  | `output_sam_calculate_flag_se_match` | 633 |  |
|  | `output_sam_calculate_flag_pe_map` | 645 |  |
|  | `output_sam_print_match_cigar` | 668 |  |
|  | `output_sam_print_opt_field_tag_X4` | 715 |  |
|  | `output_sam_print_opt_field_tag_X4_pe` | 724 |  |
|  | `output_sam_print_opt_field_tag_X5` | 735 |  |
|  | `output_sam_print_opt_field_tag_X5_pe` | 758 |  |
|  | `output_sam_print_opt_field_tag_XA_match` | 804 |  |
|  | `output_sam_print_opt_field_tag_XA_se` | 832 |  |
|  | `output_sam_print_opt_field_tag_XA_pe_end1` | 850 |  |
|  | `output_sam_print_opt_field_tag_XA_pe_end2` | 868 |  |
|  | `output_sam_print_qname` | 886 |  |
|  | `output_sam_print_seq__qualities` | 900 |  |
|  | `output_sam_print_core_fields_pe` | 952 |  |
|  | `output_sam_print_core_fields_se` | 1035 |  |
|  | `output_sam_print_optional_fields_se` | 1086 |  |
|  | `output_sam_print_optional_fields_pe` | 1125 |  |
|  | `output_sam_single_end_matches` | 1169 |  |
|  | `output_sam_paired_end_matches_unpaired` | 1231 |  |
|  | `output_sam_paired_end_matches_paired` | 1292 |  |
|  | `output_sam_paired_end_matches` | 1368 |  |

### `src/mapper/mapper.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `mapper_se_thread` | 144 |  |
|  | `mapper_pe_thread` | 221 |  |
|  | `mapper_run` | 312 |  |
|  | `mapper_se_run` | 399 |  |
|  | `mapper_pe_run` | 402 |  |

### `src/mapper/mapper_io.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `mapper_load_index` | 89 |  |
|  | `mapper_io_handler_new_se` | 102 |  |
|  | `mapper_io_handler_new_pe` | 126 |  |
|  | `mapper_io_handler_delete` | 166 |  |
|  | `mapper_se_reload_buffer` | 180 |  |
|  | `mapper_se_parse_sequence` | 184 |  |
|  | `mapper_pe_reload_buffer` | 189 |  |
|  | `mapper_pe_parse_sequence` | 235 |  |
|  | `mapper_read_sequence` | 261 |  |
|  | `mapper_read_paired_sequence` | 295 |  |
|  | `mapper_io_handler_output_matches` | 328 |  |
|  | `mapper_io_handler_output_paired_matches` | 361 |  |

### `src/mapper/mapper_parameters.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `mapper_parameters_set_defaults_io` | 32 |  |
|  | `mapper_parameters_set_defaults_system` | 58 |  |
|  | `mapper_parameters_set_defaults_cuda` | 65 |  |
|  | `mapper_parameters_set_defaults_misc` | 89 |  |
|  | `mapper_parameters_set_defaults` | 98 |  |

### `src/mapper/mapper_profile.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `mapper_profile_print_mapper_se` | 166 |  |
|  | `mapper_profile_print_mapper_pe` | 170 |  |

### `src/mapper/mapper_stats.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `mapper_stats_new` | 45 |  |
|  | `mapper_stats_clear` | 53 |  |
|  | `mapper_stats_delete` | 64 |  |
|  | `mapper_stats_template_init` | 70 |  |
|  | `mapper_stats_template_length_sample` | 91 |  |
|  | `mapper_stats_template_length_get_num_samples` | 103 |  |
|  | `mapper_stats_template_length_get_ci_min_samples` | 106 |  |
|  | `mapper_stats_template_length_is_reliable` | 116 |  |
|  | `mapper_stats_template_length_get_mean` | 122 |  |
|  | `mapper_stats_template_length_get_stddev` | 125 |  |
|  | `mapper_stats_template_length_get_expected_max` | 128 |  |
|  | `mapper_stats_template_length_get_expected_min` | 134 |  |
|  | `mapper_stats_template_length_get_sigma_dev` | 140 |  |

### `src/profiler/profiler_timer.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `TIMER_START` | 29 |  |
|  | `TIMER_STOP` | 33 |  |
|  | `TIMER_PAUSE` | 37 |  |
|  | `TIMER_CONTINUE` | 41 |  |
|  | `TIMER_RESET` | 44 |  |
|  | `TIMER_RESTART` | 47 |  |
|  | `TIMER_GET_TOTAL_NS` | 51 |  |

### `src/profiler/profiler_vtune.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `PROFILE_VTUNE_START` | 42 |  |
|  | `PROFILE_VTUNE_STOP` | 43 |  |

### `src/stats/report_stats.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `init_mapping_stats` | 30 |  |
|  | `merge_ihash` | 46 |  |
|  | `merge_mapping_stats` | 63 |  |
|  | `update_counts` | 88 |  |
|  | `update_distance_counts` | 107 |  |
|  | `update_conversion_counts` | 130 |  |
|  | `get_read_type` | 144 |  |
|  | `collect_se_mapping_stats` | 152 |  |
|  | `collect_pe_mapping_stats` | 177 |  |
|  | `output_json_uint_element` | 257 |  |
|  | `output_json_uint_array` | 265 |  |
|  | `output_mapping_stats` | 295 |  |

### `src/stats/stats_counter.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `COUNTER_RESET` | 30 |  |
|  | `COUNTER_ADD` | 33 |  |
|  | `COUNTER_GET_TOTAL` | 53 |  |
|  | `COUNTER_GET_NUM_SAMPLES` | 56 |  |
|  | `COUNTER_GET_MEAN` | 65 |  |
|  | `COUNTER_GET_VARIANCE` | 68 |  |
|  | `COUNTER_GET_STDDEV` | 71 |  |

### `src/system/commons.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `integer_num_ciphers` | 97 |  |
|  | `integer_proportion` | 118 |  |
|  | `integer_to_ascii` | 134 |  |
|  | `gem_log2` | 181 |  |
|  | `gem_loge` | 193 |  |
|  | `gem_rand_IID` | 206 |  |
|  | `system_get_num_processors` | 273 |  |
|  | `system_get_cwd` | 277 |  |
|  | `system_get_hostname` | 282 |  |
|  | `system_get_user_name` | 287 |  |
|  | `system_print_info` | 293 |  |
|  | `system_get_time` | 320 |  |

### `src/system/errors.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `gem_print_stack_trace` | 50 |  |
|  | `gem_handle_error_signals` | 76 |  |

### `src/system/gruntime.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `gruntime_init` | 33 |  |
|  | `gruntime_destroy` | 45 |  |

### `src/system/gthread.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `gem_thread_register_id` | 33 |  |
|  | `gem_thread_cleanup` | 51 |  |

### `src/system/mm_allocator.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `mm_allocator_segment_allocate` | 59 |  |
|  | `mm_allocator_segment_clear` | 75 |  |
|  | `mm_allocator_segment_free` | 82 |  |
|  | `mm_allocator_new_segment` | 91 |  |
|  | `mm_allocator_new` | 102 |  |
|  | `mm_allocator_clear` | 125 |  |
|  | `mm_allocator_delete` | 139 |  |
|  | `mm_allocator_get_segment` | 160 |  |
|  | `mm_allocator_get_num_segments` | 165 |  |
|  | `mm_allocator_allocate_reserve` | 203 |  |
|  | `mm_allocator_allocate` | 231 |  |
|  | `mm_allocator_free_recycle_segment` | 341 |  |
|  | `mm_allocator_free_contiguous_segment` | 359 |  |
|  | `mm_allocator_free_contiguous_memory` | 372 |  |
|  | `mm_allocator_free_malloc` | 381 |  |
|  | `mm_allocator_free` | 415 |  |
|  | `mm_allocator_push_memory_state` | 443 |  |
|  | `mm_allocator_pop_memory_state` | 464 |  |

### `src/system/mm_slab.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `mm_slab_add_new_segment` | 46 |  |
|  | `mm_slab_new_` | 85 |  |
|  | `mm_slab_delete` | 113 |  |
|  | `mm_slab_request` | 128 |  |
|  | `mm_slab_return` | 140 |  |
|  | `mm_slab_get_slab_size` | 147 |  |

### `src/system/report.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `gem_error_get_report_function` | 41 |  |
|  | `gem_error_set_report_function` | 44 |  |
|  | `gem_error_get_stream` | 47 |  |
|  | `gem_log_get_stream` | 53 |  |
|  | `gem_is_mute_error_stream` | 64 |  |
|  | `gem_is_mute_log_stream` | 67 |  |
|  | `vtfprintf` | 77 |  |
|  | `tfprintf` | 104 |  |
|  | `tab_fprintf` | 127 |  |
|  | `fprintf_tabs` | 172 |  |
|  | `tab_global_print` | 181 |  |
|  | `tab_global_inc` | 184 |  |
|  | `tab_global_dec` | 190 |  |
|  | `ticker_count_reset` | 227 |  |
|  | `ticker_percentage_finish` | 252 |  |
|  | `ticker_count_finish` | 266 |  |
|  | `ticker_update` | 280 |  |
|  | `ticker_update_mutex` | 301 |  |
|  | `ticker_finish` | 306 |  |
|  | `ticker_add_process_label` | 322 |  |
|  | `ticker_add_finish_label` | 326 |  |
|  | `ticker_mutex_enable` | 339 |  |
|  | `ticker_mutex_cleanup` | 342 |  |

### `src/tools/gem-mapper.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `gem_mapper_open_input_file` | 45 |  |
|  | `gem_mapper_open_input` | 79 |  |
|  | `gem_mapper_close_input` | 96 |  |
|  | `gem_mapper_open_output` | 104 |  |
|  | `gem_mapper_close_output` | 131 |  |
|  | `gem_mapper_print_profile` | 134 |  |
|  | `main` | 182 |  |

### `src/tools/interface/mapper_arguments.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `parse_integer_system` | 136 |  |
|  | `gem_mapper_print_usage` | 158 |  |
|  | `gem_mapper_parameters_presets` | 165 |  |
|  | `gem_mapper_parameters_check` | 181 |  |
|  | `gem_mapper_parse_arguments_io` | 239 |  |
|  | `gem_mapper_parse_arguments_single_end` | 370 |  |
|  | `gem_mapper_parse_arguments_paired_end` | 590 |  |
|  | `gem_mapper_parse_arguments_bisulfite` | 738 |  |
|  | `gem_mapper_parse_arguments_alignment_score` | 784 |  |
|  | `gem_mapper_parse_arguments_mapq_score` | 855 |  |
|  | `gem_mapper_parse_arguments_reporting` | 882 |  |
|  | `gem_mapper_parse_arguments_output_format` | 904 |  |
|  | `gem_mapper_parse_arguments_system` | 947 |  |
|  | `gem_mapper_parse_arguments_gpu` | 971 |  |
|  | `gem_mapper_parse_arguments_misc` | 1053 |  |
|  | `gem_mapper_parse_arguments` | 1139 |  |

### `src/utils/ihash.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `ihash_new` | 44 |  |
|  | `ihash_clear` | 50 |  |
|  | `ihash_delete` | 61 |  |
|  | `ihash_get_ihash_element` | 68 |  |
|  | `ihash_insert_element` | 73 |  |
|  | `ihash_get_element` | 91 |  |
|  | `ihash_cmp_keys` | 115 |  |
|  | `ihash_sort_by_key` | 124 |  |

### `src/utils/options_menu.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `options_get_num_options` | 26 |  |
|  | `options_adaptor_getopt` | 31 |  |
|  | `options_adaptor_getopt_short` | 48 |  |
|  | `options_fprint_menu` | 69 |  |

### `src/utils/string_buffer.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `string_init` | 35 |  |
|  | `string_resize` | 50 |  |
|  | `string_clear` | 71 |  |
|  | `string_destroy` | 75 |  |
|  | `string_set_buffer` | 85 |  |
|  | `string_get_buffer` | 93 |  |
|  | `string_set_length` | 96 |  |
|  | `string_get_length` | 101 |  |
|  | `string_append_char` | 113 |  |
|  | `string_append_eos` | 120 |  |
|  | `string_right_append_buffer` | 166 |  |
|  | `string_right_append_string` | 176 |  |
|  | `string_is_null` | 230 |  |
|  | `string_copy` | 284 |  |
|  | `gem_strncpy` | 340 |  |
|  | `gem_strcmp` | 352 |  |
|  | `gem_streq` | 358 |  |
|  | `gem_strcaseeq` | 361 |  |
|  | `gem_strcat` | 373 |  |
|  | `gem_strlen` | 381 |  |

### `src/utils/vector.c`

| Done | Function | Line | Test |
|---|---|---|---|
|  | `vector_new_` | 36 |  |
|  | `vector_reserve` | 46 |  |
|  | `vector_delete` | 78 |  |

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
