/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_verify_candidates.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "approximate_search/approximate_search_verify_candidates.h"
#include "approximate_search/approximate_search_control.h"
#include "filtering/filtering_candidates_process.h"
#include "filtering/filtering_candidates_verify.h"
#include "filtering/filtering_candidates_verify_buffered.h"
#include "filtering/filtering_candidates_extend.h"
#include "filtering/filtering_candidates_align.h"

/*
 * Debug
 */
#define DEBUG_REGION_SCHEDULE_PRINT GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Benchmark
 */
#ifdef CUDA_BENCHMARK_GENERATE_VERIFY_CANDIDATES
FILE* benchmark_verify_candidates = NULL;
#endif

/*
 * Verify Candidates
 */
void approximate_search_verify_candidates(
    approximate_search_t* const search,
    matches_t* const matches) {
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,"[GEM]>ASM::Verify Candidates\n");
    tab_global_inc();
  }
  // Parameters
  pattern_t* const pattern = &search->pattern;
  filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
  // Verify Candidates
  filtering_candidates_verify_candidates(filtering_candidates,pattern);
  filtering_candidates_align_candidates(filtering_candidates,
      pattern,search->emulated_rc_search,false,false,matches);
  search->processing_state = asearch_processing_state_candidates_verified;
  // Adjust max-differences
  asearch_control_adjust_max_differences_using_strata(search,matches);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); }
}
/*
 * Verify Extend Candidate
 */
uint64_t approximate_search_verify_extend_candidate(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const candidate_pattern,
    const match_trace_t* const extended_match,
    mapper_stats_t* const mapper_stats,
    paired_matches_t* const paired_matches,
    const sequence_end_t candidate_end) {
  return filtering_candidates_extend_match(filtering_candidates,
      candidate_pattern,extended_match,paired_matches,candidate_end,mapper_stats);
}
/*
 * Verify Candidates Buffered
 */
void approximate_search_verify_candidates_buffered_print_benchmark(approximate_search_t* const search);
void approximate_search_verify_candidates_buffered_copy(
    approximate_search_t* const search,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  // Add to GPU-Buffer
  filtering_candidates_verify_buffered_add(
      search->filtering_candidates,&search->pattern,
      gpu_buffer_align_bpm,&search->gpu_buffer_align_offset,
      &search->gpu_filtering_regions,&search->gpu_num_filtering_regions);
  // BENCHMARK
#ifdef CUDA_BENCHMARK_GENERATE_VERIFY_CANDIDATES
  approximate_search_verify_candidates_buffered_print_benchmark(search);
#endif
}
void approximate_search_verify_candidates_buffered_retrieve(
    approximate_search_t* const search,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    matches_t* const matches) {
  // Retrieve
  filtering_candidates_verify_buffered_retrieve(
      search->filtering_candidates,&search->pattern,
      gpu_buffer_align_bpm,search->gpu_buffer_align_offset,
      search->gpu_filtering_regions,search->gpu_num_filtering_regions);
  // Align
  filtering_candidates_align_candidates(search->filtering_candidates,
      &search->pattern,search->emulated_rc_search,false,false,matches);
  // Update state
  search->processing_state = asearch_processing_state_candidates_verified;
}
/*
 * Display/Benchmark
 */
void approximate_search_verify_candidates_buffered_print_benchmark(approximate_search_t* const search) {
#ifdef CUDA_BENCHMARK_GENERATE_VERIFY_CANDIDATES
  // Prepare benchmark file
  if (benchmark_verify_candidates==NULL) {
    benchmark_verify_candidates = fopen("gem3.verify.candidates.benchmark","w+");
  }
  // Print all candidates' tiles
  filtering_candidates_verify_buffered_print_benchmark(
      benchmark_verify_candidates,search->filtering_candidates,&search->pattern);
#endif
}
