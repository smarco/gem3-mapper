/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_verify_candidates.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "approximate_search_verify_candidates.h"
#include "approximate_search_control.h"
#include "filtering_candidates_process.h"
#include "filtering_candidates_verify.h"
#include "filtering_candidates_verify_buffered.h"
#include "filtering_candidates_align.h"

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
GEM_INLINE void approximate_search_verify_candidates(approximate_search_t* const search,matches_t* const matches) {
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,"[GEM]>ASM::Verify Candidates\n");
    tab_global_inc();
  }
  // Parameters
  const as_parameters_t* const actual_parameters = search->as_parameters;
  const search_parameters_t* const parameters = actual_parameters->search_parameters;
  pattern_t* const pattern = &search->pattern;
  filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
  // Verify Candidates
  filtering_candidates_verify_candidates(filtering_candidates,search->archive,
      search->text_collection,pattern,actual_parameters,matches,search->mm_stack);
  filtering_candidates_align_candidates(
      filtering_candidates,search->archive->text,search->archive->locator,search->text_collection,
      pattern,search->emulated_rc_search,actual_parameters,false,matches,search->mm_stack);
  search->processing_state = asearch_processing_state_candidates_verified;
  // Adjust max-differences
  asearch_control_adjust_max_differences_using_strata(search,matches);
  // Update MCS (maximum complete stratum)
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  search->max_matches_reached = num_matches >= parameters->search_max_matches;
  if (search->max_matches_reached) approximate_search_update_mcs(search,0);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); }
}
/*
 * Verify Candidates Buffered
 */
GEM_INLINE void approximate_search_verify_candidates_buffered_copy(
    approximate_search_t* const search,gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  // Add to GPU-Buffer
  search->gpu_buffer_align_offset = gpu_buffer_align_bpm_get_num_candidates(gpu_buffer_align_bpm);
  search->gpu_buffer_align_total = filtering_candidates_verify_buffered_add(
      search->filtering_candidates,&search->pattern,gpu_buffer_align_bpm);
  // BENCHMARK
#ifdef CUDA_BENCHMARK_GENERATE_VERIFY_CANDIDATES
  // approximate_search_region_profile_buffered_print_benchmark(search);
#endif
}
GEM_INLINE void approximate_search_verify_candidates_buffered_retrieve(
    approximate_search_t* const search,gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    matches_t* const matches) {
  // Buffer Offsets
  const uint64_t candidate_offset_begin = search->gpu_buffer_align_offset;
  const uint64_t candidate_offset_end = search->gpu_buffer_align_offset+search->gpu_buffer_align_total;
  // Retrieve
  const uint64_t num_accepted_regions = filtering_candidates_verify_buffered_retrieve(
      search->filtering_candidates,search->archive->text,search->text_collection,&search->pattern,
      gpu_buffer_align_bpm,candidate_offset_begin,candidate_offset_end,search->mm_stack);
  if (num_accepted_regions > 0) {
    // Realign
    PROFILE_START(GP_FC_REALIGN_BPM_BUFFER_CANDIDATE_REGIONS,PROFILE_LEVEL);
    filtering_candidates_align_candidates(
        search->filtering_candidates,search->archive->text,search->archive->locator,
        search->text_collection,&search->pattern,search->emulated_rc_search,
        search->as_parameters,true,matches,search->mm_stack);
    PROFILE_STOP(GP_FC_REALIGN_BPM_BUFFER_CANDIDATE_REGIONS,PROFILE_LEVEL);
  }
  // Update state
  search->processing_state = asearch_processing_state_candidates_verified;
}
