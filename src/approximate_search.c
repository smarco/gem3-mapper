/*
 * PROJECT: GEMMapper
 * FILE: approximate_search.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "approximate_search.h"
#include "approximate_search_filtering.h"
#include "approximate_search_neighborhood.h"

/*
 * Setup
 */
GEM_INLINE void approximate_search_init(
    approximate_search_t* const search,archive_t* const archive,
    as_parameters_t* const as_parameters,const bool emulated_rc_search) {
  // Index Structures & Parameters
  search->archive = archive;
  search->as_parameters = as_parameters;
  search->emulated_rc_search = emulated_rc_search;
}
GEM_INLINE void approximate_search_configure(
    approximate_search_t* const search,filtering_candidates_t* const filtering_candidates,
    text_collection_t* text_collection,interval_set_t* const interval_set,mm_stack_t* const mm_stack) {
  // Set Auxiliary Structures (external)
  search->filtering_candidates = filtering_candidates; // Filtering Candidates
  search->text_collection = text_collection;           // Text-Collection
  search->interval_set = interval_set;                 // Interval Set
  search->mm_stack = mm_stack;                         // Set MM
}
GEM_INLINE void approximate_search_reset(approximate_search_t* const search) {
  // Reset Approximate Search State
  search->search_state = asearch_begin;
  search->verify_candidates = true;
  search->stop_before_neighborhood_search = false;
  const uint64_t max_search_error = search->as_parameters->max_search_error_nominal;
  const uint64_t max_effective_filtering_error = search->pattern.max_effective_filtering_error;
  search->max_differences = MIN(max_search_error,max_effective_filtering_error);
  search->max_complete_stratum = ALL;
  search->max_matches_reached = false;
  // Prepare filtering candidates
  filtering_candidates_clear(search->filtering_candidates);
  // Prepare region profile
  if (search->max_differences > 0) {
    region_profile_new(&search->region_profile,search->pattern.key_length,search->mm_stack);
  }
}
GEM_INLINE void approximate_search_destroy(approximate_search_t* const search) { /* NOP */ }
/*
 * Accessors
 */
GEM_INLINE uint64_t approximate_search_get_num_filtering_candidates(const approximate_search_t* const search) {
  if (search->search_state == asearch_exact_matches) {
    return search->hi_exact_matches - search->lo_exact_matches;
  } else {
    const filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
    return filtering_candidates_get_num_candidate_regions(filtering_candidates);
  }
}
GEM_INLINE uint64_t approximate_search_get_num_exact_filtering_candidates(const approximate_search_t* const search) {
  return (search->search_state == asearch_exact_matches) ? search->hi_exact_matches - search->lo_exact_matches : 0;
}
GEM_INLINE void approximate_search_update_mcs(approximate_search_t* const search,const uint64_t max_complete_stratum) {
  search->max_complete_stratum = max_complete_stratum;
}
/*
 * Approximate String Matching using the FM-index
 */
GEM_INLINE void approximate_search(approximate_search_t* const search,matches_t* const matches) {
  PROF_START(GP_AS_MAIN);
  /*
   * Select mapping strategy
   */
  const search_parameters_t* const parameters = search->as_parameters->search_parameters;
  switch (parameters->mapping_mode) {
    case mapping_adaptive_filtering_fast:
    case mapping_adaptive_filtering_match:
    case mapping_adaptive_filtering_complete:
      approximate_search_filtering_adaptive(search,matches); // Adaptive & incremental mapping
      break;
    case mapping_fixed_filtering_complete:
      approximate_search_filtering_complete(search,matches); // Filtering Complete
      break;
    case mapping_neighborhood_search:
      approximate_search_neighborhood_search(search,matches); // Brute-force mapping
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  // Set matches-MCS
  if (matches) matches->max_complete_stratum = MIN(matches->max_complete_stratum,search->max_complete_stratum);
  PROF_STOP(GP_AS_MAIN);
}
/*
 * Approximate String Matching using the FM-index (Verification)
 */
GEM_INLINE void approximate_search_verify(approximate_search_t* const search,matches_t* const matches) {
  // Verify
  const uint64_t num_accepted_regions = filtering_candidates_verify_candidates(
      search->filtering_candidates,search->archive,search->text_collection,
      &search->pattern,search->as_parameters,matches,search->mm_stack);
  if (num_accepted_regions > 0) {
    // Realign
    filtering_candidates_align_candidates(search->filtering_candidates,search->archive->text,
        search->text_collection,&search->pattern,search->emulated_rc_search,
        search->as_parameters,true,matches,search->mm_stack);
  }
  // Update state
  if (search->search_state==asearch_verify_candidates) {
    search->search_state = asearch_candidates_verified;
  }
}
GEM_INLINE void approximate_search_verify_using_bpm_buffer(
    approximate_search_t* const search,
    matches_t* const matches,bpm_gpu_buffer_t* const bpm_gpu_buffer,
    const uint64_t candidate_offset_begin,const uint64_t candidate_offset_end) {
  // Retrieve
  const uint64_t num_accepted_regions = filtering_candidates_bpm_buffer_retrieve(
      search->filtering_candidates,search->archive->text,search->text_collection,
      &search->pattern,bpm_gpu_buffer,candidate_offset_begin,candidate_offset_end,search->mm_stack);
  if (num_accepted_regions > 0) {
    // Realign
    PROF_START(GP_FC_REALIGN_BPM_BUFFER_CANDIDATE_REGIONS);
    filtering_candidates_align_candidates(search->filtering_candidates,search->archive->text,
        search->text_collection,&search->pattern,search->emulated_rc_search,
        search->as_parameters,true,matches,search->mm_stack);
    PROF_STOP(GP_FC_REALIGN_BPM_BUFFER_CANDIDATE_REGIONS);
  }
  // Update state
  if (search->search_state==asearch_verify_candidates) {
    search->search_state = asearch_candidates_verified;
  }
}
GEM_INLINE void approximate_search_hold_verification_candidates(approximate_search_t* const search) {
  filtering_candidates_set_all_regions_pending(search->filtering_candidates);
}
GEM_INLINE void approximate_search_release_verification_candidates(approximate_search_t* const search) {
  filtering_candidates_set_all_regions_unverified(search->filtering_candidates);
  if (search->search_state==asearch_candidates_verified) {
    search->search_state = asearch_verify_candidates;
  }
}
