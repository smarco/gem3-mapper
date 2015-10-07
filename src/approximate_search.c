/*
 * PROJECT: GEMMapper
 * FILE: approximate_search.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "approximate_search.h"
#include "approximate_search_filtering_adaptive.h"
#include "approximate_search_filtering_complete.h"
#include "approximate_search_neighborhood.h"

/*
 * Approximate Search State
 */
const char* approximate_search_state_label[] =
{
    [0]  = "begin",
    [1]  = "no_regions",
    [2]  = "exact_matches",
    [3]  = "exact_filtering_adaptive",
    [4]  = "verify_candidates",
    [5]  = "candidates_verified",
    [6]  = "exact_filtering_boost",
    [7]  = "inexact_filtering",
    [8]  = "neighborhood",
    [9]  = "end",
    [10] = "read_recovery",
    [11] = "unbounded_alignment",
    [12] = "probe_candidates",
};

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
  search->stop_before_neighborhood_search = false;
  const uint64_t max_complete_error = search->as_parameters->complete_search_error_nominal; // FIXME + search->pattern.num_low_quality_bases;
  search->max_complete_error = MIN(max_complete_error,search->pattern.max_effective_filtering_error);
  search->max_complete_stratum = ALL;
  search->max_matches_reached = false;
  // Prepare filtering candidates
  filtering_candidates_clear(search->filtering_candidates);
  // Prepare region profile
  if (search->max_complete_error > 0) {
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
 * Modifiers
 */
GEM_INLINE void approximate_search_hold_verification_candidates(approximate_search_t* const search) {
  filtering_candidates_set_all_regions_pending(search->filtering_candidates);
}
GEM_INLINE void approximate_search_release_verification_candidates(approximate_search_t* const search) {
  filtering_candidates_set_all_regions_unverified(search->filtering_candidates);
  if (search->search_state==asearch_candidates_verified) {
    search->search_state = asearch_verify_candidates;
  }
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
    case mapping_adaptive_filtering_thorough:
      approximate_search_filtering_adaptive(search,matches); // Adaptive & incremental mapping
      break;
    case mapping_adaptive_filtering_complete:
    case mapping_fixed_filtering_complete:
      approximate_search_filtering_complete(search,matches); // Filtering Complete
      break;
    case mapping_neighborhood_search:
      approximate_search_neighborhood_search(search,matches); // Brute-force mapping
      break;
    case mapping_region_profile_fixed:
      approximate_search_filtering_adaptive_generate_regions(search);
      matches->max_complete_stratum = 0; PROF_STOP(GP_AS_MAIN); return;
      break;
    case mapping_test:
      approximate_search_test(search,matches);
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
 * Display
 */
GEM_INLINE void approximate_search_print(FILE* const stream,approximate_search_t* const search) {
  tab_fprintf(stream,"[GEM]>ApproximateSearch\n");
  tab_global_inc();
  tab_fprintf(stream,"=> Search.State %s\n",approximate_search_state_label[search->search_state]);
  tab_fprintf(stream,"=> Max.complete.error %lu\n",search->max_complete_error);
  tab_fprintf(stream,"=> MCS %lu\n",search->max_complete_stratum);
  tab_fprintf(stream,"=> Max.matches.reached %lu\n",search->max_matches_reached);
  tab_fprintf(stream,"=> Exact.lo %lu\n",search->lo_exact_matches);
  tab_fprintf(stream,"=> Exact.hi %lu\n",search->hi_exact_matches);
  tab_fprintf(stream,"=> Region.Profile\n");
  tab_global_inc();
  region_profile_print(stream,&search->region_profile,false);
  tab_global_dec();
  tab_global_dec();
}
