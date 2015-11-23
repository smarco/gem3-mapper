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
#include "filtering_candidates.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Approximate Search State/Stage
 */
const char* asearch_processing_state_label[] =
{
    // Begin
    [asearch_processing_state_begin]  = "begin",
    // Region Profile
    [asearch_processing_state_region_partitioned]  = "region-partitioned",
    [asearch_processing_state_region_profiled]  = "region-profiled",
    [asearch_processing_state_no_regions]  = "no-regions",
    [asearch_processing_state_exact_matches]  = "exact_matches",
    // Verify Candidates
    [asearch_processing_state_candidates_processed]  = "candidates-processed",
    [asearch_processing_state_candidates_verified]  = "candidates-verified",
};
const char* asearch_stage_label[] =
{
    [asearch_stage_begin]  = "begin",
    [asearch_stage_read_recovery]  = "read-recovery",
    [asearch_stage_filtering_adaptive]  = "filtering-adaptive",
    [asearch_stage_filtering_boost]  = "filtering-boost",
    [asearch_stage_inexact_filtering]  = "inexact-filtering",
    [asearch_stage_unbounded_alignment]  = "unbounded-alignment",
    [asearch_stage_neighborhood]  = "neighborhood",
    [asearch_stage_end]  = "end",
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
GEM_INLINE void approximate_search_reset(approximate_search_t* const search) {
  // Reset Approximate Search State
  search->search_stage = asearch_stage_begin;
  search->processing_state = asearch_processing_state_begin;
  search->stop_before_neighborhood_search = false;
  const uint64_t max_complete_error = search->as_parameters->complete_search_error_nominal; // FIXME + search->pattern.num_low_quality_bases;
  search->max_complete_error = MIN(max_complete_error,search->pattern.max_effective_filtering_error);
  search->max_complete_stratum = ALL;
  search->max_matches_reached = false;
  // Prepare region profile
  if (search->max_complete_error > 0) {
    region_profile_new(&search->region_profile,search->pattern.key_length,search->mm_stack);
  }
}
GEM_INLINE void approximate_search_destroy(approximate_search_t* const search) {
  /* NOP */
}
/*
 * Memory Injection (Support Data Structures)
 */
GEM_INLINE void approximate_search_inject_mm_stack(
    approximate_search_t* const search,mm_stack_t* const mm_stack) {
  search->mm_stack = mm_stack;                         // Set MM
}
GEM_INLINE void approximate_search_inject_interval_set(
    approximate_search_t* const search,interval_set_t* const interval_set) {
  search->interval_set = interval_set;                 // Interval Set
}
GEM_INLINE void approximate_search_inject_text_collection(
    approximate_search_t* const search,text_collection_t* const text_collection) {
  search->text_collection = text_collection;           // Text-Collection
}
GEM_INLINE void approximate_search_inject_filtering_candidates(
    approximate_search_t* const search,filtering_candidates_t* const filtering_candidates) {
  search->filtering_candidates = filtering_candidates; // Filtering Candidates
}
/*
 * Accessors
 */
GEM_INLINE uint64_t approximate_search_get_num_exact_filtering_candidates(const approximate_search_t* const search) {
  return (search->processing_state == asearch_processing_state_exact_matches) ?
      search->hi_exact_matches - search->lo_exact_matches : 0;
}
GEM_INLINE void approximate_search_update_mcs(approximate_search_t* const search,const uint64_t max_complete_stratum) {
  search->max_complete_stratum = max_complete_stratum;
}
GEM_INLINE uint64_t approximate_search_get_num_regions_profile(const approximate_search_t* const search) {
  const region_profile_t* const region_profile = &search->region_profile;
  return region_profile->num_filtering_regions;
}
GEM_INLINE uint64_t approximate_search_get_num_decode_candidates(const approximate_search_t* const search) {
  const filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
  return filtering_candidates_get_num_candidate_positions(filtering_candidates);
}
GEM_INLINE uint64_t approximate_search_get_num_verify_candidates(const approximate_search_t* const search) {
  if (search->processing_state == asearch_processing_state_exact_matches) {
    return search->hi_exact_matches - search->lo_exact_matches;
  } else {
    const filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
    return filtering_candidates_get_num_candidate_regions(filtering_candidates);
  }
}
/*
 * Modifiers
 */
GEM_INLINE void approximate_search_hold_verification_candidates(approximate_search_t* const search) {
  filtering_candidates_set_all_regions_pending(search->filtering_candidates);
}
GEM_INLINE void approximate_search_release_verification_candidates(approximate_search_t* const search) {
  filtering_candidates_set_all_regions_unverified(search->filtering_candidates);
  if (search->processing_state==asearch_processing_state_candidates_verified) {
    search->processing_state = asearch_processing_state_candidates_processed;
  }
}
/*
 * Approximate String Matching using the FM-index
 */
GEM_INLINE void approximate_search(approximate_search_t* const search,matches_t* const matches) {
  PROFILE_START(GP_AS_MAIN,PROFILE_LEVEL);
  /*
   * Select mapping strategy
   */
  const search_parameters_t* const parameters = search->as_parameters->search_parameters;
  switch (parameters->mapping_mode) {
    case mapping_adaptive_filtering_fast:
    case mapping_adaptive_filtering_thorough:
      approximate_search_filtering_adaptive(search,matches); // Adaptive mapping
      break;
    case mapping_adaptive_filtering_complete:
    case mapping_fixed_filtering_complete:
      approximate_search_filtering_complete(search,matches); // Filtering complete mapping
      break;
    case mapping_neighborhood_search:
      approximate_search_neighborhood_search(search,matches); // Brute-force mapping
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  PROFILE_STOP(GP_AS_MAIN,PROFILE_LEVEL);
}
/*
 * Display
 */
GEM_INLINE void approximate_search_print(FILE* const stream,approximate_search_t* const search) {
  tab_fprintf(stream,"[GEM]>ApproximateSearch\n");
  tab_global_inc();
  tab_fprintf(stream,"=> Search.Stage %s\n",asearch_stage_label[search->search_stage]);
  tab_fprintf(stream,"  => Search.State %s\n",asearch_processing_state_label[search->processing_state]);
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
