/*
 * PROJECT: GEMMapper
 * FILE: approximate_search.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "approximate_search/approximate_search.h"
#include "approximate_search/approximate_search_filtering_adaptive.h"
#include "approximate_search/approximate_search_filtering_complete.h"
#include "approximate_search/approximate_search_neighborhood.h"
#include "filtering/filtering_candidates.h"

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
    // Verify Candidates
    [asearch_processing_state_candidates_processed]  = "candidates-processed",
    [asearch_processing_state_candidates_verified]  = "candidates-verified",
};
const char* asearch_stage_label[] =
{
    [asearch_stage_begin]  = "begin",
    [asearch_stage_read_recovery]  = "read-recovery",
    [asearch_stage_filtering_adaptive]  = "filtering-adaptive",
    [asearch_stage_inexact_filtering]  = "inexact-filtering",
    [asearch_stage_local_alignment]  = "local-alignment",
    [asearch_stage_neighborhood]  = "neighborhood",
    [asearch_stage_end]  = "end",
};
/*
 * Setup
 */
void approximate_search_init(
    approximate_search_t* const search,
    archive_t* const archive,
    search_parameters_t* const search_parameters,
    const bool emulated_rc_search) {
  // Index Structures & Parameters
  search->archive = archive;
  search->search_parameters = search_parameters;
  search->emulated_rc_search = emulated_rc_search;
}
void approximate_search_reset(approximate_search_t* const search) {
  // Reset Approximate Search State
  search->search_stage = asearch_stage_begin;
  search->processing_state = asearch_processing_state_begin;
  search->stop_before_neighborhood_search = false;
  const uint64_t max_complete_error = search->search_parameters->complete_search_error_nominal;
  search->max_complete_error = MIN(max_complete_error,search->pattern.max_effective_filtering_error);
  search->max_complete_stratum = ALL;
  // Prepare region profile
  const uint64_t key_length = search->pattern.key_length;
  region_profile_new(&search->region_profile,key_length,search->mm_stack);
  // Reset metrics
  const double proper_length = fm_index_get_proper_length(search->archive->fm_index);
  const search_parameters_t* const search_parameters = search->search_parameters;
  const int32_t swg_match_score = search_parameters->swg_penalties.generic_match_score;
  approximate_search_metrics_init(&search->metrics,proper_length,key_length,swg_match_score);
}
void approximate_search_destroy(approximate_search_t* const search) {
  /* NOP */
}
/*
 * Memory Injection (Support Data Structures)
 */
void approximate_search_inject_mm_stack(
    approximate_search_t* const search,
    mm_stack_t* const mm_stack) {
  search->mm_stack = mm_stack;
}
void approximate_search_inject_interval_set(
    approximate_search_t* const search,
    interval_set_t* const interval_set) {
  search->interval_set = interval_set;
}
void approximate_search_inject_text_collection(
    approximate_search_t* const search,
    text_collection_t* const text_collection) {
  search->text_collection = text_collection;
}
void approximate_search_inject_filtering_candidates(
    approximate_search_t* const search,
    filtering_candidates_t* const filtering_candidates,
    text_collection_t* const text_collection,
    mm_stack_t* const mm_stack) {
  search->filtering_candidates = filtering_candidates;
  filtering_candidates_inject_search(
      search->filtering_candidates,search->archive,search->search_parameters);
  filtering_candidates_inject_mm_stack(search->filtering_candidates,mm_stack);
  filtering_candidates_inject_text_collection(search->filtering_candidates,text_collection);
}
/*
 * Accessors
 */
void approximate_search_update_mcs(
    approximate_search_t* const search,
    const uint64_t max_complete_stratum) {
  search->max_complete_stratum = max_complete_stratum;
}
uint64_t approximate_search_get_num_regions_profile(const approximate_search_t* const search) {
  const region_profile_t* const region_profile = &search->region_profile;
  return region_profile->num_filtering_regions;
}
uint64_t approximate_search_get_num_decode_candidates(const approximate_search_t* const search) {
  const region_profile_t* const region_profile = &search->region_profile;
  return region_profile->total_candidates;
}
uint64_t approximate_search_get_num_verify_candidates(const approximate_search_t* const search) {
  return filtering_candidates_get_num_candidate_regions(search->filtering_candidates);
}
/*
 * Approximate String Matching using the FM-index
 */
void approximate_search(approximate_search_t* const search,matches_t* const matches) {
  PROFILE_START(GP_AS_MAIN,PROFILE_LEVEL);
  /*
   * Select mapping strategy
   */
  switch (search->search_parameters->mapping_mode) {
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
void approximate_search_print(FILE* const stream,approximate_search_t* const search) {
  tab_fprintf(stream,"[GEM]>ApproximateSearch\n");
  tab_global_inc();
  tab_fprintf(stream,"=> Search.Stage %s\n",asearch_stage_label[search->search_stage]);
  tab_fprintf(stream,"  => Search.State %s\n",asearch_processing_state_label[search->processing_state]);
  tab_fprintf(stream,"=> Max.complete.error %lu\n",search->max_complete_error);
  tab_fprintf(stream,"=> MCS %lu\n",search->max_complete_stratum);
  tab_fprintf(stream,"=> Region.Profile\n");
  tab_global_inc();
  region_profile_print(stream,&search->region_profile,false);
  tab_global_dec();
  tab_global_dec();
}
