/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   Approximate-String-Matching (ASM) main module.
 *   Dispatch the search depending on the search-approach selected
 *   and provides data structures for the search
 */

#include "approximate_search/approximate_search.h"
#include "approximate_search/approximate_search_filtering_adaptive.h"
#include "approximate_search/approximate_search_neighborhood.h"
#include "approximate_search/approximate_search_hybrid.h"
#include "filtering/candidates/filtering_candidates.h"

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
    [asearch_stage_filtering_adaptive]  = "filtering-adaptive",
    [asearch_stage_filtering_adaptive_finished]  = "filtering-adaptive-finished",
    [asearch_stage_neighborhood]  = "neighborhood",
    [asearch_stage_local_alignment]  = "local-alignment",
    [asearch_stage_end]  = "end",
};

/*
 * Setup
 */
void approximate_search_init(
    approximate_search_t* const search,
    archive_t* const archive,
    search_parameters_t* const search_parameters,
    filtering_candidates_t* const filtering_candidates,
    nsearch_schedule_t* const nsearch_schedule,
    mm_allocator_t* const mm_allocator) {
  // Archive & Search Parameters
  search->archive = archive;
  search->search_parameters = search_parameters;
  // Filtering Candidates
  search->filtering_candidates = filtering_candidates;
  if (filtering_candidates!=NULL) {
    filtering_candidates_inject_handlers(
        filtering_candidates,archive,search_parameters,mm_allocator);
  }
  // Filtering Candidates Buffered
  filtering_candidates_buffered_inject_handlers(
      &search->filtering_candidates_buffered,mm_allocator);
  // Region Profile
  region_profile_inject_mm(&search->region_profile,mm_allocator);
  search->num_limited_exact_matches = 0;
  // Nsearch
  search->nsearch_schedule = nsearch_schedule;
  nsearch_schedule_inject_mm(search->nsearch_schedule,mm_allocator);
  // MM
  search->mm_allocator = mm_allocator;
}
void approximate_search_destroy(approximate_search_t* const search) {
  // Pattern
  pattern_destroy(&search->pattern,search->mm_allocator);
  // Region Profile
  region_profile_destroy(&search->region_profile);
}
/*
 * Prepare
 */
void approximate_search_prepare(
    approximate_search_t* const search,
    const bool run_length_pattern,
    sequence_t* const sequence) {
  // Pattern
  pattern_init(
      &search->pattern,sequence,&search->do_quality_search,
      search->search_parameters,run_length_pattern,search->bisulfite_conversion,search->mm_allocator);
  // Region profile
  const uint64_t key_length = search->pattern.key_length;
  region_profile_init(&search->region_profile,key_length);
  // Reset Approximate Search State
  search->search_stage = asearch_stage_begin;
  search->processing_state = asearch_processing_state_begin;
  const uint64_t max_complete_error = search->search_parameters->complete_search_error_nominal;
  search->max_search_error = MIN(max_complete_error,search->pattern.max_effective_filtering_error);
}
/*
 * Accessors
 */
uint64_t approximate_search_get_num_regions_profile(const approximate_search_t* const search) {
  const region_profile_t* const region_profile = &search->region_profile;
  return region_profile->num_filtering_regions;
}
uint64_t approximate_search_get_num_canonical_regions_profile(const approximate_search_t* const search) {
  const region_profile_t* const region_profile = &search->region_profile;
  return region_profile->num_filtering_regions;
}
uint64_t approximate_search_get_num_decode_candidates(const approximate_search_t* const search) {
  const region_profile_t* const region_profile = &search->region_profile;
  return region_profile->total_candidates;
}
uint64_t approximate_search_get_num_filtering_candidates(const approximate_search_t* const search) {
  return filtering_candidates_get_num_regions(search->filtering_candidates);
}
uint64_t approximate_search_get_num_filtering_candidates_buffered(
    const approximate_search_t* const search) {
  return filtering_candidates_buffered_get_num_regions(&search->filtering_candidates_buffered);
}
uint64_t approximate_search_get_num_filtering_canonical_candidates_buffered(
    const approximate_search_t* const search) {
  return filtering_candidates_buffered_get_num_canonical_regions(&search->filtering_candidates_buffered);
}
uint64_t approximate_search_get_num_filtering_candidate_buffered_tiles_length(
    const approximate_search_t* const search) {
  return filtering_candidates_buffered_get_region_tiles_length(&search->filtering_candidates_buffered);
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
      approximate_search_filtering_adaptive(search,matches); // Adaptive mapping
      break;
    case mapping_neighborhood_search_brute_force:
      approximate_search_neighborhood_search_brute_force(search,matches); // Brute-force mapping
      break;
    case mapping_neighborhood_search_partition:
      approximate_search_neighborhood_search_partition(search,matches); // NS-partition mapping
      break;
    case mapping_hybrid_sensitive:
    case mapping_hybrid_complete:
      approximate_search_hybrid(search,matches);
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
  tab_fprintf(stream,"=> Max.complete.error %lu\n",search->max_search_error);
  tab_fprintf(stream,"=> Region.Profile\n");
  tab_global_inc();
  region_profile_print(stream,&search->region_profile,false);
  tab_global_dec();
  tab_global_dec();
}
