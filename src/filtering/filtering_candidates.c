/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering/filtering_candidates.h"
#include "filtering/filtering_region.h"
#include "filtering/filtering_region_verify.h"
#include "filtering/filtering_region_align.h"
#include "align/align.h"
#include "matches/matches_classify.h"

/*
 * Debug
 */
#define DEBUG_FILTERING_CANDIDATES  GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Constants
 */
#define REGIONS_BUFFER_INIT                      100
#define CANDIDATE_POSITIONS_INIT                 1000

/*
 * Setup
 */
void filtering_candidates_init(filtering_candidates_t* const filtering_candidates) {
  // Candidates
  filtering_candidates->filtering_positions = vector_new(CANDIDATE_POSITIONS_INIT,filtering_position_t);
  filtering_candidates->filtering_regions = vector_new(CANDIDATE_POSITIONS_INIT,filtering_region_t);
  filtering_candidates->discarded_regions = vector_new(CANDIDATE_POSITIONS_INIT,filtering_region_t);
  filtering_candidates->verified_regions = vector_new(CANDIDATE_POSITIONS_INIT,verified_region_t);
  // Cache
  filtering_region_cache_init(&filtering_candidates->filtering_region_cache);
}
void filtering_candidates_clear(filtering_candidates_t* const filtering_candidates) {
  // Candidates
  vector_clear(filtering_candidates->filtering_positions);
  vector_clear(filtering_candidates->filtering_regions);
  vector_clear(filtering_candidates->discarded_regions);
  vector_clear(filtering_candidates->verified_regions);
}
void filtering_candidates_destroy(filtering_candidates_t* const filtering_candidates) {
  // Candidates
  vector_delete(filtering_candidates->filtering_positions);
  vector_delete(filtering_candidates->filtering_regions);
  vector_delete(filtering_candidates->discarded_regions);
  vector_delete(filtering_candidates->verified_regions);
  // Cache
  filtering_region_cache_destroy(&filtering_candidates->filtering_region_cache);
}
/*
 * Memory Injection (Support Data Structures)
 */
void filtering_candidates_inject_search(
    filtering_candidates_t* const filtering_candidates,
    archive_t* const archive,
    search_parameters_t* const search_parameters) {
  filtering_candidates->archive = archive;
  filtering_candidates->search_parameters = search_parameters;
}
void filtering_candidates_inject_mm_stack(
    filtering_candidates_t* const filtering_candidates,
    mm_stack_t* const mm_stack) {
  filtering_candidates->mm_stack = mm_stack;
}
void filtering_candidates_inject_text_collection(
    filtering_candidates_t* const filtering_candidates,
    text_collection_t* const text_collection) {
  filtering_candidates->text_collection = text_collection;
}
/*
 * Accessors
 */
uint64_t filtering_candidates_get_num_candidate_positions(
    const filtering_candidates_t* const filtering_candidates) {
  return vector_get_used(filtering_candidates->filtering_positions);
}
uint64_t filtering_candidates_get_num_candidate_regions(
    const filtering_candidates_t* const filtering_candidates) {
  return vector_get_used(filtering_candidates->filtering_regions);
}
uint64_t filtering_candidates_count_candidate_regions(
    filtering_candidates_t* const filtering_candidates_end,
    const filtering_region_status_t filtering_region_status) {
  uint64_t count = 0;
  VECTOR_ITERATE(filtering_candidates_end->filtering_regions,filtering_region,n,filtering_region_t) {
    if (filtering_region->status == filtering_region_status) ++count;
  }
  return count;
}
/*
 * Adding candidate positions
 */
void filtering_candidates_add_read_interval(
    filtering_candidates_t* const filtering_candidates,
    search_parameters_t* const search_parameters,
    const uint64_t interval_lo,
    const uint64_t interval_hi,
    const uint64_t key_length,
    const uint64_t align_distance) {
  // Check total candidates
  const uint64_t total_candidates = interval_hi-interval_lo;
  if (gem_expect_false(total_candidates==0)) return;
  // Compute number of matches to add
  select_parameters_t* const select_parameters = &search_parameters->select_parameters_align;
  uint64_t pending_candidates;
  if (align_distance==0 && select_parameters->min_reported_strata_nominal==0) {
    pending_candidates = MIN(select_parameters->max_reported_matches,total_candidates);
  } else {
    pending_candidates = total_candidates;
  }
  // Store candidate positions
  vector_t* const filtering_positions = filtering_candidates->filtering_positions;
  vector_reserve_additional(filtering_positions,pending_candidates);
  filtering_position_t* filtering_position = vector_get_free_elm(filtering_positions,filtering_position_t);
  const uint64_t interval_top = interval_lo + pending_candidates;
  uint64_t index_position;
  for (index_position=interval_lo;index_position<interval_top;++index_position) {
    filtering_position->source_region_begin = 0;
    filtering_position->source_region_end = key_length;
    filtering_position->source_region_error = align_distance;
    filtering_position->region_index_position = index_position;
    filtering_position->align_distance = align_distance;
    ++filtering_position;
  }
  vector_update_used(filtering_positions,filtering_position);
}
void filtering_candidates_add_region_interval(
    filtering_candidates_t* const filtering_candidates,
    const uint64_t interval_lo,
    const uint64_t interval_hi,
    const uint64_t region_begin_pos,
    const uint64_t region_end_pos,
    const uint64_t region_errors) {
  // Check total candidates
  const uint64_t total_candidates = interval_hi-interval_lo;
  if (gem_expect_false(total_candidates==0)) return;
  // Store candidate positions
  vector_t* const filtering_positions = filtering_candidates->filtering_positions;
  vector_reserve_additional(filtering_positions,total_candidates);
  filtering_position_t* filtering_position = vector_get_free_elm(filtering_positions,filtering_position_t);
  uint64_t index_position;
  for (index_position=interval_lo;index_position<interval_hi;++index_position) {
    filtering_position->source_region_begin = region_begin_pos;
    filtering_position->source_region_end = region_end_pos;
    filtering_position->source_region_error = region_errors;
    filtering_position->region_index_position = index_position;
    filtering_position->align_distance = ALIGN_DISTANCE_INF;
    ++filtering_position;
  }
  vector_update_used(filtering_positions,filtering_position);
}
void filtering_candidates_add_region_interval_set(
    filtering_candidates_t* const filtering_candidates,
    interval_set_t* const interval_set,
    const uint64_t region_begin_pos,
    const uint64_t region_end_pos) {
  INTERVAL_SET_ITERATE(interval_set,interval) {
    filtering_candidates_add_region_interval(
        filtering_candidates,interval->lo,interval->hi,
        region_begin_pos,region_end_pos,interval->distance);
  }
}
void filtering_candidates_add_region_interval_set_thresholded(
    filtering_candidates_t* const filtering_candidates,
    interval_set_t* const interval_set,
    const uint64_t region_begin_pos,
    const uint64_t region_end_pos,
    const uint64_t max_error) {
  INTERVAL_SET_ITERATE(interval_set,interval) {
    if (interval->distance <= max_error) {
      filtering_candidates_add_region_interval(
          filtering_candidates,interval->lo,interval->hi,
          region_begin_pos,region_end_pos,interval->distance);
    }
  }
}
/*
 * Sorting
 */
int filtering_position_cmp_position(const filtering_position_t* const a,const filtering_position_t* const b) {
  const int cmp = a->text_begin_position - b->text_begin_position;
  if (cmp!=0) return cmp;
  return a->text_end_position - b->text_end_position;
}
int filtering_region_cmp_sort_align_distance(const filtering_region_t* const a,const filtering_region_t* const b) {
  return a->region_alignment.distance_min_bound - b->region_alignment.distance_min_bound;
}
int filtering_region_cmp_sort_scaffold_coverage(const filtering_region_t* const a,const filtering_region_t* const b) {
  return b->match_scaffold.scaffolding_coverage - a->match_scaffold.scaffolding_coverage;
}
int verified_region_cmp_position(const verified_region_t* const a,const verified_region_t* const b) {
  return a->begin_position - b->begin_position;
}
void filtering_positions_sort_positions(vector_t* const filtering_positions) {
  void* array = vector_get_mem(filtering_positions,filtering_position_t);
  const size_t count = vector_get_used(filtering_positions);
  qsort(array,count,sizeof(filtering_position_t),(int (*)(const void *,const void *))filtering_position_cmp_position);
}
void filtering_regions_sort_align_distance(vector_t* const filtering_regions) {
  void* array = vector_get_mem(filtering_regions,filtering_region_t);
  const size_t count = vector_get_used(filtering_regions);
  qsort(array,count,sizeof(filtering_region_t),(int (*)(const void *,const void *))filtering_region_cmp_sort_align_distance);
}
void filtering_regions_sort_scaffold_coverage(vector_t* const filtering_regions) {
  void* array = vector_get_mem(filtering_regions,filtering_region_t);
  const size_t count = vector_get_used(filtering_regions);
  qsort(array,count,sizeof(filtering_region_t),(int (*)(const void *,const void *))filtering_region_cmp_sort_scaffold_coverage);
}
void verified_regions_sort_positions(vector_t* const verified_regions) {
  void* array = vector_get_mem(verified_regions,verified_region_t);
  const size_t count = vector_get_used(verified_regions);
  qsort(array,count,sizeof(verified_region_t),(int (*)(const void *,const void *))verified_region_cmp_position);
}
/*
 * Display
 */
void filtering_candidates_print_regions_by_status(
    FILE* const stream,
    vector_t* const filtering_regions,
    const filtering_region_status_t status,
    const text_collection_t* const text_collection,
    const bool print_matching_regions) {
  uint64_t i, total_printed = 0;
  const uint64_t num_regions = vector_get_used(filtering_regions);
  filtering_region_t* const fregion = vector_get_mem(filtering_regions,filtering_region_t);
  // Count
  for (i=0;i<num_regions;++i) {
    if (fregion[i].status!=status) continue;
    ++total_printed;
  }
  if (total_printed == 0) return;
  tab_fprintf(stream,"  => Regions.%s  (%"PRIu64")\n",filtering_region_status_label[status],total_printed);
  // Print
  tab_global_inc();
  for (i=0;i<num_regions;++i) {
    if (fregion[i].status!=status) continue;
    filtering_region_print(stream,fregion+i,text_collection,print_matching_regions);
  }
  tab_global_dec();
}
void filtering_candidates_print_regions(
    FILE* const stream,
    filtering_candidates_t* const filtering_candidates,
    const bool print_matching_regions) {
  tab_fprintf(stream,"[GEM]>Filtering.Regions\n");
  text_collection_t* const text_collection = filtering_candidates->text_collection;
  vector_t* const filtering_regions = filtering_candidates->filtering_regions;
  vector_t* const discarded_regions = filtering_candidates->discarded_regions;
  filtering_candidates_print_regions_by_status(
      stream,filtering_regions,filtering_region_pending,text_collection,print_matching_regions);
  filtering_candidates_print_regions_by_status(
      stream,filtering_regions,filtering_region_unverified,text_collection,print_matching_regions);
  filtering_candidates_print_regions_by_status(
      stream,discarded_regions,filtering_region_verified_discarded,text_collection,print_matching_regions);
  filtering_candidates_print_regions_by_status(
      stream,filtering_regions,filtering_region_accepted,text_collection,print_matching_regions);
  filtering_candidates_print_regions_by_status(
      stream,discarded_regions,filtering_region_accepted_subdominant,text_collection,print_matching_regions);
  filtering_candidates_print_regions_by_status(
      stream,filtering_regions,filtering_region_aligned,text_collection,print_matching_regions);
  const uint64_t total_regions =
      vector_get_used(filtering_candidates->filtering_regions) +
      vector_get_used(filtering_candidates->discarded_regions);
  if (total_regions > 0) tab_fprintf(stream,"  => Total.Regions %"PRIu64"\n",total_regions);
}
