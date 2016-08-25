/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates_align.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering/filtering_candidates_align.h"
#include "filtering/filtering_candidates_process.h"
#include "filtering/filtering_region_align.h"

#define FILTERING_CANDIDATES_ALIGN_CACHE
#define FILTERING_CANDIDATES_ALIGN_SELECT_PRUNE

/*
 * Debug
 */
#define DEBUG_FILTERING_CANDIDATES  GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Filtering Candidate distance bound
 */
bool filtering_candidates_align_is_subdominant(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    matches_t* const matches) {
  // Parameters
  search_parameters_t* const search_parameters = filtering_candidates->search_parameters;
  select_parameters_t* const select_parameters = &search_parameters->select_parameters_align;
  const alignment_model_t alignment_model = search_parameters->alignment_model;
  swg_penalties_t* const swg_penalties = &search_parameters->swg_penalties;
  region_alignment_t* const region_alignment = &filtering_region->region_alignment;
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  // Basic cases
  const uint64_t min_reported_strata = select_parameters->min_reported_strata_nominal;
  const uint64_t max_reported_matches = select_parameters->max_reported_matches;
  if (num_matches == 0 || num_matches < max_reported_matches || min_reported_strata > 0) return false;
  if (alignment_model != alignment_model_gap_affine) return false;
  // Bounded Cases (Only pays off to align matches that can be include within user report limits)
  // The candidate needs to have a expected max-score than the current max
  const uint64_t candidate_edit_distance_bound = region_alignment->distance_min_bound;
  const uint64_t candidate_max_score_bound = align_swg_score_compute_max_score_bound(
      swg_penalties,candidate_edit_distance_bound,pattern->key_length);
  match_trace_t** const match_traces = matches_get_match_traces(matches);
  return candidate_max_score_bound <= match_traces[max_reported_matches-1]->swg_score;
}
/*
 * Filtering Candidates Cache
 */
bool filtering_candidates_align_search_filtering_region_cache(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const region,
    const uint64_t run_length,
    matches_t* const matches,
    match_trace_t* const match_trace) {
  // Search the cache
  text_collection_t* const text_collection = filtering_candidates->text_collection;
  match_trace_t* const match_trace_cache =
      filtering_region_transient_cache_search(
          &filtering_candidates->filtering_region_cache,region,text_collection);
  if (match_trace_cache==NULL) return false;
  // Clone the match-trace found in the cache
  filtering_region_align_clone(match_trace_cache,match_trace,region,run_length);
  return true;
}
/*
 * Filtering Candidates Region Align
 */
bool filtering_candidates_align_region(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const region,
    pattern_t* const pattern,
    const bool local_alignment,
    const bool extended_match,
    matches_t* const matches) {
  // Parameters
  archive_t* const archive = filtering_candidates->archive;
  locator_t* const locator = archive->locator;
  archive_text_t* const archive_text = archive->text;
  text_collection_t* const text_collection = filtering_candidates->text_collection;
  mm_stack_t* const mm_stack = filtering_candidates->mm_stack;
  // Retrieve Candidate (if needed)
  filtering_region_retrieve_text(region,pattern,archive_text,text_collection,mm_stack);
  // Search Cache (Before jumping into aligning the region)
  match_trace_t match_trace;
  bool match_trace_aligned = !extended_match &&
      filtering_candidates_align_search_filtering_region_cache(
          filtering_candidates,region,pattern->run_length,matches,&match_trace);
  // Align the region
  if (!match_trace_aligned) {
    match_trace_aligned = filtering_region_align(filtering_candidates,
        region,pattern,local_alignment,matches,&match_trace);
    if (!match_trace_aligned) return false; // Not aligned or subdominant
  }
  // Add to matches
  const bool set_local_match_aside = (!local_alignment && !extended_match);
  if (set_local_match_aside && match_trace.type == match_type_local) {
    // Add Local Alignment (Pending)
    matches_add_local_match_pending(matches,&match_trace);
    return false; // Return (not added)
  } else {
    // Add Global Alignment
    bool match_replaced;
    match_trace_t* const match_trace_added =
        matches_add_match_trace(matches,locator,&match_trace,&match_replaced);
    if (match_trace_added==NULL) return false;
    if (extended_match) match_trace_added->type = match_type_extended;
    filtering_region_transient_cache_add(
        &filtering_candidates->filtering_region_cache,region,match_trace_added);
    return !match_replaced; // Return (Repeated?)
  }
}
/*
 * Filtering Candidates (Re)Alignment
 */
uint64_t filtering_candidates_align_candidates(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    const bool extended_match,
    const bool local_alignment,
    matches_t* const matches) {
  // DEBUG
  gem_cond_debug_block(DEBUG_FILTERING_CANDIDATES) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Filtering.Candidates (align_acepted_regions)\n");
    tab_global_inc();
  }
  // Hint to matches
  const uint64_t num_filtering_regions = vector_get_used(filtering_candidates->filtering_regions);
  if (num_filtering_regions==0) return 0;
  PROFILE_START(GP_FC_REALIGN_CANDIDATE_REGIONS,PROFILE_LEVEL);
  // Prepare Candidate Vectors
  filtering_regions_sort_align_distance(filtering_candidates->filtering_regions); // Sort wrt align_distance
  vector_reserve_additional(filtering_candidates->discarded_regions,num_filtering_regions);
  filtering_region_t* regions_discarded = vector_get_free_elm(filtering_candidates->discarded_regions,filtering_region_t);
  filtering_region_t* regions_in = vector_get_mem(filtering_candidates->filtering_regions,filtering_region_t);
  filtering_region_t* regions_out = regions_in;
  // Clear cache
  filtering_region_cache_clear(&filtering_candidates->filtering_region_cache);
  // Traverse all accepted candidates (text-space)
  uint64_t n, num_accepted_regions = 0, accepted_candidates = 0;
  for (n=0;n<num_filtering_regions;++n,++regions_in) {
    // Skip other regions
    if (regions_in->status != filtering_region_accepted) {
      *regions_out = *regions_in;
      ++regions_out;
      continue;
    }
    ++accepted_candidates;
    // Check if candidate is subdominant (check distance bounds)
    bool candidate_subdominant = false;
    if (!extended_match) {
      candidate_subdominant =
          filtering_candidates_align_is_subdominant(
              filtering_candidates,regions_in,pattern,matches);
    }
    if (candidate_subdominant) {
      PROF_INC_COUNTER(GP_FC_SELECT_PRUNE_HIT);
      matches_metrics_set_limited_candidates(&matches->metrics,true);
      *regions_discarded = *regions_in;
      regions_discarded->status = filtering_region_accepted_subdominant;
      ++regions_discarded;
      continue;
    }
    // Align Region
    const bool accepted_region =
        filtering_candidates_align_region(
            filtering_candidates,regions_in,pattern,
            local_alignment,extended_match,matches);
    if (accepted_region) ++num_accepted_regions;
  }
  // Update number of accepted regions (from verification)
  matches_metrics_add_accepted_candidates(&matches->metrics,accepted_candidates);
  // Update used
  vector_update_used(filtering_candidates->filtering_regions,regions_out);
  vector_update_used(filtering_candidates->discarded_regions,regions_discarded);
  // DEBUG
  gem_cond_debug_block(DEBUG_FILTERING_CANDIDATES) {
    tab_global_dec();
    tab_fprintf(gem_log_get_stream(),"[GEM]>Filtering.Candidates (align_acepted_regions - aftermath)\n");
    tab_global_inc();
    filtering_candidates_print_regions(gem_log_get_stream(),filtering_candidates,false);
    tab_global_dec();
  }
  // Return total accepted regions
  PROFILE_STOP(GP_FC_REALIGN_CANDIDATE_REGIONS,PROFILE_LEVEL);
  return num_accepted_regions;
}
