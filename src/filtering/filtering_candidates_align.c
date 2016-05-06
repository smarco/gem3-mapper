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
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  region_alignment_t* const region_alignment = &filtering_region->region_alignment;
  // Basic cases
  if (num_matches == 0) return false;
  if (num_matches < select_parameters->max_reported_matches) return false;
  // Bounded Cases (Only pays off to align matches that can be include within user report limits)
  switch (alignment_model) { // (Select alignment model)
    case alignment_model_hamming:
    case alignment_model_levenshtein: {
      const uint64_t candidate_min_distance_bound = region_alignment->distance_min_bound;
      match_trace_t* const last_ranked_match_trace = matches_get_ranked_match_trace(matches,select_parameters);
      // Need a candidate expected to have less distance than the current max
      return candidate_min_distance_bound >= last_ranked_match_trace->edit_distance;
    }
    case alignment_model_gap_affine: {
      const uint64_t candidate_edit_distance_bound = region_alignment->distance_min_bound;
      const uint64_t candidate_max_score_bound = align_swg_score_compute_max_score_bound(
          swg_penalties,candidate_edit_distance_bound,pattern->key_length);
      match_trace_t* const last_ranked_match_trace = matches_get_ranked_match_trace(matches,select_parameters);
      // Need a candidate expected to have better score than the current max
      return candidate_max_score_bound <= last_ranked_match_trace->swg_score;
    }
    default:
      return false;
  }
  return false;
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
  match_trace_t* const match_trace_cache = filtering_region_transient_cache_search(
      &filtering_candidates->filtering_region_cache,region,text_collection,matches);
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
    const bool emulated_rc_search,
    const bool local_alignment,
    const bool extended_match,
    matches_t* const matches) {
  // Parameters
  // Parameters
  archive_t* const archive = filtering_candidates->archive;
  locator_t* const locator = archive->locator;
  archive_text_t* const archive_text = archive->text;
  text_collection_t* const text_collection = filtering_candidates->text_collection;
  mm_stack_t* const mm_stack = filtering_candidates->mm_stack;
  search_parameters_t* const search_parameters = filtering_candidates->search_parameters;
  const alignment_model_t alignment_model = search_parameters->alignment_model;
  select_parameters_t* const select_parameters_align = &search_parameters->select_parameters_align;
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
        region,pattern,emulated_rc_search,local_alignment,matches,&match_trace);
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
    bool match_added, match_replaced;
    match_trace_t* match_trace_added;
    matches_add_match_trace__preserve_rank(
        matches,locator,&match_trace,select_parameters_align,
        alignment_model,&match_trace_added,&match_added,&match_replaced,mm_stack);
    if (match_added || match_replaced) {
      if (extended_match) match_trace_added->type = match_type_extended;
      filtering_region_transient_cache_add(
          &filtering_candidates->filtering_region_cache,
          region,match_trace_added->match_trace_offset,mm_stack);
    }
    // Return (Repeated?)
    return match_added;
  }
}
/*
 * Filtering Candidates (Re)Alignment
 */
uint64_t filtering_candidates_align_candidates(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    const bool emulated_rc_search,
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
  matches_hint_allocate_match_trace(matches,num_filtering_regions);
  // Prepare Candidate Vectors
  filtering_regions_sort_align_distance(filtering_candidates->filtering_regions); // Sort wrt align_distance
  vector_reserve_additional(filtering_candidates->discarded_regions,num_filtering_regions);
  filtering_region_t* regions_discarded = vector_get_free_elm(filtering_candidates->discarded_regions,filtering_region_t);
  filtering_region_t* regions_in = vector_get_mem(filtering_candidates->filtering_regions,filtering_region_t);
  filtering_region_t* regions_out = regions_in;
  // Clear cache
  filtering_region_cache_clear(&filtering_candidates->filtering_region_cache);
  // Traverse all accepted candidates (text-space)
  uint64_t n, num_accepted_regions = 0;
  for (n=0;n<num_filtering_regions;++n,++regions_in) {
    // Skip other regions
    if (regions_in->status != filtering_region_accepted) {
      *regions_out = *regions_in;
      ++regions_out;
      continue;
    }
    // Check if candidate is subdominant (check distance bounds)
    const bool candidate_subdominant = !extended_match &&
        filtering_candidates_align_is_subdominant(filtering_candidates,regions_in,pattern,matches);
    if (candidate_subdominant) {
      PROF_INC_COUNTER(GP_FC_SELECT_PRUNE_HIT);
      *regions_discarded = *regions_in;
      regions_discarded->status = filtering_region_accepted_subdominant;
      ++regions_discarded;
      continue;
    }
    // Align Region
    const bool accepted_region = filtering_candidates_align_region(filtering_candidates,
        regions_in,pattern,emulated_rc_search,local_alignment,extended_match,matches);
    if (accepted_region) ++num_accepted_regions;
  }
  // Update used
  matches_metrics_add_accepted_candidates(&matches->metrics,num_filtering_regions);
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
