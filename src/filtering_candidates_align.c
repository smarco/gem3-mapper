/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates_align.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering_candidates_align.h"
#include "filtering_candidates_process.h"
#include "filtering_region_align.h"

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
 * Candidate distance bound
 */
bool filtering_candidates_align_is_behond_distance_limits(
    filtering_region_t* const filtering_region,select_parameters_t* const select_parameters,
    const alignment_model_t alignment_model,swg_penalties_t* const swg_penalties,
    const bool approximated_distance,const uint64_t key_length,
    matches_t* const matches) {
#ifdef FILTERING_CANDIDATES_ALIGN_SELECT_PRUNE
  // Parameters
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  // Basic cases
  if (num_matches == 0) return false;
  if (num_matches < select_parameters->max_reported_matches) return false;
  // Bounded Cases (Only pays off to align matches that can be include within user report limits)
  switch (alignment_model) { // (Select alignment model)
    case alignment_model_hamming:
    case alignment_model_levenshtein: {
      const uint64_t candidate_min_distance_bound = (approximated_distance) ?
          filtering_region->align_distance_min_bound : filtering_region->align_distance;
      match_trace_t* const last_ranked_match_trace = matches_get_ranked_match_trace(matches,select_parameters);
      // Need a candidate expected to have less distance than the current max
      return candidate_min_distance_bound >= last_ranked_match_trace->edit_distance;
    }
    case alignment_model_gap_affine: {
      const uint64_t candidate_edit_distance_bound = (approximated_distance) ?
          filtering_region->align_distance_min_bound : filtering_region->align_distance;
      const uint64_t candidate_max_score_bound =
          align_swg_score_compute_max_score_bound(swg_penalties,candidate_edit_distance_bound,key_length);
      match_trace_t* const last_ranked_match_trace = matches_get_ranked_match_trace(matches,select_parameters);
      // Need a candidate expected to have better score than the current max
      return candidate_max_score_bound <= last_ranked_match_trace->swg_score;
    }
    default:
      return false;
  }
#endif
  return false;
}
bool filtering_candidates_align_is_subdominant(
    filtering_region_t* const filtering_region,archive_text_t* const archive_text,
    const as_parameters_t* const as_parameters,select_parameters_t* const select_parameters,
    pattern_t* const pattern,const bool approximated_distance,
    const bool emulated_rc_search,matches_t* const matches,
    text_collection_t* const text_collection,mm_stack_t* const mm_stack) {
  // Parameters
  search_parameters_t* const search_parameters = as_parameters->search_parameters;
  const alignment_model_t alignment_model = search_parameters->alignment_model;
  swg_penalties_t* const swg_penalties = &search_parameters->swg_penalties;
  const uint64_t key_length = pattern->key_length;
  // Check distance limits
  const bool behond_distance_limits =
      filtering_candidates_align_is_behond_distance_limits(
          filtering_region,select_parameters,alignment_model,
          swg_penalties,approximated_distance,key_length,matches);
  if (approximated_distance && !behond_distance_limits) {
    // Candidate is going to be (re)aligned using a distance bound
    if (alignment_model==alignment_model_gap_affine) {
      // Adjust distance bound by scaffolding
//      const uint64_t distance_bound = filtering_region->align_distance_min_bound;
      filtering_region_align_adjust_distance_by_scaffolding(filtering_region,archive_text,
          as_parameters,pattern,matches,text_collection,mm_stack);
//      gem_slog(">Filtering.Candidates.Align.AdjustBound\tDistance-Bound=%lu Adjusted=%lu\n",
//          distance_bound,filtering_region->align_distance_min_bound);
      // Re-evaluate
      const uint64_t max_error = pattern->max_effective_filtering_error;
      if (filtering_region->align_distance_min_bound > max_error) {
        return true;
      }
      return filtering_candidates_align_is_behond_distance_limits(filtering_region,
          select_parameters,alignment_model,swg_penalties,true,key_length,matches);
    }
  }
  return behond_distance_limits;
}
bool filtering_candidates_align_search_filtering_region_cache(
    filtering_candidates_t* const filtering_candidates,filtering_region_t* const region,
    text_collection_t* const text_collection,match_trace_t* const match_trace,
    matches_t* const matches) {
#ifdef FILTERING_CANDIDATES_ALIGN_CACHE
  // Search the cache
  match_trace_t* const match_trace_cache =
      filtering_region_transient_cache_search(
          &filtering_candidates->filtering_region_cache,region,text_collection,matches);
  if (match_trace_cache==NULL) return false;
  // Clone the match-trace found in the cache
  filtering_region_align_clone(match_trace_cache,match_trace,region,text_collection);
  return true;
#else
  return false;
#endif
}
bool filtering_candidates_align_region(
    filtering_candidates_t* const filtering_candidates,filtering_region_t* const region,
    archive_text_t* const archive_text,const locator_t* const locator,
    text_collection_t* const text_collection,pattern_t* const pattern,
    const bool emulated_rc_search,const as_parameters_t* const as_parameters,
    const bool approximated_distance,const bool align_always,
    matches_t* const matches,mm_stack_t* const mm_stack) {
  // Parameters
  search_parameters_t* const search_parameters = as_parameters->search_parameters;
  const alignment_model_t alignment_model = search_parameters->alignment_model;
  select_parameters_t* const select_parameters_align = &search_parameters->select_parameters_align;
  // Retrieve Candidate (if needed)
  if (region->text_trace_offset == UINT64_MAX) {
    const uint64_t text_length = region->end_position-region->begin_position;
    region->text_trace_offset = archive_text_retrieve(archive_text,
        text_collection,region->begin_position,text_length,false,mm_stack);
  }
  // Align the region (First Search Cache)
  match_trace_t match_trace;
  bool match_trace_aligned = !align_always &&
      filtering_candidates_align_search_filtering_region_cache(
          filtering_candidates,region,text_collection,&match_trace,matches);
  if (!match_trace_aligned) {
    // Align the region
    match_trace_aligned = filtering_region_align(region,archive_text,text_collection,
        as_parameters,emulated_rc_search,pattern,matches,&match_trace,mm_stack);
    if (!match_trace_aligned) return false; // Not aligned or subdominant
  }
  // Add to matches
  bool match_added, match_replaced;
  match_trace_t* match_trace_added;
  matches_add_match_trace__preserve_rank(
      matches,locator,true,&match_trace,select_parameters_align,
      alignment_model,&match_trace_added,&match_added,&match_replaced,mm_stack);
#ifdef FILTERING_CANDIDATES_ALIGN_CACHE
  if (match_added || match_replaced) {
    filtering_region_transient_cache_add(
        &filtering_candidates->filtering_region_cache,
        region,match_trace_added->match_trace_offset,mm_stack);
  }
#endif
  return match_added; // Repeated
}
/*
 * Candidate Alignment
 */
uint64_t filtering_candidates_align_accepted_regions(
    filtering_candidates_t* const filtering_candidates,
    archive_text_t* const archive_text,const locator_t* const locator,
    text_collection_t* const text_collection,pattern_t* const pattern,
    const bool emulated_rc_search,const as_parameters_t* const as_parameters,
    const bool approximated_distance,const bool align_always,
    matches_t* const matches,mm_stack_t* const mm_stack) {
  // DEBUG
  gem_cond_debug_block(DEBUG_FILTERING_CANDIDATES) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Filtering.Candidates (align_acepted_regions)\n");
    tab_global_inc();
  }
  // Parameters
  search_parameters_t* const search_parameters = as_parameters->search_parameters;
  select_parameters_t* const select_parameters = &search_parameters->select_parameters_align;
  // Prepare Filtering Regions (Sort wrt align_distance)
  filtering_regions_sort_align_distance(filtering_candidates->filtering_regions);
  // Prepare Candidate Vectors
  const uint64_t num_filtering_regions = vector_get_used(filtering_candidates->filtering_regions);
  vector_reserve_additional(filtering_candidates->discarded_regions,num_filtering_regions);
  filtering_region_t* regions_discarded = vector_get_free_elm(filtering_candidates->discarded_regions,filtering_region_t);
  filtering_region_t* regions_in = vector_get_mem(filtering_candidates->filtering_regions,filtering_region_t);
  filtering_region_t* regions_out = regions_in;
  filtering_region_cache_clear(&filtering_candidates->filtering_region_cache); // Clear cache
  // Traverse all accepted candidates (text-space)
  uint64_t n, num_accepted_regions = 0;
  for (n=0;n<num_filtering_regions;++n,++regions_in) {
    if (regions_in->status != filtering_region_accepted) { // Skip other regions
      *regions_out = *regions_in;
      ++regions_out;
    } else {
      // Check if candidate is subdominant (check distance bounds)
      const bool candidate_subdominant = !align_always &&
          filtering_candidates_align_is_subdominant(
            regions_in,archive_text,as_parameters,select_parameters,pattern,
            approximated_distance,emulated_rc_search,matches,text_collection,mm_stack);
      if (candidate_subdominant) {
        PROF_INC_COUNTER(GP_FC_SELECT_PRUNE_HIT);
        *regions_discarded = *regions_in;
        regions_discarded->status = filtering_region_accepted_subdominant;
        ++regions_discarded;
      } else {
        // Align Region
        const bool region_aligned = filtering_candidates_align_region(
            filtering_candidates,regions_in,archive_text,locator,text_collection,pattern,
            emulated_rc_search,as_parameters,approximated_distance,align_always,matches,mm_stack);
        if (region_aligned) {
          ++num_accepted_regions;
        } else {
          // Not aligned (probably because is too distant; sub-dominant)
          *regions_discarded = *regions_in;
          ++regions_discarded;
        }
      }
    }
  }
  // Update used
  matches_metrics_add_accepted_candidates(&matches->metrics,num_filtering_regions);
  vector_update_used(filtering_candidates->discarded_regions,regions_discarded);
  vector_update_used(filtering_candidates->filtering_regions,regions_out);
  // DEBUG
  gem_cond_debug_block(DEBUG_FILTERING_CANDIDATES) {
    tab_global_dec();
    tab_fprintf(gem_log_get_stream(),"[GEM]>Filtering.Candidates (align_acepted_regions - aftermath)\n");
    tab_global_inc();
    filtering_candidates_print_regions(gem_log_get_stream(),filtering_candidates,text_collection,false,false);
    tab_global_dec();
  }
  // Return total accepted regions
  return num_accepted_regions;
}
uint64_t filtering_candidates_unbounded_align_regions_by_status(
    vector_t* const filtering_regions,const filtering_region_status_t region_status,
    archive_text_t* const archive_text,const locator_t* const locator,
    text_collection_t* const text_collection,pattern_t* const pattern,
    const bool emulated_rc_search,const as_parameters_t* const as_parameters,
    matches_t* const matches,mm_stack_t* const mm_stack) {
  // Parameters
  search_parameters_t* const search_parameters = as_parameters->search_parameters;
  const uint64_t num_regions = vector_get_used(filtering_regions);
  filtering_region_t* const regions = vector_get_mem(filtering_regions,filtering_region_t);
  // Iterate over all regions
  uint64_t n, num_unbounded_alignments = 0;
  for (n=0;n<num_regions;++n) {
    filtering_region_t* const filtering_region = regions + n;
    if (filtering_region->status != region_status) continue;
    // Retrieve Candidate (if needed)
    if (filtering_region->text_trace_offset == UINT64_MAX) {
      if (filtering_region->align_distance>0 && search_parameters->alignment_model!=alignment_model_none) {
        const uint64_t text_length = filtering_region->end_position-filtering_region->begin_position;
        filtering_region->text_trace_offset = archive_text_retrieve(archive_text,text_collection,
            filtering_region->begin_position,text_length,false,mm_stack);
      }
    }
    // Align the region
    match_trace_t match_trace;
    if (filtering_region_align_unbounded(filtering_region,archive_text,text_collection,
        as_parameters,emulated_rc_search,pattern,matches,&match_trace,mm_stack)) {
      filtering_region->status = filtering_region_aligned_unbounded;
      if (matches_add_match_trace(matches,locator,true,&match_trace,mm_stack)) {
        ++num_unbounded_alignments;
      }
    }
  }
  // Return number of unbounded alignments
  return num_unbounded_alignments;
}
uint64_t filtering_candidates_unbounded_align_regions(
    filtering_candidates_t* const filtering_candidates,
    archive_text_t* const archive_text,const locator_t* const locator,
    text_collection_t* const text_collection,pattern_t* const pattern,
    const bool emulated_rc_search,const as_parameters_t* const as_parameters,
    matches_t* const matches,mm_stack_t* const mm_stack) {
  uint64_t total_unbounded_alignments = 0;
  // 1. Align unbounded aligned-subdominant
  total_unbounded_alignments += filtering_candidates_unbounded_align_regions_by_status(
      filtering_candidates->discarded_regions,filtering_region_aligned_subdominant,
      archive_text,locator,text_collection,pattern,emulated_rc_search,as_parameters,matches,mm_stack);
  // 2. Align unbounded accepted-subdominant
  if (total_unbounded_alignments == 0) {
    total_unbounded_alignments += filtering_candidates_unbounded_align_regions_by_status(
        filtering_candidates->discarded_regions,filtering_region_accepted_subdominant,
        archive_text,locator,text_collection,pattern,emulated_rc_search,as_parameters,matches,mm_stack);
  }
  // 3. Align unbounded verify-discarded
  if (total_unbounded_alignments == 0) {
    total_unbounded_alignments += filtering_candidates_unbounded_align_regions_by_status(
        filtering_candidates->discarded_regions,filtering_region_verified_discarded,
        archive_text,locator,text_collection,pattern,emulated_rc_search,as_parameters,matches,mm_stack);
  }
  // DEBUG
  gem_cond_debug_block(DEBUG_FILTERING_CANDIDATES) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Filtering.Candidates (unbounded_align)\n");
    tab_global_inc();
    filtering_candidates_print_regions(gem_log_get_stream(),filtering_candidates,text_collection,false,false);
    tab_global_dec();
  }
  // Clear discarded vector
  vector_clear(filtering_candidates->discarded_regions);
  // Return total unbounded-alignments
  return total_unbounded_alignments;
}
uint64_t filtering_candidates_align_candidates(
    filtering_candidates_t* const filtering_candidates,
    archive_text_t* const archive_text,const locator_t* const locator,
    text_collection_t* const text_collection,pattern_t* const pattern,
    const bool emulated_rc_search,const as_parameters_t* const as_parameters,
    const bool approximated_distance,const bool align_always,
    matches_t* const matches,mm_stack_t* const mm_stack) {
  // Hint to matches
  const uint64_t num_filtering_regions = vector_get_used(filtering_candidates->filtering_regions);
  if (num_filtering_regions==0) return 0;
  PROFILE_START(GP_FC_REALIGN_CANDIDATE_REGIONS,PROFILE_LEVEL);
  matches_hint_allocate_match_trace(matches,num_filtering_regions);
  // Realign
  const uint64_t aligned_regions = filtering_candidates_align_accepted_regions(
      filtering_candidates,archive_text,locator,text_collection,pattern,emulated_rc_search,
      as_parameters,approximated_distance,align_always,matches,mm_stack);
  PROFILE_STOP(GP_FC_REALIGN_CANDIDATE_REGIONS,PROFILE_LEVEL);
  // Return number of aligned regions
  return aligned_regions;
}
/*
 * Search for unbounded-alignments
 */
uint64_t filtering_candidates_align_unbounded(
    filtering_candidates_t* const filtering_candidates,
    archive_text_t* const archive_text,const locator_t* const locator,
    text_collection_t* const text_collection,pattern_t* const pattern,
    const bool emulated_rc_search,const as_parameters_t* const as_parameters,
    matches_t* const matches,mm_stack_t* const mm_stack) {
  // Check number of filtering regions
  const uint64_t discarded_candidates = vector_get_used(filtering_candidates->discarded_regions);
  if (discarded_candidates==0) return 0;
  PROFILE_START(GP_FC_UNBOUNDED_ALIGNMENT,PROFILE_LEVEL);
  // Unbounded-align all discarded candidates
  matches_hint_allocate_match_trace(matches,discarded_candidates);
  const uint64_t aligned_regions = filtering_candidates_unbounded_align_regions(
      filtering_candidates,archive_text,locator,text_collection,pattern,
      emulated_rc_search,as_parameters,matches,mm_stack);
  PROFILE_STOP(GP_FC_UNBOUNDED_ALIGNMENT,PROFILE_LEVEL);
  return aligned_regions;
}
