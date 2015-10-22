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

/*
 * Debug
 */
#define DEBUG_FILTERING_CANDIDATES  GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PMED
//#define FILTERING_CANDIDATES_CACHE

/*
 * Candidate Alignment
 */
GEM_INLINE uint64_t filtering_candidates_align_accepted_regions(
    filtering_candidates_t* const filtering_candidates,
    archive_text_t* const archive_text,const locator_t* const locator,
    text_collection_t* const text_collection,pattern_t* const pattern,
    const bool emulated_rc_search,const as_parameters_t* const as_parameters,
    const bool approximated_distance,matches_t* const matches,mm_stack_t* const mm_stack) {
  // DEBUG
  gem_cond_debug_block(DEBUG_FILTERING_CANDIDATES) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Filtering.Candidates (align_acepted_regions)\n");
    tab_global_inc();
  }
  // Parameters
  search_parameters_t* const search_parameters = as_parameters->search_parameters;
  const uint64_t max_error_after_best = as_parameters->alignment_max_error_after_best_nominal;
  // Sort filtering regions (wrt align_distance)
  filtering_regions_sort_align_distance(filtering_candidates->filtering_regions);
  // Prepare candidate vectors
  const uint64_t num_filtering_regions = vector_get_used(filtering_candidates->filtering_regions);
  vector_reserve_additional(filtering_candidates->discarded_regions,num_filtering_regions);
  filtering_region_t* regions_discarded = vector_get_free_elm(filtering_candidates->discarded_regions,filtering_region_t);
  filtering_region_t* regions_in = vector_get_mem(filtering_candidates->filtering_regions,filtering_region_t);
  filtering_region_t* regions_out = regions_in;
#ifdef FILTERING_CANDIDATES_CACHE
  filtering_region_t* last_aligned_region = NULL;
#endif
  // Traverse all accepted candidates (text-space)
  uint64_t n, num_accepted_regions = 0;
  for (n=0;n<num_filtering_regions;++n,++regions_in) {
    if (regions_in->status != filtering_region_accepted) { // Skip other regions
      *regions_out = *regions_in;
      ++regions_out;
    } else {
      // Check filtering strata-after-best
      const uint64_t min_edit_distance = matches_metrics_get_min_edit_distance(&matches->metrics);
      const uint64_t distance_bound = (approximated_distance) ?
          regions_in->align_distance_min_bound : regions_in->align_distance;
      if (distance_bound > min_edit_distance + max_error_after_best) {
        *regions_discarded = *regions_in;
        regions_discarded->status = filtering_region_accepted_subdominant;
        matches_metrics_inc_subdominant_candidates(&matches->metrics);
        ++regions_discarded;
      } else {
        // Retrieve Candidate (if needed)
        if (regions_in->text_trace_offset == UINT64_MAX) {
          if (regions_in->align_distance>0 && search_parameters->alignment_model!=alignment_model_none) {
            const uint64_t text_length = regions_in->end_position-regions_in->begin_position;
            regions_in->text_trace_offset = archive_text_retrieve(archive_text,text_collection,
                regions_in->begin_position,text_length,false,mm_stack);
          }
        }
        // Search the cache (of previous aligned regions)
        match_trace_t match_trace;
#ifdef FILTERING_CANDIDATES_CACHE
        if (last_aligned_region!=NULL && last_aligned_region->align_distance==regions_in->align_distance) {
          // Search cache
          filtering_region_cache_compute_footprint(regions_in,text_collection);
          match_trace_t* const cached_match_trace = filtering_region_cache_search(
              &filtering_candidates->filtering_region_cache,regions_in,text_collection);
          if (cached_match_trace!=NULL) {
            filtering_region_align_clone(cached_match_trace,&match_trace,regions_in,text_collection);
            // Add to matches
            if (matches_add_match_trace(matches,&match_trace,true,locator,mm_stack)) ++num_accepted_regions;
            // Next region
            continue;
          }
        }
#endif
        // Align the region
        if (filtering_region_align(regions_in,archive_text,text_collection,
            as_parameters,emulated_rc_search,pattern,matches,&match_trace,mm_stack)) {
          // Add to matches
          if (matches_add_match_trace(matches,&match_trace,true,locator,mm_stack)) {
            ++num_accepted_regions;
#ifdef FILTERING_CANDIDATES_CACHE
            last_aligned_region = regions_in;
            filtering_region_cache_compute_footprint(regions_in,text_collection); // FIXME Cond.
            filtering_region_cache_add(&filtering_candidates->filtering_region_cache,regions_in,&match_trace,mm_stack);
#endif
          }
        } else {
          *regions_discarded = *regions_in;
          ++regions_discarded;
        }
      }
    }
  }
  // Update used
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
GEM_INLINE uint64_t filtering_candidates_unbounded_align_regions_by_status(
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
      if (filtering_region->status==filtering_region_accepted_subdominant) {
        matches_metrics_dec_subdominant_candidates(&matches->metrics);
      }
      filtering_region->status = filtering_region_aligned_unbounded;
      if (matches_add_match_trace(matches,&match_trace,true,locator,mm_stack)) ++num_unbounded_alignments;
    }
  }
  // Return number of unbounded alignments
  return num_unbounded_alignments;
}
GEM_INLINE uint64_t filtering_candidates_unbounded_align_regions(
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
GEM_INLINE uint64_t filtering_candidates_align_candidates(
    filtering_candidates_t* const filtering_candidates,
    archive_text_t* const archive_text,const locator_t* const locator,
    text_collection_t* const text_collection,pattern_t* const pattern,
    const bool emulated_rc_search,const as_parameters_t* const as_parameters,
    const bool approximated_distance,matches_t* const matches,mm_stack_t* const mm_stack) {
  // Hint to matches
  const uint64_t num_filtering_regions = vector_get_used(filtering_candidates->filtering_regions);
  if (num_filtering_regions==0) return 0;
  PROFILE_START(GP_FC_REALIGN_CANDIDATE_REGIONS,PROFILE_LEVEL);
  matches_hint_allocate_match_trace(matches,num_filtering_regions);
  // Realign
  const uint64_t aligned_regions = filtering_candidates_align_accepted_regions(
      filtering_candidates,archive_text,locator,text_collection,pattern,emulated_rc_search,
      as_parameters,approximated_distance,matches,mm_stack);
  PROFILE_STOP(GP_FC_REALIGN_CANDIDATE_REGIONS,PROFILE_LEVEL);
  // Return number of aligned regions
  return aligned_regions;
}
/*
 * Search for unbounded-alignments
 */
GEM_INLINE uint64_t filtering_candidates_align_unbounded(
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
