/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering_candidates.h"
#include "archive_text.h"
#include "filtering_region.h"
#include "match_align.h"
#include "matches_classify.h"

/*
 * Debug
 */
#define DEBUG_FILTERING_REGIONS            GEM_DEEP_DEBUG
#define DEBUG_VERIFY_REGIONS               GEM_DEEP_DEBUG
#define DEBUG_ACCEPTED_REGIONS             GEM_DEEP_DEBUG
#define DEBUG_REGIONS_MATCHING             false // GEM_DEEP_DEBUG

/*
 * Mode
 */
#define SPLIT_EXTENSION_WINDOW

/*
 * Constants
 */
#define REGIONS_BUFFER_INIT                      100
#define CANDIDATE_POSITIONS_INIT                 1000
#define DECODE_NUM_POSITIONS_PREFETCHED          10
#define RETRIEVE_SAMPLE_NUM_POSITIONS_PREFETCHED 10

/*
 * Candidate Position
 */
typedef struct {
  // Source Region
  uint64_t source_region_offset;
  locator_interval_t* locator_interval;
  // Decode data
  uint64_t decode_distance;
  uint64_t decode_sampled_pos;
  // Region location
  uint64_t region_index_position;
  uint64_t region_text_position;
  // Final position
  uint64_t begin_position;       // Region effective begin position (adjusted to error boundaries)
  uint64_t end_position;         // Region effective end position (adjusted to error boundaries)
  uint64_t base_position_offset; // Offset to base filtering position (Begin position without error boundary correction)
} filtering_position_t;
/*
 * Setup
 */
GEM_INLINE void filtering_candidates_init(filtering_candidates_t* const filtering_candidates) {
  // Region Buffer
  filtering_candidates->regions_buffer = vector_new(REGIONS_BUFFER_INIT,region_search_t);
  // Candidates
  filtering_candidates->filtering_positions = vector_new(CANDIDATE_POSITIONS_INIT,filtering_position_t);
  filtering_candidates->filtering_regions = vector_new(CANDIDATE_POSITIONS_INIT,filtering_region_t);
  filtering_candidates->discarded_regions = vector_new(CANDIDATE_POSITIONS_INIT,filtering_region_t);
  filtering_candidates->verified_regions = vector_new(CANDIDATE_POSITIONS_INIT,verified_region_t);
  filtering_candidates->filtering_regions_locator = NULL;
}
GEM_INLINE void filtering_candidates_clear(filtering_candidates_t* const filtering_candidates) {
  // Region Buffer
  vector_clear(filtering_candidates->regions_buffer);
  // Candidates
  vector_clear(filtering_candidates->filtering_positions);
  vector_clear(filtering_candidates->filtering_regions);
  vector_clear(filtering_candidates->discarded_regions);
  vector_clear(filtering_candidates->verified_regions);
  if (filtering_candidates->filtering_regions_locator) vector_clear(filtering_candidates->filtering_regions_locator);
}
GEM_INLINE void filtering_candidates_destroy(filtering_candidates_t* const filtering_candidates) {
  vector_delete(filtering_candidates->regions_buffer);
  vector_delete(filtering_candidates->filtering_positions);
  vector_delete(filtering_candidates->filtering_regions);
  vector_delete(filtering_candidates->discarded_regions);
  vector_delete(filtering_candidates->verified_regions);
  if (filtering_candidates->filtering_regions_locator) vector_delete(filtering_candidates->filtering_regions_locator);
}
/*
 * Accessors
 */
GEM_INLINE uint64_t filtering_candidates_get_num_candidate_regions(const filtering_candidates_t* const filtering_candidates) {
  return vector_get_used(filtering_candidates->filtering_regions);
}
GEM_INLINE uint64_t filtering_candidates_count_candidate_regions(
    filtering_candidates_t* const filtering_candidates_end,const filtering_region_status_t filtering_region_status) {
  uint64_t count = 0;
  VECTOR_ITERATE(filtering_candidates_end->filtering_regions,filtering_region,n,filtering_region_t) {
    if (filtering_region->status == filtering_region_status) ++count;
  }
  return count;
}
GEM_INLINE void filtering_candidates_set_all_regions_pending(filtering_candidates_t* const filtering_candidates) {
  // Set all filtering regions as pending
  VECTOR_ITERATE(filtering_candidates->filtering_regions,filtering_region,n,filtering_region_t) {
    filtering_region->status = filtering_region_pending;
  }
}
GEM_INLINE void filtering_candidates_set_all_regions_unverified(filtering_candidates_t* const filtering_candidates) {
  // Set all filtering regions as pending
  VECTOR_ITERATE(filtering_candidates->filtering_regions,filtering_region,n,filtering_region_t) {
    filtering_region->status = filtering_region_unverified;
  }
}
/*
 * Adding candidate positions
 */
GEM_INLINE uint64_t filtering_candidates_add_region(
    filtering_candidates_t* const filtering_candidates,
    const uint64_t region_begin_pos,const uint64_t region_end_pos,
    const uint64_t region_errors) {
  // Store region
  region_search_t* region;
  vector_alloc_new(filtering_candidates->regions_buffer,region_search_t,region);
  region->begin = region_begin_pos;
  region->end = region_end_pos;
  region->degree = region_errors;
  return vector_get_used(filtering_candidates->regions_buffer)-1;
}
GEM_INLINE void filtering_candidates_add_interval(
    filtering_candidates_t* const filtering_candidates,
    const uint64_t interval_lo,const uint64_t interval_hi,
    const uint64_t region_begin_pos,const uint64_t region_end_pos,
    const uint64_t region_errors,mm_stack_t* const mm_stack) {
  const uint64_t num_candidates = interval_hi-interval_lo;
  if (gem_expect_false(num_candidates==0)) return;
  // Store region
  const uint64_t region_offset = filtering_candidates_add_region(
      filtering_candidates,region_begin_pos,region_end_pos,region_errors);
  // Store candidate positions (index-space)
  vector_reserve_additional(filtering_candidates->filtering_positions,num_candidates);
  filtering_position_t* position_index = vector_get_free_elm(filtering_candidates->filtering_positions,filtering_position_t);
  uint64_t index_position;
  for (index_position=interval_lo;index_position<interval_hi;++index_position) {
    position_index->source_region_offset = region_offset;
    position_index->region_index_position = index_position;
    ++position_index;
  }
  vector_update_used(filtering_candidates->filtering_positions,position_index);
}
GEM_INLINE void filtering_candidates_add_interval_set(
    filtering_candidates_t* const filtering_candidates,interval_set_t* const interval_set,
    const uint64_t region_begin_pos,const uint64_t region_end_pos,mm_stack_t* const mm_stack) {
  INTERVAL_SET_ITERATE(interval_set,interval) {
    filtering_candidates_add_interval(filtering_candidates,interval->lo,interval->hi,
        region_begin_pos,region_end_pos,interval->distance,mm_stack);
  }
}
GEM_INLINE void filtering_candidates_add_interval_set_thresholded(
    filtering_candidates_t* const filtering_candidates,interval_set_t* const interval_set,
    const uint64_t region_begin_pos,const uint64_t region_end_pos,
    const uint64_t max_error,mm_stack_t* const mm_stack) {
  INTERVAL_SET_ITERATE(interval_set,interval) {
    if (interval->distance <= max_error) {
      filtering_candidates_add_interval(filtering_candidates,interval->lo,interval->hi,
          region_begin_pos,region_end_pos,interval->distance,mm_stack);
    }
  }
}
/*
 * Sorting
 */
int filtering_position_cmp_position(const filtering_position_t* const a,const filtering_position_t* const b) {
  return a->begin_position - b->begin_position;
}
int filtering_region_cmp_sort_align_distance(const filtering_region_t* const a,const filtering_region_t* const b) {
  return a->align_distance - b->align_distance;
}
int verified_region_cmp_position(const verified_region_t* const a,const verified_region_t* const b) {
  return a->begin_position - b->begin_position;
}
GEM_INLINE void filtering_positions_sort_positions(vector_t* const filtering_positions) {
  void* array = vector_get_mem(filtering_positions,filtering_position_t);
  const size_t count = vector_get_used(filtering_positions);
  qsort(array,count,sizeof(filtering_position_t),(int (*)(const void *,const void *))filtering_position_cmp_position);
}
GEM_INLINE void filtering_regions_sort_align_distance(vector_t* const filtering_regions) {
  void* array = vector_get_mem(filtering_regions,filtering_region_t);
  const size_t count = vector_get_used(filtering_regions);
  qsort(array,count,sizeof(filtering_region_t),(int (*)(const void *,const void *))filtering_region_cmp_sort_align_distance);
}
GEM_INLINE void verified_regions_sort_positions(vector_t* const verified_regions) {
  void* array = vector_get_mem(verified_regions,verified_region_t);
  const size_t count = vector_get_used(verified_regions);
  qsort(array,count,sizeof(verified_region_t),(int (*)(const void *,const void *))verified_region_cmp_position);
}
/*
 * Candidate Alignment
 */
GEM_INLINE uint64_t filtering_candidates_align_accepted_regions(
    filtering_candidates_t* const filtering_candidates,
    archive_text_t* const archive_text,const locator_t* const locator,
    text_collection_t* const text_collection,pattern_t* const pattern,
    const bool emulated_rc_search,const as_parameters_t* const as_parameters,
    const bool approximated_distance,matches_t* const matches,mm_stack_t* const mm_stack) {
  // Parameters
  const uint64_t max_delta = as_parameters->alignment_max_error_after_best_nominal;
  // Sort filtering regions (wrt align_distance)
  filtering_regions_sort_align_distance(filtering_candidates->filtering_regions);
  // Prepare candidate vectors
  const uint64_t num_filtering_regions = vector_get_used(filtering_candidates->filtering_regions);
  vector_reserve_additional(filtering_candidates->discarded_regions,num_filtering_regions);
  filtering_region_t* regions_discarded = vector_get_free_elm(filtering_candidates->discarded_regions,filtering_region_t);
  filtering_region_t* regions_in = vector_get_mem(filtering_candidates->filtering_regions,filtering_region_t);
  filtering_region_t* regions_out = regions_in;
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
      if (distance_bound > min_edit_distance+max_delta) { // FIXME
        gem_cond_debug_block(DEBUG_ACCEPTED_REGIONS) {
          gem_slog("[GEM]>Accepted.Regions\n");
          gem_slog("  => Region [%lu,%lu) discarded (Edit-bound=%lu,Edit-min=%lu,Margin=%lu)\n",
              regions_in->begin_position,regions_in->end_position,
              distance_bound,min_edit_distance,max_delta);
        }
        *regions_discarded = *regions_in;
        regions_discarded->status = filtering_region_accepted_subdominant;
        matches_metrics_inc_subdominant_candidates(&matches->metrics);
        ++regions_discarded;
      } else {
        // DEBUG
        gem_cond_debug_block(DEBUG_ACCEPTED_REGIONS) {
          gem_slog("[GEM]>Accepted.Regions\n");
          gem_slog("  => Region [%lu,%lu) (dist=%lu,matching-regions=%lu)\n",
              regions_in->begin_position,regions_in->end_position,
              regions_in->align_distance,regions_in->match_scaffold.num_scaffold_regions);
        }
        // Align the region
        match_trace_t match_trace;
        if (filtering_region_align(regions_in,archive_text,text_collection,
            as_parameters,emulated_rc_search,pattern,matches,&match_trace,mm_stack)) {
          // Add to matches
          if (matches_add_match_trace(matches,&match_trace,true,locator,mm_stack)) ++num_accepted_regions;
        } else {
          *regions_discarded = *regions_in;
          ++regions_discarded;
        }
      }
    }
  }
  // Update vectors
  vector_update_used(filtering_candidates->discarded_regions,regions_discarded);
  vector_update_used(filtering_candidates->filtering_regions,regions_out);
  // Return total accepted regions
  return num_accepted_regions;
}
GEM_INLINE uint64_t filtering_candidates_unbounded_align_discarded_regions(
    filtering_candidates_t* const filtering_candidates,
    archive_text_t* const archive_text,const locator_t* const locator,
    text_collection_t* const text_collection,pattern_t* const pattern,
    const bool emulated_rc_search,const as_parameters_t* const as_parameters,
    matches_t* const matches,mm_stack_t* const mm_stack) {
  // Prepare candidate vectors
  const uint64_t num_discarded_regions = vector_get_used(filtering_candidates->discarded_regions);
  filtering_region_t* regions_discarded = vector_get_mem(filtering_candidates->discarded_regions,filtering_region_t);
  // Traverse all accepted candidates (text-space)
  uint64_t n, num_unbounded_alignments = 0;
  for (n=0;n<num_discarded_regions;++n,++regions_discarded) {
    // Align the region
    match_trace_t match_trace;
    if (filtering_region_align_unbounded(regions_discarded,archive_text,text_collection,
        as_parameters,emulated_rc_search,pattern,matches,&match_trace,mm_stack)) {
      // Add to matches
      if (regions_discarded->status == filtering_region_accepted_subdominant) {
        matches_metrics_dec_subdominant_candidates(&matches->metrics);
      }
      if (matches_add_match_trace(matches,&match_trace,true,locator,mm_stack)) ++num_unbounded_alignments;
    }
  }
  // Clear discarded vector
  vector_clear(filtering_candidates->discarded_regions);
  // Return total unbounded-alignments
  return num_unbounded_alignments;
}
/*
 * Candidate Verification
 */
GEM_INLINE uint64_t filtering_candidates_verify_filtering_regions(
    filtering_candidates_t* const filtering_candidates,text_collection_t* const candidates_collection,
    const pattern_t* const pattern,const as_parameters_t* const as_parameters,
    matches_t* const matches) {
  // Parameters
  search_parameters_t* const search_parameters = as_parameters->search_parameters;
  // Traverse all regions (text-space)
  const uint64_t num_filtering_regions = vector_get_used(filtering_candidates->filtering_regions);
  vector_reserve_additional(filtering_candidates->verified_regions,num_filtering_regions);
  vector_reserve_additional(filtering_candidates->discarded_regions,num_filtering_regions);
  verified_region_t* regions_verified = vector_get_free_elm(filtering_candidates->verified_regions,verified_region_t);
  filtering_region_t* regions_discarded = vector_get_free_elm(filtering_candidates->discarded_regions,filtering_region_t);
  filtering_region_t* regions_in = vector_get_mem(filtering_candidates->filtering_regions,filtering_region_t);
  filtering_region_t* regions_out = regions_in;
  uint64_t n, num_regions_accepted = 0;
  gem_cond_debug_block(DEBUG_VERIFY_REGIONS) { gem_slog("[GEM]>Verify.Regions\n"); }
  for (n=0;n<num_filtering_regions;++n,++regions_in) {
    // Check region status (Skip other than unverified)
    if (regions_in->status!=filtering_region_unverified) {
      *regions_out = *regions_in;
      ++regions_out;
    } else {
      // Verify region
      if (filtering_region_verify(regions_in,candidates_collection,search_parameters,pattern)) {
        *regions_out = *regions_in;
        ++regions_out;
        ++num_regions_accepted;
      } else {
        *regions_discarded = *regions_in;
        ++regions_discarded;
      }
      // Add to verify regions
      regions_verified->begin_position = regions_in->begin_position;
      regions_verified->end_position = regions_in->end_position;
      ++regions_verified;
    }
    // DEBUG
    gem_cond_debug_block(DEBUG_VERIFY_REGIONS) {
      gem_slog("  => Region [%lu,%lu) (distance=%lu,matching-regions=%lu) Status=%s\n",
          regions_in->begin_position,regions_in->end_position,
          regions_in->align_distance,regions_in->match_scaffold.num_scaffold_regions,
          regions_in->status==filtering_region_discarded ? "Discarded" :
          (regions_in->status==filtering_region_accepted ? "Accepted" : "?") );
    }
  }
  // Return
  vector_update_used(filtering_candidates->filtering_regions,regions_out);
  vector_update_used(filtering_candidates->discarded_regions,regions_discarded);
  vector_update_used(filtering_candidates->verified_regions,regions_verified);
  return num_regions_accepted;
}
GEM_INLINE uint64_t filtering_candidates_verify_filtering_regions_multiple_hits(
    filtering_candidates_t* const filtering_candidates,text_collection_t* const text_collection,
    const pattern_t* const pattern,const as_parameters_t* const as_parameters,
    matches_t* const matches) {
  // Parameters
  search_parameters_t* const search_parameters = as_parameters->search_parameters;
  // Traverse all regions (text-space)
  const uint64_t num_filtering_regions = vector_get_used(filtering_candidates->filtering_regions);
  vector_reserve_additional(filtering_candidates->verified_regions,num_filtering_regions);
  vector_reserve_additional(filtering_candidates->discarded_regions,num_filtering_regions);
  verified_region_t* regions_verified = vector_get_free_elm(filtering_candidates->verified_regions,verified_region_t);
  filtering_region_t* regions_discarded = vector_get_free_elm(filtering_candidates->discarded_regions,filtering_region_t);
  filtering_region_t* regions_in = vector_get_mem(filtering_candidates->filtering_regions,filtering_region_t);
  uint64_t n, total_regions_accepted = 0, num_regions_accepted;
  gem_cond_debug_block(DEBUG_VERIFY_REGIONS) { gem_slog("[GEM]>Verify.Regions\n"); }
  for (n=0;n<num_filtering_regions;++n,++regions_in) {
    // Check region status (Skip other than unverified)
    if (regions_in->status!=filtering_region_unverified) {
      filtering_region_t* regions_out;
      vector_alloc_new(filtering_candidates->filtering_regions,filtering_region_t,regions_out);
      *regions_out = *regions_in;
    } else {
      // Verify region
      num_regions_accepted = filtering_region_verify_multiple_hits(filtering_candidates->filtering_regions,
          regions_in,text_collection,search_parameters,pattern);
      if (num_regions_accepted > 0) {
        total_regions_accepted += num_regions_accepted;
      } else {
        *regions_discarded = *regions_in;
        ++regions_discarded;
      }
      // Add to verify regions
      regions_verified->begin_position = regions_in->begin_position;
      regions_verified->end_position = regions_in->end_position;
      ++regions_verified;
    }
    // DEBUG
    gem_cond_debug_block(DEBUG_VERIFY_REGIONS) {
      if (num_regions_accepted > 0) {
        uint64_t i;
        for (i=total_regions_accepted-num_regions_accepted;i<total_regions_accepted;++i) {
          filtering_region_t* regions_it = vector_get_elm(filtering_candidates->filtering_regions,i,filtering_region_t);
          gem_slog("  => Region [%lu,%lu) (distance=%lu,matching-regions=%lu) Status=Accepted\n",
              regions_it->begin_position,regions_it->end_position,
              regions_it->align_distance,regions_it->match_scaffold.num_scaffold_regions);
        }
      } else {
        gem_slog("  => Region [%lu,%lu) Status=Discarded\n",regions_in->begin_position,regions_in->end_position);
      }
      num_regions_accepted = 0;
    }
  }
  // Return
  vector_update_used(filtering_candidates->discarded_regions,regions_discarded);
  vector_update_used(filtering_candidates->verified_regions,regions_verified);
  return total_regions_accepted;
}
/*
 * Retrieve all candidates(text) from the index
 */
GEM_INLINE void filtering_candidates_retrieve_filtering_regions(
    filtering_candidates_t* const filtering_candidates,archive_text_t* const archive_text,
    text_collection_t* const text_collection,mm_stack_t* const mm_stack) {
  // Traverse all candidates (text-space)
  VECTOR_ITERATE(filtering_candidates->filtering_regions,filtering_region,candidate_pos,filtering_region_t) {
    // Retrieve only those marked as unverified to be verified
    if (filtering_region->status!=filtering_region_unverified) continue;
    // Retrieve only not retrieved yet
    if (filtering_region->text_trace_offset!=UINT64_MAX) continue;
    // Retrieve text(s)
    const uint64_t text_length = filtering_region->end_position - filtering_region->begin_position;
    filtering_region->text_trace_offset =
        archive_text_retrieve(archive_text,text_collection,
            filtering_region->begin_position,text_length,false,mm_stack);
  }
}
/*
 * Compose filtering regions
 */
GEM_INLINE void filtering_candidates_compose_matching_regions(
    filtering_candidates_t* const filtering_candidates,const uint64_t first_candidate_idx,
    const uint64_t last_candidate_idx,const bool compose_region_chaining,mm_stack_t* const mm_stack) {
  // Fetch candidate-positions
  filtering_position_t* const candidate_positions = vector_get_mem(filtering_candidates->filtering_positions,filtering_position_t);
  // Allow new matching candidate-region
  filtering_region_t* filtering_region;
  vector_alloc_new(filtering_candidates->filtering_regions,filtering_region_t,filtering_region);
  filtering_region->status = filtering_region_unverified; // Newly created region (unverified)
  filtering_region->text_trace_offset = UINT64_MAX;
  filtering_region->begin_position = candidate_positions[first_candidate_idx].begin_position;
  filtering_region->end_position = candidate_positions[last_candidate_idx].end_position;
  filtering_region->base_position_offset = candidate_positions[first_candidate_idx].base_position_offset;
  PROF_ADD_COUNTER(GP_CANDIDATE_REGION_LENGTH,filtering_region->end_position-filtering_region->begin_position);
  match_scaffold_t* const matches_scaffold = &filtering_region->match_scaffold;
  match_scaffold_init(matches_scaffold);
  // Compose regions matching
  if (compose_region_chaining) {
    const uint64_t num_regions_matching = last_candidate_idx-first_candidate_idx+1;
    matches_scaffold->scaffold_regions = mm_stack_calloc(mm_stack,num_regions_matching,region_matching_t,false);
    matches_scaffold->num_scaffold_regions = num_regions_matching;
    matches_scaffold->scaffolding_coverage = 0;
    uint64_t i;
    for (i=0;i<num_regions_matching;++i) {
      region_matching_t* const region_matching = matches_scaffold->scaffold_regions + i;
      filtering_position_t* const candidate_position = candidate_positions + first_candidate_idx + i;
      region_search_t* const source_region = vector_get_elm(
          filtering_candidates->regions_buffer,candidate_position->source_region_offset,region_search_t);
      // Region error
      region_matching->matching_type = (source_region->degree==0) ? region_matching_exact : region_matching_approximate;
      region_matching->error = source_region->degree;
      region_matching->cigar_length = 0;
      // Read coordinates
      region_matching->key_begin = source_region->begin;
      region_matching->key_end = source_region->end;
      // Text coordinates (relative to the effective begin position)
      const uint64_t region_length = region_matching->key_end - region_matching->key_begin;
      region_matching->text_begin = candidate_position->region_text_position - filtering_region->begin_position;
      region_matching->text_end = region_matching->text_begin + region_length;
      // Coverage
      matches_scaffold->scaffolding_coverage += region_length;
    }
  }
  // DEBUG
  gem_cond_debug_block(DEBUG_FILTERING_REGIONS) {
    filtering_region_t* filtering_region = vector_get_last_elm(filtering_candidates->filtering_regions,filtering_region_t);
    gem_slog("  Region-found => [%lu,%lu)=%lu bases (matching-regions=%lu)\n",
        filtering_region->begin_position,filtering_region->end_position,
        filtering_region->end_position-filtering_region->begin_position,
        filtering_region->match_scaffold.num_scaffold_regions);
  }
}
GEM_INLINE uint64_t filtering_candidates_compose_filtering_regions(
    filtering_candidates_t* const filtering_candidates,const uint64_t key_length,
    const uint64_t max_delta_difference,const bool compose_region_chaining,
    mm_stack_t* const mm_stack) {
  // Sort candidate positions (text-space)
  filtering_positions_sort_positions(filtering_candidates->filtering_positions);
  // Sort verified regions
  verified_regions_sort_positions(filtering_candidates->verified_regions);
  const uint64_t num_verified_regions = vector_get_used(filtering_candidates->verified_regions);
  verified_region_t* const verified_region = vector_get_mem(filtering_candidates->verified_regions,verified_region_t);
  // Traverse positions and eliminate duplicates
  const uint64_t num_candidate_positions = vector_get_used(filtering_candidates->filtering_positions);
  filtering_position_t* const candidate_positions = vector_get_mem(filtering_candidates->filtering_positions,filtering_position_t);
  uint64_t candidate_idx = 0, verified_region_idx = 0;
  gem_cond_debug_block(DEBUG_FILTERING_REGIONS) { gem_slog("[GEM]>Composing.Regions\n"); }
  while (candidate_idx < num_candidate_positions) {
    // Determine the positions belonging to the same region
    const uint64_t region_begin_position = candidate_positions[candidate_idx].begin_position;
    uint64_t last_position = region_begin_position;
    uint64_t group_idx = candidate_idx + 1;
    while (group_idx < num_candidate_positions) {
      const uint64_t position = candidate_positions[group_idx].begin_position;
      const uint64_t delta = position - last_position;
      if (delta > max_delta_difference) break; // Doesn't belong to the group. Stop!
      ++group_idx; // Next
    }
    // Check region against verified regions
    bool is_already_verified = false;
    while (verified_region_idx < num_verified_regions && verified_region[verified_region_idx].end_position <= region_begin_position) {
      ++verified_region_idx;
    }
    if (verified_region_idx < num_verified_regions) {
      const uint64_t verified_begin_position = verified_region[verified_region_idx].begin_position;
      const uint64_t verified_end_position = verified_region[verified_region_idx].end_position;
      const uint64_t region_end_position = candidate_positions[group_idx-1].end_position;
      is_already_verified = (verified_begin_position <= region_begin_position && region_end_position <= verified_end_position);
      gem_debug_block() {
        if (is_already_verified) {
          PROF_INC_COUNTER(GP_CANDIDATE_REGIONS_DUPLICATED);
        }
      }
    }
    if (!is_already_verified) {
      // Create a region candidate with the positions from [candidate_idx] to [group_idx-1]
      filtering_candidates_compose_matching_regions(filtering_candidates,
          candidate_idx,group_idx-1,compose_region_chaining,mm_stack);
    }
    // Next group
    const uint64_t num_regions_matching = group_idx-candidate_idx;
    candidate_idx += num_regions_matching;
  }
  // Clear candidate positions
  vector_clear(filtering_candidates->filtering_positions);
  // DEBUG
  gem_cond_debug_block(DEBUG_REGIONS_MATCHING) {
    filtering_candidates_print_matching_regions(stderr,filtering_candidates);
  }
  // Return number of filtering regions generated
  return vector_get_used(filtering_candidates->filtering_regions);
}
/*
 * Filtering adjustment of the position wrt region/seed on which the candidate is based
 */
GEM_INLINE void filtering_candidates_adjust_filtering_position(
    filtering_position_t* const filtering_position,archive_text_t* const archive_text,
    const uint64_t begin_offset,const uint64_t end_offset,const uint64_t boundary_error) {
  // Decode Sampled RL-Position (RL-text encoded)
  uint64_t region_text_position = filtering_position->region_text_position;
  uint64_t rl_increased_error = 0;
  if (archive_text->run_length) {
    const uint64_t sampled_idx = region_text_position / SAMPLED_RL_SAMPLING_RATE;
    rl_increased_error = (region_text_position % SAMPLED_RL_SAMPLING_RATE) * SAMPLED_RL_MAX_RUN_LENGTH;
    region_text_position = sampled_rl_get_sample(archive_text->sampled_rl,sampled_idx);
  }
  // Adjust Position (project to beginning of the candidate)
  const locator_interval_t* const locator_interval = filtering_position->locator_interval;
  const uint64_t position_offset = begin_offset + rl_increased_error;
  uint64_t base_begin_position = (region_text_position > position_offset) ? (region_text_position - position_offset) : 0;
  // Adjust Begin Position
  uint64_t begin_position;
  if (base_begin_position < locator_interval->begin_position) { // Adjust by locator-interval
    base_begin_position = locator_interval->begin_position; // Possible trim at the beginning
    begin_position = locator_interval->begin_position;
  } else {
    begin_position = (base_begin_position > boundary_error) ? base_begin_position-boundary_error : 0;
    if (begin_position < locator_interval->begin_position) { // Adjust by locator-interval
      begin_position = locator_interval->begin_position;
    }
  }
  // Adjust End Position
  uint64_t end_position = region_text_position + end_offset + boundary_error;
  if (end_position >= locator_interval->end_position) { // Adjust by locator-interval
    end_position = locator_interval->end_position; // Possible trim at the end
  }
  filtering_position->begin_position = begin_position;
  filtering_position->end_position = end_position;
  filtering_position->base_position_offset = base_begin_position-begin_position;
}
GEM_INLINE void filtering_candidates_compute_extension_region(
    const locator_t* const locator,const bool extension_onward,
    const uint64_t extended_begin_position,const uint64_t extended_effective_length,
    const uint64_t candidate_key_length,const uint64_t max_filtering_error,
    const uint64_t max_template_size,uint64_t* const candidate_begin_position,
    uint64_t* const candidate_end_position) {
  locator_interval_t* const locator_interval = locator_lookup_interval(locator,extended_begin_position);
  if (extension_onward) {
    *candidate_begin_position = BOUNDED_SUBTRACTION(extended_begin_position,max_filtering_error,locator_interval->begin_position);
    const uint64_t end_offset = extended_effective_length + max_template_size + candidate_key_length + max_filtering_error;
    *candidate_end_position = BOUNDED_ADDITION(extended_begin_position,end_offset,locator_interval->end_position);
  } else {
    *candidate_end_position = extended_begin_position + extended_effective_length + max_filtering_error;
    if (*candidate_end_position > locator_interval->end_position) {
      *candidate_end_position = locator_interval->end_position;
    }
    const uint64_t begin_offset = max_template_size + candidate_key_length + max_filtering_error;
    *candidate_begin_position = BOUNDED_SUBTRACTION(extended_begin_position,begin_offset,locator_interval->begin_position);
  }
}
/*
 * Decode of all candidate positions (index-space -> text-space)
 */
GEM_INLINE void filtering_candidates_decode_filtering_positions(
    filtering_candidates_t* const filtering_candidates,
    const locator_t* const locator,archive_text_t* const archive_text,
    const fm_index_t* const fm_index,const uint64_t key_length,const uint64_t boundary_error) {
  // Traverse all candidate positions in index-space
  vector_t* const candidate_text_positions = filtering_candidates->filtering_positions;
  VECTOR_ITERATE(candidate_text_positions,filtering_position,n,filtering_position_t) {
    // Lookup Position
    filtering_position->region_text_position = fm_index_lookup(fm_index,filtering_position->region_index_position);
    // Locate Position
    filtering_position->locator_interval = locator_lookup_interval(locator,filtering_position->region_text_position);
    // Adjust Position
    region_search_t* const source_region = vector_get_elm(filtering_candidates->regions_buffer,
        filtering_position->source_region_offset,region_search_t);
    filtering_candidates_adjust_filtering_position(filtering_position,archive_text,
        source_region->begin,key_length-source_region->begin,boundary_error);
  }
}
/*
 * Batch decode of all candidate positions (index-space -> text-space)
 *   (All the steps (CSA-lookup, rankQueries) are performed with prefetch-loops)
 */
typedef struct {
  uint64_t vector_rank;
  uint64_t index_position;
  uint64_t distance;
  uint64_t used_slot;
  bwt_block_locator_t bwt_block_locator;
} fc_batch_decode_candidate;
GEM_INLINE void filtering_candidates_decode_filtering_position_batch_prefetched(
    filtering_candidates_t* const filtering_candidates,const locator_t* const locator,
    archive_text_t* const archive_text,const fm_index_t* const fm_index,
    const uint64_t key_length,const uint64_t boundary_error) {
  // Init
  const bwt_t* const bwt = fm_index->bwt;
  const sampled_sa_t* const sampled_sa = fm_index->sampled_sa;
  vector_t* const candidate_text_positions = filtering_candidates->filtering_positions;
  const uint64_t num_candidate_text_positions = vector_get_used(candidate_text_positions);
  filtering_position_t* const filtering_position = vector_get_mem(candidate_text_positions,filtering_position_t);
  // Initial fill batch
  fc_batch_decode_candidate batch[DECODE_NUM_POSITIONS_PREFETCHED]; // Batch Decode
  bool is_sampled;
  uint64_t current_position, LF, i;
  for (i=0,current_position=0;i<DECODE_NUM_POSITIONS_PREFETCHED && current_position<num_candidate_text_positions;++current_position) {
    LF = bwt_LF(bwt,filtering_position[current_position].region_index_position,&is_sampled);
    if (!is_sampled) {
      batch[i].index_position = LF;
      batch[i].vector_rank = current_position;
      batch[i].distance = 0;
      batch[i].used_slot = true;
      ++i;
    } else {
      filtering_position[current_position].decode_sampled_pos = LF;
      filtering_position[current_position].decode_distance = 0; PROF_ADD_COUNTER(GP_FMIDX_LOOKUP_DIST,0);
    }
  }
  const bool full_filled_batch = (i==DECODE_NUM_POSITIONS_PREFETCHED);
  for (;i<DECODE_NUM_POSITIONS_PREFETCHED;++i) {
    batch[i].used_slot = false;
  }
  // Full-prefetch loop for sampled-LF
  if (full_filled_batch) {
    while (current_position<num_candidate_text_positions) {
      for (i=0;i<DECODE_NUM_POSITIONS_PREFETCHED;++i) {
        bwt_prefetch(bwt,batch[i].index_position,&(batch[i].bwt_block_locator));
      }
      for (i=0;i<DECODE_NUM_POSITIONS_PREFETCHED;++i) {
        ++(batch[i].distance);
        batch[i].index_position = bwt_prefetched_LF(bwt,batch[i].index_position,&is_sampled,&(batch[i].bwt_block_locator));
        if (is_sampled) {
          filtering_position[batch[i].vector_rank].decode_sampled_pos = batch[i].index_position;
          filtering_position[batch[i].vector_rank].decode_distance = batch[i].distance;
          PROF_ADD_COUNTER(GP_FMIDX_LOOKUP_DIST,batch[i].distance);
          batch[i].used_slot = false;
          // Select new candidate to decode
          while (current_position < num_candidate_text_positions) {
            LF = bwt_LF(bwt,filtering_position[current_position].region_index_position,&is_sampled);
            if (!is_sampled) break;
            filtering_position[current_position].decode_sampled_pos = LF;
            filtering_position[current_position].decode_distance = 0; PROF_ADD_COUNTER(GP_FMIDX_LOOKUP_DIST,0);
            ++current_position;
          }
          if (current_position < num_candidate_text_positions) {
            batch[i].index_position = LF;
            batch[i].vector_rank = current_position;
            batch[i].distance = 0;
            batch[i].used_slot = true;
            ++current_position;
          }
        }
      }
    }
  }
  // Solve remaining queries
  for (i=0;i<DECODE_NUM_POSITIONS_PREFETCHED;++i) {
    if (batch[i].used_slot) {
      do {
        ++(batch[i].distance);
        batch[i].index_position = bwt_LF(bwt,batch[i].index_position,&is_sampled);
      } while (!is_sampled);
      filtering_position[batch[i].vector_rank].decode_sampled_pos = batch[i].index_position;
      filtering_position[batch[i].vector_rank].decode_distance = batch[i].distance;
      PROF_ADD_COUNTER(GP_FMIDX_LOOKUP_DIST,batch[i].distance);
    }
  }
  // Prefetch SA-retrieve samples
  const uint64_t bwt_length = fm_index_get_length(fm_index);
  uint64_t num_left_positions = num_candidate_text_positions;
  current_position = 0;
  while (num_left_positions > 0) {
    const uint64_t batch_size = MIN(num_left_positions,RETRIEVE_SAMPLE_NUM_POSITIONS_PREFETCHED);
    const uint64_t batch_top = current_position+batch_size;
    for (i=current_position;i<batch_top;++i) {
      sampled_sa_prefetch_sample(sampled_sa,filtering_position[i].decode_sampled_pos);
    }
    for (i=current_position;i<batch_top;++i) {
      filtering_position[i].region_text_position =
          (sampled_sa_get_sample(sampled_sa,filtering_position[i].decode_sampled_pos) +
              filtering_position[i].decode_distance) % bwt_length;
    }
    current_position = batch_top;
    num_left_positions -= batch_size;
  }
  // Adjust decoded position to the beginning of the read
  for (current_position=0;current_position<num_candidate_text_positions;++current_position) {
    // Locate Position
    filtering_position[current_position].locator_interval =
        locator_lookup_interval(locator,filtering_position[current_position].region_text_position);
    // Adjust Position
    region_search_t* const source_region = vector_get_elm(filtering_candidates->regions_buffer,
        filtering_position[current_position].source_region_offset,region_search_t);
    filtering_candidates_adjust_filtering_position(
        filtering_position+current_position,archive_text,
        source_region->begin,key_length-source_region->begin,boundary_error);
  }
}
/*
 * Processing & Verification
 */
GEM_INLINE uint64_t filtering_candidates_process_candidates(
    filtering_candidates_t* const filtering_candidates,
    archive_t* const archive,const pattern_t* const pattern,
    const as_parameters_t* const as_parameters,
    const bool compose_region_chaining,mm_stack_t* const mm_stack) {
  // Check non-empty pending candidates set
  uint64_t pending_candidates = vector_get_used(filtering_candidates->filtering_positions);
  PROF_ADD_COUNTER(GP_CANDIDATE_POSITIONS,pending_candidates);
  if (pending_candidates==0) return 0; // Nothing to do

  // Batch decode+adjust of all positions of the candidates (cip_begin_position = decoded(candidate_region_index_position))
  PROF_START(GP_FC_PROCESS_CANDIDATES);
  PROF_START(GP_FC_DECODE_POSITIONS);
  const uint64_t key_length = pattern->key_length;
  if (pending_candidates < DECODE_NUM_POSITIONS_PREFETCHED) {
    filtering_candidates_decode_filtering_positions(
        filtering_candidates,archive->locator,archive->text,
        archive->fm_index,key_length,pattern->max_effective_bandwidth);
  } else {
    filtering_candidates_decode_filtering_position_batch_prefetched(
        filtering_candidates,archive->locator,archive->text,
        archive->fm_index,key_length,pattern->max_effective_bandwidth);
  }
  PROF_STOP(GP_FC_DECODE_POSITIONS);

  // Compose matching regions into candidate regions (also filter out duplicated positions or already checked)
  PROF_START(GP_FC_COMPOSE_REGIONS);
  search_parameters_t* const search_parameters = as_parameters->search_parameters;
  pending_candidates = filtering_candidates_compose_filtering_regions(filtering_candidates,key_length,
      pattern->max_effective_bandwidth,compose_region_chaining && search_parameters->alignment_scaffolding,mm_stack);
  PROF_STOP(GP_FC_COMPOSE_REGIONS);
  PROF_ADD_COUNTER(GP_CANDIDATE_REGIONS,pending_candidates);

  // Return total candidate regions
  PROF_STOP(GP_FC_PROCESS_CANDIDATES);
  return pending_candidates;
}
GEM_INLINE uint64_t filtering_candidates_verify_candidates(
    filtering_candidates_t* const filtering_candidates,archive_t* const archive,
    text_collection_t* const text_collection,const pattern_t* const pattern,
    const as_parameters_t* const as_parameters,
    matches_t* const matches,mm_stack_t* const mm_stack) {
  // Check number of filtering regions
  uint64_t pending_candidates = vector_get_used(filtering_candidates->filtering_regions);
  if (pending_candidates==0) return 0;
  PROF_START(GP_FC_VERIFICATION);
  // Retrieve text-candidates
  PROF_START(GP_FC_RETRIEVE_CANDIDATE_REGIONS);
  filtering_candidates_retrieve_filtering_regions(filtering_candidates,archive->text,text_collection,mm_stack);
  PROF_STOP(GP_FC_RETRIEVE_CANDIDATE_REGIONS);
  // Verify candidates
  PROF_START(GP_FC_VERIFY_CANDIDATE_REGIONS);
  pending_candidates = filtering_candidates_verify_filtering_regions(
      filtering_candidates,text_collection,pattern,as_parameters,matches);
  PROF_STOP(GP_FC_VERIFY_CANDIDATE_REGIONS);
  PROF_STOP(GP_FC_VERIFICATION);
  return pending_candidates;
}
GEM_INLINE uint64_t filtering_candidates_align_candidates(
    filtering_candidates_t* const filtering_candidates,
    archive_text_t* const archive_text,const locator_t* const locator,
    text_collection_t* const text_collection,pattern_t* const pattern,
    const bool emulated_rc_search,const as_parameters_t* const as_parameters,
    const bool approximated_distance,matches_t* const matches,mm_stack_t* const mm_stack) {
  PROF_START(GP_FC_REALIGN_CANDIDATE_REGIONS);
  // Hint to matches
  const uint64_t num_filtering_regions = vector_get_used(filtering_candidates->filtering_regions);
  if (num_filtering_regions==0) return 0;
  matches_hint_allocate_match_trace(matches,num_filtering_regions);
  // Realign
  const uint64_t aligned_regions = filtering_candidates_align_accepted_regions(
      filtering_candidates,archive_text,locator,text_collection,pattern,emulated_rc_search,
      as_parameters,approximated_distance,matches,mm_stack);
  PROF_STOP(GP_FC_REALIGN_CANDIDATE_REGIONS);
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
  PROF_START(GP_FC_UNBOUNDED_ALIGNMENT);
  // Retrieve text-candidates
  PROF_START(GP_FC_RETRIEVE_CANDIDATE_REGIONS);
  filtering_candidates_retrieve_filtering_regions(filtering_candidates,archive_text,text_collection,mm_stack);
  PROF_STOP(GP_FC_RETRIEVE_CANDIDATE_REGIONS);
  // Unbounded-align all discarded candidates
  matches_hint_allocate_match_trace(matches,discarded_candidates);
  const uint64_t aligned_regions = filtering_candidates_unbounded_align_discarded_regions(
      filtering_candidates,archive_text,locator,text_collection,pattern,
      emulated_rc_search,as_parameters,matches,mm_stack);
  PROF_STOP(GP_FC_UNBOUNDED_ALIGNMENT);
  return aligned_regions;
}
/*
 * BPM-Buffer API (Verification)
 */
GEM_INLINE uint64_t filtering_candidates_bpm_buffer_add(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,bpm_gpu_buffer_t* const bpm_gpu_buffer) {
  // Check number of pending regions
  const uint64_t pending_candidates = filtering_candidates_get_num_candidate_regions(filtering_candidates);
  if (pending_candidates==0) return 0;
  // Add the pattern(chunks) to the buffer (add new queries)
  bpm_gpu_buffer_put_pattern(bpm_gpu_buffer,pattern);
  // Fetch pattern dimensions
  const uint64_t max_error = pattern->max_effective_filtering_error;
  const uint64_t pattern_length = pattern->key_length;
  const uint64_t num_pattern_chunks = pattern->bpm_pattern.gpu_num_chunks;
  const uint64_t pattern_tile_tall = pattern->bpm_pattern.gpu_entries_per_chunk*BPM_GPU_PATTERN_ENTRY_LENGTH;
  // Traverse all candidates (text-space) & add them to the buffer
  const filtering_region_t* candidate_region = vector_get_mem(filtering_candidates->filtering_regions,filtering_region_t);
  uint64_t candidate_pos, pattern_chunk, total_candidates_added;
  for (candidate_pos=0,total_candidates_added=0;candidate_pos<pending_candidates;++candidate_pos,++candidate_region) {
    // Locate candidate sequence
    const uint64_t begin_position = candidate_region->begin_position;
    const uint64_t candidate_length = candidate_region->end_position - begin_position;
    // Calculate tile dimensions
    pattern_tiled_t pattern_tiled;
    const bool pattern_can_align = pattern_tiled_init(&pattern_tiled,
        pattern_length,pattern_tile_tall,candidate_length,max_error);
    if (!pattern_can_align) continue;
    // Initialize current tile variables
    for (pattern_chunk=0;pattern_chunk<num_pattern_chunks;++pattern_chunk,++total_candidates_added) {
      // BPM-GPU put candidate
      bpm_gpu_buffer_put_candidate(bpm_gpu_buffer,
          begin_position+pattern_tiled.tile_offset,pattern_tiled.tile_wide,pattern_chunk);
      // Calculate next tile
      pattern_tiled_calculate_next(&pattern_tiled);
    }
  }
  // Return the final number of chunk-candidates added to the buffer
  PROF_ADD_COUNTER(GP_BMP_TILED_NUM_TILES,total_candidates_added);
  PROF_ADD_COUNTER(GP_BMP_TILED_NUM_TILES_VERIFIED,total_candidates_added);
  return total_candidates_added;
}
GEM_INLINE uint64_t filtering_candidates_bpm_buffer_retrieve(
    filtering_candidates_t* const filtering_candidates,archive_text_t* const archive_text,
    text_collection_t* const text_collection,pattern_t* const pattern,
    bpm_gpu_buffer_t* const bpm_gpu_buffer,const uint64_t candidate_offset_begin,
    const uint64_t candidate_offset_end,mm_stack_t* const mm_stack) {
  /*
   * Retrieve filtering-regions from BPM-Buffer
   */
  PROF_START(GP_FC_RETRIEVE_BPM_BUFFER_CANDIDATE_REGIONS);
  // Clear filtering candidates
  filtering_candidates_clear(filtering_candidates);
  if (gem_expect_false(candidate_offset_begin==candidate_offset_end)) return 0;
  // Fetch Parameters
  const uint64_t key_length = pattern->key_length;
  const uint64_t max_error = pattern->max_effective_filtering_error;
  // Fetch tile dimensions
  const uint64_t num_chunks = pattern->bpm_pattern.gpu_num_chunks;
  const uint64_t pattern_tile_tall = pattern->bpm_pattern.gpu_entries_per_chunk*BPM_GPU_PATTERN_ENTRY_LENGTH;
  // Prepare filtering-regions vectors
  const uint64_t pending_candidates = (candidate_offset_end-candidate_offset_begin)/num_chunks;
  vector_reserve(filtering_candidates->verified_regions,pending_candidates,false);
  vector_reserve(filtering_candidates->filtering_regions,pending_candidates,false);
  vector_reserve(filtering_candidates->discarded_regions,pending_candidates,false);
  filtering_region_t* regions_accepted = vector_get_mem(filtering_candidates->filtering_regions,filtering_region_t);
  filtering_region_t* regions_discarded = vector_get_mem(filtering_candidates->discarded_regions,filtering_region_t);
  verified_region_t* regions_verified = vector_get_mem(filtering_candidates->verified_regions,verified_region_t);
  // Traverse all candidates (text-space) & sum-up their alignment distance
  uint64_t num_accepted_regions = 0;
  uint64_t candidate_idx=candidate_offset_begin, pattern_chunk, candidate_pos;
  gem_cond_debug_block(DEBUG_VERIFY_REGIONS) { gem_slog("[GEM]>Verify.Regions\n"); }
  for (candidate_pos=0;candidate_pos<pending_candidates;++candidate_pos) {
    // Get the accepted candidate
    uint64_t candidate_begin_position, candidate_end_position, candidate_length;
    uint32_t chunk_length;
    bpm_gpu_buffer_get_candidate(bpm_gpu_buffer,candidate_idx,&candidate_begin_position,&chunk_length);
    bpm_gpu_buffer_get_candidate(bpm_gpu_buffer,candidate_idx+(num_chunks-1),&candidate_end_position,&chunk_length);
    candidate_end_position += chunk_length;
    candidate_length = candidate_end_position - candidate_begin_position;
    // Calculate tile dimensions (Sum up the alignment distance of all the tiles)
    pattern_tiled_t pattern_tiled;
    const bool pattern_can_align = pattern_tiled_init(&pattern_tiled,key_length,pattern_tile_tall,candidate_length,max_error);
    if (!pattern_can_align) continue;
    bool unaligned_tiled = false;
    uint64_t global_distance = 0, distance_link_tiles = 0;
    for (pattern_chunk=0;pattern_chunk<num_chunks;++pattern_chunk,++candidate_idx) {
      // Retrieve alignment distance
      uint32_t tile_distance=0, tile_match_column=0;
      bpm_gpu_buffer_get_candidate_result(bpm_gpu_buffer,candidate_idx,&tile_distance,&tile_match_column);
      pattern_tiled.tile_distance = tile_distance;
      pattern_tiled.tile_match_column = tile_match_column;
      global_distance += tile_distance;
      if (tile_distance > max_error) unaligned_tiled = true;
      // Bound the joint distance by estimating the cost of connecting the path through the tiles
      distance_link_tiles += pattern_tiled_bound_matching_path(&pattern_tiled);
      // Calculate next tile
      pattern_tiled_calculate_next(&pattern_tiled);
    }
    // Check total distance & Compose the retrieved region
    if (!unaligned_tiled && global_distance <= max_error) {
      // Configure accepted candidate
      regions_accepted->status = filtering_region_accepted;
      regions_accepted->text_trace_offset = UINT64_MAX; // Not retrieved yet
      regions_accepted->begin_position = candidate_begin_position;
      regions_accepted->end_position = candidate_end_position;
      regions_accepted->base_position_offset = 0;
      // Configure regions matching (we sacrifice this information as to save memory)
      match_scaffold_init(&regions_accepted->match_scaffold);
      // Distance Bound estimation
      regions_accepted->align_distance_min_bound = (global_distance > max_error) ? max_error : global_distance;
      global_distance += distance_link_tiles;
      regions_accepted->align_distance = (global_distance > max_error) ? max_error : global_distance;
      regions_accepted->align_match_begin_column = 0;
      regions_accepted->align_match_end_column =
          BOUNDED_ADDITION(pattern_tiled.prev_tile_match_position,regions_accepted->align_distance,candidate_length-1);
      ++regions_accepted; ++num_accepted_regions;
      gem_cond_debug_block(DEBUG_VERIFY_REGIONS) {
        gem_slog("  => Region [%lu,%lu) (dist=%lu,regMatch=none) Status=OK\n",
            regions_accepted->begin_position,regions_accepted->end_position,regions_accepted->align_distance);
      }
    } else {
      // Configure discarded candidate
      regions_discarded->status = filtering_region_discarded;
      regions_discarded->text_trace_offset = UINT64_MAX; // Not retrieved yet
      regions_discarded->begin_position = candidate_begin_position;
      regions_discarded->end_position = candidate_end_position;
      // Distance Bound estimation
      regions_discarded->align_distance_min_bound = global_distance;
      regions_discarded->align_distance = global_distance + distance_link_tiles;
      // Configure regions matching (we sacrifice this information as to save memory)
      match_scaffold_init(&regions_discarded->match_scaffold);
      ++regions_discarded;
      gem_cond_debug_block(DEBUG_VERIFY_REGIONS) {
        gem_slog("  => Region [%lu,%lu) (dist=inf,regMatch=none) Status=DISCARDED\n",
            regions_discarded->begin_position,regions_discarded->end_position);
      }
    }
    // Add to verified regions
    regions_verified->begin_position = candidate_begin_position;
    regions_verified->end_position = candidate_end_position;
    ++regions_verified;
  }
  // Update Accepted/Discarded/Verified
  vector_update_used(filtering_candidates->verified_regions,regions_verified);
  vector_update_used(filtering_candidates->filtering_regions,regions_accepted);
  vector_update_used(filtering_candidates->discarded_regions,regions_discarded);
  PROF_STOP(GP_FC_RETRIEVE_BPM_BUFFER_CANDIDATE_REGIONS);
  // Return number of accepted regions
  return num_accepted_regions;
}
/*
 * Pair Extension
 */
GEM_INLINE uint64_t filtering_candidates_extend_match(
    filtering_candidates_t* const filtering_candidates,
    archive_text_t* const archive_text,const locator_t* const locator,
    text_collection_t* const text_collection,const match_trace_t* const extended_match,
    pattern_t* const candidate_pattern,const as_parameters_t* const candidate_actual_parameters,
    mapper_stats_t* const mapper_stats,paired_matches_t* const paired_matches,
    const sequence_end_t candidate_end,mm_stack_t* const mm_stack) {
  PROF_START(GP_FC_EXTEND_MATCH);
  // Parameters
  search_parameters_t* const search_parameters = candidate_actual_parameters->search_parameters;
  const uint64_t max_filtering_error = candidate_pattern->max_effective_filtering_error;
  /*
   * Retrieve text-candidate
   */
  PROF_START(GP_FC_EXTEND_RETRIEVE_CANDIDATE_REGIONS);
  // Compute region on the opposite strand
  const uint64_t extended_effective_length = extended_match->match_alignment.effective_length;
  const uint64_t extended_match_position = archive_text_get_projection(archive_text,
      extended_match->match_alignment.match_position,extended_effective_length);
  const bool search_onward = false;
  const uint64_t max_template_size = mapper_stats_template_length_get_expected_max(mapper_stats);
  uint64_t candidate_begin_position, candidate_end_position;
  filtering_candidates_compute_extension_region(
      locator,search_onward,extended_match_position,extended_effective_length,
      candidate_pattern->key_length,max_filtering_error,max_template_size,
      &candidate_begin_position,&candidate_end_position);
  const uint64_t candidate_length = candidate_end_position-candidate_begin_position;
  const uint64_t text_trace_offset = archive_text_retrieve(archive_text,text_collection,
      candidate_begin_position,candidate_length,false,mm_stack);
  PROF_STOP(GP_FC_EXTEND_RETRIEVE_CANDIDATE_REGIONS);
  /*
   * Verify candidate region (may contain multiple matches)
   */
  PROF_START(GP_FC_EXTEND_VERIFY_CANDIDATE_REGIONS);
  uint64_t candidates_found = filtering_region_verify_extension(
      filtering_candidates->filtering_regions,filtering_candidates->verified_regions,
      text_collection,text_trace_offset,candidate_begin_position,search_parameters,candidate_pattern);
  PROF_ADD_COUNTER(GP_FC_EXTEND_VERIFY_CANDIDATE_LENGTH,candidate_length);
  PROF_STOP(GP_FC_EXTEND_VERIFY_CANDIDATE_REGIONS);
  if (candidates_found==0) { PROF_STOP(GP_FC_EXTEND_MATCH); return 0; }
  /*
   * Align accepted candidates
   */
  PROF_START(GP_FC_EXTEND_REALIGN_CANDIDATE_REGIONS);
  matches_t* matches_candidate = (candidate_end==paired_end1) ? paired_matches->matches_end1 : paired_matches->matches_end2;
  matches_hint_allocate_match_trace(matches_candidate,candidates_found); // Hint to matches
  // Align
  candidates_found = filtering_candidates_align_accepted_regions(filtering_candidates,
      archive_text,locator,text_collection,candidate_pattern,false,
      candidate_actual_parameters,false,matches_candidate,mm_stack);
  PROF_STOP(GP_FC_EXTEND_REALIGN_CANDIDATE_REGIONS);
  PROF_STOP(GP_FC_EXTEND_MATCH);
  // Return number of extended-matches found
  return candidates_found;
}
#ifdef SPLIT_EXTENSION_WINDOW
GEM_INLINE void filtering_candidates_process_extension_candidates(
    filtering_candidates_t* const extended_filtering_candidates,
    filtering_candidates_t* const candidate_filtering_candidates,
    archive_t* const archive,text_collection_t* const text_collection,
    const pattern_t* const extended_pattern,const pattern_t* const candidate_pattern,
    const search_parameters_t* const search_parameters,mapper_stats_t* const mapper_stats,
    paired_matches_t* const paired_matches,mm_stack_t* const mm_stack) {
  // Parameters
  const uint64_t max_template_size = mapper_stats_template_length_get_expected_max(mapper_stats);
  // Allocate candidate filtering regions
  const uint64_t num_filtering_regions = vector_get_used(extended_filtering_candidates->filtering_regions);
  filtering_region_t* regions_extended = vector_get_mem(extended_filtering_candidates->filtering_regions,filtering_region_t);
  // Traverse extended-regions and generate candidate-regions
  // Split the extension candidates in chunks and add them to the candidate filtering-regions
  uint64_t n;
  for (n=0;n<num_filtering_regions;++n,++regions_extended) {
    // Compute candidate region boundaries
    const uint64_t extended_eff_begin_position = regions_extended->begin_position;
    uint64_t candidate_begin_position, candidate_end_position;
    filtering_candidates_compute_extension_region(
        archive->locator,true,extended_eff_begin_position,
        extended_pattern->key_length+extended_pattern->max_effective_filtering_error,
        candidate_pattern->key_length,candidate_pattern->max_effective_filtering_error,
        max_template_size,&candidate_begin_position,&candidate_end_position);
    const uint64_t candidate_length = candidate_end_position-candidate_begin_position;
    // Add candidate region
    const uint64_t max_error = candidate_pattern->max_effective_filtering_error;
    const uint64_t pattern_length = candidate_pattern->key_length;
    const uint64_t window_length = pattern_length + 2*max_error;
    const uint64_t overlap = pattern_length + max_error;
    const uint64_t begin_position = archive_text_get_projection(archive->text,candidate_begin_position,candidate_length);
    const uint64_t end_position = begin_position + candidate_length;
    uint64_t current_begin_position = begin_position, current_end_position = begin_position + window_length;
    // Init template
    filtering_region_t regions_candidate;
    regions_candidate.status = filtering_region_unverified; // Newly created region (unverified)
    match_scaffold_init(&regions_candidate.match_scaffold);
    uint64_t num_chunks_added = 0;
    while (current_end_position <= end_position) {
      // Add overlapping candidate
      if (current_begin_position != begin_position) {
        regions_candidate.begin_position = BOUNDED_SUBTRACTION(current_begin_position,overlap,begin_position);
        regions_candidate.end_position = BOUNDED_ADDITION(current_begin_position,overlap,end_position);
        vector_insert(candidate_filtering_candidates->filtering_regions,regions_candidate,filtering_region_t);
        ++num_chunks_added;
      }
      // Add new chunk
      regions_candidate.begin_position = current_begin_position;
      regions_candidate.end_position = current_end_position;
      vector_insert(candidate_filtering_candidates->filtering_regions,regions_candidate,filtering_region_t);
      ++num_chunks_added;
      // Next
      current_begin_position += window_length;
      current_end_position += window_length;
    }
    // Last chunk
    if (num_chunks_added==0 || (current_begin_position < end_position && end_position-current_begin_position > overlap)) {
      regions_candidate.begin_position = current_begin_position;
      regions_candidate.end_position = end_position;
      vector_insert(candidate_filtering_candidates->filtering_regions,regions_candidate,filtering_region_t);
    } else if (current_begin_position == begin_position) {
      filtering_region_t* const last_regions_candidate =
          vector_get_last_elm(candidate_filtering_candidates->filtering_regions,filtering_region_t);
      last_regions_candidate->end_position = end_position;
    }
  }
}
#else
GEM_INLINE void filtering_candidates_process_extension_candidates(
    filtering_candidates_t* const extended_filtering_candidates,
    filtering_candidates_t* const candidate_filtering_candidates,
    archive_t* const archive,text_collection_t* const text_collection,const pattern_t* const extended_pattern,
    const pattern_t* const candidate_pattern,const search_parameters_t* const search_parameters,
    paired_matches_t* const paired_matches,mm_stack_t* const mm_stack) {
  // Parameters
  const uint64_t min_template_size = mapper_stats_template_length_get_expected_min(mapper_stats);
  const uint64_t max_template_size = mapper_stats_template_length_get_expected_max(mapper_stats);
  // Allocate candidate filtering regions
  const uint64_t num_filtering_regions = vector_get_used(extended_filtering_candidates->filtering_regions);
  filtering_region_t* regions_extended = vector_get_mem(extended_filtering_candidates->filtering_regions,filtering_region_t);
  vector_reserve_additional(candidate_filtering_candidates->filtering_regions,num_filtering_regions);
  filtering_region_t* regions_candidate = vector_get_mem(candidate_filtering_candidates->filtering_regions,filtering_region_t);
  // Traverse extended-regions and generate candidate-regions
  uint64_t n;
  for (n=0;n<num_filtering_regions;++n,++regions_extended,++regions_candidate) {
    // Compute candidate region boundaries
    const uint64_t extended_eff_begin_position = regions_extended->begin_position;
    uint64_t candidate_begin_position, candidate_end_position;
    filtering_candidates_compute_extension_region(
        archive->locator,true,extended_eff_begin_position,
        extended_pattern->key_length+extended_pattern->max_effective_filtering_error,
        candidate_pattern->key_length,candidate_pattern->max_effective_filtering_error,
        min_template_size,max_template_size,&candidate_begin_position,&candidate_end_position);
    const uint64_t candidate_length = candidate_end_position-candidate_begin_position;
    // Add candidate region
    regions_candidate->status = filtering_region_unverified; // Newly created region (unverified)
    match_scaffold_init(&regions_candidate->match_scaffold);
    regions_candidate->begin_position =
        archive_text_get_projection(archive->text,candidate_begin_position,candidate_length);
    regions_candidate->end_position = regions_candidate->begin_position+candidate_length;
    regions_candidate->num_regions_matching = 0;
    regions_candidate->regions_matching = NULL;
    regions_candidate->coverage = 0;
  }
  vector_add_used(candidate_filtering_candidates->filtering_regions,num_filtering_regions);
}
#endif
/*
 * Display
 */
GEM_INLINE void filtering_candidates_print_matching_regions(
    FILE* const stream,filtering_candidates_t* const filtering_candidates) {
  int64_t i;
  fprintf(stream,"[GEM]>Matching.Regions\n");
  fprintf(stream,"  => Initial.Regions\n");
  const uint64_t num_regions = vector_get_used(filtering_candidates->regions_buffer);
  region_search_t* const regions = vector_get_mem(filtering_candidates->regions_buffer,region_search_t);
  for (i=num_regions-1;i>=0;--i) {
    fprintf(stream,"    #%lu -> [%lu,%lu) \n",num_regions-i-1,regions[i].begin,regions[i].end);
  }
  fprintf(stream,"  => Matching.Regions\n");
  const uint64_t num_candidate_regions = vector_get_used(filtering_candidates->filtering_regions);
  filtering_region_t* const filtering_region = vector_get_mem(filtering_candidates->filtering_regions,filtering_region_t);
  for (i=0;i<num_candidate_regions;++i) {
    match_scaffold_print(stream,NULL,&filtering_region->match_scaffold);
  }
}


