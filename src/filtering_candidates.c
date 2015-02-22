/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering_candidates.h"
#include "matches_align.h"
#include "archive_text_retrieve.h"
#include "filtering_region.h"

/*
 * Debug
 */
#define DEBUG_FILTERING_REGIONS            false
#define DEBUG_VERIFY_REGIONS               false
#define DEBUG_ACCEPTED_REGIONS             false
#define DEBUG_REGIONS_MATCHING             false

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
  uint64_t begin_position;
  uint64_t effective_begin_position;
  uint64_t effective_end_position;
} filtering_position_t;
/*
 * Setup
 */
GEM_INLINE void filtering_candidates_init(filtering_candidates_t* const filtering_candidates) {
  // Filtering Status
  filtering_candidates->total_candidates_accepted = 0;
  // Region Buffer
  filtering_candidates->regions_buffer = vector_new(REGIONS_BUFFER_INIT,region_search_t);
  // Candidates
  filtering_candidates->filtering_positions = vector_new(CANDIDATE_POSITIONS_INIT,filtering_position_t);
  filtering_candidates->verified_positions = ihash_new();
  filtering_candidates->filtering_regions = vector_new(CANDIDATE_POSITIONS_INIT,filtering_region_t);
  filtering_candidates->verified_regions = vector_new(CANDIDATE_POSITIONS_INIT,verified_region_t);
}
GEM_INLINE void filtering_candidates_clear(filtering_candidates_t* const filtering_candidates) {
  // Filtering Status
  filtering_candidates->total_candidates_accepted = 0;
  // Region Buffer
  vector_clear(filtering_candidates->regions_buffer);
  // Candidates
  vector_clear(filtering_candidates->filtering_positions);
  ihash_clear(filtering_candidates->verified_positions);
  vector_clear(filtering_candidates->filtering_regions);
  vector_clear(filtering_candidates->verified_regions);
}
GEM_INLINE void filtering_candidates_destroy(filtering_candidates_t* const filtering_candidates) {
  vector_delete(filtering_candidates->regions_buffer);
  vector_delete(filtering_candidates->filtering_positions);
  ihash_delete(filtering_candidates->verified_positions);
  vector_delete(filtering_candidates->filtering_regions);
  vector_delete(filtering_candidates->verified_regions);
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
/*
 * Adding candidate positions
 */
GEM_INLINE uint64_t filtering_candidates_add_region(
    filtering_candidates_t* const filtering_candidates,const uint64_t region_start_pos,
    const uint64_t region_end_pos,const uint64_t region_errors) {
  // Store region
  region_search_t* region;
  vector_alloc_new(filtering_candidates->regions_buffer,region_search_t,region);
  region->start = region_start_pos;
  region->end = region_end_pos;
  region->degree = region_errors;
  return vector_get_used(filtering_candidates->regions_buffer)-1;
}
GEM_INLINE void filtering_candidates_index_candidate_position(
    filtering_candidates_t* const filtering_candidates,
    const uint64_t index_position,mm_stack_t* const mm_stack) {
  ihash_insert(filtering_candidates->verified_positions,index_position,1);
}
GEM_INLINE bool filtering_candidates_lookup_candidate_position(
    filtering_candidates_t* const filtering_candidates,const uint64_t index_position) {
  uint64_t* const verified_index_position =
      ihash_get(filtering_candidates->verified_positions,index_position,uint64_t);
  return (verified_index_position!=NULL);
}
GEM_INLINE void filtering_candidates_add_interval(
    filtering_candidates_t* const filtering_candidates,
    const uint64_t interval_lo,const uint64_t interval_hi,const uint64_t region_start_pos,
    const uint64_t region_end_pos,const uint64_t region_errors,mm_stack_t* const mm_stack) {
  const uint64_t num_candidates = interval_hi-interval_lo;
  if (gem_expect_false(num_candidates==0)) return;
  // Store region
  const uint64_t region_offset = filtering_candidates_add_region(
      filtering_candidates,region_start_pos,region_end_pos,region_errors);
  // Store candidate positions (index-space)
  vector_reserve_additional(filtering_candidates->filtering_positions,num_candidates);
  filtering_position_t* position_index = vector_get_free_elm(filtering_candidates->filtering_positions,filtering_position_t);
  uint64_t index_position;
  for (index_position=interval_lo;index_position<interval_hi;++index_position) {
    if (filtering_candidates_lookup_candidate_position(filtering_candidates,index_position)) {
      PROF_INC_COUNTER(GP_CANDIDATE_POSITIONS_DUPLICATED);
      continue;
    }
    filtering_candidates_index_candidate_position(filtering_candidates,index_position,mm_stack);
    position_index->source_region_offset = region_offset;
    position_index->region_index_position = index_position;
    ++position_index;
  }
  vector_update_used(filtering_candidates->filtering_positions,position_index);
}
GEM_INLINE void filtering_candidates_add_interval_set(
    filtering_candidates_t* const filtering_candidates,interval_set_t* const interval_set,
    const uint64_t region_start_pos,const uint64_t region_end_pos,mm_stack_t* const mm_stack) {
  INTERVAL_SET_ITERATE(interval_set,interval) {
    filtering_candidates_add_interval(filtering_candidates,interval->lo,interval->hi,
        region_start_pos,region_end_pos,interval->distance,mm_stack);
  }
}
GEM_INLINE void filtering_candidates_add_interval_set_thresholded(
    filtering_candidates_t* const filtering_candidates,interval_set_t* const interval_set,
    const uint64_t region_start_pos,const uint64_t region_end_pos,const uint64_t max_error,
    mm_stack_t* const mm_stack) {
  INTERVAL_SET_ITERATE(interval_set,interval) {
    if (interval->distance <= max_error) {
      filtering_candidates_add_interval(filtering_candidates,interval->lo,interval->hi,
          region_start_pos,region_end_pos,interval->distance,mm_stack);
    }
  }
}
/*
 * Sorting
 */
int filtering_position_cmp_position(const filtering_position_t* const a,const filtering_position_t* const b) {
  return a->begin_position - b->begin_position;
}
int filtering_region_cmp_align_distance(const filtering_region_t* const a,const filtering_region_t* const b) {
  return a->align_distance - b->align_distance;
}
int filtering_region_cmp_position(const filtering_region_t* const a,const filtering_region_t* const b) {
  return a->begin_position - b->begin_position;
}
int verified_region_cmp_position(const verified_region_t* const a,const verified_region_t* const b) {
  return a->effective_begin_position - b->effective_begin_position;
}
GEM_INLINE void filtering_positions_sort_positions(vector_t* const filtering_positions) {
  void* array = vector_get_mem(filtering_positions,filtering_position_t);
  const size_t count = vector_get_used(filtering_positions);
  qsort(array,count,sizeof(filtering_position_t),(int (*)(const void *,const void *))filtering_position_cmp_position);
}
GEM_INLINE void filtering_regions_align_distance(vector_t* const filtering_regions) {
  void* array = vector_get_mem(filtering_regions,filtering_region_t);
  const size_t count = vector_get_used(filtering_regions);
  qsort(array,count,sizeof(filtering_region_t),(int (*)(const void *,const void *))filtering_region_cmp_align_distance);
}
GEM_INLINE void filtering_regions_sort_position(vector_t* const filtering_regions) {
  void* array = vector_get_mem(filtering_regions,filtering_region_t);
  const size_t count = vector_get_used(filtering_regions);
  qsort(array,count,sizeof(filtering_region_t),(int (*)(const void *,const void *))filtering_region_cmp_position);
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
    text_collection_t* const candidates_collection,
    const pattern_t* const pattern,const strand_t search_strand,
    const search_actual_parameters_t* const search_actual_parameters,
    const bool extended_matches,const uint64_t extension_candidate_length,
    const bool approximated_distance,matches_t* const matches,mm_stack_t* const mm_stack) {
  // Parameters
  search_parameters_t* const search_parameters = search_actual_parameters->search_parameters;
  const uint64_t max_delta = search_actual_parameters->max_filtering_strata_after_best_nominal;
  // Sort filtering regions (wrt align_distance)
  filtering_regions_align_distance(filtering_candidates->filtering_regions);
  // Traverse all accepted candidates (text-space)
  const uint64_t num_filtering_regions = vector_get_used(filtering_candidates->filtering_regions);
  vector_reserve_additional(filtering_candidates->verified_regions,num_filtering_regions);
  filtering_region_t* regions_in = vector_get_mem(filtering_candidates->filtering_regions,filtering_region_t);
  filtering_region_t* regions_out = regions_in;
  verified_region_t* verified_regions = vector_get_free_elm(filtering_candidates->verified_regions,verified_region_t);
  uint64_t n, num_accepted_regions = 0;
  gem_cond_debug_block(DEBUG_ACCEPTED_REGIONS) { gem_slog("[GEM]>Accepted.Regions\n"); }
  for (n=0;n<num_filtering_regions;++n,++regions_in) {
    if (regions_in->status != filtering_region_accepted) { // Skip other regions
      *regions_out = *regions_in;
      ++regions_out;
    } else {
      // Check filtering strata-after-best
      const uint64_t min_counters_value = matches_counters_get_min(matches);
      if (approximated_distance && regions_in->align_distance_min_bound > min_counters_value+max_delta) {
        verified_regions->status = filtering_region_accepted_subdominant;
      } else if (!approximated_distance && regions_in->align_distance > min_counters_value+max_delta) {
        verified_regions->status = filtering_region_accepted_subdominant;
      } else {
        // DEBUG
        gem_cond_debug_block(DEBUG_ACCEPTED_REGIONS) {
          gem_slog("  => Region [%lu,%lu) (dist=%lu,regMatch=%lu)\n",
              regions_in->effective_begin_position,regions_in->effective_end_position,
              regions_in->align_distance,regions_in->num_regions_matching);
        }
        // Align the region
        match_trace_t match_trace;
        if (filtering_region_align(regions_in,candidates_collection,
              search_parameters,search_strand,pattern,matches,&match_trace,mm_stack)) {
          ++num_accepted_regions;
          // Correct match position
          if (extended_matches && search_strand==Reverse) {
            match_trace.index_position = 2*regions_in->begin_position -
                match_trace.index_position + extension_candidate_length  - match_trace.effective_length;
          }
          #ifdef GEM_DEBUG
          match_trace.regions_matching = regions_in->regions_matching;
          match_trace.num_regions_matching = regions_in->num_regions_matching;
          #endif
          matches_add_match_trace_t(matches,&match_trace,true,mm_stack);
          verified_regions->status = filtering_region_aligned;
        } else {
          verified_regions->status = filtering_region_aligned_subdominant;
        }
      }
      // Store as verified
      verified_regions->effective_begin_position = regions_in->effective_begin_position;
      verified_regions->effective_end_position = regions_in->effective_end_position;
      ++verified_regions;
    }
  }
  // Update vectors
  vector_update_used(filtering_candidates->verified_regions,verified_regions);
  vector_update_used(filtering_candidates->filtering_regions,regions_out);
  // Return total accepted regions
  filtering_candidates->total_candidates_accepted += num_accepted_regions;
  return num_accepted_regions;
}
/*
 * Candidate Verification
 */
GEM_INLINE uint64_t filtering_candidates_verify_filtering_regions(
    filtering_candidates_t* const filtering_candidates,text_collection_t* const candidates_collection,
    const pattern_t* const pattern,const strand_t search_strand,
    const search_actual_parameters_t* const search_actual_parameters,matches_t* const matches) {
  // Parameters
  search_parameters_t* const search_parameters = search_actual_parameters->search_parameters;
  const uint64_t max_search_matches = search_actual_parameters->search_parameters->max_search_matches;
  // Traverse all regions (text-space)
  uint64_t total_accepted_regions = filtering_candidates->total_candidates_accepted;
  const uint64_t num_filtering_regions = vector_get_used(filtering_candidates->filtering_regions);
  vector_reserve_additional(filtering_candidates->verified_regions,num_filtering_regions);
  filtering_region_t* regions_in = vector_get_mem(filtering_candidates->filtering_regions,filtering_region_t);
  filtering_region_t* regions_out = regions_in;
  verified_region_t* regions_discarded = vector_get_free_elm(filtering_candidates->verified_regions,verified_region_t);
  uint64_t n, num_regions_accepted = 0;
  gem_cond_debug_block(DEBUG_VERIFY_REGIONS) { gem_slog("[GEM]>Verify.Regions\n"); }
  for (n=0;n<num_filtering_regions;++n,++regions_in) {
    // Verify region
    if (regions_in->status==filtering_region_pending) { // Skip pending regions
      *regions_out = *regions_in;
      ++regions_out;
    } else if (filtering_region_verify(regions_in,candidates_collection,search_parameters,pattern)) {
      *regions_out = *regions_in;
      ++regions_out; ++num_regions_accepted;
      if (gem_expect_false(++total_accepted_regions >= max_search_matches)) break; // Threshold matches accepted
    } else {
      regions_discarded->effective_begin_position = regions_in->effective_begin_position;
      regions_discarded->effective_end_position = regions_in->effective_end_position;
      ++regions_discarded;
    }
    // DEBUG
    gem_cond_debug_block(DEBUG_VERIFY_REGIONS) {
      gem_slog("  => Region [%lu,%lu) (dist=%lu,regMatch=%lu) Status=%s\n",
          regions_in->effective_begin_position,regions_in->effective_end_position,
          regions_in->align_distance,regions_in->num_regions_matching,
          regions_in->status==filtering_region_discarded ? "Discarded" :
          (regions_in->status==filtering_region_accepted ? "Accepted" : "?") );
    }
  }
  // Return
  vector_update_used(filtering_candidates->filtering_regions,regions_out);
  vector_update_used(filtering_candidates->verified_regions,regions_discarded);
  return num_regions_accepted;
}
/*
 * Retrieve all candidates(text) from the index
 */
GEM_INLINE void filtering_candidates_retrieve_filtering_regions(
    filtering_candidates_t* const filtering_candidates,
    text_collection_t* const candidates_collection,
    const archive_t* const archive,mm_stack_t* const mm_stack) {
  // Traverse all candidates (text-space)
  VECTOR_ITERATE(filtering_candidates->filtering_regions,filtering_region,candidate_pos,filtering_region_t) {
    // Retrieve text(s)
    const uint64_t text_length = filtering_region->effective_end_position - filtering_region->effective_begin_position;
    archive_text_retrieve(archive,candidates_collection,filtering_region->effective_begin_position,
        text_length,false,&filtering_region->text_trace_offset,mm_stack);
  }
}
/*
 * Compose filtering regions
 */
GEM_INLINE void filtering_candidates_compose_matching_regions(
    filtering_candidates_t* const filtering_candidates,
    const uint64_t first_candidate_idx,const uint64_t last_candidate_idx,mm_stack_t* const mm_stack) {
  // Fetch candidate-positions
  filtering_position_t* const candidate_positions = vector_get_mem(filtering_candidates->filtering_positions,filtering_position_t);
  const uint64_t num_regions_matching = last_candidate_idx-first_candidate_idx+1;
  // Allow new matching candidate-region
  filtering_region_t* filtering_region;
  vector_alloc_new(filtering_candidates->filtering_regions,filtering_region_t,filtering_region);
  filtering_region->status = filtering_region_unverified; // Newly created region (unverified)
  filtering_region->begin_position = candidate_positions[first_candidate_idx].begin_position;
  filtering_region->effective_begin_position = candidate_positions[first_candidate_idx].effective_begin_position;
  filtering_region->effective_end_position = candidate_positions[last_candidate_idx].effective_end_position;
  filtering_region->regions_matching = mm_stack_calloc(mm_stack,num_regions_matching,region_matching_t,false);
  filtering_region->num_regions_matching = num_regions_matching;
  PROF_ADD_COUNTER(GP_CANDIDATE_REGION_LENGTH,filtering_region->effective_end_position-filtering_region->effective_begin_position);
  uint64_t i;
  for (i=0;i<num_regions_matching;++i) {
    region_matching_t* const region_matching = filtering_region->regions_matching + i;
    filtering_position_t* const candidate_position = candidate_positions + first_candidate_idx + i;
    region_search_t* const source_region = vector_get_elm(
        filtering_candidates->regions_buffer,candidate_position->source_region_offset,region_search_t);
    // Region error
    region_matching->error = source_region->degree;
    // Read coordinates
    region_matching->key_begin = source_region->end;
    region_matching->key_end = source_region->start;
    // Text coordinates (relative to the effective begin position)
    const uint64_t region_length = region_matching->key_end - region_matching->key_begin;
    region_matching->text_begin = candidate_position->region_text_position - filtering_region->effective_begin_position;
    region_matching->text_end = region_matching->text_begin + region_length;
  }
}
GEM_INLINE uint64_t filtering_candidates_compose_filtering_regions(
    filtering_candidates_t* const filtering_candidates,const dna_text_t* const enc_text,
    const uint64_t key_length,const uint64_t max_delta_difference,mm_stack_t* const mm_stack) {
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
      // Next
      last_position = position;
      ++group_idx;
    }
    // Check region against verified regions
    bool is_already_verified = false;
    while (verified_region_idx < num_verified_regions &&
        verified_region[verified_region_idx].effective_end_position <= region_begin_position) {
      ++verified_region_idx;
    }
    if (verified_region_idx < num_verified_regions) {
      const uint64_t verified_begin_position = verified_region[verified_region_idx].effective_begin_position;
      const uint64_t verified_end_position = verified_region[verified_region_idx].effective_end_position;
      const uint64_t region_end_position = candidate_positions[group_idx-1].effective_end_position;
      is_already_verified = (verified_begin_position <= region_begin_position && region_end_position <= verified_end_position);
      gem_debug_block() {
        if (is_already_verified) {
          PROF_INC_COUNTER(GP_CANDIDATE_REGIONS_DUPLICATED);
        }
      }
    }
    if (!is_already_verified) {
      // Create a region candidate with the positions from [candidate_idx] to [group_idx-1]
      filtering_candidates_compose_matching_regions(filtering_candidates,candidate_idx,group_idx-1,mm_stack);
      gem_cond_debug_block(DEBUG_FILTERING_REGIONS) {
        filtering_region_t* filtering_region = vector_get_last_elm(filtering_candidates->filtering_regions,filtering_region_t);
        gem_slog("  Region-found => [%lu,%lu)=%lu bases (regMatch=%lu)\n",
            filtering_region->effective_begin_position,filtering_region->effective_end_position,
            filtering_region->effective_end_position-filtering_region->effective_begin_position,
            filtering_region->num_regions_matching);
      }
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

//  // Add to verified positions // TODO
//  uint64_t candidate_idx = 0, verified_idx = 0;
//  const uint64_t num_verified_positions = vector_get_used(filtering_candidates->verified_positions);
//  const uint64_t* const verified_positions = vector_get_mem(filtering_candidates->verified_positions,uint64_t);
//  // Merge verified candidate positions with accepted positions
//  filtering_candidates_add_verified_positions(filtering_candidates,num_accepted_positions);

  // Return number of filtering regions generated
  return vector_get_used(filtering_candidates->filtering_regions);
}
/*
 * Filtering adjustment of the position wrt region/seed on which the candidate is based
 */
GEM_INLINE void filtering_candidates_adjust_filtering_position(
    filtering_position_t* const filtering_position,const uint64_t begin_offset,
    const uint64_t end_offset,const uint64_t boundary_error) {
  // Adjust Position
  const locator_interval_t* const locator_interval = filtering_position->locator_interval;
  uint64_t begin_position = (filtering_position->region_text_position > begin_offset) ?
      (filtering_position->region_text_position - begin_offset) : 0;
  uint64_t effective_begin_position;
  if (begin_position < locator_interval->begin_position) { // Adjust by locator-interval
    begin_position = locator_interval->begin_position; // Possible trim at the beginning
    effective_begin_position = locator_interval->begin_position;
  } else {
    effective_begin_position = (begin_position > boundary_error) ? begin_position-boundary_error : 0;
    if (effective_begin_position < locator_interval->begin_position) { // Adjust by locator-interval
      effective_begin_position = locator_interval->begin_position;
    }
  }
  uint64_t effective_end_position = filtering_position->region_text_position + end_offset + boundary_error;
  if (effective_end_position >= locator_interval->end_position) { // Adjust by locator-interval
    effective_end_position = locator_interval->end_position; // Possible trim at the end
  }
  filtering_position->begin_position = begin_position;
  filtering_position->effective_begin_position = effective_begin_position;
  filtering_position->effective_end_position = effective_end_position;
}
/*
 * Decode of all candidate positions (index-space -> text-space)
 */
GEM_INLINE void filtering_candidates_decode_filtering_positions(
    filtering_candidates_t* const filtering_candidates,const locator_t* const locator,
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
    filtering_candidates_adjust_filtering_position(filtering_position,
        source_region->end,key_length-source_region->end,boundary_error);
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
    filtering_candidates_t* const filtering_candidates,
    const locator_t* const locator,const fm_index_t* const fm_index,
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
    filtering_candidates_adjust_filtering_position(filtering_position+current_position,
        source_region->end,key_length-source_region->end,boundary_error);
  }
}
/*
 * Processing & Verification
 */
GEM_INLINE uint64_t filtering_candidates_process_candidates(
    filtering_candidates_t* const filtering_candidates,
    archive_t* const archive,const pattern_t* const pattern,
    mm_stack_t* const mm_stack) {
  // Check non-empty pending candidates set
  uint64_t pending_candidates = vector_get_used(filtering_candidates->filtering_positions);
  PROF_ADD_COUNTER(GP_CANDIDATE_POSITIONS,pending_candidates);
  if (pending_candidates==0) return 0; // Nothing to do

  // Batch decode+adjust of all positions of the candidates (cip_begin_position = decoded(candidate_region_index_position))
  PROF_START(GP_FC_PROCESS_CANDIDATES);
  PROF_START(GP_FC_DECODE_POSITIONS);
  const uint64_t key_length = pattern->key_length;
  if (pending_candidates < DECODE_NUM_POSITIONS_PREFETCHED) {
    filtering_candidates_decode_filtering_positions(filtering_candidates,
        archive->locator,archive->fm_index,key_length,pattern->max_effective_bandwidth);
  } else {
    filtering_candidates_decode_filtering_position_batch_prefetched(filtering_candidates,
        archive->locator,archive->fm_index,key_length,pattern->max_effective_bandwidth);
  }
  PROF_STOP(GP_FC_DECODE_POSITIONS);

  // Compose matching regions into candidate regions (also filter out duplicated positions or already checked)
  PROF_START(GP_FC_COMPOSE_REGIONS);
  pending_candidates = filtering_candidates_compose_filtering_regions(
      filtering_candidates,archive->enc_text,key_length,pattern->max_effective_bandwidth,mm_stack);
  PROF_STOP(GP_FC_COMPOSE_REGIONS);
  PROF_ADD_COUNTER(GP_CANDIDATE_REGIONS,pending_candidates);

  // Return total candidate regions
  PROF_STOP(GP_FC_PROCESS_CANDIDATES);
  return pending_candidates;
}
GEM_INLINE uint64_t filtering_candidates_verify_candidates(
    filtering_candidates_t* const filtering_candidates,
    archive_t* const archive,text_collection_t* const text_collection,
    const pattern_t* const pattern,const strand_t search_strand,
    const search_actual_parameters_t* const search_actual_parameters,
    matches_t* const matches,mm_stack_t* const mm_stack) {
  // Check number of filtering regions
  uint64_t pending_candidates = vector_get_used(filtering_candidates->filtering_regions);
  if (pending_candidates==0) return 0;

  // Retrieve text-candidates
  PROF_START(GP_FC_VERIFICATION);
  PROF_START(GP_FC_RETRIEVE_CANDIDATE_REGIONS);
  filtering_candidates_retrieve_filtering_regions(filtering_candidates,text_collection,archive,mm_stack);
  PROF_STOP(GP_FC_RETRIEVE_CANDIDATE_REGIONS);

  // Verify candidates
  PROF_START(GP_FC_VERIFY_CANDIDATE_REGIONS);
  pending_candidates = filtering_candidates_verify_filtering_regions(
      filtering_candidates,text_collection,pattern,search_strand,search_actual_parameters,matches);
  PROF_STOP(GP_FC_VERIFY_CANDIDATE_REGIONS);
  if (pending_candidates==0) { PROF_STOP(GP_FC_VERIFICATION); return 0; }

  // Align accepted candidates
  PROF_START(GP_FC_REALIGN_CANDIDATE_REGIONS);
  matches_hint_allocate_match_trace(matches,pending_candidates); // Hint to matches
  const uint64_t accepted_regions = filtering_candidates_align_accepted_regions(filtering_candidates,
      text_collection,pattern,search_strand,search_actual_parameters,false,0,false,matches,mm_stack);
  PROF_STOP(GP_FC_REALIGN_CANDIDATE_REGIONS);

  PROF_STOP(GP_FC_VERIFICATION);
  return accepted_regions;
}
/*
 * BPM-Buffer API (Verification)
 */
GEM_INLINE uint64_t filtering_candidates_bpm_buffer_add(
    filtering_candidates_t* const filtering_candidates,
    archive_t* const archive,pattern_t* const pattern,const strand_t search_strand,
    const search_actual_parameters_t* const search_actual_parameters,
    bpm_gpu_buffer_t* const bpm_gpu_buffer,mm_stack_t* const mm_stack) {
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
    const uint64_t begin_position = candidate_region->effective_begin_position;
    const uint64_t candidate_length = candidate_region->effective_end_position - begin_position;
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
GEM_INLINE void filtering_candidates_bpm_buffer_align(
    filtering_candidates_t* const filtering_candidates,
    archive_t* const archive,text_collection_t* const text_collection,
    pattern_t* const pattern,const strand_t search_strand,const search_actual_parameters_t* const search_actual_parameters,
    bpm_gpu_buffer_t* const bpm_gpu_buffer,const uint64_t candidate_offset_begin,const uint64_t candidate_offset_end,
    matches_t* const matches,mm_stack_t* const mm_stack) {
  // Count total candidates
  const uint64_t pending_candidates = filtering_candidates_get_num_candidate_regions(filtering_candidates);
  if (gem_expect_false(pending_candidates==0)) return;
  /*
   * Retrieve filtering-regions from BPM-Buffer
   */
  PROF_START(GP_FC_RETRIEVE_BPM_BUFFER_CANDIDATE_REGIONS);
  // Hint to matches
  matches_hint_allocate_match_trace(matches,pending_candidates);
  // Fetch Parameters
  const uint64_t key_length = pattern->key_length;
  const uint64_t max_error = pattern->max_effective_filtering_error;
  // Fetch tile dimensions
  const uint64_t num_chunks = pattern->bpm_pattern.gpu_num_chunks;
  const uint64_t pattern_tile_tall = pattern->bpm_pattern.gpu_entries_per_chunk*BPM_GPU_PATTERN_ENTRY_LENGTH;
  // Prepare filtering-regions vectors
  vector_reserve_additional(filtering_candidates->verified_regions,pending_candidates);
  filtering_region_t* candidate_regions = vector_get_mem(filtering_candidates->filtering_regions,filtering_region_t);
  filtering_region_t* accepted_regions = candidate_regions;
  verified_region_t* verified_regions = vector_get_mem(filtering_candidates->verified_regions,verified_region_t);
  // Traverse all candidates (text-space) & sum-up their alignment distance
  uint64_t candidate_idx=candidate_offset_begin, pattern_chunk, candidate_pos;
  gem_cond_debug_block(DEBUG_VERIFY_REGIONS) { gem_slog("[GEM]>Verify.Regions\n"); }
  for (candidate_pos=0;candidate_pos<pending_candidates;++candidate_pos,++candidate_regions) {
    // Get the accepted candidate
    const uint64_t begin_position = candidate_regions->effective_begin_position;
    const uint64_t candidate_length = candidate_regions->effective_end_position - begin_position;
    // Calculate tile dimensions
    pattern_tiled_t pattern_tiled;
    const bool pattern_can_align = pattern_tiled_init(&pattern_tiled,key_length,pattern_tile_tall,candidate_length,max_error);
    if (!pattern_can_align) continue;
    // Sum up the alignment distance of all the tiles
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
    // Check total distance
    if (!unaligned_tiled && global_distance <= max_error) {
      // Allocate text-trace
      const uint64_t text_trace_offset = text_collection_new_trace(text_collection);
      text_trace_t* const text_trace = text_collection_get_trace(text_collection,text_trace_offset);
      text_trace->text = dna_text_retrieve_sequence(archive->enc_text,begin_position,candidate_length,mm_stack);
      text_trace->length = candidate_length;
      // Configure accepted candidate
      accepted_regions->status = filtering_region_accepted;
      accepted_regions->text_trace_offset = text_trace_offset;
      accepted_regions->begin_position = begin_position;
      accepted_regions->effective_begin_position = candidate_regions->effective_begin_position;
      accepted_regions->effective_end_position = candidate_regions->effective_end_position;
      // Configure regions matching
      accepted_regions->num_regions_matching = candidate_regions->num_regions_matching;
      accepted_regions->regions_matching = candidate_regions->regions_matching;
      accepted_regions->coverage = candidate_regions->coverage;
      // Distance Bound estimation
      accepted_regions->align_distance_min_bound = (global_distance > max_error) ? max_error : global_distance;
      global_distance += distance_link_tiles;
      accepted_regions->align_distance = (global_distance > max_error) ? max_error : global_distance;
      accepted_regions->align_match_begin_column = 0;
      accepted_regions->align_match_end_column =
          BOUNDED_ADDITION(pattern_tiled.prev_tile_match_position,accepted_regions->align_distance,candidate_length-1);
      ++accepted_regions;
      gem_cond_debug_block(DEBUG_VERIFY_REGIONS) {
        gem_slog("  => Region [%lu,%lu) (dist=%lu,regMatch=%lu) Status=OK\n",
            accepted_regions->effective_begin_position,accepted_regions->effective_end_position,
            accepted_regions->align_distance,accepted_regions->num_regions_matching);
      }
    } else {
      verified_regions->status = filtering_region_discarded;
      verified_regions->effective_begin_position = candidate_regions->effective_begin_position;
      verified_regions->effective_end_position = candidate_regions->effective_end_position;
      ++verified_regions;
      gem_cond_debug_block(DEBUG_VERIFY_REGIONS) {
        gem_slog("  => Region [%lu,%lu) (dist=inf,regMatch=%lu) Status=DISCARDED\n",
            candidate_regions->effective_begin_position,candidate_regions->effective_end_position,
            accepted_regions->num_regions_matching);
      }
    }
  }
  // Update Accepted/Discarded
  vector_update_used(filtering_candidates->filtering_regions,accepted_regions);
  vector_update_used(filtering_candidates->verified_regions,verified_regions);
  PROF_STOP(GP_FC_RETRIEVE_BPM_BUFFER_CANDIDATE_REGIONS);
  /*
   * Realign
   */
  PROF_START(GP_FC_REALIGN_CANDIDATE_REGIONS);
  filtering_candidates_align_accepted_regions(
      filtering_candidates,text_collection,pattern,search_strand,
      search_actual_parameters,false,0,true,matches,mm_stack);
  PROF_STOP(GP_FC_REALIGN_CANDIDATE_REGIONS);
}
/*
 * Paired Verification
 */
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
GEM_INLINE void filtering_candidates_paired_regions_filtering(
    filtering_candidates_t* const filtering_candidates_end1,filtering_candidates_t* const filtering_candidates_end2,
    const uint64_t min_template_length,const uint64_t max_template_length,const bool absolute_distance) {
  // Locate compatible filtering regions
  const uint64_t num_candidates_end1 = vector_get_used(filtering_candidates_end1->filtering_regions);
  filtering_region_t* candidate_end1 = vector_get_mem(filtering_candidates_end1->filtering_regions,filtering_region_t);
  const uint64_t num_candidates_end2 = vector_get_used(filtering_candidates_end2->filtering_regions);
  filtering_region_t* candidate_end2 = vector_get_mem(filtering_candidates_end2->filtering_regions,filtering_region_t);
  uint64_t pos_candidates_end1 = 0, pos_candidates_end2 = 0;
  while (pos_candidates_end1 < num_candidates_end1 && pos_candidates_end2 < num_candidates_end2) {
    // Use an auxiliary sentinel to check compatible filtering regions
    filtering_region_t* candidate_aux = candidate_end2;
    uint64_t pos_candidates_aux = pos_candidates_end2;
    while (pos_candidates_aux < num_candidates_end2) {
      // Calculate template length
      const uint64_t template_length = paired_match_get_template_observed_length(
          candidate_end1->effective_begin_position,candidate_end1->effective_end_position,
          candidate_end2->effective_begin_position,candidate_end2->effective_end_position);
      // Check relative position
      if (candidate_end2->effective_begin_position < candidate_end1->effective_begin_position) {
        if (absolute_distance && min_template_length <= template_length && template_length <= max_template_length) {
          // Accepted
          candidate_end1->status = filtering_region_unverified;
          candidate_end2->status = filtering_region_unverified;
        } else {
          ++candidate_end2; ++pos_candidates_end2; // Increase base-end2 sentinel
        }
      } else {
        if (min_template_length <= template_length && template_length <= max_template_length) {
          // Accepted
          candidate_end1->status = filtering_region_unverified;
          candidate_end2->status = filtering_region_unverified;
        } else if (template_length > max_template_length) {
          break; // Passed the limit nothing is compatible
        }
      }
      // Increase auxiliary sentinel
      ++candidate_aux; ++pos_candidates_aux;
    }
    // Increase end1 sentinel
    ++candidate_end1; ++pos_candidates_end1;
  }
}
GEM_INLINE uint64_t filtering_candidates_extend_match(
    filtering_candidates_t* const filtering_candidates,
    archive_t* const archive,text_collection_t* const text_collection,
    const match_trace_t* const extended_match,const pattern_t* const candidate_pattern,
    const strand_t candidate_search_strand,const bool search_onward,
    const search_actual_parameters_t* const candidate_actual_parameters,
    paired_matches_t* const paired_matches,const sequence_end_t candidate_end,
    mm_stack_t* const mm_stack) {
  PROF_START(GP_FC_EXTEND_MATCH);
  // Parameters
  search_parameters_t* const search_parameters = candidate_actual_parameters->search_parameters;
  const bool search_reverse_strand = candidate_search_strand==Reverse;
  const uint64_t max_filtering_error = candidate_pattern->max_effective_filtering_error;
  /*
   * Retrieve text-candidate
   */
  PROF_START(GP_FC_EXTEND_RETRIEVE_CANDIDATE_REGIONS);
  const uint64_t max_expected_template_size =
      (COUNTER_GET_MEAN(&paired_matches->unique_template_size) +
      (2.0*COUNTER_GET_STDDEV(&paired_matches->unique_template_size)));
  locator_interval_t* const locator_interval = locator_lookup_interval(archive->locator,extended_match->index_position);
  uint64_t text_trace_offset;
  uint64_t index_begin_position, index_end_position;
  if (search_onward) {
    index_begin_position = extended_match->index_position;
    index_end_position = extended_match->index_position + extended_match->effective_length +
        max_expected_template_size + candidate_pattern->key_length + max_filtering_error;
    if (index_end_position > locator_interval->end_position) index_end_position = locator_interval->end_position;
  } else {
    index_end_position = extended_match->index_position + extended_match->effective_length;
    const uint64_t offset = max_expected_template_size + candidate_pattern->key_length + max_filtering_error;
    index_begin_position = (extended_match->index_position > offset) ? extended_match->index_position - offset : 0 ;
    if (index_begin_position < locator_interval->begin_position) index_begin_position = locator_interval->begin_position;
  }
  const uint64_t candidate_length = index_end_position-index_begin_position;
  archive_text_retrieve(archive,text_collection,index_begin_position,
      candidate_length,search_reverse_strand,&text_trace_offset,mm_stack);
  PROF_STOP(GP_FC_EXTEND_RETRIEVE_CANDIDATE_REGIONS);
  /*
   * Verify candidate region (may contain multiple matches)
   */
  PROF_START(GP_FC_EXTEND_VERIFY_CANDIDATE_REGIONS);
  uint64_t candidates_found = filtering_region_verify_extension(
      filtering_candidates->filtering_regions,text_collection,
      text_trace_offset,index_begin_position,search_parameters,candidate_pattern);
  PROF_STOP(GP_FC_EXTEND_VERIFY_CANDIDATE_REGIONS);
  if (candidates_found==0) { PROF_STOP(GP_FC_EXTEND_MATCH); return 0; }
  /*
   * Align accepted candidates
   */
  PROF_START(GP_FC_EXTEND_REALIGN_CANDIDATE_REGIONS);
  matches_t* matches_candidate = (candidate_end==paired_end1) ? paired_matches->matches_end1 : paired_matches->matches_end2;
  matches_hint_allocate_match_trace(matches_candidate,candidates_found); // Hint to matches
  // Align
  candidates_found = filtering_candidates_align_accepted_regions(
      filtering_candidates,text_collection,candidate_pattern,candidate_search_strand,
      candidate_actual_parameters,true,candidate_length,false,matches_candidate,mm_stack);
  PROF_STOP(GP_FC_EXTEND_REALIGN_CANDIDATE_REGIONS);
  PROF_STOP(GP_FC_EXTEND_MATCH);
  // Return number of extended-matches found
  return candidates_found;
}
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
    fprintf(stream,"    #%lu -> [%lu,%lu) \n",num_regions-i-1,regions[i].end,regions[i].start);
  }
  fprintf(stream,"  => Matching.Regions\n");
  const uint64_t num_candidate_regions = vector_get_used(filtering_candidates->filtering_regions);
  filtering_region_t* const candidate_regions = vector_get_mem(filtering_candidates->filtering_regions,filtering_region_t);
  for (i=0;i<num_candidate_regions;++i) {
    filtering_region_print_matching_regions(stream,
        candidate_regions[i].regions_matching,candidate_regions[i].num_regions_matching,
        candidate_regions[i].effective_begin_position,candidate_regions[i].effective_end_position,i);
  }
}


