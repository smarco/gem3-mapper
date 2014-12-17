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
  region_search_t* source_region;
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
 * Verified Region
 */
typedef struct {
  uint64_t position;
  uint64_t length;
} verified_region_t;
/*
 * Setup
 */
GEM_INLINE void filtering_candidates_init(filtering_candidates_t* const filtering_candidates) {
  filtering_candidates->filtering_positions = vector_new(CANDIDATE_POSITIONS_INIT,filtering_position_t);
  filtering_candidates->filtering_regions = vector_new(CANDIDATE_POSITIONS_INIT,filtering_region_t);
  filtering_candidates->verified_regions = vector_new(CANDIDATE_POSITIONS_INIT,verified_region_t);
  filtering_candidates->num_candidates_accepted = 0;
  filtering_candidates->regions_buffer = vector_new(REGIONS_BUFFER_INIT,region_search_t);
}
GEM_INLINE void filtering_candidates_clear(filtering_candidates_t* const filtering_candidates) {
  vector_clear(filtering_candidates->filtering_positions);
  vector_clear(filtering_candidates->filtering_regions);
  vector_clear(filtering_candidates->verified_regions);
  filtering_candidates->num_candidates_accepted = 0;
  vector_clear(filtering_candidates->regions_buffer);
}
GEM_INLINE void filtering_candidates_destroy(filtering_candidates_t* const filtering_candidates) {
  vector_delete(filtering_candidates->filtering_positions);
  vector_delete(filtering_candidates->filtering_regions);
  vector_delete(filtering_candidates->verified_regions);
  vector_delete(filtering_candidates->regions_buffer);
}
/*
 * Accessors
 */
GEM_INLINE uint64_t filtering_candidates_get_num_candidate_regions(const filtering_candidates_t* const filtering_candidates) {
  return vector_get_used(filtering_candidates->filtering_regions);
}
/*
 * Adding candidate positions
 */
GEM_INLINE void filtering_candidates_add_interval(
    filtering_candidates_t* const filtering_candidates,const uint64_t interval_lo,const uint64_t interval_hi,
    const uint64_t region_start_pos,const uint64_t region_end_pos,const uint64_t region_errors) {
  const uint64_t num_candidates = interval_hi-interval_lo;
  if (gem_expect_false(num_candidates==0)) return;
  // Store region
  region_search_t* region;
  vector_alloc_new(filtering_candidates->regions_buffer,region_search_t,region);
  region->start = region_start_pos;
  region->end = region_end_pos;
  region->degree = region_errors;
  // Store candidate positions (index-space)
  vector_reserve_additional(filtering_candidates->filtering_positions,num_candidates);
  filtering_position_t* position_index = vector_get_free_elm(filtering_candidates->filtering_positions,filtering_position_t);
  uint64_t index_position;
  for (index_position=interval_lo;index_position<interval_hi;++index_position) {
    position_index->source_region = region;
    position_index->region_index_position = index_position;
    ++position_index;
  }
  vector_add_used(filtering_candidates->filtering_positions,num_candidates);
}
GEM_INLINE void filtering_candidates_add_interval_set(
    filtering_candidates_t* const filtering_candidates,interval_set_t* const interval_set,
    const uint64_t region_start_pos,const uint64_t region_end_pos) {
  GEM_NOT_IMPLEMENTED(); // TODO
}
GEM_INLINE void filtering_candidates_add_interval_set_thresholded(
    filtering_candidates_t* const filtering_candidates,interval_set_t* const interval_set,
    const uint64_t region_start_pos,const uint64_t region_end_pos,const uint64_t max_error) {
  GEM_NOT_IMPLEMENTED(); // TODO
}
//  uint64_t it;
//  interval_t* result_interval = (interval_t*)vector_get_mem(result_vector) + init_int;
//  for (it=init_int; it<end_int; ++it, ++result_interval) {
//    if (result_interval->misms <= num_misms) {
//      ADD_INTERVAL_TO_FILTER_QUERIES(buffer_queries, result_interval->lo,
//        result_interval->hi, start_pos, end_pos, result_interval->misms);
//    }
//  }
// }
/*
 * Sorting
 */
int filtering_position_cmp_position(const filtering_position_t* const a,const filtering_position_t* const b) {
  return a->begin_position - b->begin_position;
}
GEM_INLINE void filtering_candidates_sort_filtering_positions(const filtering_candidates_t* const filtering_candidates) {
  // Sort filtering positions (filtering_position_t) wrt their position in the text
  void* array = vector_get_mem(filtering_candidates->filtering_positions,filtering_position_t);
  const size_t count = vector_get_used(filtering_candidates->filtering_positions);
  qsort(array,count,sizeof(filtering_position_t),(int (*)(const void *,const void *))filtering_position_cmp_position);
}
/*
 * Candidate Alignment
 */
GEM_INLINE void filtering_candidates_align_accepted_regions(
    filtering_candidates_t* const filtering_candidates,text_collection_t* const candidates_collection,
    const pattern_t* const pattern,const strand_t search_strand,
    const search_actual_parameters_t* const search_actual_parameters,
    matches_t* const matches,mm_stack_t* const mm_stack) {
  // Model
  const alignment_model_t alignment_model = search_actual_parameters->search_parameters->alignment_model;
  // Pattern
  const uint8_t* const key = pattern->key;
  const uint64_t key_length = pattern->key_length;
  search_parameters_t* const search_parameters = search_actual_parameters->search_parameters;
  const bool* const allowed_enc = search_parameters->allowed_enc;
  // Traverse all accepted candidates (text-space)
  VECTOR_ITERATE(filtering_candidates->filtering_regions,filtering_region,candidate_pos,filtering_region_t) {
    if (filtering_region->status == filtering_region_accepted) {
      // Align region
      // TODO measure bound for lev-distance estimation and then implement no-align subdomint-acepted
      filtering_region_align(filtering_region,candidates_collection,
          alignment_model,allowed_enc,&search_parameters->swg_penalties,
          search_strand,pattern,key,key_length,matches,mm_stack);
    }
  }
}
/*
 * Candidate Verification
 */
GEM_INLINE uint64_t filtering_candidates_verify_filtering_regions(
    filtering_candidates_t* const filtering_candidates,text_collection_t* const candidates_collection,
    const pattern_t* const pattern,const strand_t search_strand,
    const search_actual_parameters_t* const search_actual_parameters,matches_t* const matches) {
  // Model
  const alignment_model_t alignment_model = search_actual_parameters->search_parameters->alignment_model;
  // Pattern
  const uint8_t* const key = pattern->key;
  const uint64_t key_length = pattern->key_length;
  const bool* const allowed_enc = search_actual_parameters->search_parameters->allowed_enc;
  // Matching Constraints
  const uint64_t max_search_matches = search_actual_parameters->search_parameters->max_search_matches;
  const uint64_t max_effective_filtering_error = pattern->max_effective_filtering_error;
  uint64_t total_accepted_regions = filtering_candidates->num_candidates_accepted;
  // Traverse all regions (text-space)
  VECTOR_ITERATE(filtering_candidates->filtering_regions,filtering_region,n,filtering_region_t) {
    // Check matches accepted
    if (gem_expect_false(total_accepted_regions > max_search_matches)) break;
    // Verify the region
    const bool accepted = filtering_region_verify(
        filtering_region,candidates_collection,alignment_model,
        allowed_enc,pattern,key,key_length,max_effective_filtering_error);
    if (accepted) ++total_accepted_regions;
  }
  // Update
  const uint64_t accepted_regions = total_accepted_regions-filtering_candidates->num_candidates_accepted;
  filtering_candidates->num_candidates_accepted = total_accepted_regions;
  // Return
  return accepted_regions;
}
/*
 * Retrieve all candidates(text) from the index
 */
GEM_INLINE void filtering_candidates_retrieve_filtering_regions(
    filtering_candidates_t* const filtering_candidates,text_collection_t* const candidates_collection,
    const locator_t* const locator,const dna_text_t* const enc_text,mm_stack_t* const mm_stack) {
  // Traverse all candidates (text-space)
  VECTOR_ITERATE(filtering_candidates->filtering_regions,filtering_region,candidate_pos,filtering_region_t) {
    // Retrieve text(s)
    const uint64_t text_length = filtering_region->effective_end_position - filtering_region->effective_begin_position;
    archive_text_retrieve(locator,NULL,enc_text,candidates_collection,
        filtering_region->effective_begin_position,text_length,&filtering_region->text_trace_offset,mm_stack);
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
  filtering_region->status = filtering_region_new; // Newly created region
  filtering_region->begin_position = candidate_positions[first_candidate_idx].begin_position;
  filtering_region->effective_begin_position = candidate_positions[first_candidate_idx].effective_begin_position;
  filtering_region->effective_end_position = candidate_positions[last_candidate_idx].effective_end_position;
  filtering_region->regions_matching = mm_stack_calloc(mm_stack,num_regions_matching,region_matching_t,false);
  filtering_region->num_regions_matching = num_regions_matching;
  uint64_t i;
  for (i=0;i<num_regions_matching;++i) {
    region_matching_t* const region_matching = filtering_region->regions_matching + i;
    filtering_position_t* const candidate_position = candidate_positions + first_candidate_idx + i;
    // Region error
    region_matching->error = candidate_position->source_region->degree;
    // Read coordinates
    region_matching->key_begin = candidate_position->source_region->end;
    region_matching->key_end = candidate_position->source_region->start;
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
  filtering_candidates_sort_filtering_positions(filtering_candidates);
  // Traverse positions and eliminate duplicates
  const uint64_t num_candidate_positions = vector_get_used(filtering_candidates->filtering_positions);
  filtering_position_t* const candidate_positions = vector_get_mem(filtering_candidates->filtering_positions,filtering_position_t);
  uint64_t candidate_idx = 0;
  while (candidate_idx < num_candidate_positions) {
    // Determine the positions belonging to the same region
    uint64_t last_position = candidate_positions[candidate_idx].begin_position;
    uint64_t group_idx = candidate_idx + 1;
    while (group_idx < num_candidate_positions) {
      const uint64_t position = candidate_positions[group_idx].begin_position;
      const uint64_t delta = position - last_position;
      if (delta > max_delta_difference) break; // Doesn't belong to the group. Stop!
      // Next
      last_position = position;
      ++group_idx;
    }
    // Create a region candidate with the positions from [candidate_idx] to [group_idx-1]
    filtering_candidates_compose_matching_regions(filtering_candidates,candidate_idx,group_idx-1,mm_stack);
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
    const locator_t* const locator,const fm_index_t* const fm_index,
    vector_t* const candidate_text_positions,const uint64_t key_length,const uint64_t boundary_error) {
  // Traverse all candidate positions in index-space
  VECTOR_ITERATE(candidate_text_positions,filtering_position,n,filtering_position_t) {
    // Lookup Position
    filtering_position->region_text_position = fm_index_lookup(fm_index,filtering_position->region_index_position);
    // Locate Position
    filtering_position->locator_interval = locator_lookup_interval(locator,filtering_position->region_text_position);
    // Adjust Position
    filtering_candidates_adjust_filtering_position(filtering_position,
        filtering_position->source_region->end,key_length-filtering_position->source_region->end,boundary_error);
  }
}
/*
 * Batch decode of all candidate positions (index-space -> text-space)
 *   (All the steps (CSA-lookup, rankQueries) are performed with prefetch-loops)
 */
//typedef struct {
//  uint64_t vector_rank;
//  uint64_t index_position;
//  uint64_t distance;
//  uint64_t used_slot;
//  bwt_block_locator_t bwt_block_locator;
//} fc_batch_decode_candidate;
//GEM_INLINE void filtering_candidates_decode_filtering_position_batch_prefetched(
//    const locator_t* const locator,const fm_index_t* const fm_index,
//    vector_t* const candidate_text_positions,const uint64_t key_length,const uint64_t boundary_error) {
//  // Init
//  const uint64_t bwt_length = fm_index_get_length(fm_index);
//  const bwt_t* const bwt = fm_index->bwt;
//  const sampled_sa_t* const sampled_sa = fm_index->sampled_sa;
//  candidate_position_t* const candidates = vector_get_mem(candidate_text_positions,candidate_position_t);
//  // Batch Decode
//  fc_batch_decode_candidate batch[DECODE_NUM_POSITIONS_PREFETCHED];
//  const uint64_t num_candidate_text_positions = vector_get_used(candidate_text_positions);
//  // Initial fill batch
//  uint64_t current_position=0, i;
//  for (i=0;i<DECODE_NUM_POSITIONS_PREFETCHED && current_position<num_candidate_text_positions;++current_position) {
//    if (!sampled_sa_is_sampled(sampled_sa,candidates[current_position].candidate_region_index_position)) {
//      batch[i].index_position = candidates[current_position].candidate_region_index_position;
//      batch[i].vector_rank = current_position;
//      batch[i].distance = 0;
//      batch[i].used_slot = true;
//      ++i;
//    }
//  }
//  const bool full_filled_batch = (i==DECODE_NUM_POSITIONS_PREFETCHED);
//  for (;i<DECODE_NUM_POSITIONS_PREFETCHED;++i) {
//    batch[i].used_slot = false;
//  }
//  // Full-prefetch loop for sampled-LF
//  if (full_filled_batch==DECODE_NUM_POSITIONS_PREFETCHED) {
//    while (current_position<num_candidate_text_positions) {
//      for (i=0;i<DECODE_NUM_POSITIONS_PREFETCHED;++i) {
//        bwt_prefetch(bwt,batch[i].index_position,&(batch[i].bwt_block_locator));
//      }
//      for (i=0;i<DECODE_NUM_POSITIONS_PREFETCHED;++i) {
//        batch[i].index_position = bwt_prefetched_LF(bwt,batch[i].index_position,&(batch[i].bwt_block_locator));
//        ++(batch[i].distance);
//        if (sampled_sa_is_sampled(sampled_sa,batch[i].index_position)) {
//          candidates[batch[i].vector_rank].candidate_decode_sampled_pos = batch[i].index_position;
//          candidates[batch[i].vector_rank].candidate_decode_distance = batch[i].distance;
//          batch[i].used_slot = false;
//          // Select new candidate to decode
//          while (current_position < num_candidate_text_positions &&
//                 !sampled_sa_is_sampled(sampled_sa,candidates[current_position].candidate_region_index_position)) {
//            ++current_position;
//          }
//          if (current_position < num_candidate_text_positions) {
//            batch[i].index_position = candidates[current_position].candidate_region_index_position;
//            batch[i].vector_rank = current_position;
//            batch[i].distance = 0;
//            batch[i].used_slot = true;
//          }
//        }
//      }
//    }
//  }
//  // Solve remaining queries
//  for (i=0;i<DECODE_NUM_POSITIONS_PREFETCHED;++i) {
//    if (batch[i].used_slot) {
//      do {
//        batch[i].index_position = bwt_LF(bwt,batch[i].index_position);
//        ++(batch[i].distance);
//      } while (!sampled_sa_is_sampled(sampled_sa,batch[i].index_position));
//      candidates[batch[i].vector_rank].candidate_decode_sampled_pos = batch[i].index_position;
//      candidates[batch[i].vector_rank].candidate_decode_distance = batch[i].distance;
//    }
//  }
//  // Prefetch SA-retrieve samples
//  uint64_t num_left_positions = num_candidate_text_positions;
//  current_position = 0;
//  while (num_left_positions < num_candidate_text_positions) {
//    const uint64_t batch_size = MIN(num_left_positions,RETRIEVE_SAMPLE_NUM_POSITIONS_PREFETCHED);
//    const uint64_t batch_top = current_position+batch_size;
//    for (i=current_position;i<batch_top;++i) {
//      sampled_sa_prefetch_sample(sampled_sa,candidates[i].candidate_decode_sampled_pos);
//    }
//    for (i=current_position;i<batch_top;++i) {
//      candidates[i].candidate_region_text_position =
//          (sampled_sa_get_sample(sampled_sa,candidates[i].candidate_decode_sampled_pos) + candidates[i].candidate_decode_distance) % bwt_length;
//    }
//    current_position = batch_top;
//    num_left_positions -= batch_size;
//  }
//  // Adjust decoded position to the beginning of the read
//  candidate_position_t* candidate = vector_get_mem(candidate_text_positions,candidate_position_t);
//  for (current_position=0;current_position<num_candidate_text_positions;++current_position,++candidate) {
//    // Locate Position
//    candidates->locator_interval = locator_lookup_interval(locator,candidates->candidate_region_text_position);
//    // Adjust Position
//    filtering_candidates_adjust_filtering_position(candidates,
//        candidates->candidate_region->end,key_length-candidates->candidate_region->end,boundary_error);
//  }
//}
/*
 * Processing & Verification
 */
GEM_INLINE uint64_t filtering_candidates_process_candidates(
    filtering_candidates_t* const filtering_candidates,
    const locator_t* const locator,const fm_index_t* const fm_index,
    const dna_text_t* const enc_text,const pattern_t* const pattern,
    const search_actual_parameters_t* const search_actual_parameters,mm_stack_t* const mm_stack) {
  // Check non-empty pending candidates set
  uint64_t pending_candidates = vector_get_used(filtering_candidates->filtering_positions);
  PROF_ADD_COUNTER(GP_CANDIDATE_POSITIONS,pending_candidates);
  if (pending_candidates==0) return 0; // Nothing to do

  // Batch decode+adjust of all positions of the candidates (cip_begin_position = decoded(candidate_region_index_position))
  PROF_START(GP_FC_PROCESS_CANDIDATES);
  PROF_START(GP_FC_DECODE_POSITIONS);
  const uint64_t key_length = pattern->key_length;
  const uint64_t boundary_error = pattern->max_effective_filtering_error;
//  if (pending_candidates < DECODE_NUM_POSITIONS_PREFETCHED) {
    filtering_candidates_decode_filtering_positions(
        locator,fm_index,filtering_candidates->filtering_positions,key_length,boundary_error);
//  } else {  // TODO Enable batch decode
//    filtering_candidates_decode_filtering_position_batch_prefetched(
//        locator,fm_index,filtering_candidates->candidate_text_positions,key_length,boundary_error);
//  }
  PROF_STOP(GP_FC_DECODE_POSITIONS);

  // Compose matching regions into candidate regions (also filter out duplicated positions or already checked)
  PROF_START(GP_FC_COMPOSE_REGIONS);
  pending_candidates = filtering_candidates_compose_filtering_regions(
      filtering_candidates,enc_text,key_length,boundary_error,mm_stack);
  PROF_STOP(GP_FC_COMPOSE_REGIONS);
  PROF_ADD_COUNTER(GP_CANDIDATE_REGIONS,pending_candidates);

  // Return total candidate regions
  PROF_STOP(GP_FC_PROCESS_CANDIDATES);
  return pending_candidates;
}
GEM_INLINE uint64_t filtering_candidates_verify_candidates(
    filtering_candidates_t* const filtering_candidates,text_collection_t* const text_collection,
    const locator_t* const locator,const fm_index_t* const fm_index,
    const dna_text_t* const enc_text,const pattern_t* const pattern,const strand_t search_strand,
    const search_actual_parameters_t* const search_actual_parameters,
    matches_t* const matches,mm_stack_t* const mm_stack) {
  // Process candidates
  uint64_t pending_candidates = filtering_candidates_process_candidates(filtering_candidates,
      locator,fm_index,enc_text,pattern,search_actual_parameters,mm_stack);
  if (pending_candidates==0) return 0;

  // Retrieve text-candidates
  PROF_START(GP_FC_VERIFICATION);
  PROF_START(GP_FC_RETRIEVE_CANDIDATE_REGIONS);
  filtering_candidates_retrieve_filtering_regions(filtering_candidates,text_collection,locator,enc_text,mm_stack);
  PROF_STOP(GP_FC_RETRIEVE_CANDIDATE_REGIONS);

  // Verify candidates
  PROF_START(GP_FC_VERIFY_CANDIDATE_REGIONS);
  pending_candidates = filtering_candidates_verify_filtering_regions(
      filtering_candidates,text_collection,pattern,search_strand,search_actual_parameters,matches);
  PROF_STOP(GP_FC_VERIFY_CANDIDATE_REGIONS);
  if (pending_candidates==0) { PROF_STOP(GP_FC_VERIFICATION); return 0; }

  // Align accepted candidates
  PROF_START(GP_FC_REALIGN_CANDIDATE_REGIONS);
  matches_hint_add_match_trace(matches,pending_candidates); // Hint to matches
  filtering_candidates_align_accepted_regions(filtering_candidates,text_collection,
      pattern,search_strand,search_actual_parameters,matches,mm_stack);
  PROF_STOP(GP_FC_REALIGN_CANDIDATE_REGIONS);

  PROF_STOP(GP_FC_VERIFICATION);
  return pending_candidates;
}
/*
 * BPM-Buffer API (Verification)
 */
GEM_INLINE uint64_t filtering_candidates_bpm_buffer_add(
    filtering_candidates_t* const filtering_candidates,
    const locator_t* const locator,const fm_index_t* const fm_index,
    const dna_text_t* const enc_text,pattern_t* const pattern,const strand_t search_strand,
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
    filtering_candidates_t* const filtering_candidates,const text_collection_t* const text_collection,
    const dna_text_t* const enc_text,pattern_t* const pattern,const strand_t search_strand,
    const search_actual_parameters_t* const search_actual_parameters,
    bpm_gpu_buffer_t* const bpm_gpu_buffer,const uint64_t candidate_offset_begin,const uint64_t candidate_offset_end,
    matches_t* const matches,mm_stack_t* const mm_stack) {
  // Count total candidates
  const uint64_t pending_candidates = filtering_candidates_get_num_candidate_regions(filtering_candidates);
  if (gem_expect_false(pending_candidates==0)) return;
  PROF_START(GP_FC_REALIGN_CANDIDATE_REGIONS);
  // Hint to matches
  matches_hint_add_match_trace(matches,pending_candidates);
  // Fetch Parameters
  const uint8_t* const key = pattern->key;
  const uint64_t key_length = pattern->key_length;
  const uint64_t max_error = pattern->max_effective_filtering_error;
  search_parameters_t* const search_parameters = search_actual_parameters->search_parameters;
  // Fetch tile dimensions
  const uint64_t num_chunks = pattern->bpm_pattern.gpu_num_chunks;
  const uint64_t pattern_tile_tall = pattern->bpm_pattern.gpu_entries_per_chunk*BPM_GPU_PATTERN_ENTRY_LENGTH;
  // Traverse all candidates (text-space) & sum-up their alignment distance
  filtering_region_t* candidate_region = vector_get_mem(filtering_candidates->filtering_regions,filtering_region_t);
  uint64_t candidate_idx=candidate_offset_begin, pattern_chunk, candidate_pos;
  for (candidate_pos=0;candidate_pos<pending_candidates;++candidate_pos,++candidate_region) {
    // Get the accepted candidate
    const uint64_t begin_position = candidate_region->effective_begin_position;
    const uint64_t candidate_length = candidate_region->effective_end_position - begin_position;
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
      text_trace->text = dna_text_retrieve_sequence(enc_text,begin_position,candidate_length,mm_stack);
      text_trace->length = candidate_length;
      // Configure accepted candidate
      candidate_region->text_trace_offset = text_trace_offset;
      candidate_region->begin_position = begin_position;
      global_distance += distance_link_tiles;

      // FIXME Wrong Bound estimation
      candidate_region->align_distance = max_error;
      candidate_region->align_match_column = candidate_length-1;
//      candidate_region->align_distance = (global_distance > max_error) ? max_error : global_distance;
//      candidate_region->align_match_column =
//          (candidate_region->align_distance==0) ? pattern_tiled.prev_tile_match_position : candidate_length-1;
      filtering_region_align(candidate_region,text_collection,
          search_parameters->alignment_model,search_parameters->allowed_enc,
          &search_parameters->swg_penalties,search_strand,pattern,key,key_length,matches,mm_stack);
    }
  }
  PROF_STOP(GP_FC_REALIGN_CANDIDATE_REGIONS);
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
    filtering_region_print_matching_regions(stream,candidate_regions+i,i);
  }
}


