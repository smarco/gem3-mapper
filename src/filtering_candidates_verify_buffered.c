/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates_bpm_buffer.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering_candidates_verify_buffered.h"
#include "align_bpm_pattern.h"
#include "align_bpm_distance.h"

/*
 * Debug
 */
#define DEBUG_FILTERING_CANDIDATES  GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PMED


/*
 * BPM-Buffer API (Verification)
 */
uint64_t filtering_candidates_verify_buffered_add(
    filtering_candidates_t* const filtering_candidates,pattern_t* const pattern,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  // Check number of pending regions
  const uint64_t pending_candidates = filtering_candidates_get_num_candidate_regions(filtering_candidates);
  if (pending_candidates==0) return 0;
  // Add the pattern to the buffer (add new queries)
  gpu_buffer_align_bpm_add_pattern(gpu_buffer_align_bpm,&pattern->bpm_pattern);
  gpu_buffer_align_bpm_record_query_length(gpu_buffer_align_bpm,pattern->key_length);
  // Fetch pattern dimensions
  const uint64_t max_error = pattern->max_effective_filtering_error;
  const uint64_t gpu_pattern_length = pattern->key_length;
  const uint64_t gpu_pattern_num_tiles = pattern->bpm_pattern.gpu_num_tiles;
  const uint64_t gpu_pattern_tile_tall = pattern->bpm_pattern.gpu_entries_per_tile*gpu_bpm_pattern_get_entry_length();
  // Traverse all candidates (text-space) & add them to the buffer
  const filtering_region_t* candidate_region = vector_get_mem(filtering_candidates->filtering_regions,filtering_region_t);
  uint64_t candidate_pos, tile_offset, total_candidates_added;
  for (candidate_pos=0,total_candidates_added=0;candidate_pos<pending_candidates;++candidate_pos,++candidate_region) {
    // Locate candidate sequence
    const uint64_t begin_position = candidate_region->begin_position;
    const uint64_t candidate_length = candidate_region->end_position - begin_position;
    gem_slog("> Candidate #%lu (%lu nt)\n",candidate_pos,candidate_length);
    // Calculate tile dimensions
    pattern_tiled_t pattern_tiled;
    const bool pattern_can_align = pattern_tiled_init(&pattern_tiled,
        gpu_pattern_length,gpu_pattern_tile_tall,candidate_length,max_error);
    if (!pattern_can_align) continue;
    // Initialize current tile variables
    for (tile_offset=0;tile_offset<gpu_pattern_num_tiles;++tile_offset,++total_candidates_added) {
      // BPM-GPU put candidate
      gpu_buffer_align_bpm_add_candidate(gpu_buffer_align_bpm,tile_offset,
          begin_position+pattern_tiled.tile_offset,pattern_tiled.tile_wide);
      gem_slog("  => Tile #%lu (Offsets=%lu-%lu Tall=%lu Wide=%lu)\n",
          tile_offset,pattern_tiled.tile_offset,
          pattern_tiled.tile_offset+pattern_tiled.tile_wide,
          pattern_tiled.tile_tall,pattern_tiled.tile_wide);
      // Calculate next tile
      pattern_tiled_calculate_next(&pattern_tiled);
    }
  }
  gpu_buffer_align_bpm_record_candidates_per_query(gpu_buffer_align_bpm,total_candidates_added);
  // Return the final number of candidate-tiles added to the buffer
  PROF_ADD_COUNTER(GP_BMP_TILED_NUM_TILES,total_candidates_added);
  PROF_ADD_COUNTER(GP_BMP_TILED_NUM_TILES_VERIFIED,total_candidates_added);
  return total_candidates_added;
}
void filtering_candidates_verify_buffered_get_candidate(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,const uint64_t candidate_offset,
    const uint64_t num_tiles,uint64_t* const candidate_begin_position,
    uint64_t* const candidate_end_position,uint64_t* const candidate_length) {
  uint32_t tile_length;
  gpu_buffer_align_bpm_get_candidate(gpu_buffer_align_bpm,
      candidate_offset,candidate_begin_position,&tile_length);
  gpu_buffer_align_bpm_get_candidate(gpu_buffer_align_bpm,
      candidate_offset+(num_tiles-1),candidate_end_position,&tile_length);
  *candidate_end_position += tile_length;
  *candidate_length = *candidate_end_position - *candidate_begin_position;
}
bool filtering_candidates_verify_buffered_get_result(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    const uint64_t candidate_offset,const uint64_t num_tiles,
    bpm_pattern_t* const bpm_pattern,const uint64_t key_length,
    const uint64_t pattern_tile_tall,const uint64_t candidate_length,
    const uint64_t max_error,uint64_t* const global_distance,
    uint64_t* const distance_link_tiles,uint64_t* const last_tile_match_position,
    text_collection_t* const text_collection,mm_stack_t* const mm_stack) {
  // Calculate tile dimensions (Sum up the alignment distance of all the tiles)
  pattern_tiled_t pattern_tiled;
  const bool pattern_can_align = pattern_tiled_init(
      &pattern_tiled,key_length,pattern_tile_tall,candidate_length,max_error);
  if (!pattern_can_align) return false;
  // Put together all tiles from the candidate (join results)
  *global_distance = 0;
  *distance_link_tiles = 0;
  bool unaligned_tiled = false;
  uint64_t pattern_tile;
  for (pattern_tile=0;pattern_tile<num_tiles;++pattern_tile) {
    // Retrieve alignment distance
    uint32_t tile_distance=0, tile_match_column=0;
    gpu_buffer_align_bpm_get_result(gpu_buffer_align_bpm,
        candidate_offset+pattern_tile,&tile_distance,&tile_match_column);
    pattern_tiled.tile_distance = tile_distance;
    pattern_tiled.tile_match_column = tile_match_column;
    *global_distance += tile_distance;
    if (tile_distance > max_error) unaligned_tiled = true; // return false;
    // Bound the joint distance by estimating the cost of connecting the path through the tiles
    *distance_link_tiles += pattern_tiled_bound_matching_path(&pattern_tiled);
    // Calculate next tile
    pattern_tiled_calculate_next(&pattern_tiled);
    // DEBUG
#ifdef CUDA_CHECK_BUFFERED_VERIFY_CANDIDATES
    mm_stack_push_state(mm_stack);
    // Get Candidate & Pattern
    uint64_t candidate_text_position;
    uint32_t candidate_length;
    bpm_pattern_t bpm_pattern_tile;
    gpu_buffer_align_bpm_retrieve_pattern(gpu_buffer_align_bpm,
        candidate_offset+pattern_tile,&bpm_pattern_tile,mm_stack);
    gpu_buffer_align_bpm_get_candidate(gpu_buffer_align_bpm,
        candidate_offset+pattern_tile,&candidate_text_position,&candidate_length);
    // Get Candidate Text
    const uint64_t text_trace_offset = archive_text_retrieve(
        gpu_buffer_align_bpm->archive_text,text_collection,
        candidate_text_position,candidate_length,false,mm_stack); // Retrieve text(s)
    const text_trace_t* const text_trace = text_collection_get_trace(text_collection,text_trace_offset);
    const uint8_t* const text = text_trace->text; // Candidate
    uint64_t i, uncalled_bases_text = 0;
    for (i=0;i<candidate_length;++i) {
      if (text[i]==ENC_DNA_CHAR_N) ++uncalled_bases_text;
    }
    // Align BPM & Set result
    uint64_t check_tile_match_end_column, check_tile_distance;
    bpm_compute_edit_distance_cutoff(&bpm_pattern_tile,text,candidate_length,
        &check_tile_match_end_column,&check_tile_distance,bpm_pattern_tile.pattern_length,false);
    if (tile_distance!=check_tile_distance || tile_match_column!=check_tile_match_end_column) {
      gem_error_msg("Filtering.Candidates.Verify.Buffered. Check verify candidate "
          "(Distance:%d!=%lu) (MatchPos:%d!=%lu) (Text.Uncalled.bases=%lu)",
          tile_distance,check_tile_distance,tile_match_column,
          check_tile_match_end_column,uncalled_bases_text);
    }
    //      // Whole read // TODO + STATS
    //      uint64_t match_end_column, match_distance;
    //      bpm_compute_edit_distance_cutoff(global_pattern,text,candidate_length,
    //          &match_end_column,&match_distance,global_pattern->pattern_length,false);
    mm_stack_pop_state(mm_stack,false);
#endif
  }
  *last_tile_match_position = pattern_tiled.prev_tile_match_position;
  return !unaligned_tiled;
}
uint64_t filtering_candidates_verify_buffered_retrieve(
    filtering_candidates_t* const filtering_candidates,archive_text_t* const archive_text,
    text_collection_t* const text_collection,pattern_t* const pattern,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,const uint64_t candidate_offset_begin,
    const uint64_t candidate_offset_end,mm_stack_t* const mm_stack) {
  /*
   * Retrieve filtering-regions from BPM-Buffer
   */
  if (gem_expect_false(candidate_offset_begin==candidate_offset_end)) return 0;
  PROFILE_START(GP_FC_RETRIEVE_BPM_BUFFER_CANDIDATE_REGIONS,PROFILE_LEVEL);
  // Fetch Parameters
  const uint64_t max_error = pattern->max_effective_filtering_error;
  const uint64_t gpu_pattern_length = pattern->key_length;
  const uint64_t gpu_pattern_num_tiles = pattern->bpm_pattern.gpu_num_tiles;
  const uint64_t gpu_pattern_tile_tall = pattern->bpm_pattern.gpu_entries_per_tile*gpu_bpm_pattern_get_entry_length();
  // Prepare filtering-regions vectors
  const uint64_t pending_candidates = (candidate_offset_end-candidate_offset_begin)/gpu_pattern_num_tiles;
  vector_reserve(filtering_candidates->verified_regions,pending_candidates,false);
  vector_reserve(filtering_candidates->filtering_regions,pending_candidates,false);
  vector_reserve(filtering_candidates->discarded_regions,pending_candidates,false);
  filtering_region_t* regions_accepted = vector_get_mem(filtering_candidates->filtering_regions,filtering_region_t);
  filtering_region_t* regions_discarded = vector_get_mem(filtering_candidates->discarded_regions,filtering_region_t);
  verified_region_t* regions_verified = vector_get_mem(filtering_candidates->verified_regions,verified_region_t);
  // Traverse all candidates (text-space) & sum-up their alignment distance
  uint64_t num_accepted_regions = 0;
  uint64_t candidate_idx=candidate_offset_begin, candidate_pos;
  for (candidate_pos=0;candidate_pos<pending_candidates;++candidate_pos) {
    // Get the candidate dimensions
    uint64_t candidate_begin_position, candidate_end_position, candidate_length;
    filtering_candidates_verify_buffered_get_candidate(
        gpu_buffer_align_bpm,candidate_idx,gpu_pattern_num_tiles,
        &candidate_begin_position,&candidate_end_position,&candidate_length);
    // Get the candidate result
    uint64_t global_distance = 0, distance_link_tiles = 0, last_tile_match_position = 0;
    const bool aligned = filtering_candidates_verify_buffered_get_result(
        gpu_buffer_align_bpm,candidate_idx,gpu_pattern_num_tiles,&pattern->bpm_pattern,
        gpu_pattern_length,gpu_pattern_tile_tall,candidate_length,max_error,
        &global_distance,&distance_link_tiles,&last_tile_match_position,
        text_collection,mm_stack);
    candidate_idx += gpu_pattern_num_tiles; // Next
    // Check total distance & Compose the retrieved region
    if (aligned && global_distance <= max_error) {
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
      regions_accepted->align_match_end_column =
          BOUNDED_ADDITION(last_tile_match_position,regions_accepted->align_distance,candidate_length-1);
      ++regions_accepted;
      ++num_accepted_regions;
    } else {
      // Configure discarded candidate
      regions_discarded->status = filtering_region_verified_discarded;
      regions_discarded->text_trace_offset = UINT64_MAX; // Not retrieved yet
      regions_discarded->begin_position = candidate_begin_position;
      regions_discarded->end_position = candidate_end_position;
      // Distance Bound estimation
      regions_discarded->align_distance_min_bound = global_distance;
      regions_discarded->align_distance = global_distance + distance_link_tiles;
      // Configure regions matching (we sacrifice this information as to save memory)
      match_scaffold_init(&regions_discarded->match_scaffold);
      ++regions_discarded;
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
  PROFILE_STOP(GP_FC_RETRIEVE_BPM_BUFFER_CANDIDATE_REGIONS,PROFILE_LEVEL);
  // DEBUG
  gem_cond_debug_block(DEBUG_FILTERING_CANDIDATES) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Filtering.Candidates (verify_regions_BPM_buffer)\n");
    tab_global_inc();
    filtering_candidates_print_regions(gem_log_get_stream(),filtering_candidates,text_collection,false,false);
    tab_global_dec();
  }
  // Return number of accepted regions
  return num_accepted_regions;
}

