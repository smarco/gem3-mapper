/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates_verify_buffered.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering/filtering_candidates_verify_buffered.h"
#include "align/align_bpm_pattern.h"
#include "align/align_bpm_distance.h"
#include "align/align.h"

/*
 * Debug
 */
#define DEBUG_FILTERING_CANDIDATES  GEM_DEEP_DEBUG
#define DEBUG_TILE_DIMENSIONS       false

#define VERIFY_BUFFERED_KMER_FILTER_LENGTH_THRESHOLD  256

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Kmer Filter
 */
bool filtering_candidates_verify_buffered_kmer_filter(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern) {
  // Parameters
  archive_text_t* const archive_text = archive_text;
  text_collection_t* const text_collection = text_collection;
  mm_stack_t* const mm_stack = filtering_candidates->mm_stack;
  // Retrieve text-candidate
  filtering_region_retrieve_text(filtering_region,pattern,archive_text,text_collection,mm_stack);
  const text_trace_t* const text_trace =
      text_collection_get_trace(text_collection,filtering_region->text_trace_offset);
  const uint8_t* const text = text_trace->text; // Candidate
  const uint64_t eff_text_length = filtering_region->text_end_position - filtering_region->text_begin_position;
  // Generalized Kmer-Counting filter
  const uint64_t test_positive = kmer_counting_filter(&pattern->kmer_counting,text,eff_text_length);
  // Return
  return test_positive!=ALIGN_DISTANCE_INF;
}
/*
 * BPM-Buffered Add (Candidates Verification)
 */
uint64_t filtering_candidates_verify_buffered_add(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  // Check number of pending regions
  const uint64_t pending_candidates = filtering_candidates_get_num_candidate_regions(filtering_candidates);
  if (pending_candidates==0) return 0;
  // Add the pattern to the buffer (add new queries)
  gpu_buffer_align_bpm_add_pattern(gpu_buffer_align_bpm,pattern->bpm_pattern);
  gpu_buffer_align_bpm_record_query_length(gpu_buffer_align_bpm,pattern->key_length);
  // Fetch pattern dimensions
  const uint64_t max_error = pattern->max_effective_filtering_error;
  const uint64_t gpu_pattern_length = pattern->key_length;
  const uint64_t gpu_pattern_num_tiles = pattern->bpm_pattern->gpu_num_tiles;
  const uint64_t gpu_pattern_tile_tall = pattern->bpm_pattern->gpu_entries_per_tile*gpu_bpm_pattern_get_entry_length();
  // Traverse all candidates (text-space) & add them to the buffer
  filtering_region_t* filtering_region = vector_get_mem(filtering_candidates->filtering_regions,filtering_region_t);
  uint64_t candidate_pos, tile_offset, total_candidates_added;
  for (candidate_pos=0,total_candidates_added=0;candidate_pos<pending_candidates;++candidate_pos,++filtering_region) {
    // Kmer filtering
    if (!filtering_region->key_trimmed && pattern->key_length > VERIFY_BUFFERED_KMER_FILTER_LENGTH_THRESHOLD) {
      const bool accept_candidate =
          filtering_candidates_verify_buffered_kmer_filter(filtering_candidates,filtering_region,pattern);
      if (!accept_candidate) {
        PROF_INC_COUNTER(GP_FC_KMER_COUNTER_FILTER_DISCARDED);
        continue;
      }
    }
    PROF_INC_COUNTER(GP_FC_KMER_COUNTER_FILTER_ACCEPTED);
    // Locate candidate sequence
    const uint64_t begin_position = filtering_region->text_begin_position;
    const uint64_t candidate_length = filtering_region->text_end_position - begin_position;
    if (DEBUG_TILE_DIMENSIONS) gem_slog("> Candidate #%lu (%lu nt)\n",candidate_pos,candidate_length);
    // Calculate tile dimensions
    pattern_tiled_t pattern_tiled;
    pattern_tiled_init(&pattern_tiled,
        gpu_pattern_length,gpu_pattern_tile_tall,candidate_length,max_error);
    // Initialize current tile variables
    for (tile_offset=0;tile_offset<gpu_pattern_num_tiles;++tile_offset,++total_candidates_added) {
      // BPM-GPU put candidate
      gpu_buffer_align_bpm_add_candidate(gpu_buffer_align_bpm,tile_offset,
          begin_position+pattern_tiled.tile_offset,pattern_tiled.tile_wide);
      if (DEBUG_TILE_DIMENSIONS) {
        gem_slog("  => Tile #%lu (Offsets=%lu-%lu Tall=%lu Wide=%lu)\n",
            tile_offset,pattern_tiled.tile_offset,
            pattern_tiled.tile_offset+pattern_tiled.tile_wide,
            pattern_tiled.tile_tall,pattern_tiled.tile_wide);
      }
      // Calculate next tile
      pattern_tiled_calculate_next(&pattern_tiled);
    }
  }
  gpu_buffer_align_bpm_record_candidates_per_query(gpu_buffer_align_bpm,total_candidates_added);
  // Return the final number of candidate-tiles added to the buffer
  PROF_ADD_COUNTER(GP_BMP_DISTANCE_NUM_TILES,total_candidates_added);
  PROF_ADD_COUNTER(GP_BMP_DISTANCE_NUM_TILES_VERIFIED,total_candidates_added);
  return total_candidates_added;
}
/*
 * BPM-Buffered Retrieve Checkers (Candidates Verification)
 */
void filtering_candidates_verify_buffered_check_tile_distance(
    filtering_candidates_t* const filtering_candidates,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    const uint64_t candidate_offset,
    const uint64_t pattern_tile,
    const uint32_t tile_distance,
    const uint32_t tile_match_column) {
  // Parameters
  text_collection_t* const text_collection = filtering_candidates->text_collection;
  mm_stack_t* const mm_stack = filtering_candidates->mm_stack;
  // Push
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
  const uint64_t text_trace_offset = archive_text_retrieve_collection(
      gpu_buffer_align_bpm->archive_text,text_collection,
      candidate_text_position,candidate_length,false,false,mm_stack); // Retrieve text(s)
  const text_trace_t* const text_trace = text_collection_get_trace(text_collection,text_trace_offset);
  const uint8_t* const text = text_trace->text; // Candidate
  uint64_t i, uncalled_bases_text = 0;
  for (i=0;i<candidate_length;++i) {
    if (text[i]==ENC_DNA_CHAR_N) ++uncalled_bases_text;
  }
  // Align BPM & Set result
  uint64_t check_tile_match_end_column, check_tile_distance;
  bpm_compute_edit_distance(&bpm_pattern_tile,text,candidate_length,&check_tile_distance,
      &check_tile_match_end_column,bpm_pattern_tile.pattern_length,false);
  if (tile_distance!=check_tile_distance || tile_match_column!=check_tile_match_end_column) {
    if (uncalled_bases_text == 0) {
      gem_error_msg("Filtering.Candidates.Verify.Buffered. Check verify candidate "
          "(Distance:%d!=%lu) (MatchPos:%d!=%lu) (Text.Uncalled.bases=%lu)",
          tile_distance,check_tile_distance,tile_match_column,
          check_tile_match_end_column,uncalled_bases_text);
    }
  }
  // Pop
  mm_stack_pop_state(mm_stack);
}
void filtering_candidates_verify_buffered_check_global_distance(
    filtering_candidates_t* const filtering_candidates,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    const uint64_t candidate_offset,
    const uint64_t num_tiles,
    const uint64_t global_distance,
    const uint64_t distance_link_tiles,
    bpm_pattern_t* const bpm_pattern) {
  // Parameters
  text_collection_t* const text_collection = filtering_candidates->text_collection;
  mm_stack_t* const mm_stack = filtering_candidates->mm_stack;
  // Push
  mm_stack_push_state(mm_stack);
  // Get Candidate Position
  uint64_t first_candidate_text_position, last_candidate_text_position;
  uint32_t first_candidate_length, last_candidate_length;
  gpu_buffer_align_bpm_get_candidate(gpu_buffer_align_bpm,
      candidate_offset,&first_candidate_text_position,&first_candidate_length);
  gpu_buffer_align_bpm_get_candidate(gpu_buffer_align_bpm,
      candidate_offset+num_tiles-1,&last_candidate_text_position,&last_candidate_length);
  // Get Whole-Candidate Text
  const uint64_t candidate_length = last_candidate_text_position+last_candidate_length-first_candidate_text_position;
  const uint64_t text_trace_offset = archive_text_retrieve_collection(
      gpu_buffer_align_bpm->archive_text,text_collection,
      first_candidate_text_position,candidate_length,false,false,mm_stack); // Retrieve text(s)
  const text_trace_t* const text_trace = text_collection_get_trace(text_collection,text_trace_offset);
  const uint8_t* const text = text_trace->text; // Candidate
  // Check Whole-Read
  uint64_t match_end_column, match_distance;
  bpm_compute_edit_distance(bpm_pattern,text,candidate_length,
      &match_distance,&match_end_column,bpm_pattern->pattern_length,false);
//  if (!(global_distance <= match_distance && match_distance <= global_distance+distance_link_tiles)) {
    gem_slog(">FC.Verify.Candidate.Buffered.Distance\t"
        "Whole.Read=%lu\tTileWise={bound=%lu,estimated=%lu}\tDiff=%lu\n",
        match_distance,global_distance,global_distance+distance_link_tiles,ABS(match_distance-global_distance));
//  }
  PROF_ADD_COUNTER(GP_FC_VERIFY_CANDIDATES_BUFFERED_DDIFF,ABS(match_distance-global_distance));
  // Pop
  mm_stack_pop_state(mm_stack);
}
/*
 * BPM-Buffered Retrieve (Candidates Verification)
 */
void filtering_candidates_verify_buffered_get_candidate(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    const uint64_t candidate_offset,
    const uint64_t num_tiles,
    uint64_t* const candidate_begin_position,
    uint64_t* const candidate_end_position,
    uint64_t* const candidate_length) {
  uint32_t tile_length;
  gpu_buffer_align_bpm_get_candidate(gpu_buffer_align_bpm,
      candidate_offset,candidate_begin_position,&tile_length);
  gpu_buffer_align_bpm_get_candidate(gpu_buffer_align_bpm,
      candidate_offset+(num_tiles-1),candidate_end_position,&tile_length);
  *candidate_end_position += tile_length;
  *candidate_length = *candidate_end_position - *candidate_begin_position;
}
void filtering_candidates_verify_buffered_get_result(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    bpm_pattern_t* const bpm_pattern,
    const uint64_t candidate_offset,
    const uint64_t num_tiles,
    const uint64_t pattern_length,
    const uint64_t pattern_tile_tall,
    const uint64_t candidate_length,
    const uint64_t max_error,
    region_alignment_t* const region_alignment,
    text_collection_t* const text_collection,
    mm_stack_t* const mm_stack) {
  // Calculate tile dimensions (Sum up the alignment distance of all the tiles)
  pattern_tiled_t pattern_tiled;
  pattern_tiled_init(&pattern_tiled,
      pattern_length,pattern_tile_tall,candidate_length,max_error);
  // Init region-alignment
  region_alignment->num_tiles = num_tiles;
  region_alignment_tile_t* const alignment_tiles =
      mm_stack_calloc(mm_stack,num_tiles,region_alignment_tile_t,false);
  region_alignment->alignment_tiles = alignment_tiles;
  // Put together all tiles from the candidate (join results)
  uint64_t tile_pos, global_distance=0, distance_link_tiles=0;
  for (tile_pos=0;tile_pos<num_tiles;++tile_pos) {
    // Retrieve alignment distance
    uint32_t tile_distance=0, tile_match_column=0;
    gpu_buffer_align_bpm_get_result(gpu_buffer_align_bpm,
        candidate_offset+tile_pos,&tile_distance,&tile_match_column);
    pattern_tiled.tile_distance = tile_distance;
    pattern_tiled.tile_match_column = tile_match_column;
    // Compute distance
    global_distance += pattern_tiled.tile_distance;
    distance_link_tiles += pattern_tiled_bound_matching_path(&pattern_tiled);
    alignment_tiles[tile_pos].match_distance = pattern_tiled.tile_distance;
    const uint64_t tile_end_offset = pattern_tiled.tile_match_column+1;
    alignment_tiles[tile_pos].text_end_offset = pattern_tiled.tile_offset + tile_end_offset;
    const uint64_t tile_begin_offset = BOUNDED_SUBTRACTION(tile_end_offset,
        pattern_tiled.tile_tall+pattern_tiled.tile_distance,0);
    alignment_tiles[tile_pos].text_begin_offset = pattern_tiled.tile_offset + tile_begin_offset;
    // Calculate next tile
    pattern_tiled_calculate_next(&pattern_tiled);
    // DEBUG
    #ifdef CUDA_CHECK_BUFFERED_VERIFY_CANDIDATES
    filtering_candidates_verify_buffered_check_tile_distance(
        gpu_buffer_align_bpm,candidate_offset,pattern_tile,
        tile_distance,tile_match_column,text_collection,mm_stack);
    #endif
  }
  // DEBUG
  #ifdef CUDA_CHECK_BUFFERED_VERIFY_CANDIDATES
  filtering_candidates_verify_buffered_check_global_distance(
      gpu_buffer_align_bpm,candidate_offset,num_tiles,*global_distance,
      *distance_link_tiles,bpm_pattern,text_collection,mm_stack);
  #endif
  // Setup alignment result
  region_alignment->distance_min_bound = global_distance;
  global_distance += distance_link_tiles;
  region_alignment->distance_max_bound = MIN(global_distance,max_error);
}
uint64_t filtering_candidates_verify_buffered_retrieve(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    const uint64_t candidate_offset_begin,
    const uint64_t candidate_offset_end) {
  if (gem_expect_false(candidate_offset_begin==candidate_offset_end)) return 0;
  PROFILE_START(GP_FC_VERIFY_CANDIDATES_BUFFERED,PROFILE_LEVEL);
  // Parameters
  text_collection_t* const text_collection = filtering_candidates->text_collection;
  mm_stack_t* const mm_stack = filtering_candidates->mm_stack;
  const uint64_t max_error = pattern->max_effective_filtering_error;
  const uint64_t gpu_pattern_length = pattern->key_length;
  const uint64_t gpu_pattern_num_tiles = pattern->bpm_pattern->gpu_num_tiles;
  const uint64_t gpu_pattern_tile_tall = pattern->bpm_pattern->gpu_entries_per_tile*gpu_bpm_pattern_get_entry_length();
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
    region_alignment_t region_alignment;
    filtering_candidates_verify_buffered_get_result(
        gpu_buffer_align_bpm,pattern->bpm_pattern,candidate_idx,gpu_pattern_num_tiles,
        gpu_pattern_length,gpu_pattern_tile_tall,candidate_length,max_error,
        &region_alignment,text_collection,mm_stack);
    candidate_idx += gpu_pattern_num_tiles; // Next
    // Check total distance & Compose the retrieved region
    if (region_alignment.distance_min_bound <= max_error) {
      // Configure accepted candidate
      regions_accepted->status = filtering_region_accepted;
      regions_accepted->text_trace_offset = UINT64_MAX; // Not retrieved yet
      regions_accepted->text_begin_position = candidate_begin_position;
      regions_accepted->text_end_position = candidate_end_position;
      // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO Store & Retrive location data
      regions_accepted->text_base_begin_offset = 0;
      regions_accepted->text_base_end_offset = candidate_end_position-candidate_begin_position;
      regions_accepted->key_trim_left = 0;
      regions_accepted->key_trim_right = 0;
      // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
      regions_accepted->region_alignment = region_alignment;
      match_scaffold_init(&regions_accepted->match_scaffold); // We sacrifice this information as to save memory
      // Next
      ++regions_accepted;
      ++num_accepted_regions;
      PROF_INC_COUNTER(GP_ACCEPTED_REGIONS);
    } else {
      // Configure discarded candidate
      regions_discarded->status = filtering_region_verified_discarded;
      regions_discarded->text_trace_offset = UINT64_MAX; // Not retrieved yet
      regions_discarded->text_begin_position = candidate_begin_position;
      regions_discarded->text_end_position = candidate_end_position;
      regions_discarded->region_alignment = region_alignment;
      match_scaffold_init(&regions_discarded->match_scaffold); // We sacrifice this information as to save memory
      // Next
      ++regions_discarded;
      PROF_INC_COUNTER(GP_DISCARDED_REGIONS);
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
  PROFILE_STOP(GP_FC_VERIFY_CANDIDATES_BUFFERED,PROFILE_LEVEL);
  // DEBUG
  gem_cond_debug_block(DEBUG_FILTERING_CANDIDATES) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Filtering.Candidates (verify_regions_BPM_buffer)\n");
    tab_global_inc();
    filtering_candidates_print_regions(gem_log_get_stream(),filtering_candidates,false,false);
    tab_global_dec();
  }
  // Return number of accepted regions
  return num_accepted_regions;
}
/*
 * Display/Benchmark
 */
void filtering_candidates_verify_buffered_print_benchmark_tile(
    FILE* const stream,filtering_region_t* const filtering_region,
    pattern_t* const pattern,bpm_pattern_t* const bpm_pattern_tiles,
    const uint64_t max_error,const uint64_t tile_pos,
    const uint64_t tile_tall,const bool print_tile,
    text_collection_t* const text_collection) {
  // Parameters
  const uint64_t begin_position = filtering_region->text_begin_position;
  const uint64_t candidate_length = filtering_region->text_end_position - begin_position;
  const text_trace_t* const text_trace = text_collection_get_trace(text_collection,filtering_region->text_trace_offset);
  const uint8_t* const text = text_trace->text;
  // Calculate tile dimensions
  pattern_tiled_t pattern_tiled;
  pattern_tiled_init(&pattern_tiled,
      pattern->key_length,tile_tall,candidate_length,max_error);
  // Initialize current tile variables
  uint64_t tile_offset, total_candidates_added, acc_tall, i;
  for (tile_offset=0,acc_tall=0;tile_offset<=tile_pos;++tile_offset,++total_candidates_added) {
    if (tile_offset==tile_pos) {
      // Print tile key
      if (print_tile) {
        const uint8_t* const key = pattern->key;
        for (i=0;i<pattern_tiled.tile_tall;++i) {
          fprintf(stream,"%c",dna_decode(key[acc_tall+i]));
        }
      }
      // (Re)Align the tile
      bpm_compute_edit_distance(bpm_pattern_tiles+tile_offset,
          text+pattern_tiled.tile_offset,pattern_tiled.tile_wide,
          &pattern_tiled.tile_distance,&pattern_tiled.tile_match_column,
          pattern_tiled.tile_tall,false);
      // Print tile candidate (position, wide, distance)
      fprintf(stream,"\tchrX:+:%lu:%lu:%ld",begin_position,pattern_tiled.tile_wide,
          pattern_tiled.tile_distance==ALIGN_DISTANCE_INF ? -1 : pattern_tiled.tile_distance);
    }
    // Calculate next tile
    acc_tall += pattern_tiled.tile_tall;
    pattern_tiled_calculate_next(&pattern_tiled);
  }
}
void filtering_candidates_verify_buffered_print_benchmark(
    FILE* const stream,
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern) {
  // Parameters
  archive_text_t* const archive_text = filtering_candidates->archive->text;
  mm_stack_t* const mm_stack = filtering_candidates->mm_stack;
  // Allocate text collection
  text_collection_t text_collection;
  text_collection_init(&text_collection);
  // Check number of pending regions
  const uint64_t pending_candidates = filtering_candidates_get_num_candidate_regions(filtering_candidates);
  if (pending_candidates==0) return;
  // Fetch pattern dimensions
  const uint64_t max_error = pattern->max_effective_filtering_error;
  const uint64_t gpu_pattern_num_tiles = pattern->bpm_pattern->gpu_num_tiles;
  const uint64_t gpu_pattern_tile_tall = pattern->bpm_pattern->gpu_entries_per_tile*gpu_bpm_pattern_get_entry_length();
  // Compile pattern tiles
  bpm_pattern_t* const bpm_pattern = pattern->bpm_pattern;
  bpm_pattern_t* const bpm_pattern_tiles = bpm_pattern_compile_tiles(
      bpm_pattern,2*bpm_pattern->gpu_entries_per_tile,max_error,mm_stack);
  // Traverse all candidates (text-space)
  uint64_t tile_pos, candidate_pos;
  for (tile_pos=0;tile_pos<gpu_pattern_num_tiles;++tile_pos) {
    filtering_region_t* candidate_region = vector_get_mem(filtering_candidates->filtering_regions,filtering_region_t);
    for (candidate_pos=0;candidate_pos<pending_candidates;++candidate_pos,++candidate_region) {
      // Retrieve Candidate (if needed)
      filtering_region_retrieve_text(candidate_region,pattern,archive_text,&text_collection,mm_stack);
      // Print candidate tile
      filtering_candidates_verify_buffered_print_benchmark_tile(
          stream,candidate_region,pattern,bpm_pattern_tiles,max_error,tile_pos,
          gpu_pattern_tile_tall,candidate_pos==0,&text_collection);
    }
    fprintf(stream,"\n");
  }
  // Restore & Free
  filtering_region_t* candidate_region = vector_get_mem(filtering_candidates->filtering_regions,filtering_region_t);
  for (candidate_pos=0;candidate_pos<pending_candidates;++candidate_pos,++candidate_region) {
    candidate_region->text_trace_offset = UINT64_MAX;
  }
  text_collection_destroy(&text_collection);
}
