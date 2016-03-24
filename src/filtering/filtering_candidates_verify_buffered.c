/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates_verify_buffered.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering/filtering_candidates_verify_buffered.h"
#include "filtering/filtering_region_verify.h"
#include "align/align_bpm_pattern.h"
#include "align/align_bpm_distance.h"
#include "align/align.h"

/*
 * Debug
 */
#define DEBUG_FILTERING_CANDIDATES GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Constants
 */
#define VERIFY_BUFFERED_KMER_FILTER_LENGTH_THRESHOLD  256

/*
 * Kmer Filter
 */
void filtering_candidates_verify_buffered_kmer_filter(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    text_collection_t* const text_collection) {
  // Parameters
  archive_text_t* const archive_text = filtering_candidates->archive->text;
  mm_stack_t* const mm_stack = filtering_candidates->mm_stack;
  // Compile kmer-filter
  if (!pattern->kmer_counting.enabled) {
    kmer_counting_compile(&pattern->kmer_counting,pattern->key,
        pattern->key_length,pattern->max_effective_filtering_error,mm_stack);
  }
  // Retrieve text-candidate
  filtering_region_retrieve_text(filtering_region,pattern,archive_text,text_collection,mm_stack);
  const text_trace_t* const text_trace =
      text_collection_get_trace(text_collection,filtering_region->text_trace_offset);
  const uint8_t* const text = text_trace->text; // Candidate
  const uint64_t eff_text_length = filtering_region->text_end_position - filtering_region->text_begin_position;
  // Generalized Kmer-Counting filter
  const uint64_t test_positive = kmer_counting_filter(&pattern->kmer_counting,text,eff_text_length);
  if (test_positive==ALIGN_DISTANCE_INF) {
    filtering_region->region_alignment.distance_min_bound = ALIGN_DISTANCE_INF;
    PROF_INC_COUNTER(GP_FC_KMER_COUNTER_FILTER_DISCARDED);
  } else {
    PROF_INC_COUNTER(GP_FC_KMER_COUNTER_FILTER_ACCEPTED);
  }
}
/*
 * Add filtering region to buffer
 */
void filtering_candidates_verify_buffered_load_region(
    filtering_region_t* const filtering_region,
    filtering_region_buffered_t* const filtering_region_buffered,
    pattern_t* const pattern) {
  /* Source Region Offset */
  filtering_region->text_source_region_offset = filtering_region_buffered->text_source_region_offset;
  filtering_region->key_source_region_offset = filtering_region_buffered->key_source_region_offset;
  /* Text */
  filtering_region->text_trace_offset = UINT64_MAX; // Not retrieved yet
  filtering_region->text_begin_position = filtering_region_buffered->text_begin_position;
  filtering_region->text_end_position = filtering_region_buffered->text_end_position;
  /* Key */
  filtering_region_compute_key_trims(filtering_region,pattern);
  /* Alignment */
  filtering_region->max_error = pattern->max_effective_filtering_error;
  filtering_region->max_bandwidth = pattern->max_effective_bandwidth;
  filtering_region->region_alignment = filtering_region_buffered->region_alignment;
  match_scaffold_init(&filtering_region->match_scaffold); // We sacrifice this information as to save memory
  filtering_region->match_scaffold.scaffold_regions = filtering_region_buffered->scaffold_regions;
  filtering_region->match_scaffold.num_scaffold_regions = filtering_region_buffered->num_scaffold_regions;
}
void filtering_candidates_verify_buffered_store_region(
    filtering_region_buffered_t* const filtering_region_buffered,
    filtering_region_t* const filtering_region) {
  // Source Region Offset
  filtering_region_buffered->text_source_region_offset = filtering_region->text_source_region_offset;
  filtering_region_buffered->key_source_region_offset = filtering_region->key_source_region_offset;
  // Text
  filtering_region_buffered->text_begin_position = filtering_region->text_begin_position;
  filtering_region_buffered->text_end_position = filtering_region->text_end_position;
  // Alignment
  filtering_region_buffered->region_alignment = filtering_region->region_alignment;
  filtering_region_buffered->scaffold_regions = filtering_region->match_scaffold.scaffold_regions;
  filtering_region_buffered->num_scaffold_regions = filtering_region->match_scaffold.num_scaffold_regions;
}
/*
 * BPM-Buffered Add (Candidates Verification)
 */
void filtering_candidates_verify_buffered_add(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    uint64_t* const gpu_buffer_align_offset,
    filtering_region_buffered_t** const filtering_region_buffered,
    uint64_t* const gpu_num_filtering_regions) {
  // Check number of pending filtering-regions
  const uint64_t num_filtering_regions = filtering_candidates_get_num_candidate_regions(filtering_candidates);
  if (num_filtering_regions==0) {
    *filtering_region_buffered = NULL;
    *gpu_num_filtering_regions = 0;
    *gpu_buffer_align_offset = 0;
    return;
  }
  // Allocate buffered filtering regions
  mm_stack_t* const mm_stack = gpu_buffer_align_bpm->mm_stack;
  filtering_region_buffered_t* const filtering_region_buffer =
      mm_stack_calloc(mm_stack,num_filtering_regions,filtering_region_buffered_t,false);
  *filtering_region_buffered = filtering_region_buffer;
  *gpu_num_filtering_regions = num_filtering_regions;
  *gpu_buffer_align_offset = gpu_buffer_align_bpm_get_num_candidates(gpu_buffer_align_bpm);
  // Add the pattern to the buffer (add new queries)
  gpu_buffer_align_bpm_add_pattern(gpu_buffer_align_bpm,pattern->bpm_pattern,pattern->bpm_pattern_tiles);
  gpu_buffer_align_bpm_record_query_length(gpu_buffer_align_bpm,pattern->key_length);
  // Traverse all candidates (text-space)
  bpm_pattern_t* const bpm_pattern = pattern->bpm_pattern;
  bpm_pattern_t* const bpm_pattern_tiles = pattern->bpm_pattern_tiles;
  filtering_region_t* filtering_region = vector_get_mem(filtering_candidates->filtering_regions,filtering_region_t);
  uint64_t candidate_pos, total_candidates_added = 0;
  for (candidate_pos=0;candidate_pos<num_filtering_regions;++candidate_pos,++filtering_region) {
    // Filter out key-trimmed regions
    if (filtering_region->key_trimmed) {
      filtering_candidates_verify_buffered_store_region(
          filtering_region_buffer+candidate_pos,filtering_region);
      continue; // Next
    }
    // Prepare Alignment-Tiles
    filtering_region_alignment_prepare(filtering_region,bpm_pattern,bpm_pattern_tiles,mm_stack);
    // Kmer filtering
    if (pattern->key_length > VERIFY_BUFFERED_KMER_FILTER_LENGTH_THRESHOLD) {
      filtering_candidates_verify_buffered_kmer_filter(filtering_candidates,
          filtering_region,pattern,gpu_buffer_align_bpm->text_collection);
      if (filtering_region->region_alignment.distance_min_bound==ALIGN_DISTANCE_INF) {
        filtering_candidates_verify_buffered_store_region(
            filtering_region_buffer+candidate_pos,filtering_region);
        continue; // Next
      }
    }
    // BPM-GPU put all candidates (tiles)
    const uint64_t num_pattern_tiles = bpm_pattern_tiles->num_pattern_tiles;
    region_alignment_t* const region_alignment = &filtering_region->region_alignment;
    region_alignment_tile_t* const alignment_tiles = region_alignment->alignment_tiles;
    uint64_t tile_pos;
    for (tile_pos=0;tile_pos<num_pattern_tiles;++tile_pos) {
      region_alignment_tile_t* const alignment_tile = alignment_tiles + tile_pos;
      const uint64_t candidate_text_position = filtering_region->text_begin_position + alignment_tile->text_begin_offset;
      const uint64_t candidate_length = alignment_tile->text_end_offset-alignment_tile->text_begin_offset;
      gpu_buffer_align_bpm_add_candidate(gpu_buffer_align_bpm,tile_pos,candidate_text_position,candidate_length);
    }
    total_candidates_added += num_pattern_tiles;
    PROF_ADD_COUNTER(GP_ASSW_VERIFY_CANDIDATES_TILES_COPIED,num_pattern_tiles);
    // Add the filtering region to the buffer
    filtering_candidates_verify_buffered_store_region(filtering_region_buffer+candidate_pos,filtering_region);
  }
  gpu_buffer_align_bpm_record_candidates_per_tile(gpu_buffer_align_bpm,num_filtering_regions);
  PROF_ADD_COUNTER(GP_BMP_DISTANCE_NUM_TILES,total_candidates_added);
  PROF_ADD_COUNTER(GP_BMP_DISTANCE_NUM_TILES_VERIFIED,total_candidates_added);
}
/*
 * BPM-Buffered Retrieve Checkers (Candidates Verification)
 */
void filtering_candidates_verify_buffered_check_tile_distance(
    filtering_candidates_t* const filtering_candidates,
    bpm_pattern_t* const bpm_pattern_tile,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    const uint64_t candidate_idx,
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
  gpu_buffer_align_bpm_get_candidate(gpu_buffer_align_bpm,
      candidate_idx,&candidate_text_position,&candidate_length);
  // Get Candidate Text
  const uint64_t text_trace_offset = archive_text_retrieve_collection(gpu_buffer_align_bpm->archive_text,
      text_collection,candidate_text_position,candidate_length,false,false,mm_stack); // Retrieve text(s)
  const text_trace_t* const text_trace = text_collection_get_trace(text_collection,text_trace_offset);
  const uint8_t* const text = text_trace->text; // Candidate
  uint64_t i, uncalled_bases_text = 0;
  for (i=0;i<candidate_length;++i) {
    if (text[i]==ENC_DNA_CHAR_N) ++uncalled_bases_text;
  }
  // Align BPM & Set result
  uint64_t check_tile_match_end_column, check_tile_distance;
  bpm_compute_edit_distance(bpm_pattern_tile,text,candidate_length,&check_tile_distance,
      &check_tile_match_end_column,bpm_pattern_tile->pattern_length,false);
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
    filtering_region_buffered_t* const filtering_region,
    bpm_pattern_t* const bpm_pattern,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    const uint64_t global_distance) {
  // Parameters
  text_collection_t* const text_collection = filtering_candidates->text_collection;
  mm_stack_t* const mm_stack = filtering_candidates->mm_stack;
  // Push
  mm_stack_push_state(mm_stack);
  // Retrieve text
  const uint64_t candidate_position = filtering_region->text_begin_position;
  const uint64_t candidate_length = filtering_region->text_end_position-filtering_region->text_begin_position;
  const uint64_t text_trace_offset = archive_text_retrieve_collection(
      gpu_buffer_align_bpm->archive_text,text_collection,
      candidate_position,candidate_length,false,false,mm_stack);
  const text_trace_t* const text_trace = text_collection_get_trace(text_collection,text_trace_offset);
  // Check Whole-Read
  uint64_t match_end_column, match_distance;
  bpm_compute_edit_distance(bpm_pattern,text_trace->text,text_trace->text_length,
      &match_distance,&match_end_column,bpm_pattern->pattern_length,false);
  //  if (!(global_distance <= match_distance && match_distance <= global_distance+distance_link_tiles)) {
  //  gem_slog(">FC.Verify.Candidate.Buffered.Distance\t"
  //      "Whole.Read=%lu\tTileWise={bound=%lu}\tDiff=%lu\n",
  //      match_distance,global_distance,ABS(match_distance-global_distance));
  //  }
  PROF_ADD_COUNTER(GP_FC_VERIFY_CANDIDATES_BUFFERED_DDIFF,ABS(match_distance-global_distance));
  // Pop
  mm_stack_pop_state(mm_stack);
}
/*
 * Retrieve filtering region from the buffer
 */
void filtering_candidates_verify_buffered_retrieve_region_alignment(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_buffered_t* const region_buffered,
    region_alignment_t* const region_alignment,
    bpm_pattern_t* const bpm_pattern,
    bpm_pattern_t* const bpm_pattern_tiles,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    uint64_t candidate_idx) {
  // Traverse all tiles
  region_alignment_tile_t* const alignment_tiles = region_alignment->alignment_tiles;
  const uint64_t num_tiles = region_alignment->num_tiles;
  uint64_t tile_pos, global_distance=0;
  for (tile_pos=0;tile_pos<num_tiles;++tile_pos) {
    bpm_pattern_t* const bpm_pattern_tile = bpm_pattern_tiles + tile_pos;
    region_alignment_tile_t* const alignment_tile = alignment_tiles + tile_pos;
    // Retrieve alignment distance
    uint32_t tile_distance=0, tile_match_column=0;
    gpu_buffer_align_bpm_get_result(gpu_buffer_align_bpm,candidate_idx,&tile_distance,&tile_match_column);
    const uint64_t tile_offset = alignment_tile->text_begin_offset;
    const uint64_t tile_end_offset = tile_match_column+1;
    const uint64_t tile_tall = bpm_pattern_tile->pattern_length;
    const uint64_t tile_begin_offset = BOUNDED_SUBTRACTION(tile_end_offset,tile_tall+tile_distance,0);
    alignment_tile->match_distance = tile_distance;
    alignment_tile->text_end_offset = tile_offset + tile_end_offset;
    alignment_tile->text_begin_offset = tile_offset + tile_begin_offset;
    global_distance += tile_distance;
    // DEBUG
    #ifdef CUDA_CHECK_BUFFERED_VERIFY_CANDIDATES
    filtering_candidates_verify_buffered_check_tile_distance(filtering_candidates,
        bpm_pattern_tile,gpu_buffer_align_bpm,candidate_idx,tile_distance,tile_match_column);
    #endif
    // Next
    ++candidate_idx;
  }
  region_alignment->distance_min_bound = global_distance;
  PROF_ADD_COUNTER(GP_ASSW_VERIFY_CANDIDATES_TILES_RETRIVED,num_tiles);
  // DEBUG
  #ifdef CUDA_CHECK_BUFFERED_VERIFY_CANDIDATES
  filtering_candidates_verify_buffered_check_global_distance(filtering_candidates,
      region_buffered,bpm_pattern,gpu_buffer_align_bpm,global_distance);
  #endif
}
void filtering_candidates_verify_buffered_retrieve(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    const uint64_t candidate_offset_begin,
    filtering_region_buffered_t* const filtering_region_buffered,
    uint64_t const num_filtering_regions) {
  PROFILE_START(GP_FC_VERIFY_CANDIDATES_BUFFERED,PROFILE_LEVEL);
  // Parameters
  const uint64_t key_length = pattern->key_length;
  const uint64_t max_error = pattern->max_effective_filtering_error;
  bpm_pattern_t* const bpm_pattern = pattern->bpm_pattern;
  bpm_pattern_t* const bpm_pattern_tiles = pattern->bpm_pattern_tiles;
  // Prepare filtering-regions vectors
  vector_clear(filtering_candidates->verified_regions);
  vector_clear(filtering_candidates->filtering_regions);
  vector_clear(filtering_candidates->discarded_regions);
  if (num_filtering_regions==0) return; // No filtering regions
  vector_reserve(filtering_candidates->filtering_regions,num_filtering_regions,false);
  vector_reserve(filtering_candidates->discarded_regions,num_filtering_regions,false);
  filtering_region_t* regions_accepted = vector_get_mem(filtering_candidates->filtering_regions,filtering_region_t);
  filtering_region_t* regions_discarded = vector_get_mem(filtering_candidates->discarded_regions,filtering_region_t);
  // Traverse all filtering regions buffered
  uint64_t region_pos, candidate_idx = candidate_offset_begin;
  for (region_pos=0;region_pos<num_filtering_regions;++region_pos) {
    filtering_region_buffered_t* const region_buffered = filtering_region_buffered + region_pos;
    region_alignment_t* const region_alignment = &region_buffered->region_alignment;
    // Retrieve trimmed-key region
    const uint64_t text_length = region_buffered->text_end_position - region_buffered->text_begin_position;
    if (key_length > text_length) {
      filtering_region_t trimmed_region;
      filtering_candidates_verify_buffered_load_region(&trimmed_region,region_buffered,pattern);
      if (filtering_region_verify(filtering_candidates,&trimmed_region,pattern,false)) {
        *regions_accepted = trimmed_region;
        ++regions_accepted;
        PROF_INC_COUNTER(GP_ACCEPTED_REGIONS);
      } else {
        *regions_discarded = trimmed_region;
        ++regions_discarded;
        PROF_INC_COUNTER(GP_DISCARDED_REGIONS);
      }
      continue; // Next
    }
    // Retrieve already discarded region
    if (region_alignment->distance_min_bound==ALIGN_DISTANCE_INF) {
      filtering_candidates_verify_buffered_load_region(regions_discarded,region_buffered,pattern);
      regions_discarded->status = filtering_region_verified_discarded;
      ++regions_discarded;
      continue; // Next
    }
    // Retrieve & compose verified region
    filtering_candidates_verify_buffered_retrieve_region_alignment(
        filtering_candidates,region_buffered,region_alignment,
        bpm_pattern,bpm_pattern_tiles,gpu_buffer_align_bpm,candidate_idx);
    candidate_idx += region_alignment->num_tiles; // Skip tiles
    if (region_alignment->distance_min_bound <= max_error) {
      filtering_candidates_verify_buffered_load_region(regions_accepted,region_buffered,pattern);
      regions_accepted->status = filtering_region_accepted; // Accepted candidate
      ++regions_accepted;
      PROF_INC_COUNTER(GP_ACCEPTED_REGIONS);
    } else {
      region_alignment->distance_min_bound = ALIGN_DISTANCE_INF; // To force CPU/GPU same
      filtering_candidates_verify_buffered_load_region(regions_discarded,region_buffered,pattern);
      regions_discarded->status = filtering_region_verified_discarded; // Discarded candidate
      ++regions_discarded;
      PROF_INC_COUNTER(GP_DISCARDED_REGIONS);
    }
  }
  // Update Accepted/Discarded/Verified
  vector_update_used(filtering_candidates->filtering_regions,regions_accepted);
  vector_update_used(filtering_candidates->discarded_regions,regions_discarded);
  PROFILE_STOP(GP_FC_VERIFY_CANDIDATES_BUFFERED,PROFILE_LEVEL);
  // DEBUG
  gem_cond_debug_block(DEBUG_FILTERING_CANDIDATES) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Filtering.Candidates (verify_regions_BPM_buffer)\n");
    tab_global_inc();
    filtering_candidates_print_regions(gem_log_get_stream(),filtering_candidates,false);
    tab_global_dec();
  }
}
///*
// * Display/Benchmark
// */
//void filtering_candidates_verify_buffered_print_benchmark_tile(
//    FILE* const stream,
//    filtering_region_t* const filtering_region,
//    pattern_t* const pattern,
//    bpm_pattern_t* const bpm_pattern_tiles,
//    const uint64_t max_error,
//    const uint64_t tile_pos,
//    const uint64_t tile_tall,
//    const bool print_tile,
//    text_collection_t* const text_collection) {
//  // Parameters
//  const uint64_t begin_position = filtering_region->text_begin_position;
//  const uint64_t candidate_length = filtering_region->text_end_position - begin_position;
//  const text_trace_t* const text_trace = text_collection_get_trace(text_collection,filtering_region->text_trace_offset);
//  const uint8_t* const text = text_trace->text;
//  // Calculate tile dimensions
//  pattern_tiled_t pattern_tiled;
//  pattern_tiled_init(&pattern_tiled,
//      pattern->key_length,tile_tall,candidate_length,max_error);
//  // Initialize current tile variables
//  uint64_t tile_offset, total_candidates_added, acc_tall, i;
//  for (tile_offset=0,acc_tall=0;tile_offset<=tile_pos;++tile_offset,++total_candidates_added) {
//    if (tile_offset==tile_pos) {
//      // Print tile key
//      if (print_tile) {
//        const uint8_t* const key = pattern->key;
//        for (i=0;i<pattern_tiled.tile_tall;++i) {
//          fprintf(stream,"%c",dna_decode(key[acc_tall+i]));
//        }
//      }
//      // (Re)Align the tile
//      bpm_compute_edit_distance(bpm_pattern_tiles+tile_offset,
//          text+pattern_tiled.tile_offset,pattern_tiled.tile_wide,
//          &pattern_tiled.tile_distance,&pattern_tiled.tile_match_column,
//          pattern_tiled.tile_tall,false);
//      // Print tile candidate (position, wide, distance)
//      fprintf(stream,"\tchrX:+:%lu:%lu:%ld",begin_position,pattern_tiled.tile_wide,
//          pattern_tiled.tile_distance==ALIGN_DISTANCE_INF ? -1 : pattern_tiled.tile_distance);
//    }
//    // Calculate next tile
//    acc_tall += pattern_tiled.tile_tall;
//    pattern_tiled_calculate_next(&pattern_tiled);
//  }
//}
//void filtering_candidates_verify_buffered_print_benchmark(
//    FILE* const stream,
//    filtering_candidates_t* const filtering_candidates,
//    pattern_t* const pattern) {
//  // Parameters
//  archive_text_t* const archive_text = filtering_candidates->archive->text;
//  mm_stack_t* const mm_stack = filtering_candidates->mm_stack;
//  // Allocate text collection
//  text_collection_t text_collection;
//  text_collection_init(&text_collection);
//  // Check number of pending regions
//  const uint64_t pending_candidates = filtering_candidates_get_num_candidate_regions(filtering_candidates);
//  if (pending_candidates==0) return;
//  // Traverse all candidates (text-space)
//  const uint64_t max_error = pattern->max_effective_filtering_error;
//  bpm_pattern_t* const bpm_pattern = pattern->bpm_pattern;
//  bpm_pattern_t* const bpm_pattern_tiles = pattern->bpm_pattern_tiles;
//  uint64_t tile_pos, candidate_pos;
//  for (tile_pos=0;tile_pos<gpu_pattern_num_tiles;++tile_pos) {
//    filtering_region_t* candidate_region = vector_get_mem(filtering_candidates->filtering_regions,filtering_region_t);
//    for (candidate_pos=0;candidate_pos<pending_candidates;++candidate_pos,++candidate_region) {
//      // Retrieve Candidate (if needed)
//      filtering_region_retrieve_text(candidate_region,pattern,archive_text,&text_collection,mm_stack);
//      // Print candidate tile
//      filtering_candidates_verify_buffered_print_benchmark_tile(
//          stream,candidate_region,pattern,bpm_pattern_tiles,max_error,tile_pos,
//          gpu_pattern_tile_tall,candidate_pos==0,&text_collection);
//    }
//    fprintf(stream,"\n");
//  }
//  // Restore & Free
//  filtering_region_t* candidate_region = vector_get_mem(filtering_candidates->filtering_regions,filtering_region_t);
//  for (candidate_pos=0;candidate_pos<pending_candidates;++candidate_pos,++candidate_region) {
//    candidate_region->text_trace_offset = UINT64_MAX;
//  }
//  text_collection_destroy(&text_collection);
//}
