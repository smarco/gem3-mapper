/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   Filtering candidates module provides functions to verify filtering-regions
 *   against its corresponding region of text in the index and compute the
 *   distance of the alignment between both
 *   This "buffered" module operates in batches of filtering-regions and
 *   makes use of GPU-buffers to offload the verification/alignment of
 *   regions to a GPU
 */

#include "align/alignment.h"
#include "filtering/candidates/filtering_candidates_verify_buffered.h"
#include "filtering/region/filtering_region_verify.h"
#include "align/align_bpm_pattern.h"
#include "align/align_bpm_distance.h"

/*
 * Debug
 */
#define DEBUG_FILTERING_CANDIDATES GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Load/Store Buffered-Region
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
  filtering_region->alignment = filtering_region_buffered->alignment;
  match_scaffold_init(&filtering_region->match_scaffold);
  filtering_region->match_scaffold.alignment_regions = filtering_region_buffered->alignment_regions;
  filtering_region->match_scaffold.num_alignment_regions = filtering_region_buffered->num_alignment_regions;
  filtering_region->match_scaffold.scaffolding_coverage = filtering_region_buffered->scaffold_coverage;
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
  filtering_region_buffered->alignment = filtering_region->alignment;
  filtering_region_buffered->alignment_regions = filtering_region->match_scaffold.alignment_regions;
  filtering_region_buffered->num_alignment_regions = filtering_region->match_scaffold.num_alignment_regions;
  filtering_region_buffered->scaffold_coverage = filtering_region->match_scaffold.scaffolding_coverage;
}
/*
 * BPM-Buffered Add (Candidates Verification)
 */
uint64_t filtering_candidates_verify_buffered_add_region(
    filtering_candidates_t* const filtering_candidates,
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    const uint64_t candidate_pos,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  // Parameters
  archive_text_t* const archive_text = filtering_candidates->archive->text;
  text_collection_t* const text_collection = &filtering_candidates->text_collection;
  filtering_region_buffered_t* const filtering_region_buffer = filtering_candidates_buffered->regions_buffered;
  const uint64_t max_error = pattern->max_effective_filtering_error;
  // Prepare Alignment-Tiles
  alignment_t* const alignment = &filtering_region->alignment;
  alignment_filters_t* const filters = &pattern->alignment_filters;
  const uint64_t text_length = filtering_region->text_end_position - filtering_region->text_begin_position;
  filtering_candidates_init_alignment(filtering_candidates,
      alignment,pattern,text_length,max_error,false);
  // BPM-GPU put all candidates (tiles)
  const bool use_kmer_filter = (pattern->key_length > BPM_PREFERED_TILE_LENGTH);
  alignment_tile_t* const alignment_tiles = alignment->alignment_tiles;
  alignment_filters_tile_t* const filters_tiles = filters->tiles;
  const uint64_t num_tiles = alignment->num_tiles;
  uint64_t tile_pos, total_candidates_added=0;
  // Check pre-filters
  text_trace_t* text_trace = NULL;
  if (use_kmer_filter) {
    // Retrieve text-candidate
    filtering_region_retrieve_text(filtering_region,pattern,archive_text,text_collection);
    text_trace = text_collection_get_trace(text_collection,filtering_region->text_trace_offset);
    // Traverse all tiles
    for (tile_pos=0;tile_pos<num_tiles;++tile_pos) {
      // Fetch tiles
      alignment_tile_t* const alignment_tile = alignment_tiles + tile_pos;
      alignment_filters_tile_t* const filters_tile = filters_tiles + tile_pos;
      // K-mer filter
      if (alignment_tile->distance==ALIGN_DISTANCE_UNKNOWN) {
        alignment_tile->distance =
            alignment_verify_levenshtein_kmer_filter(
                alignment_tile,filters_tile,pattern->key,
                text_trace->text_padded,filters->mm_stack);
        if (alignment_tile->distance==ALIGN_DISTANCE_INF) {
          alignment->distance_min_bound = ALIGN_DISTANCE_INF; // Skip candidate
          break;
        }
      }
    }
  }
  // Add candidate to GPU BPM-buffer
  if (alignment->distance_min_bound!=ALIGN_DISTANCE_INF) {
    for (tile_pos=0;tile_pos<num_tiles;++tile_pos) {
      // Fetch tiles
      alignment_tile_t* const alignment_tile = alignment_tiles + tile_pos;
      // Add tile to GPU BPM-buffer
      const uint64_t candidate_text_position = filtering_region->text_begin_position + alignment_tile->text_begin_offset;
      const uint64_t candidate_length = alignment_tile->text_end_offset-alignment_tile->text_begin_offset;
      gpu_buffer_align_bpm_add_candidate(gpu_buffer_align_bpm,tile_pos,candidate_text_position,candidate_length);
      ++total_candidates_added;
    }
  }
  PROF_ADD_COUNTER(GP_ASSW_VERIFY_CANDIDATES_TILES_COPIED,num_tiles);
  // Add the filtering region to the buffer
  filtering_candidates_verify_buffered_store_region(filtering_region_buffer+candidate_pos,filtering_region);
  // Return total candidates added
  return total_candidates_added;
}
void filtering_candidates_verify_buffered_add(
    filtering_candidates_t* const filtering_candidates,
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    pattern_t* const pattern,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    uint64_t* const gpu_buffer_align_offset) {
  // Check number of pending filtering-regions
  const uint64_t num_filtering_regions = filtering_candidates_get_num_regions(filtering_candidates);
  if (num_filtering_regions==0) {
    filtering_candidates_buffered_clear(filtering_candidates_buffered);
    *gpu_buffer_align_offset = 0;
    return;
  }
  // Allocate buffered filtering regions
  filtering_candidates_buffered_allocate_regions(
      filtering_candidates,filtering_candidates_buffered,num_filtering_regions);
  filtering_region_buffered_t* const filtering_region_buffer = filtering_candidates_buffered->regions_buffered;
  *gpu_buffer_align_offset = gpu_buffer_align_bpm_get_num_candidates(gpu_buffer_align_bpm);
  // Add the pattern to the buffer (add new queries)
  gpu_buffer_align_bpm_add_pattern(gpu_buffer_align_bpm,pattern);
  gpu_buffer_align_bpm_record_query_length(gpu_buffer_align_bpm,pattern->key_length);
  // Traverse all candidates (text-space)
  filtering_region_t** const regions_in = filtering_candidates_get_regions(filtering_candidates);
  uint64_t candidate_pos, total_candidates_added = 0;
  for (candidate_pos=0;candidate_pos<num_filtering_regions;++candidate_pos) {
    // Filtering region
    filtering_region_t* const filtering_region = regions_in[candidate_pos];
    // Filter out exact-matches & key-trimmed regions
    if (filtering_region->alignment.distance_min_bound==0 || filtering_region->key_trimmed) {
      filtering_candidates_verify_buffered_store_region(
          filtering_region_buffer+candidate_pos,filtering_region);
      continue; // Next
    }
    // Add filtering-region
    total_candidates_added += filtering_candidates_verify_buffered_add_region(
        filtering_candidates,filtering_candidates_buffered,
        filtering_region,pattern,candidate_pos,gpu_buffer_align_bpm);
  }
  gpu_buffer_align_bpm_record_candidates_per_tile(gpu_buffer_align_bpm,total_candidates_added);
  PROF_ADD_COUNTER(GP_CANDIDATE_TILES,total_candidates_added);
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
  text_collection_t* const text_collection = &filtering_candidates->text_collection;
  archive_text_t* const archive_text = filtering_candidates->archive->text;
  mm_stack_t* const mm_stack = filtering_candidates->mm->mm_general;
  // Push
  mm_stack_push_state(mm_stack);
  // Get Candidate & Pattern
  uint64_t candidate_text_position;
  uint32_t candidate_length;
  gpu_buffer_align_bpm_get_candidate(gpu_buffer_align_bpm,
      candidate_idx,&candidate_text_position,&candidate_length);
  // Get Candidate Text
  const uint64_t text_trace_offset = archive_text_retrieve_collection(
      archive_text,text_collection,candidate_text_position,candidate_length,false,false); // Retrieve text(s)
  const text_trace_t* const text_trace = text_collection_get_trace(text_collection,text_trace_offset);
  const uint8_t* const text = text_trace->text; // Candidate
  uint64_t i, uncalled_bases_text = 0;
  for (i=0;i<candidate_length;++i) {
    if (text[i]==ENC_DNA_CHAR_N) ++uncalled_bases_text;
  }
  // Align BPM & Set result
  uint64_t check_tile_match_end_column, check_tile_distance;
  bpm_compute_edit_distance(
      bpm_pattern_tile,text,candidate_length,&check_tile_distance,
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
    const uint64_t global_distance) {
  // Parameters
  text_collection_t* const text_collection = &filtering_candidates->text_collection;
  archive_text_t* const archive_text = filtering_candidates->archive->text;
  // Retrieve text
  const uint64_t candidate_position = filtering_region->text_begin_position;
  const uint64_t candidate_length = filtering_region->text_end_position-filtering_region->text_begin_position;
  const uint64_t text_trace_offset = archive_text_retrieve_collection(
      archive_text,text_collection,candidate_position,candidate_length,false,false);
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
  PROF_ADD_COUNTER(GP_FC_VERIFY_CANDIDATES_BUFFERED_DDIFF,ABS((int64_t)match_distance-(int64_t)global_distance));
}
/*
 * Compute filtering-region alignment
 */
void filtering_candidates_verify_buffered_compute_alignment(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_buffered_t* const region_buffered,
    alignment_t* const alignment,
    pattern_t* const pattern) {
  // Parameters
  text_collection_t* const text_collection = &filtering_candidates->text_collection;
  archive_text_t* const archive_text = filtering_candidates->archive->text;
  mm_stack_t* const mm_stack = filtering_candidates->mm->mm_general;
  // Push
  mm_stack_push_state(mm_stack);
  // Retrieve text
  const uint64_t candidate_position = region_buffered->text_begin_position;
  const uint64_t candidate_length = region_buffered->text_end_position-region_buffered->text_begin_position;
  const uint64_t text_trace_offset = archive_text_retrieve_collection(
      archive_text,text_collection,candidate_position,candidate_length,false,false);
  text_trace_t* const text_trace = text_collection_get_trace(text_collection,text_trace_offset);
  // Myers's BPM algorithm [EditFilter]
  alignment_verify_levenshtein(
      alignment,&pattern->alignment_filters,
      pattern->key,text_trace->text_padded,
      pattern->max_effective_filtering_error);
  // Pop
  mm_stack_pop_state(mm_stack);
}
/*
 * Retrieve filtering-region alignment from the buffer
 */
uint64_t filtering_candidates_verify_buffered_retrieve_alignment(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_buffered_t* const region_buffered,
    alignment_t* const alignment,
    pattern_t* const pattern,
    const uint64_t max_error,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    const uint64_t candidate_tile_idx) {
  // Parameters
  const uint64_t num_tiles = alignment->num_tiles;
  alignment_tile_t* const alignment_tiles = alignment->alignment_tiles;
  // Traverse all tiles
  uint64_t tile_pos;
  alignment->distance_min_bound = 0;
  for (tile_pos=0;tile_pos<num_tiles;++tile_pos) {
    alignment_tile_t* const alignment_tile = alignment_tiles + tile_pos;
    alignment_filters_tile_t* const filters_tile = pattern->alignment_filters.tiles + tile_pos;
    if (alignment_tile->distance!=ALIGN_DISTANCE_UNKNOWN) {
      if (alignment_tile->distance==ALIGN_DISTANCE_INF) {
        alignment->distance_min_bound = ALIGN_DISTANCE_INF;
      } else {
        alignment->distance_min_bound += alignment_tile->distance;
      }
    } else {
      // Retrieve alignment distance/column
      uint32_t tile_distance=0, tile_match_column=0;
      gpu_buffer_align_bpm_get_result(gpu_buffer_align_bpm,
          candidate_tile_idx+tile_pos,&tile_distance,&tile_match_column);
      if (tile_distance > filters_tile->max_error) { // As CPU version
        alignment_tile->distance = ALIGN_DISTANCE_INF;
        alignment->distance_min_bound = ALIGN_DISTANCE_INF;
      }
      // Offsets
      const uint64_t tile_offset = alignment_tile->text_begin_offset;
      const uint64_t tile_end_offset = tile_match_column+1;
      const uint64_t tile_tall = filters_tile->tile_length;
      const uint64_t tile_begin_offset = BOUNDED_SUBTRACTION(tile_end_offset,tile_tall+tile_distance,0);
      alignment_tile->distance = tile_distance;
      alignment_tile->text_end_offset = tile_offset + tile_end_offset;
      alignment_tile->text_begin_offset = tile_offset + tile_begin_offset;
      // DEBUG
      #ifdef CUDA_CHECK_BUFFERED_VERIFY_CANDIDATES
      filtering_candidates_verify_buffered_check_tile_distance(
          filtering_candidates,filters_tile->bpm_pattern_tile,gpu_buffer_align_bpm,
          candidate_tile_idx+tile_pos,tile_distance,tile_match_column);
      #endif
      // Check global distance
      if (alignment->distance_min_bound != ALIGN_DISTANCE_INF) {
        alignment->distance_min_bound += alignment_tile->distance;
      }
    }
  }
  // Check global distance
  if (alignment->distance_min_bound != ALIGN_DISTANCE_INF &&
      alignment->distance_min_bound > max_error) {
    alignment->distance_min_bound = ALIGN_DISTANCE_INF;
  }
  PROF_ADD_COUNTER(GP_ASSW_VERIFY_CANDIDATES_TILES_RETRIVED,num_tiles);
  // DEBUG
  #ifdef CUDA_CHECK_BUFFERED_VERIFY_CANDIDATES
  filtering_candidates_verify_buffered_check_global_distance(
      filtering_candidates,region_buffered,
      pattern->alignment_filters.bpm_pattern,
      alignment->distance_min_bound);
  #endif
  // Return
  return num_tiles;
}
void filtering_candidates_verify_buffered_retrieve(
    filtering_candidates_t* const filtering_candidates,
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    pattern_t* const pattern,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    const uint64_t gpu_buffer_align_offset) {
  PROFILE_START(GP_FC_VERIFY_CANDIDATES_BUFFERED,PROFILE_LEVEL);
  // Check
  const uint64_t num_filtering_regions = filtering_candidates_buffered->num_regions;
  if (num_filtering_regions==0) return; // No filtering regions
  // Parameters
  const uint64_t key_length = pattern->key_length;
  const uint64_t max_error = pattern->max_effective_filtering_error;
  filtering_region_buffered_t* const filtering_region_buffer = filtering_candidates_buffered->regions_buffered;
  // Traverse all filtering regions buffered
  uint64_t region_pos, candidate_tile_idx = gpu_buffer_align_offset;
  for (region_pos=0;region_pos<num_filtering_regions;++region_pos) {
    // Retrieve region
    filtering_region_buffered_t* const region_buffered = filtering_region_buffer + region_pos;
    alignment_t* const alignment = &region_buffered->alignment;
    // Detect exact-matches
    if (region_buffered->alignment.distance_min_bound == 0) {
      filtering_region_t* const region_accepted = filtering_candidates_allocate_region(filtering_candidates);
      filtering_candidates_verify_buffered_load_region(region_accepted,region_buffered,pattern);
      region_accepted->status = filtering_region_accepted; // Accepted candidate
      PROF_INC_COUNTER(GP_ACCEPTED_REGIONS);
      continue; // Next
    }
    // Detect trimmed matches
    const uint64_t text_length = region_buffered->text_end_position - region_buffered->text_begin_position;
    if (key_length > text_length) {
      filtering_region_t trimmed_region;
      filtering_candidates_verify_buffered_load_region(&trimmed_region,region_buffered,pattern);
      if (filtering_region_verify(filtering_candidates,&trimmed_region,pattern,false)) {
        filtering_region_t* const region_accepted =
            filtering_candidates_allocate_region(filtering_candidates);
        *region_accepted = trimmed_region;
        PROF_INC_COUNTER(GP_ACCEPTED_REGIONS);
      } else {
        filtering_region_t* const region_discarded =
            filtering_candidates_allocate_discarded_region(filtering_candidates);
        *region_discarded = trimmed_region;
        PROF_INC_COUNTER(GP_DISCARDED_REGIONS);
      }
      continue; // Next
    }
    // Detect already discarded region
    if (alignment->distance_min_bound==ALIGN_DISTANCE_INF) {
      filtering_region_t* const region_discarded =
          filtering_candidates_allocate_discarded_region(filtering_candidates);
      filtering_candidates_verify_buffered_load_region(region_discarded,region_buffered,pattern);
      region_discarded->status = filtering_region_verified_discarded;
      continue; // Next
    }
    // Retrieve & compose verified region
    if (gpu_buffer_align_bpm->align_bpm_enabled) {
      candidate_tile_idx += filtering_candidates_verify_buffered_retrieve_alignment(
          filtering_candidates,region_buffered,alignment,pattern,
          max_error,gpu_buffer_align_bpm,candidate_tile_idx);
    } else {
      filtering_candidates_verify_buffered_compute_alignment(
          filtering_candidates,region_buffered,alignment,pattern);
    }
    if (alignment->distance_min_bound <= max_error) {
      filtering_region_t* const region_accepted =
          filtering_candidates_allocate_region(filtering_candidates);
      filtering_candidates_verify_buffered_load_region(region_accepted,region_buffered,pattern);
      region_accepted->status = filtering_region_accepted; // Accepted candidate
      PROF_INC_COUNTER(GP_ACCEPTED_REGIONS);
    } else {
      filtering_region_t* const region_discarded =
          filtering_candidates_allocate_discarded_region(filtering_candidates);
      alignment->distance_min_bound = ALIGN_DISTANCE_INF; // To force CPU/GPU same
      filtering_candidates_verify_buffered_load_region(region_discarded,region_buffered,pattern);
      region_discarded->status = filtering_region_verified_discarded; // Discarded candidate
      PROF_INC_COUNTER(GP_DISCARDED_REGIONS);
    }
  }
  // DEBUG
  PROFILE_STOP(GP_FC_VERIFY_CANDIDATES_BUFFERED,PROFILE_LEVEL);
  gem_cond_debug_block(DEBUG_FILTERING_CANDIDATES) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Filtering.Candidates (verify_regions_BPM_buffer)\n");
    tab_global_inc();
    filtering_candidates_print_regions(gem_log_get_stream(),filtering_candidates,false);
    tab_global_dec();
  }
}
