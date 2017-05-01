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

#include "filtering/candidates/filtering_candidates_buffered_bpm_align.h"
#include "filtering/candidates/filtering_candidates_classify.h"
#include "filtering/region/filtering_region_verify.h"
#include "align/alignment.h"
#include "align/align_bpm.h"
#include "align/align_bpm_pattern.h"
#include "align/align_bpm_distance.h"
#include "matches/scaffold/match_scaffold_levenshtein.h"
#include "matches/scaffold/match_scaffold_chain.h"

/*
 * Debug
 */
#define DEBUG_FILTERING_CANDIDATES GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PMED


/*
 * BPM-Buffered Sort
 */
int64_t filtering_region_cmp_buffered_regions(filtering_region_t** a,filtering_region_t** b) {
  // Check state
  if ((*a)->status != filtering_region_accepted && (*b)->status == filtering_region_accepted) return  1;
  if ((*a)->status == filtering_region_accepted && (*b)->status != filtering_region_accepted) return -1;
  // Compare distance
  const int64_t cmp = (int64_t)(*a)->alignment.distance_min_bound - (int64_t)(*b)->alignment.distance_min_bound;
  if (cmp) return cmp;
  // Compare position (Helps stability)
  return (int64_t)(*a)->text_begin_position - (int64_t)(*b)->text_begin_position;
}
#define VECTOR_SORT_NAME                 buffered_regions
#define VECTOR_SORT_TYPE                 filtering_region_t*
#define VECTOR_SORT_CMP(a,b)             filtering_region_cmp_buffered_regions(a,b)
#include "utils/vector_sort.h"
void filtering_candidates_buffered_bpm_align_sort_regions(
    filtering_candidates_buffered_t* const filtering_candidates_buffered) {
  buffer_sort_buffered_regions(
      filtering_candidates_buffered->regions,
      filtering_candidates_buffered->num_regions);
}
/*
 * BPM-Buffered Select
 */
void filtering_candidates_buffered_bpm_align_discard_subdominant_edit(
    filtering_candidates_t* const filtering_candidates,
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    pattern_t* const pattern) {
  // Parameters
  const uint64_t num_regions = filtering_candidates_buffered->num_regions;
  // Traverse all regions and discard subdominant
  uint64_t i, best_distance = UINT64_MAX;
  for (i=0;i<num_regions;++i) {
    // Fetch region
    filtering_region_t* const filtering_region = filtering_candidates_buffered->regions[i];
    if (filtering_region->status != filtering_region_accepted) break;
    // Check subdominant
    const uint64_t edit_distance_bound = filtering_region->alignment.distance_min_bound;
    const bool is_subdominant =
        filtering_candidates_classify_subdominant_region_edit(
            filtering_candidates,i,edit_distance_bound,best_distance);
    if (is_subdominant) break;
    // Accept & update distance
    best_distance = edit_distance_bound;
  }
  // Set remaining regions as subdominant
  for (;i<num_regions;++i) {
    filtering_candidates_buffered->regions[i]->status = filtering_region_accepted_subdominant;
  }
}
void filtering_candidates_buffered_bpm_align_discard_subdominant_swg(
    filtering_candidates_t* const filtering_candidates,
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    pattern_t* const pattern) {
  // Parameters
  search_parameters_t* const search_parameters = filtering_candidates->search_parameters;
  const uint64_t num_regions = filtering_candidates_buffered->num_regions;
  swg_penalties_t* const swg_penalties = &search_parameters->swg_penalties;
  // Traverse all regions and discard subdominant
  uint64_t i;
  int32_t best_score = 0;
  for (i=0;i<num_regions;++i) {
    // Fetch region
    filtering_region_t* const filtering_region = filtering_candidates_buffered->regions[i];
    if (filtering_region->status != filtering_region_accepted) break;
    // Check subdominant
    const uint64_t edit_distance_bound = filtering_region->alignment.distance_min_bound;
    const int32_t max_score_bound = align_swg_score_compute_max_score_bound(
        swg_penalties,edit_distance_bound,pattern->key_length);
    const bool is_subdominant =
        filtering_candidates_classify_subdominant_region_swg(
            filtering_candidates,i,max_score_bound,best_score);
    if (is_subdominant) break;
    // Accept & update distance
    const int32_t min_score_bound = align_swg_score_compute_min_score_bound(
        swg_penalties,edit_distance_bound,pattern->key_length);
    best_score = MIN(best_score,min_score_bound);
  }
  // Set remaining regions as subdominant
  for (;i<num_regions;++i) {
    filtering_candidates_buffered->regions[i]->status = filtering_region_accepted_subdominant;
  }
}
void filtering_candidates_buffered_bpm_align_discard_subdominant(
    filtering_candidates_t* const filtering_candidates,
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    pattern_t* const pattern) {
  // Select alignment model
  search_parameters_t* const search_parameters = filtering_candidates->search_parameters;
  switch (search_parameters->match_alignment_model) {
    case match_alignment_model_hamming:
    case match_alignment_model_levenshtein: {
      filtering_candidates_buffered_bpm_align_discard_subdominant_edit(
          filtering_candidates,filtering_candidates_buffered,pattern);
      break;
    }
    case match_alignment_model_gap_affine: {
      filtering_candidates_buffered_bpm_align_discard_subdominant_swg(
          filtering_candidates,filtering_candidates_buffered,pattern);
      break;
    }
    default:
      break;
  }
}
/*
 * BPM-Buffered Prepare Scaffold
 */
bool filtering_candidates_buffered_bpm_align_prepare_scaffold(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern) {
  // Parameters
  search_parameters_t* const search_parameters = filtering_candidates->search_parameters;
  const uint64_t global_min_identity = search_parameters->alignment_global_min_identity_nominal;
  match_scaffold_t* const match_scaffold = &filtering_region->match_scaffold;
  // Scaffold the alignment
  text_trace_t* const text_trace = &filtering_region->text_trace;
  match_scaffold_region_chain(
      match_scaffold,pattern,
      text_trace,NULL);/* No need of matches as no RL-enabled */
  // Check coverage
  const uint64_t distance_min_bound = filtering_region->alignment.distance_min_bound;
  const uint64_t max_coverage_bound = BOUNDED_SUBTRACTION(pattern->key_length,distance_min_bound,0);
  if (max_coverage_bound < global_min_identity ||
      match_scaffold->scaffolding_coverage < max_coverage_bound) {
    match_scaffold->scaffold_type = scaffold_region_chain;
    return true; // Add to buffer
  } else {
    match_scaffold->scaffold_type = scaffold_levenshtein;
    return false; // Done
  }
}
/*
 * BPM-Buffered Add (Candidates Alignment)
 */
void filtering_candidates_buffered_bpm_align_add_region(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    const uint64_t candidate_pos,
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align) {
  // Parameters
  archive_text_t* const archive_text = filtering_candidates->archive->text;
  alignment_t* const alignment = &filtering_region->alignment;
  const uint64_t num_tiles = alignment->num_tiles;
  // Add candidate to GPU BPM-buffer
  uint64_t tile_pos;
  for (tile_pos=0;tile_pos<num_tiles;++tile_pos) {
    // Fetch tiles
    alignment_tile_t* const alignment_tile = alignment->alignment_tiles + tile_pos;
    pattern_tile_t* const pattern_tile = pattern->pattern_tiled.tiles + tile_pos;
    // Retrieve Candidate (if needed)
    filtering_region_retrieve_text(filtering_region,
        pattern,archive_text,filtering_candidates->mm_allocator);
    // Add tile to GPU BPM-buffer
    const uint64_t match_distance = alignment_tile->distance;
    if (match_distance!=ALIGN_DISTANCE_INF) {
      const uint64_t tile_offset = pattern_tile->tile_offset;
      text_trace_t* const text_trace = &filtering_region->text_trace;
      const bool left_gap_alignment = (text_trace->strand==Forward);
      uint8_t* const text = text_trace->text_padded + alignment_tile->text_begin_offset;
      const uint64_t text_length = alignment_tile->text_end_offset - alignment_tile->text_begin_offset;
      gpu_buffer_bpm_align_add_candidate(
          gpu_buffer_bpm_align,tile_pos,tile_offset,
          text,text_length,left_gap_alignment);
    }
  }
}
void filtering_candidates_buffered_bpm_align_add_filtering_regions(
    filtering_candidates_t* const filtering_candidates,
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    pattern_t* const pattern,
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align) {
  // Parameters
  archive_text_t* const archive_text = filtering_candidates->archive->text;
  // Traverse all regions and add those accepted
  const uint64_t num_regions = filtering_candidates_buffered->num_regions;
  uint64_t candidate_pos;
  for (candidate_pos=0;candidate_pos<num_regions;++candidate_pos) {
    // Fetch region
    filtering_region_t* const filtering_region = filtering_candidates_buffered->regions[candidate_pos];
    // Filter out exact-matching regions & not-accepted regions
    if (filtering_region->status != filtering_region_accepted) continue; // Next
    if (filtering_region->alignment.distance_min_bound==0) continue; // Next
    // Retrieve Candidate (if needed)
    filtering_region_retrieve_text(filtering_region,
        pattern,archive_text,filtering_candidates->mm_allocator);
    // Prepare scaffold
    const bool scaffold_candidate =
        filtering_candidates_buffered_bpm_align_prepare_scaffold(
            filtering_candidates,filtering_region,pattern);
    if (!scaffold_candidate) continue; // Next
    // Add filtering-region
    filtering_candidates_buffered_bpm_align_add_region(
        filtering_candidates,filtering_region,
        pattern,candidate_pos,gpu_buffer_bpm_align);
  }
}
void filtering_candidates_buffered_bpm_align_add(
    filtering_candidates_t* const filtering_candidates,
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    pattern_t* const pattern,
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    uint64_t* const gpu_buffer_align_offset) {
  // Check number of pending filtering-regions
  const uint64_t num_regions = filtering_candidates_buffered->num_regions;
  search_parameters_t* const search_parameters = filtering_candidates->search_parameters;
  if (num_regions==0 || search_parameters->alignment_force_full_swg) {
    filtering_candidates_buffered->num_regions = 0;
    return;
  }
  // Sort buffered regions
  filtering_candidates_buffered_bpm_align_sort_regions(filtering_candidates_buffered);
  // Discard subdominant regions (prepare eligible regions)
  filtering_candidates_buffered_bpm_align_discard_subdominant(
      filtering_candidates,filtering_candidates_buffered,pattern);
  // Store buffer offset
  *gpu_buffer_align_offset = gpu_buffer_bpm_align_get_num_candidates(gpu_buffer_bpm_align);
  // Add the pattern to the buffer (add new queries)
  gpu_buffer_bpm_align_add_pattern(gpu_buffer_bpm_align,pattern);
  // Add all eligible filtering regions
  filtering_candidates_buffered_bpm_align_add_filtering_regions(
      filtering_candidates,filtering_candidates_buffered,pattern,gpu_buffer_bpm_align);
  // PROF
  PROF_ADD_COUNTER(GP_CANDIDATE_TILES,gpu_buffer_bpm_align->current_candidates_added);
  PROF_ADD_COUNTER(GP_BMP_DISTANCE_NUM_TILES_VERIFIED,gpu_buffer_bpm_align->current_candidates_added);
}
/*
 * BPM-Align Buffered Retrieve Check
 */
void filtering_candidates_buffered_bpm_align_retrieve_scaffold_check(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    pattern_tile_t* const pattern_tile,
    alignment_tile_t* const alignment_tile,
    match_alignment_t* const match_alignment_buffered,
    matches_t* const matches) {
  // Paramters
  const uint64_t key_offset = pattern_tile->tile_offset;
  const uint64_t text_begin = alignment_tile->text_begin_offset;
  text_trace_t* const text_trace = &filtering_region->text_trace;
  // Check
  bpm_align_matrix_t bpm_align_matrix;
  align_bpm_compute_matrix(
      &pattern_tile->bpm_pattern_tile,
      text_trace->text_padded + text_begin,
      alignment_tile->text_end_offset - alignment_tile->text_begin_offset,
      alignment_tile->distance,&bpm_align_matrix,filtering_candidates->mm_allocator);
  match_alignment_t match_alignment_check;
  match_alignment_check.match_position = text_trace->position + text_begin;
  align_bpm_backtrace_matrix(
      &pattern_tile->bpm_pattern_tile,
      pattern->key + key_offset,
      text_trace->text_padded + text_begin,
      text_trace->strand==Forward,
      &bpm_align_matrix,
      &match_alignment_check,
      matches->cigar_vector);
  const int cmp_cigars = matches_cigar_cmp(
      matches->cigar_vector,match_alignment_buffered->cigar_offset,
      match_alignment_buffered->cigar_length,matches->cigar_vector,
      match_alignment_check.cigar_offset,match_alignment_check.cigar_length);
  if (match_alignment_buffered->match_position != match_alignment_check.match_position ||
      cmp_cigars != 0) {
    fprintf(stderr,
        "Filtering.Candidates.BPM.Align.Buffered.\n"
        "\tCheck Alignment Tile-CIGAR Failed (PosBuffered=%lu) != (PosChecked=%lu)\n",
        match_alignment_buffered->match_position,match_alignment_check.match_position);
    fprintf(stderr,"\tCIGAR.GPU: ");
    match_cigar_print(stderr,matches->cigar_vector,
        match_alignment_buffered->cigar_offset,match_alignment_buffered->cigar_length);
    fprintf(stderr,"\n\tCIGAR.CPU: ");
    match_cigar_print(stderr,matches->cigar_vector,
        match_alignment_check.cigar_offset,match_alignment_check.cigar_length);
    fprintf(stderr,"\n");
  }
}
/*
 * BPM-Align Buffered Retrieve (Candidates Scaffolding)
 */
void filtering_candidates_buffered_bpm_align_compute_scaffold(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    matches_t* const matches) {
  // Parameters
  search_parameters_t* const search_parameters = filtering_candidates->search_parameters;
  alignment_t* const alignment = &filtering_region->alignment;
  // Allocate scaffold
  match_scaffold_t* const match_scaffold = &filtering_region->match_scaffold;
  const uint64_t matching_min_length = search_parameters->alignment_scaffolding_min_matching_length_nominal;
  match_scaffold_levenshtein_allocate(match_scaffold,pattern->key_length,matching_min_length);
  // Traverse all tiles from candidate
  text_trace_t* const text_trace = &filtering_region->text_trace;
  const uint64_t num_tiles = alignment->num_tiles;
  uint64_t tile_pos;
  for (tile_pos=0;tile_pos<num_tiles;++tile_pos) {
    alignment_tile_t* const alignment_tile = alignment->alignment_tiles + tile_pos;
    pattern_tile_t* const pattern_tile = pattern->pattern_tiled.tiles + tile_pos;
    match_scaffold_levenshtein_tile(
        match_scaffold,pattern,text_trace,alignment_tile,
        pattern_tile,matching_min_length,matches);
  }
}
void filtering_candidates_buffered_bpm_align_retrieve_scaffold_tile(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    pattern_tile_t* const pattern_tile,
    alignment_tile_t* const alignment_tile,
    matches_t* const matches,
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    uint64_t* const gpu_buffer_offset) {
  // Paramters
  search_parameters_t* const search_parameters = filtering_candidates->search_parameters;
  // Check match-distance
  const uint64_t match_distance = alignment_tile->distance;
  if (match_distance!=ALIGN_DISTANCE_INF) {
    // Retrieve Alignment (CIGAR)
    const uint64_t key_offset = pattern_tile->tile_offset;
    const uint64_t text_begin = alignment_tile->text_begin_offset;
    text_trace_t* const text_trace = &filtering_region->text_trace;
    uint8_t* const text = text_trace->text_padded + text_begin;
    const uint64_t match_position = text_trace->position;
    match_alignment_t match_alignment;
    match_alignment.match_position = match_position + text_begin;
    gpu_buffer_bpm_align_retrieve_alignment(
        gpu_buffer_bpm_align,text,*gpu_buffer_offset,
        &match_alignment,matches->cigar_vector);
    // CHECK
#ifdef GPU_CHECK_BPM_ALIGN
    filtering_candidates_buffered_bpm_align_retrieve_scaffold_check(
        filtering_candidates,filtering_region,
        pattern,pattern_tile,alignment_tile,
        &match_alignment,matches);
#endif
    // Store the offset (from the beginning of the text)
    // (accounting for the text-padding offset)
    const uint64_t alignment_offset = match_alignment.match_position - match_position;
    const uint64_t text_padded_left = text_trace->text_padded_left;
    gem_fatal_check_msg(text_begin + alignment_offset < text_padded_left,
        "Scaffold levenshtein. Negative coordinates because of padding");
    // Offset wrt text (without padding)
    match_alignment.match_text_offset = text_begin + alignment_offset - text_padded_left;
//      // DEBUG
//      match_alignment_print_pretty(stderr,&match_alignment,matches->cigar_vector,
//        pattern->key + key_offset,pattern_tile->bpm_pattern_tile.pattern_length,
//        text_trace->text_padded + text_begin + alignment_offset,
//        alignment_tile->text_end_offset - match_alignment.match_text_offset,mm_allocator);
    // Add the alignment to the scaffold
    const uint64_t matching_min_length = search_parameters->alignment_scaffolding_min_matching_length_nominal;
    match_scaffold_t* const match_scaffold = &filtering_region->match_scaffold;
    match_scaffold_levenshtein_compose_alignment(
        match_scaffold,&match_alignment,
        key_offset,matching_min_length,matches);
    // Increment offset
    ++(*gpu_buffer_offset);
  }
}
void filtering_candidates_buffered_bpm_align_retrieve_scaffold(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    matches_t* const matches,
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    uint64_t* const gpu_buffer_offset) {
  // Parameters
  search_parameters_t* const search_parameters = filtering_candidates->search_parameters;
  alignment_t* const alignment = &filtering_region->alignment;
  // Allocate scaffold
  match_scaffold_t* const match_scaffold = &filtering_region->match_scaffold;
  const uint64_t matching_min_length = search_parameters->alignment_scaffolding_min_matching_length_nominal;
  match_scaffold_levenshtein_allocate(match_scaffold,pattern->key_length,matching_min_length);
  // Traverse all tiles from candidate
  const uint64_t num_tiles = alignment->num_tiles;
  uint64_t tile_pos;
  for (tile_pos=0;tile_pos<num_tiles;++tile_pos) {
    pattern_tile_t* const pattern_tile = pattern->pattern_tiled.tiles + tile_pos;
    alignment_tile_t* const alignment_tile = alignment->alignment_tiles + tile_pos;
    filtering_candidates_buffered_bpm_align_retrieve_scaffold_tile(
        filtering_candidates,filtering_region,
        pattern,pattern_tile,alignment_tile,matches,
        gpu_buffer_bpm_align,gpu_buffer_offset);
  }
//  // CHECK
//  match_scaffold_print(stderr,matches,&filtering_region->match_scaffold);
//  filtering_candidates_buffered_bpm_align_compute_scaffold(
//      filtering_candidates,filtering_region,pattern,matches);
//  match_scaffold_print(stderr,matches,&filtering_region->match_scaffold);
}
void filtering_candidates_buffered_bpm_align_retrieve_region(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    matches_t* const matches,
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    uint64_t* const gpu_buffer_offset) {
  // Discarded regions
  if (filtering_region->status != filtering_region_accepted) {
    vector_insert(filtering_candidates->discarded_regions,filtering_region,filtering_region_t*);
    return;
  }
  // Exact matching regions
  if (filtering_region->alignment.distance_min_bound==0) {
    vector_insert(filtering_candidates->filtering_regions,filtering_region,filtering_region_t*);
    return;
  }
  // Candidates filter from scaffold
  if (filtering_region->match_scaffold.scaffold_type==scaffold_levenshtein) {
    vector_insert(filtering_candidates->filtering_regions,filtering_region,filtering_region_t*);
    return;
  }
  // Retrieve CIGAR & compose scaffold
  if (gpu_buffer_bpm_align->bpm_align_enabled) {
    filtering_candidates_buffered_bpm_align_retrieve_scaffold(
        filtering_candidates,filtering_region,pattern,
        matches,gpu_buffer_bpm_align,gpu_buffer_offset);
  } else {
    filtering_candidates_buffered_bpm_align_compute_scaffold(
        filtering_candidates,filtering_region,pattern,matches);
  }
  // Chains scaffolds
  match_scaffold_chain(
      &filtering_region->match_scaffold,pattern,
      &filtering_region->text_trace,false);
  filtering_region->match_scaffold.scaffold_type = scaffold_levenshtein;
  // Add to filtering regions
  vector_insert(filtering_candidates->filtering_regions,filtering_region,filtering_region_t*);
}
void filtering_candidates_buffered_bpm_align_retrieve(
    filtering_candidates_t* const filtering_candidates,
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    pattern_t* const pattern,
    matches_t* const matches,
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    const uint64_t gpu_buffer_align_offset) {
  // Check filtering regions
  const uint64_t num_regions = filtering_candidates_buffered->num_regions;
  if (num_regions==0) return;
  // Traverse all filtering regions buffered
  uint64_t candidate_pos, gpu_buffer_offset = gpu_buffer_align_offset;
  for (candidate_pos=0;candidate_pos<num_regions;++candidate_pos) {
    // Retrieve region
    filtering_region_t* const filtering_region = filtering_candidates_buffered->regions[candidate_pos];
    filtering_candidates_buffered_bpm_align_retrieve_region(
        filtering_candidates,filtering_region,pattern,
        matches,gpu_buffer_bpm_align,&gpu_buffer_offset);
  }
}



