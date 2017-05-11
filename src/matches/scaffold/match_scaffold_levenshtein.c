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
 */

#include "align/pattern/pattern.h"
#include "align/alignment.h"
#include "matches/scaffold/match_scaffold_levenshtein.h"
#include "matches/scaffold/match_scaffold_compose.h"
#include "matches/scaffold/match_scaffold_chain.h"
#include "align/align_bpm.h"

/*
 * Debug
 */
#define DEBUG_MATCH_SCAFFOLD_EDIT   false // GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PLOW

/*
 * Setup
 */
void match_scaffold_levenshtein_allocate(
    match_scaffold_t* const match_scaffold,
    const uint64_t key_length,
    const uint64_t matching_min_length,
    mm_allocator_t* const mm_allocator) {
  // Init & allocate
  match_scaffold->num_alignment_regions = 0;
  if (match_scaffold->alignment_regions!=NULL) {
    match_scaffold_free_alignment_region(
        match_scaffold,match_scaffold->alignment_regions,mm_allocator);
  }
  const uint64_t max_alignment_regions = DIV_CEIL(key_length,matching_min_length)+1;
  match_scaffold->alignment_regions =
      match_scaffold_allocate_alignment_region(
          match_scaffold,max_alignment_regions,mm_allocator);
}
/*
 * Compose the scaffolding
 */
void match_scaffold_levenshtein_compose_alignment(
    match_scaffold_t* const match_scaffold,
    const match_alignment_t* const match_alignment,
    uint64_t key_offset,
    const uint64_t matching_min_length,
    matches_t* const matches) {
  // Traverse CIGAR elements and compose scaffolding
  const uint64_t cigar_offset = match_alignment->cigar_offset;
  const uint64_t cigar_length = match_alignment->cigar_length;
  cigar_element_t* const cigar_buffer = vector_get_elm(matches->cigar_vector,cigar_offset,cigar_element_t);
  match_alignment_region_t* last_alignment_region = NULL;
  uint64_t text_offset = match_alignment->match_text_offset;
  uint64_t i = 0;
  for (i=0;i<cigar_length;) {
    switch (cigar_buffer[i].type) {
      case cigar_match: {
        // Check matching length
        const uint64_t match_length = cigar_buffer[i].length;
        if (match_length < matching_min_length) {
          text_offset += match_length;
          key_offset += match_length;
          last_alignment_region = NULL; // Nullify last alignment-region
        } else {
          last_alignment_region =
              match_scaffold_compose_add_exact_match(
                  match_scaffold,&key_offset,&text_offset,
                  cigar_offset + i,match_length); // Add Match
        }
        ++i;
        break;
      }
      case cigar_mismatch:
        // Tolerate mismatches if well delimited between exact alignment-regions
        if (last_alignment_region!=NULL && i+1<cigar_length &&
            cigar_buffer[i+1].type==cigar_match && cigar_buffer[i+1].length >= matching_min_length) {
          match_scaffold_compose_add_mismatch(
              match_scaffold,last_alignment_region,
              &key_offset,&text_offset,cigar_buffer[i+1].length+1);
          i += 2;
        } else {
          last_alignment_region = NULL; // Nullify last alignment-region
          ++text_offset; ++key_offset;
          ++i;
        }
        break;
      case cigar_ins:
        text_offset += cigar_buffer[i].length;
        last_alignment_region = NULL; // Nullify last alignment-region
        ++i;
        break;
      case cigar_del:
        key_offset += cigar_buffer[i].length;
        last_alignment_region = NULL; // Nullify last alignment-region
        ++i;
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
}
/*
 * Levenshtein Scaffold Tiled
 */
void match_scaffold_levenshtein_tile_align_bpm(
    match_scaffold_t* const match_scaffold,
    uint8_t* const key,
    const uint64_t key_offset,
    bpm_pattern_t* const bpm_pattern_tile,
    text_trace_t* const text_trace,
    const uint64_t text_begin,
    const uint64_t text_end,
    const uint64_t max_distance,
    const uint64_t matching_min_length,
    matches_t* const matches,
    mm_allocator_t* const mm_allocator) {
  // Parameters
  uint8_t* const text = text_trace->text_padded + text_begin;
  const uint64_t text_length = text_end - text_begin;
  const bool left_gap_alignment = (text_trace->strand==Forward);
  // Fill Matrix (Pv,Mv)
  bpm_align_matrix_t bpm_align_matrix;
  align_bpm_compute_matrix(
      bpm_pattern_tile,text,text_length,
      max_distance,&bpm_align_matrix,mm_allocator);
  PROF_ADD_COUNTER(GP_MATCH_SCAFFOLD_EDIT_CELLS,bpm_pattern_tile->pattern_length*text_length);
  if (bpm_align_matrix.min_score != ALIGN_DISTANCE_INF) {
    // Backtrace and generate CIGAR
    match_alignment_t match_alignment;
    match_alignment.match_position = text_trace->position + text_begin;
    align_bpm_backtrace_matrix(
        bpm_pattern_tile,key+key_offset,text,
        left_gap_alignment,&bpm_align_matrix,
        &match_alignment,matches->cigar_vector);
    // Store the offset (from the beginning of the text)
    // (accounting for the text-padding offset)
    const uint64_t alignment_offset = match_alignment.match_position - text_trace->position;
    const uint64_t text_padding_left = text_trace->text_padded_left;
    gem_fatal_check_msg(alignment_offset < text_padding_left,
        "Scaffold levenshtein. Negative coordinates because of padding");
    // Offset wrt text (without padding)
    match_alignment.match_text_offset = alignment_offset - text_padding_left;
    //    // DEBUG
    //    alignment_print_pretty(stderr,
    //        key+key_offset,100-key_offset,
    //        text_trace->text+match_alignment.match_text_offset,
    //        text_trace->text_length-match_alignment.match_text_offset,matches->cigar_vector,
    //        match_alignment.cigar_offset,match_alignment.cigar_length,mm_allocator);
    // Add the alignment to the scaffold
    match_scaffold_levenshtein_compose_alignment(
        match_scaffold,&match_alignment,
        key_offset,matching_min_length,
        matches);
  }
}
void match_scaffold_levenshtein_tile(
    match_scaffold_t* const match_scaffold,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    alignment_tile_t* const alignment_tile,
    pattern_tile_t* const pattern_tile,
    const uint64_t matching_min_length,
    matches_t* const matches,
    mm_allocator_t* const mm_allocator) {
  // Check tile distance
  const uint64_t match_distance = alignment_tile->distance;
  //  if (match_distance==0) { // TODO
  //    match_scaffold_compose_add_exact_match(...); // Add Match
  //  } else
  if (match_distance!=ALIGN_DISTANCE_INF) {
    // Parameters Key
    uint8_t* const key = pattern->key;
    const uint64_t key_offset = pattern_tile->tile_offset;
    bpm_pattern_t* const bpm_pattern_tile = &pattern_tile->bpm_pattern_tile;
    // Parameters Text
    const uint64_t text_begin = alignment_tile->text_begin_offset;
    const uint64_t text_end = alignment_tile->text_end_offset;
    const uint64_t max_distance = MIN(match_distance,pattern_tile->max_error);
    // BPM Scaffold
    mm_allocator_push_state(mm_allocator); // Push allocator state
    match_scaffold_levenshtein_tile_align_bpm(
        match_scaffold,key,key_offset,bpm_pattern_tile,
        text_trace,text_begin,text_end,max_distance,
        matching_min_length,matches,mm_allocator);
    mm_allocator_pop_state(mm_allocator); // Pop allocator state
    // PROF
    PROF_ADD_COUNTER(GP_MATCH_SCAFFOLD_EDIT_TILES_DISTANCE_BOUND,max_distance);
    PROF_INC_COUNTER(GP_MATCH_SCAFFOLD_EDIT_TILES_ALIGN);
  } else {
    PROF_ADD_COUNTER(GP_MATCH_SCAFFOLD_EDIT_TILES_SKIPPED,1);
  }
}
/*
 * Levenshtein Scaffold
 */
bool match_scaffold_levenshtein(
    match_scaffold_t* const match_scaffold,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    alignment_t* const alignment,
    const uint64_t matching_min_length,
    matches_t* const matches,
    mm_allocator_t* const mm_allocator) {
  PROF_INC_COUNTER(GP_MATCH_SCAFFOLD_EDIT_SCAFFOLDS);
  PROFILE_START(GP_MATCH_SCAFFOLD_EDIT,PROFILE_LEVEL);
  // Allocate scaffold
  match_scaffold_levenshtein_allocate(
      match_scaffold,pattern->key_length,
      matching_min_length,mm_allocator);
  // Compute dimensions
  pattern_tiled_t* const pattern_tiled = &pattern->pattern_tiled;
  const uint64_t num_tiles = pattern_tiled->num_tiles;
  PROF_ADD_COUNTER(GP_MATCH_SCAFFOLD_EDIT_TILES_TOTAL,num_tiles);
  // Scaffold each tile
  uint64_t tile_pos;
  for (tile_pos=0;tile_pos<num_tiles;++tile_pos) {
    alignment_tile_t* const alignment_tile = alignment->alignment_tiles + tile_pos;
    pattern_tile_t* const pattern_tile = pattern->pattern_tiled.tiles + tile_pos;
    match_scaffold_levenshtein_tile(
        match_scaffold,pattern,text_trace,
        alignment_tile,pattern_tile,
        matching_min_length,matches,
        mm_allocator);
  }
  // Chains scaffolds
  match_scaffold_chain(match_scaffold,pattern,text_trace,false,mm_allocator);
  // DEBUG
  gem_cond_debug_block(DEBUG_MATCH_SCAFFOLD_EDIT) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Match.Scaffold.Levenshtein\n");
    tab_global_inc();
    match_scaffold_print(gem_log_get_stream(),matches,match_scaffold);
    tab_global_dec();
  }
  PROFILE_STOP(GP_MATCH_SCAFFOLD_EDIT,PROFILE_LEVEL);
  return match_scaffold->num_alignment_regions > 0;
}

