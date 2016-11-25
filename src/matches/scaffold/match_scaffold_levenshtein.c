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

#include "text/pattern.h"
#include "align/alignment.h"
#include "matches/scaffold/match_scaffold_levenshtein.h"
#include "matches/scaffold/match_scaffold_compose.h"
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
 * Compose the scaffolding
 */
void match_scaffold_levenshtein_compose_alignment(
    match_scaffold_t* const match_scaffold,
    matches_t* const matches,
    const match_alignment_t* const match_alignment,
    uint64_t key_offset,
    const uint64_t text_padding,
    const uint64_t matching_min_length) {
  // Traverse CIGAR elements and compose scaffolding
  const uint64_t cigar_offset = match_alignment->cigar_offset;
  const uint64_t cigar_length = match_alignment->cigar_length;
  cigar_element_t* const cigar_buffer =
      vector_get_elm(matches->cigar_vector,cigar_offset,cigar_element_t);
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
          last_alignment_region = match_scaffold_compose_add_exact_match(
              match_scaffold,&key_offset,&text_offset,cigar_offset + i,match_length); // Add Match
          gem_fatal_check_msg((int64_t)text_offset < (int64_t)text_padding,
              "Scaffold levenshtein. Negative coordinates because of padding");
        }
        ++i;
        break;
      }
      case cigar_mismatch:
        // Tolerate mismatches if well delimited between exact alignment-regions
        if (last_alignment_region!=NULL && i+1<cigar_length &&
            cigar_buffer[i+1].type==cigar_match && cigar_buffer[i+1].length >= matching_min_length) {
          match_scaffold_compose_add_mismatch(match_scaffold,last_alignment_region,
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
 * Levenshtein Scaffold Single Tile
 */
void match_scaffold_levenshtein_tiled(
    match_scaffold_t* const match_scaffold,
    match_align_input_t* const align_input,
    bpm_pattern_t* const bpm_pattern_tile,
    const uint64_t key_offset,
    const uint64_t text_begin,
    const uint64_t text_end,
    const uint64_t distance_bound,
    const uint64_t matching_min_length,
    const bool left_gap_alignment,
    matches_t* const matches,
    mm_stack_t* const mm_stack) {
  // Parameters
  uint8_t* const key = align_input->key + key_offset;
  uint8_t* const text = align_input->text_padded + text_begin;
  const uint64_t text_length = text_end - text_begin;
  const uint64_t max_distance = distance_bound;
  // Fill Matrix (Pv,Mv)
  bpm_align_matrix_t bpm_align_matrix;
  align_bpm_compute_matrix(
      bpm_pattern_tile,text,text_length,
      max_distance,&bpm_align_matrix,mm_stack);
  PROF_ADD_COUNTER(GP_MATCH_SCAFFOLD_EDIT_CELLS,bpm_pattern_tile->pattern_length*text_length);
  if (bpm_align_matrix.min_score != ALIGN_DISTANCE_INF) {
    // Backtrace and generate CIGAR
    const uint64_t match_position = align_input->text_position + text_begin;
    match_alignment_t match_alignment;
    match_alignment.match_position = match_position;
    align_bpm_backtrace_matrix(
        bpm_pattern_tile,key,text,left_gap_alignment,
        &bpm_align_matrix,&match_alignment,matches->cigar_vector);
    // Store the offset (from the beginning of the text)
    // (accounting for the text-padding offset)
    const uint64_t alignment_offset = match_alignment.match_position - match_position;
    const uint64_t text_padding = align_input->text_padding;
    match_alignment.match_text_offset = text_begin + alignment_offset - text_padding;
    //    // DEBUG
    //    match_alignment_print_pretty(stderr,&match_alignment,
    //        matches->cigar_vector,align_input->key + key_offset,bpm_pattern_tile->pattern_length,
    //        align_input->text+match_alignment.match_text_offset,
    //        text_end-match_alignment.match_text_offset,mm_stack);
    // Add the alignment to the scaffold
    match_scaffold_levenshtein_compose_alignment(
        match_scaffold,matches,&match_alignment,
        key_offset,text_padding,matching_min_length);
  }
}
/*
 * Levenshtein Scaffold Tiled
 */
bool match_scaffold_levenshtein(
    match_scaffold_t* const match_scaffold,
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    matches_t* const matches,
    mm_stack_t* const mm_stack) {
  PROF_INC_COUNTER(GP_MATCH_SCAFFOLD_EDIT_SCAFFOLDS);
  PROFILE_START(GP_MATCH_SCAFFOLD_EDIT,PROFILE_LEVEL);
  // Parameters
  const uint64_t key_length =  align_input->key_length;
  const uint64_t matching_min_length = align_parameters->scaffolding_matching_min_length;
  const uint64_t max_alignment_regions = DIV_CEIL(key_length,matching_min_length);
  const bool left_gap_alignment = align_parameters->left_gap_alignment;
  // Init Scaffold
  match_scaffold->num_alignment_regions = 0;
  match_scaffold->alignment_regions =
      mm_stack_calloc(mm_stack,max_alignment_regions,match_alignment_region_t,false);
  // Push stack state
  mm_stack_push_state(mm_stack);
  // Compute dimensions
  alignment_tile_t* const alignment_tiles = align_input->alignment->alignment_tiles;
  alignment_filters_t* const filters = align_input->alignment_filters;
  alignment_filters_tile_t* const filters_tiles = filters->tiles;
  const uint64_t num_tiles = filters->num_tiles;
  PROF_ADD_COUNTER(GP_MATCH_SCAFFOLD_EDIT_TILES_TOTAL,num_tiles);
  // Compute the scaffolding of each tile
  uint64_t tile_pos, key_offset = 0;
  for (tile_pos=0;tile_pos<num_tiles;++tile_pos) {
    // Scaffold tile
    alignment_tile_t* const alignment_tile = alignment_tiles + tile_pos;
    alignment_filters_tile_t* const filters_tile = filters_tiles + tile_pos;
    const uint64_t max_distance = filters_tile->max_error;
    const uint64_t tile_length = filters_tile->tile_length;
    PROF_ADD_COUNTER(GP_MATCH_SCAFFOLD_EDIT_TILES_SKIPPED,(alignment_tile->distance==ALIGN_DISABLED)?1:0);
    const uint64_t match_distance = alignment_tile->distance;
    if (match_distance!=ALIGN_DISABLED && match_distance!=ALIGN_DISTANCE_INF) {
      // TODO if (match_distance==0) { match_scaffold_compose_add_exact_match(...); // Add Match
      const uint64_t text_begin_offset = alignment_tile->text_begin_offset;
      const uint64_t text_end_offset = alignment_tile->text_end_offset;
      const uint64_t distance_bound = MIN(match_distance,max_distance);
      mm_stack_push_state(mm_stack); // Push stack state
      match_scaffold_levenshtein_tiled(
          match_scaffold,align_input,filters_tile->bpm_pattern_tile,
          key_offset,text_begin_offset,text_end_offset,
          distance_bound,matching_min_length,left_gap_alignment,
          matches,mm_stack);
      mm_stack_pop_state(mm_stack); // Pop stack state
      PROF_INC_COUNTER(GP_MATCH_SCAFFOLD_EDIT_TILES_ALIGN);
    }
    // Next
    key_offset += tile_length;
  }
  // Pop stack state
  mm_stack_pop_state(mm_stack);
  // DEBUG
  gem_cond_debug_block(DEBUG_MATCH_SCAFFOLD_EDIT) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Match.Scaffold.Edit.Tiled\n");
    tab_global_inc();
    match_scaffold_print(gem_log_get_stream(),matches,match_scaffold);
    tab_global_dec();
  }
  PROFILE_STOP(GP_MATCH_SCAFFOLD_EDIT,PROFILE_LEVEL);
  return match_scaffold->num_alignment_regions > 0;
}

