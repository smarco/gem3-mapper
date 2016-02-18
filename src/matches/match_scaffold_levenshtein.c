/*
 * PROJECT: GEMMapper
 * FILE: match_scaffold_levenshtein.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "matches/match_scaffold_levenshtein.h"
#include "matches/match_scaffold_compose.h"
#include "align/align.h"
#include "align/align_bpm.h"
#include "data_structures/pattern.h"

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
    uint8_t* const text,
    uint64_t key_offset,
    const uint64_t matching_min_length,
    mm_stack_t* const mm_stack) {
  // Traverse CIGAR elements and compose scaffolding
  const uint64_t cigar_offset = match_alignment->cigar_offset;
  const uint64_t cigar_length = match_alignment->cigar_length;
  cigar_element_t* const cigar_buffer =
      vector_get_elm(matches->cigar_vector,cigar_offset,cigar_element_t);
  region_matching_t* last_scaffold_region = NULL;
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
          last_scaffold_region = NULL; // Nullify last region-matching
        } else {
          last_scaffold_region = match_scaffold_compose_add_matching_exact(
              text,&key_offset,&text_offset,cigar_offset + i,match_length,match_scaffold); // Add Match
        }
        ++i;
        break;
      }
      case cigar_mismatch:
        // Tolerate mismatches if well delimited between exact matching regions
        if (last_scaffold_region!=NULL && i+1<cigar_length &&
            cigar_buffer[i+1].type==cigar_match && cigar_buffer[i+1].length >= matching_min_length) {
          match_scaffold_compose_add_mismatch(&key_offset,&text_offset,
              cigar_buffer[i+1].length+1,match_scaffold,last_scaffold_region);
          i += 2;
        } else {
          last_scaffold_region = NULL; // Nullify last region-matching
          ++text_offset; ++key_offset;
          ++i;
        }
        break;
      case cigar_ins:
        text_offset += cigar_buffer[i].length;
        last_scaffold_region = NULL; // Nullify last region-matching
        ++i;
        break;
      case cigar_del:
        key_offset += cigar_buffer[i].length;
        last_scaffold_region = NULL; // Nullify last region-matching
        ++i;
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }

  }
}
/*
 * Levenshtein Align (for scaffolding)
 *   @align_input->key
 *   @align_input->key_length
 *   @align_input->bpm_pattern
 *   @align_input->text_position
 *   @align_input->text
 *   @align_input->text_length
 *   @align_input->text_offset_begin
 *   @align_input->text_offset_end
 */
void match_scaffold_levenshtein_align(
    match_scaffold_t* const match_scaffold,
    match_align_input_t* const align_input,
    const uint64_t max_error,
    const bool left_gap_alignment,
    matches_t* const matches,
    mm_stack_t* const mm_stack) {
  // Fill Matrix (Pv,Mv)
  mm_stack_push_state(mm_stack); // Save stack state
  match_align_input_t bpm_align_input = {
      .key = align_input->key,
      .key_length = align_input->key_length,
      .bpm_pattern = align_input->bpm_pattern,
      .text = align_input->text + align_input->text_offset_begin,
      .text_length = align_input->text_offset_end - align_input->text_offset_begin,
  };
  bpm_align_matrix_t bpm_align_matrix;
  align_bpm_compute_matrix(&bpm_align_input,max_error,&bpm_align_matrix,mm_stack);
  // Set distance
  match_scaffold->match_alignment.score = bpm_align_matrix.min_score;
  if (bpm_align_matrix.min_score == ALIGN_DISTANCE_INF) {
    match_scaffold->match_alignment.cigar_length = 0;
    mm_stack_pop_state(mm_stack); // Free
    return;
  }
  // Backtrace and generate CIGAR
  match_alignment_t* const match_alignment = &match_scaffold->match_alignment;
  const uint64_t match_position = align_input->text_position + align_input->text_offset_begin;
  match_alignment->match_position = match_position; // Sentinel
  align_bpm_backtrace_matrix(&bpm_align_input,left_gap_alignment,
      &bpm_align_matrix,match_alignment,matches->cigar_vector);
  mm_stack_pop_state(mm_stack); // Free
  // Store the offset
  match_alignment->match_text_offset = align_input->text_offset_begin +
      (match_alignment->match_position - match_position);
}
/*
 * Levenshtein Scaffold Single Tile
 *   @align_input->key
 *   @align_input->text_position
 *   @align_input->text
 */
void match_scaffold_levenshtein_tiled(
    match_scaffold_t* const match_scaffold,
    match_align_input_t* const align_input,
    const uint64_t matching_min_length,
    const bool left_gap_alignment,
    bpm_pattern_t* const bpm_pattern_tile,
    const uint64_t key_offset,
    const uint64_t text_begin,
    const uint64_t text_end,
    matches_t* const matches,
    mm_stack_t* const mm_stack) {
  // Fill Matrix (Pv,Mv)
  bpm_align_matrix_t bpm_align_matrix;
  match_align_input_t bpm_align_input = {
      .key = align_input->key + key_offset,
      .bpm_pattern = bpm_pattern_tile,
      .text = align_input->text + text_begin ,
      .text_length = text_end - text_begin,
  };
  align_bpm_compute_matrix(&bpm_align_input,
      bpm_pattern_tile->pattern_length,&bpm_align_matrix,mm_stack);
  PROF_ADD_COUNTER(GP_MATCH_SCAFFOLD_EDIT_CELLS,
      bpm_pattern_tile->pattern_length*bpm_align_input.text_length);
  if (bpm_align_matrix.min_score != ALIGN_DISTANCE_INF) {
    // Backtrace and generate CIGAR
    const uint64_t match_position = align_input->text_position + text_begin;
    match_alignment_t match_alignment;
    match_alignment.match_position = match_position;
    align_bpm_backtrace_matrix(&bpm_align_input,left_gap_alignment,
        &bpm_align_matrix,&match_alignment,matches->cigar_vector);
    // Store the offset (from the beginning of the text)
    match_alignment.match_text_offset = text_begin +
        (match_alignment.match_position - match_position);
    //    // DEBUG
    //    match_alignment_print_pretty(stderr,&match_alignment,
    //        matches->cigar_vector,align_input->key + key_offset,bpm_pattern_tile->pattern_length,
    //        align_input->text+match_alignment.match_text_offset,
    //        text_end-match_alignment.match_text_offset,mm_stack);
    // Add the alignment to the scaffold
    match_scaffold_levenshtein_compose_alignment(
        match_scaffold,matches,&match_alignment,
        align_input->text,key_offset,matching_min_length,mm_stack);
  }
}
/*
 * Levenshtein Scaffold Tiled
 *   @align_input->key
 *   @align_input->key_length
 *   @align_input->key_trim_left
 *   @align_input->key_trim_right
 *   @align_input->bpm_pattern
 *   @align_input->bpm_pattern_tiles
 *   @align_input->text_position
 *   @align_input->text
 *   @align_input->text_length
 *   @align_input->text_offset_base_begin
 *   @align_input->text_offset_base_end
 *   @align_parameters->max_error
 *   @align_parameters->left_gap_alignment
 *   @align_parameters->min_matching_length
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
  const uint64_t key_trim_left = align_input->key_trim_left;
  const uint64_t key_trim_right = align_input->key_trim_right;
  const uint64_t key_length =  align_input->key_length - key_trim_left - key_trim_right;
  const uint64_t matching_min_length = align_parameters->scaffolding_matching_min_length;
  const uint64_t max_scaffold_regions = DIV_CEIL(key_length,matching_min_length);
  const bool left_gap_alignment = align_parameters->left_gap_alignment;
  // Init Scaffold
  match_scaffold->num_scaffold_regions = 0;
  match_scaffold->scaffold_regions = mm_stack_calloc(
      mm_stack,max_scaffold_regions,region_matching_t,false); // Alloc
  // Save stack state
  mm_stack_push_state(mm_stack);
  // Compute dimensions
  bpm_pattern_t* const bpm_pattern_tiles = align_input->bpm_pattern_tiles;
  const uint64_t num_tiles = bpm_pattern_tiles->num_pattern_tiles;
//  const uint64_t tile_tall = bpm_pattern_tiles->words_per_tile * UINT64_LENGTH;
  const double max_error_rate = (double)align_parameters->max_error/(double)key_length;
  const uint64_t text_length = align_input->text_length;
  const uint64_t text_offset_begin = align_input->text_offset_base_begin;
  const uint64_t text_offset_end = align_input->text_offset_base_end;
  const uint64_t text_effective_length = text_offset_end - text_offset_begin;
  const uint64_t text_extension = (text_effective_length > key_length) ? text_effective_length-key_length : 0;
  // Compute the scaffolding of each tile
  uint64_t key_offset = key_trim_left, text_begin = text_offset_begin;
  uint64_t tile_pos;
  for (tile_pos=0;tile_pos<num_tiles;++tile_pos) {
    mm_stack_push_state(mm_stack); // Save stack state
    bpm_pattern_t* const bpm_pattern_tile = bpm_pattern_tiles + tile_pos;
    // Compute text region
    const uint64_t text_tile_length = bpm_pattern_tile->pattern_length + text_extension;
    const uint64_t boundary_error = max_error_rate * (double)bpm_pattern_tile->pattern_length;
    const uint64_t eff_text_begin = BOUNDED_SUBTRACTION(text_begin,boundary_error,0);
    const uint64_t eff_text_end = BOUNDED_ADDITION(text_begin,text_tile_length+boundary_error,text_length);
    // Scaffold tile
    match_scaffold_levenshtein_tiled(match_scaffold,align_input,
        matching_min_length,left_gap_alignment,bpm_pattern_tile,
        key_offset,eff_text_begin,eff_text_end,matches,mm_stack);
    // Next
    text_begin += bpm_pattern_tile->pattern_length;
    key_offset += bpm_pattern_tile->pattern_length;
    mm_stack_pop_state(mm_stack); // Free
  }
  mm_stack_pop_state(mm_stack); // Free
  // DEBUG
  gem_cond_debug_block(DEBUG_MATCH_SCAFFOLD_EDIT) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Match.Scaffold.Edit.Tiled\n");
    tab_global_inc();
    match_scaffold_print(gem_log_get_stream(),matches,match_scaffold);
    tab_global_dec();
  }
  PROFILE_STOP(GP_MATCH_SCAFFOLD_EDIT,PROFILE_LEVEL);
  return match_scaffold->num_scaffold_regions > 0;
}


////    /*
////     * Debug OND
////     */
////    match_alignment_t ond_match_alignment;
////    ond_match_alignment.match_position = 0;
////    match_align_input_t ond_align_input;
////    ond_align_input.key = align_input->key;
////    ond_align_input.key_length = align_input->key_length;
////    ond_align_input.text = align_input->text; //  + align_input->align_match_begin_column;
////    ond_align_input.text_length = align_input->text_length; // align_input->align_match_end_column - align_input->align_match_begin_column;
////    align_ond_match(&ond_align_input,ond_align_input.text_length,
////        &ond_match_alignment,matches->cigar_vector,mm_stack);
////    // Display
////    match_alignment_print_pretty(stderr,&ond_match_alignment,
////        matches->cigar_vector,ond_align_input.key,ond_align_input.key_length,
////        ond_align_input.text,ond_align_input.text_length,mm_stack);
////    return match_scaffold->num_scaffold_regions > 0;



