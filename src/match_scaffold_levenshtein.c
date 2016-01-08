/*
 * PROJECT: GEMMapper
 * FILE: match_scaffold_levenshtein.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "match_scaffold_levenshtein.h"
#include "match_scaffold_compose.h"
#include "align.h"
#include "align_bpm.h"

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
 *   @align_input->key
 *   @align_input->text
 *   @align_parameters->min_matching_length
 *   @align_parameters->min_context_length
 */
void match_scaffold_levenshtein_compose_matching_regions(
    matches_t* const matches,match_align_input_t* const align_input,
    uint64_t key_pos,uint64_t text_pos,const uint64_t cigar_offset,
    const uint64_t cigar_length,match_align_parameters_t* const align_parameters,
    match_scaffold_t* const match_scaffold,mm_stack_t* const mm_stack) {
  // Parameters
  uint8_t* const text = align_input->text;
  const uint64_t matching_min_length = align_parameters->scaffolding_matching_min_length;
  cigar_element_t* const cigar_buffer = vector_get_elm(matches->cigar_vector,cigar_offset,cigar_element_t);
  // Init scaffolding
  match_scaffold->scaffold_regions = mm_stack_calloc(mm_stack,cigar_length,region_matching_t,false);
  match_scaffold->num_scaffold_regions = 0;
  match_scaffold->scaffolding_coverage = 0;
  // Traverse CIGAR elements and compose scaffolding
  region_matching_t* last_scaffold_region = NULL;
  uint64_t i = 0;
  for (i=0;i<cigar_length;) {
    switch (cigar_buffer[i].type) {
      case cigar_match: {
        // Check matching length
        const uint64_t match_length = cigar_buffer[i].length;
        if (match_length < matching_min_length) {
          text_pos += match_length;
          key_pos += match_length;
          last_scaffold_region = NULL; // Nullify last region-matching
        } else {
          last_scaffold_region = match_scaffold_compose_add_matching_exact(
              text,&key_pos,&text_pos,cigar_offset + i,match_length,match_scaffold); // Add Match
        }
        ++i;
        break;
      }
      case cigar_mismatch:
        // Tolerate mismatches if well delimited between exact matching regions
        if (last_scaffold_region!=NULL && i+1<cigar_length &&
            cigar_buffer[i+1].type==cigar_match && cigar_buffer[i+1].length >= matching_min_length) {
          match_scaffold_compose_add_mismatch(&key_pos,&text_pos,
              cigar_buffer[i+1].length+1,match_scaffold,last_scaffold_region);
          i += 2;
        } else {
          last_scaffold_region = NULL; // Nullify last region-matching
          ++text_pos; ++key_pos;
          ++i;
        }
        break;
      case cigar_ins:
        text_pos += cigar_buffer[i].length;
        last_scaffold_region = NULL; // Nullify last region-matching
        ++i;
        break;
      case cigar_del:
        key_pos += cigar_buffer[i].length;
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
 *   @align_parameters->left_gap_alignment
 *   @align_parameters->max_error
 */
void match_scaffold_levenshtein_align(
    matches_t* const matches,match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    match_scaffold_t* const match_scaffold,mm_stack_t* const mm_stack) {
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
  align_bpm_compute_matrix(&bpm_align_input,align_parameters->max_error,&bpm_align_matrix,mm_stack);
  // Set distance
  match_scaffold->match_alignment.score = bpm_align_matrix.min_score;
  if (bpm_align_matrix.min_score == ALIGN_DISTANCE_INF) {
    match_scaffold->match_alignment.cigar_length = 0;
    mm_stack_pop_state(mm_stack,false); // Free
    return;
  }
  // Backtrace and generate CIGAR
  match_alignment_t* const match_alignment = &match_scaffold->match_alignment;
  const uint64_t match_position = align_input->text_position + align_input->text_offset_begin;
  match_alignment->match_position = match_position; // Sentinel
  align_bpm_backtrace_matrix(&bpm_align_input,align_parameters->left_gap_alignment,
      &bpm_align_matrix,match_alignment,matches->cigar_vector);
  mm_stack_pop_state(mm_stack,false); // Free
  // Store the offset
  match_alignment->match_position =
      (match_alignment->match_position - match_position) + align_input->text_offset_begin;
}
/*
 * Levenshtein Scaffold (based on levenshtein alignment)
 *   @align_input->key
 *   @align_input->key_length
 *   @align_input->bpm_pattern
 *   @align_input->text_position
 *   @align_input->text
 *   @align_input->text_length
 *   @align_input->text_offset_begin
 *   @align_input->text_offset_end
 *   @align_parameters->left_gap_alignment
 *   @align_parameters->min_matching_length
 *   @align_parameters->min_context_length
 */
bool match_scaffold_levenshtein(
    matches_t* const matches,match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    match_scaffold_t* const match_scaffold,mm_stack_t* const mm_stack) {
  PROF_INC_COUNTER(GP_MATCH_SCAFFOLD_EDIT_SCAFFOLDS);
  PROFILE_START(GP_MATCH_SCAFFOLD_EDIT,PROFILE_LEVEL);
  // Init

//  gem_timer_t timer; TIMER_RESET(&timer); TIMER_START(&timer);

  match_scaffold->scaffold_type = scaffold_levenshtein;
  // Compute the levenshtein alignment (CIGAR)
  match_scaffold_levenshtein_align(matches,align_input,align_parameters,match_scaffold,mm_stack);

//  TIMER_STOP(&timer);
//  fprintf(stdout,"> %lu %lu %lu %lu",align_input->key_length,
//      align_input->text_length,align_input->text_offset_end-align_input->text_offset_begin,
//      TIMER_GET_TOTAL_MS(&timer));

  match_alignment_t* const match_alignment = &match_scaffold->match_alignment;
  if (match_alignment->score==ALIGN_DISTANCE_INF) {
    match_scaffold->num_scaffold_regions = 0;
    match_scaffold->scaffolding_coverage = 0;
    PROFILE_STOP(GP_MATCH_SCAFFOLD_EDIT,PROFILE_LEVEL);
    return false;
  }
  // Compose matching region of the scaffolding
  match_scaffold_levenshtein_compose_matching_regions(matches,
      align_input,0,match_alignment->match_position,match_alignment->cigar_offset,
      match_alignment->cigar_length,align_parameters,match_scaffold,mm_stack);
  PROFILE_STOP(GP_MATCH_SCAFFOLD_EDIT,PROFILE_LEVEL);
  PROF_ADD_COUNTER(GP_MATCH_SCAFFOLD_EDIT_COVERAGE,(100*match_scaffold->scaffolding_coverage)/align_input->key_length);
  // DEBUG
  gem_cond_debug_block(DEBUG_MATCH_SCAFFOLD_EDIT) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Match.Scaffold.Edit\n");
    tab_global_inc();
    match_scaffold_print(gem_log_get_stream(),matches,match_scaffold);
    tab_global_dec();
//    /*
//     * Debug
//     */
//    match_alignment_t ond_match_alignment;
//    ond_match_alignment.match_position = 0;
//    match_align_input_t ond_align_input;
//    ond_align_input.key = align_input->key;
//    ond_align_input.key_length = align_input->key_length;
//    ond_align_input.text = align_input->text; //  + align_input->align_match_begin_column;
//    ond_align_input.text_length = align_input->text_length; // align_input->align_match_end_column - align_input->align_match_begin_column;
//    align_ond_match(&ond_align_input,ond_align_input.text_length,
//        &ond_match_alignment,matches->cigar_vector,mm_stack);
//    // Display
//    match_alignment_print_pretty(stderr,&ond_match_alignment,
//        matches->cigar_vector,ond_align_input.key,ond_align_input.key_length,
//        ond_align_input.text,ond_align_input.text_length,mm_stack);
//    return match_scaffold->num_scaffold_regions > 0;
  }
  return match_scaffold->num_scaffold_regions > 0;
}
