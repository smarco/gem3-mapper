/*
 * PROJECT: GEMMapper
 * FILE: filtering_region_align_configure.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering/filtering_region_align_configure.h"

/*
 * Configure Basic Alignment
 */
void filtering_region_align_configure_exact(
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    filtering_region_t* const filtering_region,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    const bool emulated_rc_search) {
  // Parameters
  const uint64_t key_length = pattern->key_length;
  swg_penalties_t* const swg_penalties = &search_parameters->swg_penalties;
  // Align input
  align_input->key_length         = key_length;
  align_input->text_position      = filtering_region->text_begin_position;
  align_input->text_trace_offset  = filtering_region->text_trace_offset;
  align_input->region_alignment   = &filtering_region->region_alignment;
  // Align Parameters
  align_parameters->emulated_rc_search = emulated_rc_search;
  align_parameters->swg_penalties      = swg_penalties;
}
void filtering_region_align_configure_hamming(
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    filtering_region_t* const filtering_region,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    const bool emulated_rc_search) {
  // Parameters
  uint8_t* const key = pattern->key;
  const uint64_t key_length = pattern->key_length;
  bool* const allowed_enc = search_parameters->allowed_enc;
  // Align input
  align_input->key                = key;
  align_input->key_length         = key_length;
  align_input->text_trace_offset  = filtering_region->text_trace_offset;
  const uint64_t text_offset = filtering_region->text_source_region_offset - filtering_region->key_source_region_offset;
  align_input->text_position      = filtering_region->text_begin_position + text_offset; // Base position
  align_input->text               = text_trace->text;
  align_input->text_length        = text_trace->text_length;
  align_input->region_alignment   = &filtering_region->region_alignment;
  // Align Parameters
  align_parameters->emulated_rc_search = emulated_rc_search;
  align_parameters->allowed_enc        = allowed_enc;
}
void filtering_region_align_configure_levenshtein(
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    filtering_region_t* const filtering_region,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    const bool emulated_rc_search,
    const bool left_gap_alignment,
    mm_stack_t* const mm_stack) {
  // Parameters
  uint8_t* const key = pattern->key;
  const uint64_t key_length = pattern->key_length;
  // Adjust alignment boundaries
  region_alignment_t* const region_alignment = &filtering_region->region_alignment;
  const uint64_t align_distance = pattern->max_effective_filtering_error;
  // Align input
  filtering_region_bpm_pattern_select(filtering_region,pattern,
      &align_input->bpm_pattern,&align_input->bpm_pattern_tiles,mm_stack);
  align_input->key                = key;
  align_input->key_length         = key_length;
  align_input->text_trace_offset  = filtering_region->text_trace_offset;
  align_input->text_position      = filtering_region->text_begin_position;
  align_input->text               = text_trace->text;
  align_input->text_length        = text_trace->text_length;
  align_input->region_alignment   = region_alignment;
  // Align Parameters
  align_parameters->emulated_rc_search = emulated_rc_search;
  align_parameters->max_error          = align_distance;
  align_parameters->left_gap_alignment = left_gap_alignment;
  align_parameters->swg_penalties      = &search_parameters->swg_penalties;
}
/*
 * Configure SWG-based Alignment
 */
void filtering_region_align_configure_swg(
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    filtering_region_t* const filtering_region,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    const bool emulated_rc_search,
    const bool left_gap_alignment,
    const bool local_alignment,
    mm_stack_t* const mm_stack) {
  // Parameters
  const uint64_t key_length = pattern->key_length;
  // Align input
  align_input->key                      = pattern->key;
  align_input->key_length               = key_length;
  align_input->key_trim_left            = filtering_region->key_trim_left;
  align_input->key_trim_right           = filtering_region->key_trim_right;
  align_input->text_trace_offset        = filtering_region->text_trace_offset;
  align_input->text_position            = filtering_region->text_begin_position;
  align_input->text                     = text_trace->text;
  align_input->text_length              = text_trace->text_length;
  align_input->region_alignment         = &filtering_region->region_alignment;
  filtering_region_bpm_pattern_select(filtering_region,pattern,
      &align_input->bpm_pattern,&align_input->bpm_pattern_tiles,mm_stack);
  // RL-Input
  align_input->run_length               = pattern->run_length;
  align_input->rl_key_runs_acc          = pattern->rl_runs_acc;
  align_input->rl_text_runs_acc         = text_trace->rl_runs_acc;
  // Align Parameters
  align_parameters->emulated_rc_search              = emulated_rc_search;
  align_parameters->max_error                       = filtering_region->max_error;
  align_parameters->max_bandwidth                   = filtering_region->max_bandwidth;
  align_parameters->left_gap_alignment              = left_gap_alignment;
  align_parameters->global_min_identity             = search_parameters->alignment_global_min_identity_nominal;
  align_parameters->global_min_swg_threshold        = search_parameters->alignment_global_min_swg_threshold_nominal;
  align_parameters->local_alignment                 = local_alignment;
  align_parameters->local_min_identity              = search_parameters->alignment_local_min_identity_nominal;
  align_parameters->local_min_swg_threshold         = search_parameters->alignment_local_min_swg_threshold_nominal;
  align_parameters->max_aligned_gap_length          = search_parameters->alignment_max_aligned_gap_length_nominal;
  align_parameters->force_full_swg                  = search_parameters->force_full_swg;
  align_parameters->scaffolding_min_coverage        = search_parameters->alignment_scaffolding_min_coverage_nominal;
  align_parameters->scaffolding_matching_min_length = search_parameters->alignment_scaffolding_min_matching_length_nominal;
  align_parameters->allowed_enc                     = search_parameters->allowed_enc;
  align_parameters->swg_penalties                   = &search_parameters->swg_penalties;
  align_parameters->cigar_curation                  = search_parameters->cigar_curation;
  align_parameters->cigar_curation_min_end_context  = search_parameters->cigar_curation_min_end_context_nominal;
}
