/*
 * PROJECT: GEMMapper
 * FILE: filtering_region_align_configure.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering_region_align_configure.h"

/*
 * Configure Basic Alignment
 */
void filtering_region_align_configure_exact(
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters,
    filtering_region_t* const filtering_region,const as_parameters_t* const as_parameters,
    pattern_t* const pattern,const bool emulated_rc_search) {
  // Parameters
  search_parameters_t* const search_parameters = as_parameters->search_parameters;
  const uint64_t key_length = pattern->key_length;
  swg_penalties_t* const swg_penalties = &search_parameters->swg_penalties;
  // Align input
  align_input->key_length         = key_length;
  align_input->text_position      = filtering_region->begin_position;
  align_input->text_trace_offset  = filtering_region->text_trace_offset;
  align_input->text_offset_begin  = filtering_region->align_match_end_column+1 - key_length;
  align_input->text_offset_end    = filtering_region->align_match_end_column+1;
  // Align Parameters
  align_parameters->emulated_rc_search = emulated_rc_search;
  align_parameters->swg_penalties      = swg_penalties;
}
void filtering_region_align_configure_hamming(
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters,
    filtering_region_t* const filtering_region,const as_parameters_t* const as_parameters,
    pattern_t* const pattern,text_trace_t* const text_trace,const bool emulated_rc_search) {
  // Parameters
  search_parameters_t* const search_parameters = as_parameters->search_parameters;
  uint8_t* const key = pattern->key;
  const uint64_t key_length = pattern->key_length;
  bool* const allowed_enc = search_parameters->allowed_enc;
  // Align input
  align_input->key                = key;
  align_input->key_length         = key_length;
  align_input->text_trace_offset  = filtering_region->text_trace_offset;
  align_input->text_position      = filtering_region->begin_position + filtering_region->base_begin_position_offset; // Base position
  align_input->text               = text_trace->text;
  align_input->text_length        = text_trace->text_length;
  align_input->text_offset_begin  = filtering_region->align_match_begin_column;
  align_input->text_offset_end    = filtering_region->align_match_begin_column + key_length;
  // Align Parameters
  align_parameters->emulated_rc_search = emulated_rc_search;
  align_parameters->allowed_enc        = allowed_enc;
}
void filtering_region_align_configure_levenshtein(
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters,
    filtering_region_t* const filtering_region,const as_parameters_t* const as_parameters,
    pattern_t* const pattern,text_trace_t* const text_trace,
    const bool emulated_rc_search,const bool left_gap_alignment) {
  // Parameters
  search_parameters_t* const search_parameters = as_parameters->search_parameters;
  uint8_t* const key = pattern->key;
  const uint64_t key_length = pattern->key_length;
  const uint64_t text_length = text_trace->text_length;
  // Adjust alignment boundaries
  const uint64_t match_effective_range = key_length + 2*filtering_region->align_distance;
  const uint64_t match_begin_column = BOUNDED_SUBTRACTION(
      filtering_region->align_match_end_column,match_effective_range,0);
  const uint64_t match_end_column = BOUNDED_ADDITION(
      filtering_region->align_match_end_column,filtering_region->align_distance,text_length);
  PROF_ADD_COUNTER(GP_ALIGNED_REGIONS_LENGTH,(match_end_column-match_end_column));
  // Align input
  align_input->key                = key;
  align_input->key_length         = key_length;
  align_input->bpm_pattern        = pattern->bpm_pattern;
  align_input->text_trace_offset  = filtering_region->text_trace_offset;
  align_input->text_position      = filtering_region->begin_position;
  align_input->text               = text_trace->text;
  align_input->text_length        = text_trace->text_length;
  align_input->text_offset_begin  = match_begin_column;
  align_input->text_offset_end    = match_end_column;
  // Align Parameters
  align_parameters->emulated_rc_search = emulated_rc_search;
  align_parameters->max_error          = filtering_region->align_distance;
  align_parameters->left_gap_alignment = left_gap_alignment;
  align_parameters->swg_penalties      = &search_parameters->swg_penalties;
}
/*
 * Configure SWG-based Alignment
 */
void filtering_region_align_configure_scaffold(
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters,
    filtering_region_t* const filtering_region,const as_parameters_t* const as_parameters,
    pattern_t* const pattern,text_trace_t* const text_trace,
    const bool left_gap_alignment) {
  // Parameters
  uint8_t* const key = pattern->key;
  const uint64_t key_length = pattern->key_length;
  // Adjust alignment boundaries (to allow optimization)
  const uint64_t max_bandwidth = pattern->max_effective_bandwidth;
  const uint64_t match_effective_range = key_length + 2*filtering_region->align_distance + max_bandwidth;
  const uint64_t text_length = text_trace->text_length;
  const uint64_t text_offset_begin = BOUNDED_SUBTRACTION(filtering_region->align_match_end_column,match_effective_range,0);
  const uint64_t text_offset_end = BOUNDED_ADDITION(filtering_region->align_match_end_column,max_bandwidth,text_length-1);
  // Align input
  align_input->key                      = key;
  align_input->key_length               = key_length;
  align_input->key_trim_left            = filtering_region->key_trim_left;
  align_input->key_trim_right           = filtering_region->key_trim_right;
  align_input->bpm_pattern              = (!filtering_region->key_trimmed) ? pattern->bpm_pattern : filtering_region->bpm_pattern_trimmed;
  align_input->bpm_pattern_tiles        = (!filtering_region->key_trimmed) ? pattern->bpm_pattern_tiles : filtering_region->bpm_pattern_trimmed_tiles;
  align_input->text_position            = filtering_region->begin_position;
  align_input->text                     = text_trace->text;
  align_input->text_length              = text_trace->text_length;
  align_input->text_offset_begin        = text_offset_begin;
  align_input->text_offset_end          = text_offset_end;
  align_input->text_offset_base_begin   = filtering_region->base_begin_position_offset;
  align_input->text_offset_base_end     = filtering_region->base_end_position_offset;
  // Align Parameters
  align_parameters->max_error                           = pattern->max_effective_filtering_error;
  align_parameters->left_gap_alignment                  = left_gap_alignment;
  align_parameters->scaffolding_min_coverage            = as_parameters->alignment_scaffolding_min_coverage_nominal;
  align_parameters->scaffolding_matching_min_length     = as_parameters->alignment_scaffolding_min_matching_length_nominal;
}
void filtering_region_align_configure_swg(
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters,
    filtering_region_t* const filtering_region,const as_parameters_t* const as_parameters,
    pattern_t* const pattern,text_trace_t* const text_trace,
    const bool emulated_rc_search,const bool left_gap_alignment) {
  // Parameters
  search_parameters_t* const search_parameters = as_parameters->search_parameters;
  const uint64_t key_length = pattern->key_length;
  // Adjust alignment boundaries (to allow optimization)
  const uint64_t align_distance_bound = filtering_region->align_distance;
  const uint64_t max_bandwidth = pattern->max_effective_bandwidth;
  const uint64_t match_effective_range = key_length + 2*filtering_region->align_distance + max_bandwidth;
  const uint64_t text_length = text_trace->text_length;
  const uint64_t text_offset_begin = BOUNDED_SUBTRACTION(filtering_region->align_match_end_column,match_effective_range,0);
  const uint64_t text_offset_end = BOUNDED_ADDITION(filtering_region->align_match_end_column,max_bandwidth,text_length-1);
  // Align input
  align_input->key                      = pattern->key;
  align_input->run_length               = pattern->run_length;
  align_input->key_length               = key_length;
  align_input->key_trim_left            = filtering_region->key_trim_left;
  align_input->key_trim_right           = filtering_region->key_trim_right;
  align_input->bpm_pattern              = (!filtering_region->key_trimmed) ? pattern->bpm_pattern : filtering_region->bpm_pattern_trimmed;
  align_input->bpm_pattern_tiles        = (!filtering_region->key_trimmed) ? pattern->bpm_pattern_tiles : filtering_region->bpm_pattern_trimmed_tiles;
  align_input->text_trace_offset        = filtering_region->text_trace_offset;
  align_input->text_position_translated = filtering_region->begin_position_translated;
  align_input->text_position            = filtering_region->begin_position;
  align_input->text                     = text_trace->text;
  align_input->text_length              = text_trace->text_length;
  align_input->text_offset_begin        = text_offset_begin;
  align_input->text_offset_end          = text_offset_end;
  align_input->text_offset_base_begin   = filtering_region->base_begin_position_offset;
  align_input->text_offset_base_end     = filtering_region->base_end_position_offset;
  align_input->align_distance_bound     = align_distance_bound;
  // RL-Input
  align_input->run_length               = pattern->run_length;
  align_input->rl_key_runs              = pattern->rl_runs;
  align_input->rl_text_runs             = text_trace->rl_runs;
  // Align Parameters
  align_parameters->emulated_rc_search                  = emulated_rc_search;
  align_parameters->max_error                           = pattern->max_effective_filtering_error;
  align_parameters->max_bandwidth                       = max_bandwidth;
  align_parameters->left_gap_alignment                  = left_gap_alignment;
  align_parameters->global_min_identity                 = as_parameters->alignment_global_min_identity_nominal;
  align_parameters->global_min_swg_threshold            = as_parameters->alignment_global_min_swg_threshold_nominal;
  align_parameters->local_min_identity                  = as_parameters->alignment_local_min_identity_nominal;
  align_parameters->local_min_swg_threshold             = as_parameters->alignment_local_min_swg_threshold_nominal;
  align_parameters->scaffolding                         = search_parameters->alignment_scaffolding;
  align_parameters->scaffolding_min_coverage            = as_parameters->alignment_scaffolding_min_coverage_nominal;
  align_parameters->scaffolding_matching_min_length     = as_parameters->alignment_scaffolding_min_matching_length_nominal;
  align_parameters->allowed_enc                         = search_parameters->allowed_enc;
  align_parameters->swg_penalties                       = &search_parameters->swg_penalties;
  align_parameters->cigar_curation                      = search_parameters->cigar_curation;
  align_parameters->cigar_curation_min_end_context      = as_parameters->cigar_curation_min_end_context_nominal;
}
