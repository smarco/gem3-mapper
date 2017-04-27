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
 *   Filtering region module provides functions to configure all the parameters
 *   involved in a full-alignment process of a filtering-region against a text-region
 */

#include "filtering/region/filtering_region_align_configure.h"

/*
 * Configure Basic Alignment
 */
void filtering_region_align_configure_exact(
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    filtering_region_t* const filtering_region,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern) {
  // Parameters
  const uint64_t key_length = pattern->key_length;
  swg_penalties_t* const swg_penalties = &search_parameters->swg_penalties;
  // Align input
  align_input->key_length          = key_length;
  align_input->sequence_clip_left  = pattern->clip_left;
  align_input->sequence_clip_right = pattern->clip_right;
  align_input->text_position       = filtering_region->text_begin_position;
  align_input->text_trace          = &filtering_region->text_trace;
  align_input->alignment           = &filtering_region->alignment;
  // Align Parameters
  align_parameters->swg_penalties  = swg_penalties;
}
void filtering_region_align_configure_hamming(
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    filtering_region_t* const filtering_region,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern) {
  // Parameters
  uint8_t* const key = pattern->key;
  const uint64_t key_length = pattern->key_length;
  swg_penalties_t* const swg_penalties = &search_parameters->swg_penalties;
  // Align input
  align_input->key                 = key;
  align_input->key_length          = key_length;
  align_input->sequence_clip_left  = pattern->clip_left;
  align_input->sequence_clip_right = pattern->clip_right;
  align_input->pattern_tiled       = &pattern->pattern_tiled;
  align_input->text_trace          = &filtering_region->text_trace;
  align_input->text_position       = filtering_region->text_begin_position; // Base position
  align_input->text                = filtering_region->text_trace.text;
  align_input->text_length         = key_length;
  align_input->alignment           = &filtering_region->alignment;
  // Align Parameters
  align_parameters->swg_penalties  = swg_penalties;
}
void filtering_region_align_configure_levenshtein(
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    filtering_region_t* const filtering_region,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    mm_allocator_t* const mm_allocator) {
  // Parameters
  uint8_t* const key = pattern->key;
  const uint64_t key_length = pattern->key_length;
  // Adjust alignment boundaries
  alignment_t* const alignment = &filtering_region->alignment;
  const uint64_t align_distance = pattern->max_effective_filtering_error;
  // Align input
  align_input->sequence_clip_left  = pattern->clip_left;
  align_input->sequence_clip_right = pattern->clip_right;
  align_input->key                 = key;
  align_input->key_length          = key_length;
  align_input->pattern_tiled       = &pattern->pattern_tiled;
  align_input->text_trace          = &filtering_region->text_trace;
  align_input->text_position       = filtering_region->text_begin_position;
  align_input->text                = filtering_region->text_trace.text;
  align_input->text_length         = filtering_region->text_trace.text_length;
  align_input->text_padded         = filtering_region->text_trace.text_padded;
  align_input->text_padding        = filtering_region->key_trim_left;
  align_input->alignment           = alignment;
  // Align Parameters
  align_parameters->max_error          = align_distance;
  align_parameters->left_gap_alignment = (filtering_region->text_trace.strand==Forward);
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
    const bool local_alignment,
    mm_allocator_t* const mm_allocator) {
  // Parameters
  const uint64_t key_length = pattern->key_length;
  // Align input
  align_input->sequence_clip_left       = pattern->clip_left;
  align_input->sequence_clip_right      = pattern->clip_right;
  align_input->key                      = pattern->key;
  align_input->key_length               = key_length;
  align_input->key_trim_left            = filtering_region->key_trim_left;
  align_input->key_trim_right           = filtering_region->key_trim_right;
  align_input->pattern_tiled            = &pattern->pattern_tiled;
  align_input->text_trace               = &filtering_region->text_trace;
  align_input->text_position            = filtering_region->text_begin_position;
  align_input->text                     = filtering_region->text_trace.text;
  align_input->text_length              = filtering_region->text_trace.text_length;
  align_input->text_padded              = filtering_region->text_trace.text_padded;
  align_input->text_padding             = filtering_region->key_trim_left;
  align_input->alignment                = &filtering_region->alignment;
  // RL-Input
  align_input->run_length               = pattern->run_length;
  align_input->rl_key_runs_acc          = pattern->rl_runs_acc;
  align_input->rl_text_runs_acc         = filtering_region->text_trace.rl_runs_acc;
  // Align Parameters
  align_parameters->max_error                       = pattern->max_effective_filtering_error;
  align_parameters->max_bandwidth                   = pattern->max_effective_bandwidth;
  align_parameters->left_gap_alignment              = (filtering_region->text_trace.strand==Forward);
  align_parameters->global_min_identity             = search_parameters->alignment_global_min_identity_nominal;
  align_parameters->global_min_swg_threshold        = search_parameters->alignment_global_min_swg_threshold_nominal;
  align_parameters->local_alignment                 = local_alignment;
  align_parameters->local_min_identity              = search_parameters->alignment_local_min_identity_nominal;
  align_parameters->local_min_swg_threshold         = search_parameters->alignment_local_min_swg_threshold_nominal;
  align_parameters->max_aligned_gap_length          = search_parameters->alignment_max_aligned_gap_length_nominal;
  align_parameters->alignment_force_full_swg        = search_parameters->alignment_force_full_swg;
  align_parameters->scaffolding_min_coverage        = search_parameters->alignment_scaffolding_min_coverage_nominal;
  align_parameters->scaffolding_matching_min_length = search_parameters->alignment_scaffolding_min_matching_length_nominal;
  align_parameters->swg_penalties                   = &search_parameters->swg_penalties;
  align_parameters->cigar_curation                  = search_parameters->cigar_curation;
  align_parameters->cigar_curation_min_end_context  = search_parameters->cigar_curation_min_end_context_nominal;
}
