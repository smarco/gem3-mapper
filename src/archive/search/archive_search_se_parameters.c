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
 *   Archive-Search Single-End Parameters (initializers & handlers)
 */

#include "archive/search/archive_search_se_parameters.h"
#include "stats/report_stats.h"

/*
 * Macro Utils
 */
#define SEARCH_INSTANTIATE_VALUE(parameters,field_name,pattern_length) \
  parameters->field_name##_nominal = integer_proportion(parameters->field_name,pattern_length)

/*
 * Init
 */
void search_parameters_init_alignment_swg_penalties(search_parameters_t* const search_parameters) {
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_A][ENC_DNA_CHAR_A] = +1;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_A][ENC_DNA_CHAR_C] = -4;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_A][ENC_DNA_CHAR_G] = -4;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_A][ENC_DNA_CHAR_T] = -4;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_A][ENC_DNA_CHAR_N] = -4;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_C][ENC_DNA_CHAR_A] = -4;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_C][ENC_DNA_CHAR_C] = +1;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_C][ENC_DNA_CHAR_G] = -4;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_C][ENC_DNA_CHAR_T] = -4;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_C][ENC_DNA_CHAR_N] = -4;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_G][ENC_DNA_CHAR_A] = -4;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_G][ENC_DNA_CHAR_C] = -4;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_G][ENC_DNA_CHAR_G] = +1;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_G][ENC_DNA_CHAR_T] = -4;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_G][ENC_DNA_CHAR_N] = -4;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_T][ENC_DNA_CHAR_A] = -4;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_T][ENC_DNA_CHAR_C] = -4;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_T][ENC_DNA_CHAR_G] = -4;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_T][ENC_DNA_CHAR_T] = +1;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_T][ENC_DNA_CHAR_N] = -4;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_N][ENC_DNA_CHAR_A] = -4;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_N][ENC_DNA_CHAR_C] = -4;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_N][ENC_DNA_CHAR_G] = -4;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_N][ENC_DNA_CHAR_T] = -4;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_N][ENC_DNA_CHAR_N] = -4;
  search_parameters->swg_penalties.generic_match_score = 1;
  search_parameters->swg_penalties.generic_mismatch_score = -4;
  search_parameters->swg_penalties.gap_open_score = -6;
  search_parameters->swg_penalties.gap_extension_score = -1;
}
void search_parameters_init_alignment(search_parameters_t* const search_parameters) {
  /* Global Alignment */
  search_parameters->complete_search_error = 0.04;
  search_parameters->complete_strata_after_best = 1.0;
  search_parameters->alignment_max_error = 0.12;
  search_parameters->alignment_max_extension_error = 0.20;
  search_parameters->alignment_max_bandwidth = 0.20;
  search_parameters->alignment_force_full_swg = false;
  search_parameters->alignment_global_min_identity = 0.80;
  search_parameters->alignment_global_min_swg_threshold = 0.40; // 0.40*read_length*match_score
  /* Local Alignment */
  search_parameters->alignment_local = local_alignment_if_unmapped;
  search_parameters->alignment_local_min_identity = 40.0;
  search_parameters->alignment_local_min_swg_threshold = 20.0;
  search_parameters->swg_score_dropoff = 100;
  search_parameters->alignment_scaffolding_min_coverage = 0.80;
  search_parameters->alignment_scaffolding_min_matching_length = 10.0;
  search_parameters->cigar_curation = true;
  search_parameters->cigar_curation_min_end_context = 0.0;
  /* Candidate verification */
  search_parameters->candidate_verification.verification_strategy = verification_chained;
  search_parameters->candidate_verification.candidate_local_drop_off = false;
  search_parameters->candidate_verification.kmer_tiles = 2;
  search_parameters->candidate_verification.kmer_length = 6;
  // Alignment Model
  search_parameters->match_alignment_model = match_alignment_model_gap_affine;
  search_parameters_init_alignment_swg_penalties(search_parameters);
}
void select_parameters_init_mapq(search_parameters_t* const search_parameters) {
  search_parameters->mapq_model_se = mapq_model_gem;
  search_parameters->mapq_model_pe = mapq_model_gem;
}
void select_parameters_init_gpu(search_parameters_t* const search_parameters) {
  search_parameters->gpu_stage_region_profile_enabled = false;
  search_parameters->gpu_stage_decode_enabled = false;
  search_parameters->gpu_stage_kmer_filter_enabled = false;
  search_parameters->gpu_stage_bpm_distance_enabled = false;
  search_parameters->gpu_stage_bpm_align_enabled = false;
}
void search_parameters_init(search_parameters_t* const search_parameters) {
  // Mapping strategy
  search_parameters->mapping_mode = mapping_adaptive_filtering_fast;
  // Bisulfite
  search_parameters->bisulfite_read = bisulfite_inferred_C2T_G2A;
	search_parameters->control_sequences[0] = SEQUENCING_CONTROL;
	search_parameters->control_sequences[1] = UNDERCONVERSION_CONTROL;
	search_parameters->control_sequences[2] = OVERCONVERSION_CONTROL;
  search_parameters->rrbs = false;
  // Restriction restriction_sites
  search_parameters->restriction_sites = NULL;
  // Clipping
  search_parameters->clipping = clipping_uncalled;
  search_parameters->clip_left = 0;
  search_parameters->clip_right = 0;
  // Qualities
  search_parameters->qualities_format = sequence_qualities_offset_33;
  // Approximate Search
  region_profile_model_init(&search_parameters->region_profile_model);
  nsearch_parameters_init(&search_parameters->nsearch_parameters);
  // Alignment Parameters
  search_parameters_init_alignment(search_parameters);
  // Select Parameters
  select_parameters_init(&search_parameters->select_parameters);
  // MAPQ Score
  select_parameters_init_mapq(search_parameters);
  // Paired End
  search_paired_parameters_init(&search_parameters->search_paired_parameters);
  // GPU
  select_parameters_init_gpu(search_parameters);
  // Check
  search_parameters->check_type = archive_check_nothing;
}
/*
 * Accessors
 */
void search_configure_mapping_strategy(
    search_parameters_t* const search_parameters,
    const mapping_mode_t mapping_mode) {
  search_parameters->mapping_mode = mapping_mode;
}
void search_configure_quality_model(
    search_parameters_t* const search_parameters,
    const sequence_qualities_format_t qualities_format,
    const uint64_t quality_threshold) {
  search_parameters->qualities_format = qualities_format;
}
void search_configure_alignment_model(
    search_parameters_t* const search_parameters,
    const match_alignment_model_t match_alignment_model) {
  search_parameters->match_alignment_model = match_alignment_model;
}
void search_configure_alignment_match_scores(
    search_parameters_t* const search_parameters,
    const int32_t matching_score) {
  // Match
  search_parameters->swg_penalties.generic_match_score = matching_score;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_A][ENC_DNA_CHAR_A] = matching_score;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_C][ENC_DNA_CHAR_C] = matching_score;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_G][ENC_DNA_CHAR_G] = matching_score;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_T][ENC_DNA_CHAR_T] = matching_score;
}
void search_configure_alignment_mismatch_scores(
    search_parameters_t* const search_parameters,
    const int32_t mismatch_penalty) {
  // Mismatch
  const int64_t mismatch_score = -((int64_t)mismatch_penalty);
  search_parameters->swg_penalties.generic_mismatch_score = mismatch_penalty;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_A][ENC_DNA_CHAR_C] = mismatch_score;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_A][ENC_DNA_CHAR_G] = mismatch_score;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_A][ENC_DNA_CHAR_T] = mismatch_score;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_A][ENC_DNA_CHAR_N] = mismatch_score;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_C][ENC_DNA_CHAR_A] = mismatch_score;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_C][ENC_DNA_CHAR_G] = mismatch_score;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_C][ENC_DNA_CHAR_T] = mismatch_score;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_C][ENC_DNA_CHAR_N] = mismatch_score;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_G][ENC_DNA_CHAR_A] = mismatch_score;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_G][ENC_DNA_CHAR_C] = mismatch_score;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_G][ENC_DNA_CHAR_T] = mismatch_score;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_G][ENC_DNA_CHAR_N] = mismatch_score;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_T][ENC_DNA_CHAR_A] = mismatch_score;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_T][ENC_DNA_CHAR_C] = mismatch_score;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_T][ENC_DNA_CHAR_G] = mismatch_score;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_T][ENC_DNA_CHAR_N] = mismatch_score;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_N][ENC_DNA_CHAR_A] = mismatch_score;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_N][ENC_DNA_CHAR_C] = mismatch_score;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_N][ENC_DNA_CHAR_G] = mismatch_score;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_N][ENC_DNA_CHAR_T] = mismatch_score;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_N][ENC_DNA_CHAR_N] = mismatch_score;
}
void search_configure_alignment_gap_scores(
    search_parameters_t* const search_parameters,
    const uint64_t gap_open_penalty,
    const uint64_t gap_extension_penalty) {
  // Gaps
  search_parameters->swg_penalties.gap_open_score = -((int32_t)gap_open_penalty);
  search_parameters->swg_penalties.gap_extension_score = -((int32_t)gap_extension_penalty);
}
void search_instantiate_values(
    search_parameters_t* const search_parameters,
    const uint64_t pattern_length) {
  // Instantiate Search Parameters Values (Nominal search parameters; evaluated to read-length)
  const uint64_t max_swg_score = pattern_length * search_parameters->swg_penalties.generic_match_score;
  SEARCH_INSTANTIATE_VALUE(search_parameters,complete_search_error,pattern_length);
  SEARCH_INSTANTIATE_VALUE(search_parameters,complete_strata_after_best,pattern_length);
  SEARCH_INSTANTIATE_VALUE(search_parameters,alignment_max_error,pattern_length);
  SEARCH_INSTANTIATE_VALUE(search_parameters,alignment_max_extension_error,pattern_length);
  SEARCH_INSTANTIATE_VALUE(search_parameters,alignment_max_bandwidth,pattern_length);
  SEARCH_INSTANTIATE_VALUE(search_parameters,alignment_global_min_identity,pattern_length);
  SEARCH_INSTANTIATE_VALUE(search_parameters,alignment_global_min_swg_threshold,max_swg_score);
  SEARCH_INSTANTIATE_VALUE(search_parameters,alignment_local_min_identity,pattern_length);
  SEARCH_INSTANTIATE_VALUE(search_parameters,alignment_local_min_swg_threshold,max_swg_score);
  SEARCH_INSTANTIATE_VALUE(search_parameters,alignment_scaffolding_min_coverage,pattern_length);
  SEARCH_INSTANTIATE_VALUE(search_parameters,alignment_scaffolding_min_matching_length,pattern_length);
  SEARCH_INSTANTIATE_VALUE(search_parameters,cigar_curation_min_end_context,pattern_length);
  // Instantiate Select Parameters Values
  select_parameters_instantiate_values(&search_parameters->select_parameters,pattern_length);
}
/*
 * Display
 */
void search_parameters_print(
    FILE* const stream,
    search_parameters_t* const search_parameters) {
  tab_fprintf(stream,"[GEM]>AS.Parameters\n");
  tab_fprintf(stream,"=> Search.Parameters\n");
  tab_global_inc();
  tab_fprintf(stream,"=> Nominal.Parameters\n");
  tab_fprintf(stream,"  => Complete.search.error %lu\n",search_parameters->complete_search_error_nominal);
  tab_fprintf(stream,"  => Complete.strata.after.best %lu\n",search_parameters->complete_strata_after_best_nominal);
  tab_fprintf(stream,"  => Alingment.max.error %lu\n",search_parameters->alignment_max_error_nominal);
  tab_fprintf(stream,"  => Alingment.max.bandwidth %lu\n",search_parameters->alignment_max_bandwidth_nominal);
  tab_fprintf(stream,"  => Alignment.Global.min.identity %lu\n",search_parameters->alignment_global_min_identity_nominal);
  tab_fprintf(stream,"  => Alignment.Global.min.swg.threshold %d\n",search_parameters->alignment_global_min_swg_threshold_nominal);
  tab_fprintf(stream,"  => Alignment.Local.min.identity %lu\n",search_parameters->alignment_local_min_identity_nominal);
  tab_fprintf(stream,"  => Alignment.Local.min.swg.threshold %d\n",search_parameters->alignment_local_min_swg_threshold_nominal);
  tab_fprintf(stream,"  => Alignment.scaffolding.min.coverage %lu\n",search_parameters->alignment_scaffolding_min_coverage_nominal);
  tab_fprintf(stream,"  => Alignment.scaffolding.min.matching.length %lu\n",search_parameters->alignment_scaffolding_min_matching_length_nominal);
  tab_fprintf(stream,"  => Cigar.curation.min.end.context %lu\n",search_parameters->cigar_curation_min_end_context_nominal);
  tab_global_dec();
}
