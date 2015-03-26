/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_parameters.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "approximate_search_parameters.h"

/*
 * Region profile default parameters
 */
GEM_INLINE void approximate_search_initialize_replacements(search_parameters_t* const search_parameters) {
  // Reset
  memset(search_parameters->allowed_chars,0,256*sizeof(bool));
  memset(search_parameters->allowed_enc,0,DNA_EXT_RANGE*sizeof(bool));
  search_parameters->mismatch_alphabet[0] = DNA_CHAR_A;
  search_parameters->mismatch_alphabet[1] = DNA_CHAR_C;
  search_parameters->mismatch_alphabet[2] = DNA_CHAR_G;
  search_parameters->mismatch_alphabet[3] = DNA_CHAR_T;
  search_parameters->mismatch_alphabet_length = 4;
  search_parameters->allowed_chars[DNA_CHAR_A] = true;
  search_parameters->allowed_chars[DNA_CHAR_C] = true;
  search_parameters->allowed_chars[DNA_CHAR_G] = true;
  search_parameters->allowed_chars[DNA_CHAR_T] = true;
  search_parameters->allowed_chars[DNA_CHAR_N] = false;
  search_parameters->allowed_enc[ENC_DNA_CHAR_A] = true;
  search_parameters->allowed_enc[ENC_DNA_CHAR_C] = true;
  search_parameters->allowed_enc[ENC_DNA_CHAR_G] = true;
  search_parameters->allowed_enc[ENC_DNA_CHAR_T] = true;
  search_parameters->allowed_enc[ENC_DNA_CHAR_N] = false;
}
GEM_INLINE void approximate_search_parameters_init(search_parameters_t* const search_parameters) {
  /*
   * Single End
   */
  // Mapping strategy
  search_parameters->mapping_mode = mapping_adaptive_filtering_match;
  search_parameters->filtering_degree = 0;
  // Qualities
  search_parameters->quality_model = quality_model_type_gem;
  search_parameters->quality_format = qualities_ignore;
  search_parameters->quality_threshold = 26;
  // Mismatch/Indels Parameters
  search_parameters->max_search_error = 0.04;
  search_parameters->max_filtering_error = 0.20;
  search_parameters->max_filtering_strata_after_best = 0.04;
  search_parameters->max_bandwidth = 0.20;
  search_parameters->complete_strata_after_best = 0.0;
  // Matches search
  search_parameters->max_search_matches = ALL;
  // Replacements
  approximate_search_initialize_replacements(search_parameters);
  // Regions handling
  search_parameters->allow_region_chaining = true;
  search_parameters->min_matching_length = 0.20;
  search_parameters->region_scaffolding_min_length = 10;
  search_parameters->region_scaffolding_coverage_threshold = 0.80;
  // Alignment Model/Score
  search_parameters->alignment_model = alignment_model_gap_affine;
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
	// Bisulfite mode
  search_parameters->bisulfite_mode = false;
  search_parameters->bisulfite_read = bisulfite_read_inferred;
  /*
   * Paired End
   */
  /* Paired-end mode/alg */
  search_parameters->paired_end = false;
  search_parameters->paired_mapping_mode = paired_mapping_map_extension;
  search_parameters->pair_discordant_search = pair_discordant_search_only_if_no_concordant;
  search_parameters->max_extendable_candidates = 20;
  search_parameters->max_matches_per_extension = 2;
  search_parameters->min_unique_pair_samples = 1000;
  /* Template allowed length */
  search_parameters->min_template_length = 0;    // Automatically adjusted in runtime
  search_parameters->max_template_length = 2000; // Automatically adjusted in runtime
  /* Concordant Orientation */
  search_parameters->pair_discordant_search = pair_discordant_search_never;
  search_parameters->pair_orientation_FR = pair_orientation_concordant;
  search_parameters->pair_orientation_RF = pair_orientation_invalid;
  search_parameters->pair_orientation_FF = pair_orientation_invalid;
  search_parameters->pair_orientation_RR = pair_orientation_invalid;
  /* Pair allowed lay-outs */
  search_parameters->pair_layout_separate = true;
  search_parameters->pair_layout_overlap = true;
  search_parameters->pair_layout_contain = true;
  search_parameters->pair_layout_dovetail = true;
  /*
   * Internals
   */
  // Region-Minimal Scheme = (20,4,2,2)
  search_parameters->rp_minimal.region_th = 20;
  search_parameters->rp_minimal.max_steps = 4;
  search_parameters->rp_minimal.dec_factor = 2;
  search_parameters->rp_minimal.region_type_th = 2;
  // Region-Delimit Scheme = (50,10,4,2)
  search_parameters->rp_delimit.region_th = 100;
  search_parameters->rp_delimit.max_steps = 4;
  search_parameters->rp_delimit.dec_factor = 2;
  search_parameters->rp_delimit.region_type_th = 2;
  // Region-Recovery Scheme = (200,1,8,2)
  search_parameters->rp_recovery.region_th = 200;
  search_parameters->rp_recovery.max_steps = 1;
  search_parameters->rp_recovery.dec_factor = 8;
  search_parameters->rp_recovery.region_type_th = 2;
  // Filtering Thresholds
  search_parameters->filtering_region_factor = 1.0;
  search_parameters->filtering_threshold = 350;
  search_parameters->pa_filtering_threshold = 2500;
}
GEM_INLINE void approximate_search_configure_mapping_strategy(
    search_parameters_t* const search_parameters,
    const mapping_mode_t mapping_mode,const float filtering_degree) {
  search_parameters->mapping_mode = mapping_mode;
  search_parameters->filtering_degree = filtering_degree;
}
GEM_INLINE void approximate_search_configure_quality_model(
    search_parameters_t* const search_parameters,
    const quality_model_t quality_model,const quality_format_t quality_format,const uint64_t quality_threshold) {
  search_parameters->quality_model = quality_model;
  search_parameters->quality_format = quality_format;
  search_parameters->quality_threshold = quality_threshold;
}
GEM_INLINE void approximate_search_configure_error_model(
    search_parameters_t* const search_parameters,float max_search_error,
    float max_filtering_error,float max_filtering_strata_after_best,
    float max_bandwidth,float complete_strata_after_best) {
  search_parameters->max_search_error = max_search_error;
  search_parameters->max_filtering_error = max_filtering_error;
  search_parameters->max_filtering_strata_after_best = max_filtering_strata_after_best;
  search_parameters->max_bandwidth = max_bandwidth;
  search_parameters->complete_strata_after_best = complete_strata_after_best;
}
GEM_INLINE void approximate_search_configure_matches(
    search_parameters_t* const search_parameters,const uint64_t max_search_matches) {
  search_parameters->max_search_matches = max_search_matches;
}
GEM_INLINE void approximate_search_configure_replacements(
    search_parameters_t* const search_parameters,
    const char* const mismatch_alphabet,const uint64_t mismatch_alphabet_length) {
  // Reset
  approximate_search_initialize_replacements(search_parameters);
  // Filter replacements
  uint64_t i, count;
  for (i=0,count=0;i<mismatch_alphabet_length;i++) {
    if (is_dna(mismatch_alphabet[i])) {
      const char c = dna_normalized(mismatch_alphabet[i]);
      search_parameters->mismatch_alphabet[count] = c;
      search_parameters->allowed_chars[(uint8_t)c] = true;
      search_parameters->allowed_enc[dna_encode(c)] = true;
      ++count;
    }
  }
  gem_cond_fatal_error(count==0,ASP_REPLACEMENT_EMPTY);
  search_parameters->mismatch_alphabet_length = count;
}
GEM_INLINE void approximate_search_configure_region_handling(
    search_parameters_t* const search_parameters,
    const float min_matching_length,const bool allow_region_chaining,
    const float region_scaffolding_min_length,const float region_scaffolding_coverage_threshold) {
  search_parameters->allow_region_chaining = allow_region_chaining;
  search_parameters->min_matching_length = min_matching_length;
  search_parameters->region_scaffolding_min_length = region_scaffolding_min_length;
  search_parameters->region_scaffolding_coverage_threshold = region_scaffolding_coverage_threshold;
}
GEM_INLINE void approximate_search_configure_alignment_model(
    search_parameters_t* const search_parameters,const alignment_model_t alignment_model) {
  search_parameters->alignment_model = alignment_model;
}
GEM_INLINE void approximate_search_configure_alignment_match_scores(
    search_parameters_t* const search_parameters,const uint64_t matching_score) {
  // Match
  search_parameters->swg_penalties.generic_match_score = matching_score;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_A][ENC_DNA_CHAR_A] = matching_score;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_C][ENC_DNA_CHAR_C] = matching_score;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_G][ENC_DNA_CHAR_G] = matching_score;
  search_parameters->swg_penalties.matching_score[ENC_DNA_CHAR_T][ENC_DNA_CHAR_T] = matching_score;
}
GEM_INLINE void approximate_search_configure_alignment_mismatch_scores(
    search_parameters_t* const search_parameters,const uint64_t mismatch_penalty) {
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
GEM_INLINE void approximate_search_configure_alignment_gap_scores(
    search_parameters_t* const search_parameters,
    const uint64_t gap_open_penalty,const uint64_t gap_extension_penalty) {
  // Gaps
  search_parameters->swg_penalties.gap_open_score = -((int32_t)gap_open_penalty);
  search_parameters->swg_penalties.gap_extension_score = -((int32_t)gap_extension_penalty);
}
GEM_INLINE void approximate_search_instantiate_values(
    search_actual_parameters_t* const search_actual_parameters,const uint64_t pattern_length) {
  const search_parameters_t* const search_parameters = search_actual_parameters->search_parameters;
  search_actual_parameters->filtering_degree_nominal = integer_proportion(search_parameters->filtering_degree,pattern_length);
  search_actual_parameters->max_search_error_nominal = integer_proportion(search_parameters->max_search_error,pattern_length);
  search_actual_parameters->max_filtering_error_nominal = integer_proportion(search_parameters->max_filtering_error,pattern_length);
  search_actual_parameters->max_filtering_strata_after_best_nominal = integer_proportion(search_parameters->max_filtering_strata_after_best,pattern_length);
  search_actual_parameters->max_bandwidth_nominal = integer_proportion(search_parameters->max_bandwidth,pattern_length);
  search_actual_parameters->complete_strata_after_best_nominal = integer_proportion(search_parameters->complete_strata_after_best,pattern_length);
  search_actual_parameters->min_matching_length_nominal = integer_proportion(search_parameters->min_matching_length,pattern_length);
  search_actual_parameters->region_scaffolding_min_length_nominal = integer_proportion(search_parameters->region_scaffolding_min_length,pattern_length);
  search_actual_parameters->region_scaffolding_coverage_threshold_nominal = integer_proportion(search_parameters->region_scaffolding_coverage_threshold,pattern_length);
}
