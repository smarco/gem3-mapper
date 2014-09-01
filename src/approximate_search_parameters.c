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
// Probing Scheme = (500,1,8,2)
#define PRP_REGION_THRESHOLD 1000
#define PRP_MAX_STEPS 3
#define PRP_DEC_FACTOR 3
#define PRP_REGION_TYPE_THRESHOLD 2
// Loose Scheme = (20,3,3,2)
#define SRP_REGION_THRESHOLD 20
#define SRP_MAX_STEPS 4
#define SRP_DEC_FACTOR 2
#define SRP_REGION_TYPE_THRESHOLD 2
// Tight Scheme = (50,7,3,2)
#define HRP_REGION_THRESHOLD 50
#define HRP_MAX_STEPS 10
#define HRP_DEC_FACTOR 4
#define HRP_REGION_TYPE_THRESHOLD 2
// Recovery Scheme = (20,4,2,2)
#define RRP_REGION_THRESHOLD 200
#define RRP_MAX_STEPS 1
#define RRP_DEC_FACTOR 8
#define RRP_REGION_TYPE_THRESHOLD 2
// Filtering thresholds
#define FILTERING_THRESHOLD 350
#define PA_FILTERING_THRESHOLD 2500
#define FILTERING_REGION_FACTOR ((double)1.0)

GEM_INLINE void approximate_search_initialize_replacements(approximate_search_parameters_t* const search_parameters) {
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
GEM_INLINE void approximate_search_parameters_init(approximate_search_parameters_t* const search_parameters) {
  /*
   * Initialize (DEFAULTS)
   */
  // Mapping strategy
  search_parameters->mapping_mode = mapping_adaptive_filtering;
  search_parameters->fast_mapping_degree = 0;
  // Qualities
  search_parameters->quality_model = quality_model_type_gem;
  search_parameters->quality_format = qualities_ignore;
  search_parameters->quality_threshold = 26;
  // Mismatch/Indels Parameters
  search_parameters->max_search_error = 0.04;
  search_parameters->max_filtering_error = 0.2;
  search_parameters->complete_strata_after_best = 0.0;
  search_parameters->min_matching_length = 0.2;
  // Matches search
  search_parameters->max_search_matches = ALL;
  // Replacements
  approximate_search_initialize_replacements(search_parameters);
  // Soft RP
  search_parameters->srp_region_th = SRP_REGION_THRESHOLD;
  search_parameters->srp_max_steps = SRP_MAX_STEPS;
  search_parameters->srp_dec_factor = SRP_DEC_FACTOR;
  search_parameters->srp_region_type_th = SRP_REGION_TYPE_THRESHOLD;
  // Hard RP
  search_parameters->hrp_region_th = HRP_REGION_THRESHOLD;
  search_parameters->hrp_max_steps = HRP_MAX_STEPS;
  search_parameters->hrp_dec_factor = HRP_DEC_FACTOR;
  search_parameters->hrp_region_type_th = HRP_REGION_TYPE_THRESHOLD;
  // Recover Read
  search_parameters->rrp_region_th = RRP_REGION_THRESHOLD;
  search_parameters->rrp_max_steps = RRP_MAX_STEPS;
  search_parameters->rrp_dec_factor = RRP_DEC_FACTOR;
  search_parameters->rrp_region_type_th = RRP_REGION_TYPE_THRESHOLD;
  // Filtering Thresholds
  search_parameters->filtering_region_factor = FILTERING_REGION_FACTOR;
  search_parameters->filtering_threshold = FILTERING_THRESHOLD;
  search_parameters->pa_filtering_threshold = PA_FILTERING_THRESHOLD;
  // Check alignments
  search_parameters->check_matches = check_none;
}
GEM_INLINE void approximate_search_configure_mapping_strategy(
    approximate_search_parameters_t* const search_parameters,
    const mapping_mode_t mapping_mode,const float mapping_degree) {
  search_parameters->mapping_mode = mapping_mode;
  search_parameters->fast_mapping_degree = mapping_degree;
}
GEM_INLINE void approximate_search_configure_quality_model(
    approximate_search_parameters_t* const search_parameters,
    const quality_model_t quality_model,const quality_format_t quality_format,const uint64_t quality_threshold) {
  search_parameters->quality_model = quality_model;
  search_parameters->quality_format = quality_format;
  search_parameters->quality_threshold = quality_threshold;
}
GEM_INLINE void approximate_search_configure_error_model(
    approximate_search_parameters_t* const search_parameters,
    float max_search_error,float max_filtering_error,
    float complete_strata_after_best,float min_matching_length) {
  search_parameters->max_search_error = max_search_error;
  search_parameters->max_filtering_error = max_filtering_error;
  search_parameters->complete_strata_after_best = complete_strata_after_best;
  search_parameters->min_matching_length = min_matching_length;
}
GEM_INLINE void approximate_search_configure_replacements(
    approximate_search_parameters_t* const search_parameters,
    char* const mismatch_alphabet,const uint64_t mismatch_alphabet_length) {
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
GEM_INLINE void approximate_search_configure_matches(
    approximate_search_parameters_t* const search_parameters,const uint64_t max_search_matches) {
  search_parameters->max_search_matches = max_search_matches;
}
GEM_INLINE void approximate_search_instantiate_values(
    approximate_search_parameters_t* const search_parameters,const uint64_t pattern_length) {
  search_parameters->fast_mapping_degree_nominal = integer_proportion(search_parameters->fast_mapping_degree,pattern_length);
  search_parameters->max_search_error_nominal = integer_proportion(search_parameters->max_search_error,pattern_length);
  search_parameters->max_filtering_error_nominal = integer_proportion(search_parameters->max_filtering_error,pattern_length);
  search_parameters->complete_strata_after_best_nominal = integer_proportion(search_parameters->complete_strata_after_best,pattern_length);
  search_parameters->min_matching_length_nominal = integer_proportion(search_parameters->min_matching_length,pattern_length);
}
