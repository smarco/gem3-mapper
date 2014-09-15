/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_parameters.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef APPROXIMATE_SEARCH_PARAMETERS_H_
#define APPROXIMATE_SEARCH_PARAMETERS_H_

#include "essentials.h"
#include "quality_model.h"
#include "dna_text.h"

// Approximate Search Internals
typedef enum { check_none, check_pmatches_correctness, check_imatches_correctness, check_correctness, check_completness } check_matches_t;
typedef enum {
  mapping_incremental_mapping,
  mapping_adaptive_filtering,
  mapping_fixed_filtering,
  mapping_fast,
  mapping_neighborhood_search
} mapping_mode_t;
typedef struct {
  /*
   * Search parameters
   */
  /* Mapping strategy (Mapping mode + properties) */
  mapping_mode_t mapping_mode;
  float mapping_degree;
  /* Qualities */
  quality_model_t quality_model;
  quality_format_t quality_format;
  uint64_t quality_threshold;
  /* Error Model (Regulates the number of Mismatch/Indels) */
  float max_search_error;
  float max_filtering_error;
  float complete_strata_after_best;
  float min_matching_length;
  /* Matches search (Regulates the number of matches) */
  uint64_t max_search_matches;
  /* Replacements (Regulates the bases that can be replaced/mismatched) */
  char mismatch_alphabet[DNA_EXT_RANGE];
  uint64_t mismatch_alphabet_length;
  bool allowed_chars[256];
  bool allowed_enc[DNA_EXT_RANGE];
  /* Alignment Score */
  /*
   * Internals
   */
  /* Soft Region Profile Parameters */
  uint64_t srp_region_th; // Max. number of candidates allowed per region
  uint64_t srp_max_steps; // Max. number of characters to explore to improve the region
  uint64_t srp_dec_factor; // Decreasing factor per step in region exploration
  uint64_t srp_region_type_th; // Threshold to classify regions {ZERO,NON_ZERO}
  /* Hard Region Profile Parameters */
  uint64_t hrp_region_th;
  uint64_t hrp_max_steps;
  uint64_t hrp_dec_factor;
  uint64_t hrp_region_type_th;
  /* Read recovery parameters */
  uint64_t rrp_region_th;
  uint64_t rrp_max_steps;
  uint64_t rrp_dec_factor;
  uint64_t rrp_region_type_th;
  /* Filtering parameters */
  uint64_t filtering_threshold;
  float filtering_region_factor;
  uint64_t pa_filtering_threshold;
  /* Checkers */
  check_matches_t check_matches;
} search_parameters_t;
typedef struct {
  /*
   * Search parameters (Functional)
   */
  search_parameters_t* search_parameters;
  /*
   * Actual Search parameters (Evaluated to read-length)
   */
  /* Mapping strategy (Mapping mode + properties) */
  uint64_t fast_mapping_degree_nominal;
  /* Error Model (Regulates the number of Mismatch/Indels) */
  uint64_t max_search_error_nominal;           // Maximum number of error/differences while searching (edit distance)
  uint64_t max_filtering_error_nominal;        // Maximum tolerated error at filtering (verifying candidates)
  uint64_t complete_strata_after_best_nominal; // Maximum complete strata from first matching stratum
  uint64_t min_matching_length_nominal;        // Minimum mapping segment size (verifying candidates)
} search_actual_parameters_t;

/*
 * Approximate Search Parameters
 */
GEM_INLINE void approximate_search_parameters_init(search_parameters_t* const search_parameters);

GEM_INLINE void approximate_search_configure_mapping_strategy(
    search_parameters_t* const search_parameters,
    const mapping_mode_t mapping_mode,const float mapping_degree);
GEM_INLINE void approximate_search_configure_quality_model(
    search_parameters_t* const search_parameters,
    const quality_model_t quality_model,const quality_format_t quality_format,const uint64_t quality_threshold);
GEM_INLINE void approximate_search_configure_error_model(
    search_parameters_t* const search_parameters,
    float max_search_error,float max_filtering_error,
    float complete_strata_after_best,float min_matching_length);
GEM_INLINE void approximate_search_configure_replacements(
    search_parameters_t* const search_parameters,
    const char* const mismatch_alphabet,const uint64_t mismatch_alphabet_length);
GEM_INLINE void approximate_search_configure_matches(
    search_parameters_t* const search_parameters,const uint64_t max_search_matches);

GEM_INLINE void approximate_search_instantiate_values(
    search_actual_parameters_t* const search_actual_parameters,const uint64_t pattern_length);

/*
 * Error Msg
 */
#define GEM_ERROR_ASP_REPLACEMENT_EMPTY "Approximate Search Parameters. Valid replacements set cannot be empty"

#endif /* APPROXIMATE_SEARCH_PARAMETERS_H_ */

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//typedef struct {
//  ch_t *forward_key;
//  ch_t *reverse_key;
//  uint64_t key_len;
//  ch_t *mismatch_mask;
//  direction_t search_direction;
//  /* Mismatch constraints */
//  uint64_t max_distance;
//  uint64_t min_ins_size;
//  uint64_t max_ins_size;
//  uint64_t base_key_len;
//  uint64_t max_extended_matches;
//  uint64_t min_anchor_size;
//  bool check_duplicates;
//  uint64_t unique_pairing;
//  /* Replacements */
//  slch_t repls[CH_T_RANGE];
//  uint64_t repls_len;
//  bool allowed_chars[CH_T_RANGE];
//  bool allowed_repl[CH_T_RANGE];
//  uint64_t num_wildcards;
//  /* Bit vector for Myers-DP */
//  uint64_t *forward_peq;
//  uint64_t *reverse_peq;
//  /* Key tracking id */
//  uint64_t key_id;
//  /* Internals */
//  uint64_t check_alignments;
//} fm_asm_pe_extend_parameters;
