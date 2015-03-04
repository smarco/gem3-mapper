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
#include "matches_align.h"
#include "region_profile.h"

// Approximate Search Internals
typedef enum {
  mapping_adaptive_filtering_fast,     // Fast Adaptive Filtering (First filtering stage)
  mapping_adaptive_filtering_match,    // Adaptive Filtering until a match is found (Never going into full-NS)
  mapping_adaptive_filtering_complete, // Full search (From Adaptive Filtering to Pure Neighborhood Search)
  mapping_neighborhood_search,         // Pure Neighborhood Search (brute-force)
  mapping_filtering_complete           // Pure Filtering Complete
} mapping_mode_t;
typedef enum {
  paired_mapping_map_both_ends,
  paired_mapping_paired_filtering,
  paired_mapping_map_extension
} paired_mapping_mode_t;
typedef enum {
  pair_layout_separate,
  pair_layout_overlap,
  pair_layout_contain,
  pair_layout_dovetail,
  pair_layout_invalid
} pair_layout_t;
typedef enum {
  pair_orientation_concordant,
  pair_orientation_discordant,
  pair_orientation_invalid
} pair_orientation_t;
typedef enum {
  pair_discordant_search_always,
  pair_discordant_search_only_if_no_concordant,
  pair_discordant_search_never
} pair_discordant_search_t;
typedef struct {
  /*
   * Search parameters
   */
  /* Mapping strategy (Mapping mode + properties) */
  mapping_mode_t mapping_mode;
  float filtering_degree;
  /* Qualities */
  quality_model_t quality_model;
  quality_format_t quality_format;
  uint64_t quality_threshold;
  /* Error Model (Regulates the number of Mismatch/Indels) */
  float max_search_error;
  float max_filtering_error;
  float max_filtering_strata_after_best;
  float max_bandwidth;
  float complete_strata_after_best;
  /* Matches search (Regulates the number of matches) */
  uint64_t max_search_matches;
  /* Replacements (Regulates the bases that can be replaced/mismatched) */
  char mismatch_alphabet[DNA_EXT_RANGE];
  uint64_t mismatch_alphabet_length;
  bool allowed_chars[256];
  bool allowed_enc[DNA_EXT_RANGE];
  /* Regions matching */
  float min_matching_length;                   // Minimum length of the read required to match (dangling ends)
  bool allow_region_chaining;                  // Allows chaining regions to compute the alignment
  float region_scaffolding_min_length;         // Minimum length of the matching region (chaining regions)
  float region_scaffolding_coverage_threshold; // Minimum coverage not to resort to scaffolding
  /* Alignment Model/Score */
  alignment_model_t alignment_model;
  swg_penalties_t swg_penalties;
  /*
   * Paired End
   */
  /* Paired-end mode/alg */
  bool paired_end;
  paired_mapping_mode_t paired_mapping_mode;
  pair_discordant_search_t pair_discordant_search;
  uint64_t max_extendable_candidates;
  uint64_t max_matches_per_extension;
  uint64_t min_unique_pair_samples;
  /* Template allowed length */
  uint64_t min_template_length;
  uint64_t max_template_length;
  /* Pair Orientation */
  pair_orientation_t pair_orientation_FR;
  pair_orientation_t pair_orientation_RF;
  pair_orientation_t pair_orientation_FF;
  pair_orientation_t pair_orientation_RR;
  /* Pair allowed lay-outs */
  bool pair_layout_separate;
  bool pair_layout_overlap;
  bool pair_layout_contain;
  bool pair_layout_dovetail;
  /*
   * Internals
   */
  /* Region-Minimal Profile */
  region_profile_model_t rp_minimal;
  /* Region-Delimit Profile */
  region_profile_model_t rp_delimit;
  /* Region-Recovery Profile */
  region_profile_model_t rp_recovery;
  /* Filtering parameters */
  uint64_t filtering_threshold;
  float filtering_region_factor;
  uint64_t pa_filtering_threshold;
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
  uint64_t filtering_degree_nominal;
  /* Error Model (Regulates the number of Mismatch/Indels) */
  uint64_t max_search_error_nominal;                // Maximum number of error/differences while searching (edit distance)
  uint64_t max_filtering_error_nominal;             // Maximum tolerated error at filtering (verifying candidates)
  uint64_t max_filtering_strata_after_best_nominal; // Maximum distance allow from the first matching stratum (at filtering)
  uint64_t max_bandwidth_nominal;                   // Maximum band allowed
  uint64_t complete_strata_after_best_nominal;      // Maximum complete strata from the first matching stratum
  /* Regions Matching */
  uint64_t min_matching_length_nominal;                   // Minimum length of the read required to match (dangling ends)
  uint64_t region_scaffolding_min_length_nominal;         // Minimum length of the matching region (chaining regions)
  uint64_t region_scaffolding_coverage_threshold_nominal; // Minimum coverage not to resort to scaffolding
} search_actual_parameters_t;

/*
 * Approximate Search Parameters
 */
GEM_INLINE void approximate_search_parameters_init(search_parameters_t* const search_parameters);

GEM_INLINE void approximate_search_configure_mapping_strategy(
    search_parameters_t* const search_parameters,
    const mapping_mode_t mapping_mode,const float filtering_degree);
GEM_INLINE void approximate_search_configure_quality_model(
    search_parameters_t* const search_parameters,
    const quality_model_t quality_model,const quality_format_t quality_format,const uint64_t quality_threshold);
GEM_INLINE void approximate_search_configure_error_model(
    search_parameters_t* const search_parameters,float max_search_error,
    float max_filtering_error,float max_filtering_strata_after_best,
    float max_bandwidth,float complete_strata_after_best);
GEM_INLINE void approximate_search_configure_matches(
    search_parameters_t* const search_parameters,const uint64_t max_search_matches);
GEM_INLINE void approximate_search_configure_replacements(
    search_parameters_t* const search_parameters,
    const char* const mismatch_alphabet,const uint64_t mismatch_alphabet_length);
GEM_INLINE void approximate_search_configure_region_handling(
    search_parameters_t* const search_parameters,
    const float min_matching_length,const bool allow_region_chaining,
    const float region_scaffolding_min_length,const float region_scaffolding_coverage_threshold);
GEM_INLINE void approximate_search_configure_alignment_model(
    search_parameters_t* const search_parameters,const alignment_model_t alignment_model);
GEM_INLINE void approximate_search_configure_alignment_match_scores(
    search_parameters_t* const search_parameters,const uint64_t matching_score);
GEM_INLINE void approximate_search_configure_alignment_mismatch_scores(
    search_parameters_t* const search_parameters,const uint64_t mismatch_penalty);

GEM_INLINE void approximate_search_instantiate_values(
    search_actual_parameters_t* const search_actual_parameters,const uint64_t pattern_length);

/*
 * Error Msg
 */
#define GEM_ERROR_ASP_REPLACEMENT_EMPTY "Approximate Search Parameters. Valid replacements set cannot be empty"

#endif /* APPROXIMATE_SEARCH_PARAMETERS_H_ */
