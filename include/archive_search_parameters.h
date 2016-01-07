/*
 * PROJECT: GEMMapper
 * FILE: archive_search_parameters.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef ARCHIVE_SEARCH_PARAMETERS_H_
#define ARCHIVE_SEARCH_PARAMETERS_H_

#include "essentials.h"
#include "quality_model.h"
#include "region_profile.h"
#include "align_swg_score.h"
#include "archive_search_paired_parameters.h"
#include "archive_select_parameters.h"
#include "match_align.h"

/*
 * Search Modes
 */
typedef enum {
  mapping_adaptive_filtering_fast,     // Adaptive filtering until high confidence results reached
  mapping_adaptive_filtering_thorough, // Adaptive filtering thoroughly exploring filtering candidates
  mapping_adaptive_filtering_complete, // Full complete search (using neighborhood search if needed)
  mapping_fixed_filtering_complete,    // Pure Filtering Complete
  mapping_neighborhood_search,         // Pure Neighborhood Search (brute-force)
  mapping_test
} mapping_mode_t;
typedef enum {
  unbounded_alignment_never,
  unbounded_alignment_if_unmapped,
} unbounded_alignment_t;

/*
 * Bisulfite Reads-Mode
 */
typedef enum {
  bisulfite_read_inferred,
  bisulfite_read_1,
  bisulfite_read_2,
  bisulfite_read_interleaved
} bisulfite_read_t;

/*
 * Check Type
 */
typedef enum {
  archive_check_nothing,
  archive_check_correct,
  archive_check_correct__first_optimum,
  archive_check_correct__all_optimum,
  archive_check_correct__complete
} archive_check_type;

/*
 * Score Model & Sorting
 */
typedef enum {
  mapq_model_none,            // None
  mapq_model_gem,             // GEM Score (Case stratification + Logistic Regression)
  mapq_model_classify,        // GEM Classification (Towards score calibration to adjust thresholds)
  mapq_model_dump_predictors, // GEM Print predictors
} mapq_model_t;
typedef enum {
  matches_sorting_distance,
  matches_sorting_mapq
} matches_sorting_t;

/*
 * Search Parameters
 */
typedef struct {
  /* Mapping strategy (Mapping mode + properties) */
  mapping_mode_t mapping_mode;
  /* Qualities */
  quality_model_t quality_model;   // TODO Revision pending
  quality_format_t quality_format; // TODO Revision pending
  uint64_t quality_threshold;      // TODO Revision pending
  /* Replacements (Regulates the bases that can be replaced/mismatched) */
  char mismatch_alphabet[DNA_EXT_RANGE];  // TODO Revision pending
  uint64_t mismatch_alphabet_length;      // TODO Revision pending
  bool allowed_chars[256];                // TODO Revision pending
  bool allowed_enc[DNA_EXT_RANGE];        // TODO Revision pending
  /* Search error model */
  float complete_search_error;
  float complete_strata_after_best;
  float alignment_max_error;
  unbounded_alignment_t unbounded_alignment;
  float max_bandwidth;
  float alignment_min_identity;                        // Minimum identity of the alignment (applies to dangling ends too)
  bool alignment_scaffolding;
  float alignment_scaffolding_min_coverage;            // Minimum length of the matching region (chaining regions)
  float alignment_scaffolding_homopolymer_min_context; // Minimum matching bases to support a homopolymer error
  float alignment_scaffolding_min_matching_length;     // Minimum matching chunk to be considered
  bool cigar_curation;
  float cigar_curation_min_end_context;
  /* Alignment Model/Score */
  alignment_model_t alignment_model;
  swg_penalties_t swg_penalties;
  double swg_threshold;
  /* Bisulfite mode */
  bisulfite_read_t bisulfite_read;
  /* Paired-end */
  search_paired_parameters_t search_paired_parameters;
  /* Filtering parameters */
  uint64_t filtering_threshold;
  float filtering_region_factor;
  region_profile_model_t rp_lightweight;  // Region-Lightweight
  region_profile_model_t rp_heavyweight;  // Region-Heavyweight
  region_profile_model_t rp_boost;        // Region-Boost
  region_profile_model_t rp_delimit;      // Region-Delimit
  /* Select Parameters */
  select_parameters_t select_parameters_report;  // Select parameters (Used to report results)
  select_parameters_t select_parameters_align;   // Select parameters (Used to prune work internally)
  /* Check */
  archive_check_type check_type;
  /* MAPQ Score */
  mapq_model_t mapq_model;
  uint8_t mapq_threshold;
} search_parameters_t;

/*
 * Actual Search Parameters (Instantiated for the actual read-length)
 */
typedef struct {
  /* Search parameters (Generic) [Shared between several searches going in parallel] */
  search_parameters_t* search_parameters;
  /* Nominal search parameters (Evaluated to read-length) */
  uint64_t complete_search_error_nominal;
  uint64_t complete_strata_after_best_nominal;
  uint64_t alignment_max_error_nominal;
  uint64_t max_bandwidth_nominal;
  uint64_t alignment_min_identity_nominal;
  uint64_t alignment_scaffolding_min_coverage_nominal;
  uint64_t alignment_scaffolding_homopolymer_min_context_nominal;
  uint64_t alignment_scaffolding_min_matching_length_nominal;
  uint64_t cigar_curation_min_end_context_nominal;
  int64_t swg_threshold_nominal;
} as_parameters_t;

/*
 * Search Parameters
 */
void search_parameters_init(search_parameters_t* const search_parameters);

void search_configure_mapping_strategy(
    search_parameters_t* const search_parameters,const mapping_mode_t mapping_mode);
void search_configure_quality_model(
    search_parameters_t* const search_parameters,const quality_model_t quality_model,
    const quality_format_t quality_format,const uint64_t quality_threshold);
void search_configure_replacements(
    search_parameters_t* const search_parameters,
    const char* const mismatch_alphabet,const uint64_t mismatch_alphabet_length);
void search_configure_alignment_model(
    search_parameters_t* const search_parameters,const alignment_model_t alignment_model);
void search_configure_alignment_match_scores(
    search_parameters_t* const search_parameters,const uint64_t matching_score);
void search_configure_alignment_mismatch_scores(
    search_parameters_t* const search_parameters,const uint64_t mismatch_penalty);

void search_instantiate_values(
    as_parameters_t* const as_parameters,const uint64_t pattern_length);

/*
 * Display
 */
void search_parameters_print(FILE* const stream,search_parameters_t* const search_parameters);
void as_parameters_print(FILE* const stream,as_parameters_t* const parameters);

/*
 * Error Msg
 */
#define GEM_ERROR_ASP_REPLACEMENT_EMPTY "Search Parameters. Valid replacements set cannot be empty"

#endif /* ARCHIVE_SEARCH_PARAMETERS_H_ */
