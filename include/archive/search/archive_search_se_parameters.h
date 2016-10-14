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

#ifndef ARCHIVE_SEARCH_PARAMETERS_H_
#define ARCHIVE_SEARCH_PARAMETERS_H_

#include "utils/essentials.h"
#include "text/sequence_qualities_model.h"
#include "align/align_swg_score.h"
#include "archive/search/archive_search_pe_parameters.h"
#include "archive/search/archive_select_parameters.h"
#include "filtering/region_profile/region_profile.h"
#include "neighborhood_search/nsearch_parameters.h"
#include "matches/align/match_align.h"

/*
 * Search Modes
 */
typedef enum {
  mapping_adaptive_filtering_fast,         // Adaptive filtering (exact region filtering)
  mapping_adaptive_filtering_complete,     // Adaptive filtering guaranteeing complete search (forcing min regions)
  mapping_neighborhood_search_brute_force, // Neighborhood Search using brute-force
  mapping_neighborhood_search_partition,   // Neighborhood Search using partitions & bidirectional search
  mapping_hybrid_sensitive,                // Hybrid search stopping until high-quality mapping(s) reached
  mapping_hybrid_complete,                 // Hybrid Complete search (filtering & NS-search)
} mapping_mode_t;
typedef enum {
  local_alignment_never,
  local_alignment_if_unmapped,
} local_alignment_t;

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
 * Clipping
 */
typedef enum {
  clipping_disabled,
  clipping_uncalled,
  clipping_masked,
  clipping_fixed
} clipping_type;

/*
 * Search Parameters
 */
typedef struct {
  /* Mapping strategy (Mapping mode + properties) */
  mapping_mode_t mapping_mode;                         // Mapping mode
  search_paired_parameters_t search_paired_parameters; // Paired-end
  bisulfite_read_t bisulfite_read;                     // Bisulfite mode
  /* Clipping */
  clipping_type clipping;
  uint64_t clip_left;
  uint64_t clip_right;
  /* Qualities */
  sequence_qualities_model_t qualities_model;          // TODO Revision pending
  sequence_qualities_format_t qualities_format;        // TODO Revision pending
  uint64_t quality_threshold;             // TODO Revision pending
  /* Replacements (Regulates the bases that can be replaced/mismatched) */
  char mismatch_alphabet[DNA_EXT_RANGE];  // TODO Revision pending
  uint64_t mismatch_alphabet_length;      // TODO Revision pending
  bool allowed_chars[256];                // TODO Revision pending
  bool allowed_enc[DNA_EXT_RANGE];        // TODO Revision pending
  /* Search error model */
  double complete_search_error;
  double complete_strata_after_best;
  double alignment_max_error;
  double alignment_max_bandwidth;
  double alignment_max_aligned_gap_length;              // Maximum length of gap to be aligned
  double alignment_global_min_identity;                 // Alignment minimum identity to be global
  double alignment_global_min_swg_threshold;            // Alignment minimum SWG score to be global
  region_profile_model_t region_profile_model;          // Region-Profile Model (Lightweight)
  nsearch_parameters_t nsearch_parameters;              // NS parameters
  /* Local Alignment */
  local_alignment_t local_alignment;
  double alignment_local_min_identity;                  // Alignment minimum identity to be local
  double alignment_local_min_swg_threshold;             // Alignment minimum SWG score to be local
  /* Scaffolding */
  bool force_full_swg;                                  // Force full SWG-Alignment
  double alignment_scaffolding_min_coverage;            // Minimum length of the alignment-region (chaining alignment-regions)
  double alignment_scaffolding_min_matching_length;     // Minimum matching chunk to be considered
  bool cigar_curation;
  double cigar_curation_min_end_context;
  /* Alignment Model/Score */
  match_alignment_model_t match_alignment_model;
  swg_penalties_t swg_penalties;
  /* Select Parameters */
  select_parameters_t select_parameters_report;         // Select parameters (Used to report results)
  select_parameters_t select_parameters_align;          // Select parameters (Used to prune work internally)
  /* MAPQ Score */
  mapq_model_t mapq_model_se;
  mapq_model_t mapq_model_pe;
  uint8_t mapq_threshold;
  /* Nominal search parameters (Evaluated to read-length) */
  uint64_t complete_search_error_nominal;
  uint64_t complete_strata_after_best_nominal;
  uint64_t alignment_max_error_nominal;
  uint64_t alignment_max_bandwidth_nominal;
  uint64_t alignment_max_aligned_gap_length_nominal;
  uint64_t alignment_global_min_identity_nominal;
  uint64_t alignment_global_min_swg_threshold_nominal;
  uint64_t alignment_local_min_identity_nominal;
  uint64_t alignment_local_min_swg_threshold_nominal;
  uint64_t alignment_scaffolding_min_coverage_nominal;
  uint64_t alignment_scaffolding_min_matching_length_nominal;
  uint64_t cigar_curation_min_end_context_nominal;
  /* Check */
  archive_check_type check_type;
} search_parameters_t;

/*
 * Search Parameters
 */
void search_parameters_init(search_parameters_t* const search_parameters);

void search_configure_mapping_strategy(
    search_parameters_t* const search_parameters,
    const mapping_mode_t mapping_mode);
void search_configure_quality_model(
    search_parameters_t* const search_parameters,
    const sequence_qualities_model_t qualities_model,
    const sequence_qualities_format_t qualities_format,
    const uint64_t quality_threshold);
void search_configure_replacements(
    search_parameters_t* const search_parameters,
    const char* const mismatch_alphabet,
    const uint64_t mismatch_alphabet_length);
void search_configure_match_alignment_model(
    search_parameters_t* const search_parameters,
    const match_alignment_model_t match_alignment_model);
void search_configure_alignment_match_scores(
    search_parameters_t* const search_parameters,
    const int32_t matching_score);
void search_configure_alignment_mismatch_scores(
    search_parameters_t* const search_parameters,
    const int32_t mismatch_penalty);

void search_instantiate_values(
    search_parameters_t* const search_parameters,
    const uint64_t pattern_length);

/*
 * Display
 */
void search_parameters_print(
    FILE* const stream,
    search_parameters_t* const search_parameters);

/*
 * Error Msg
 */
#define GEM_ERROR_ASP_REPLACEMENT_EMPTY "Search Parameters. Valid replacements set cannot be empty"

#endif /* ARCHIVE_SEARCH_PARAMETERS_H_ */
