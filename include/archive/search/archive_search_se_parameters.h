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

/*
 * Search Modes
 */
typedef enum {
  mapping_adaptive_filtering_fast,         // Adaptive filtering (exact region filtering)
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
 * Candidate verification strategies
 */
typedef enum {
  verification_BPM = 0,
  verification_chained = 1,
} verification_strategy_t;
typedef struct {
  // Verification Strategy
  verification_strategy_t verification_strategy;
  bool candidate_local_drop_off;
  // Kmer filter
  uint64_t kmer_tiles;
  uint64_t kmer_length;
} candidate_verification_t;

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
  bisulfite_read_t bisulfite_read;                     // Bisulfite mode
  bool rrbs;
	char* control_sequences[3];
  vector_t *restriction_sites;

  /* Clipping */
  clipping_type clipping;
  uint64_t clip_left;
  uint64_t clip_right;
  /* Qualities */
  sequence_qualities_format_t qualities_format;
  /* Approximate Search */
  region_profile_model_t region_profile_model;          // Region-Profile Model (Lightweight)
  candidate_verification_t candidate_verification;      // Candidate verification strategy
  nsearch_parameters_t nsearch_parameters;              // NS parameters
  bool gpu_stage_region_profile_enabled;
  bool gpu_stage_decode_enabled;
  bool gpu_stage_kmer_filter_enabled;
  bool gpu_stage_bpm_distance_enabled;
  bool gpu_stage_bpm_align_enabled;
  /* Global alignment */
  double complete_search_error;
  double complete_strata_after_best;
  double alignment_max_error;
  double alignment_max_extension_error;
  double alignment_max_bandwidth;
  bool   alignment_force_full_swg;                      // Force full SWG-Alignment
  double alignment_global_min_identity;                 // Alignment minimum identity to be global
  double alignment_global_min_swg_threshold;            // Alignment minimum SWG score to be global
  match_alignment_model_t match_alignment_model;        // Alignment Model/Score
  swg_penalties_t swg_penalties;
  /* Local Alignment */
  local_alignment_t alignment_local;
  double alignment_local_min_identity;                  // Alignment minimum identity to be local
  double alignment_local_min_swg_threshold;             // Alignment minimum SWG score to be local
  int32_t swg_score_dropoff;                            // SWG-Alignment local dropoff score
  /* Scaffolding */
  double alignment_scaffolding_min_coverage;            // Minimum length of the alignment-region (chaining alignment-regions)
  double alignment_scaffolding_min_matching_length;     // Minimum matching chunk to be considered
  bool cigar_curation;
  double cigar_curation_min_end_context;
  /* Select Parameters */
  select_parameters_t select_parameters;                // Select parameters (Used to report results)
  /* MAPQ Score */
  mapq_model_t mapq_model_se;
  mapq_model_t mapq_model_pe;
  /* Paired-end Search */
  search_paired_parameters_t search_paired_parameters;  // Paired-end
  /* Nominal search parameters (Evaluated to read-length) */
  uint64_t complete_search_error_nominal;
  uint64_t complete_strata_after_best_nominal;
  uint64_t alignment_max_error_nominal;
  uint64_t alignment_max_extension_error_nominal;
  uint64_t alignment_max_bandwidth_nominal;
  uint64_t alignment_global_min_identity_nominal;
  int32_t alignment_global_min_swg_threshold_nominal;
  uint64_t alignment_local_min_identity_nominal;
  int32_t alignment_local_min_swg_threshold_nominal;
  uint64_t alignment_scaffolding_min_coverage_nominal;
  uint64_t alignment_scaffolding_min_matching_length_nominal;
  uint64_t cigar_curation_min_end_context_nominal;
  /* Check */
  archive_check_type check_type;
} search_parameters_t;

/*
 * Init
 */
void search_parameters_init(search_parameters_t* const search_parameters);

/*
 * Accessors
 */
void search_configure_mapping_strategy(
    search_parameters_t* const search_parameters,
    const mapping_mode_t mapping_mode);
void search_configure_quality_model(
    search_parameters_t* const search_parameters,
    const sequence_qualities_format_t qualities_format,
    const uint64_t quality_threshold);
void search_configure_alignment_model(
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

#endif /* ARCHIVE_SEARCH_PARAMETERS_H_ */
