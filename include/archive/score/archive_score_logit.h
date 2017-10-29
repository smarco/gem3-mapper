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
 */

#ifndef ARCHIVE_SCORE_LOGIT_H_
#define ARCHIVE_SCORE_LOGIT_H_

#include "utils/essentials.h"
#include "align/align_swg.h"
#include "matches/matches.h"
#include "matches/classify/matches_predictors.h"
#include "matches/classify/matches_classify.h"
#include "matches/paired_matches.h"

/*
 * Classify logit-coefficients
 */
typedef struct {
  double coeff_intercept;
  // Best Match
  double coeff_primary_edit_distance_norm;
  double coeff_primary_event_distance_norm;
  double coeff_primary_swg_score_norm;
  double coeff_primary_error_quality;
  double coeff_primary_template_size_sigma_norm;
  double coeff_mapq_end1;
  double coeff_mapq_end2;
  // Sub-dominant Match
  double coeff_subdominant_edit_distance_norm;
  double coeff_subdominant_event_distance_norm;
  double coeff_subdominant_swg_score_norm;
  double coeff_subdominant_template_size_sigma_norm;
  // Search Scope
  double coeff_candidates_subdominant;
  double coeff_candidates_accepted;
  double coeff_mcs_end1;
  double coeff_mcs_end2;
  double coeff_delta_group;
  double coeff_wdelta_group;
  // Mappability
  double coeff_sequence_length_norm;
  double coeff_avg_region_length_norm;
  double coeff_max_region_length_norm;
  double coeff_kmer_frequency;
} archive_score_logit_coeff_t;
typedef struct {
  archive_score_logit_coeff_t unique_logit_coeff;
  archive_score_logit_coeff_t mmap_logit_coeff;
} archive_score_logit_model_t;

/*
 * Compute Probabilities
 */
double archive_score_logit(
    const archive_score_logit_coeff_t* const logit_coeff,
    const matches_predictors_t* const matches_predictors,
    const matches_classification_t* const matches_classification);

/*
 * Matches Classify Probabilities (using Logistic Regression)
 */
double archive_score_logit_unique(
    const archive_score_logit_model_t* const logit_model,
    const matches_predictors_t* const matches_predictors,
    const matches_classification_t* const matches_classification);
double archive_score_logit_mmap(
    const archive_score_logit_model_t* const logit_model,
    const matches_predictors_t* const matches_predictors,
    const matches_classification_t* const matches_classification);

#endif /* ARCHIVE_SCORE_LOGIT_H_ */
