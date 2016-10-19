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

#ifndef MATCHES_CLASSIFY_LOGIT_H_
#define MATCHES_CLASSIFY_LOGIT_H_

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
  double coeff_primary_template_size_sigma_norm;
  double coeff_mapq_end1;
  double coeff_mapq_end2;
  // Sub-dominant Match
  double coeff_subdominant_edit_distance_norm;
  double coeff_subdominant_event_distance_norm;
  double coeff_subdominant_swg_score_norm;
  double coeff_subdominant_template_size_sigma_norm;
  // Search Scope
  double coeff_accepted_candidates;
  double coeff_aligned_candidates;
  double coeff_accepted_matches;
  double coeff_mcs_end1;
  double coeff_mcs_end2;
  // Mappability
  double coeff_avg_region_length_norm;
  double coeff_max_region_length_norm;
  double coeff_kmer_frequency;
} matches_classify_logit_coeff_t;
typedef struct {
  matches_classify_logit_coeff_t unique_logit_coeff;
  matches_classify_logit_coeff_t mmaps_logit_coeff;
  matches_classify_logit_coeff_t mmaps_d1_logit_coeff;
} matches_classify_logit_model_t;

/*
 * Compute Probabilities
 */
double matches_classify_logit(
    const matches_classify_logit_coeff_t* const logit_coeff,
    const matches_predictors_t* const matches_predictors);

/*
 * Matches Classify Probabilities (using Logistic Regression)
 */
double matches_classify_logit_unique(
    const matches_classify_logit_model_t* const logit_model,
    const matches_predictors_t* const matches_predictors);
double matches_classify_logit_mmaps(
    const matches_classify_logit_model_t* const logit_model,
    const matches_predictors_t* const matches_predictors);
double matches_classify_logit_mmaps_d1(
    const matches_classify_logit_model_t* const logit_model,
    const matches_predictors_t* const matches_predictors);

#endif /* MATCHES_CLASSIFY_LOGIT_H_ */
