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

#include "archive/score/archive_score_logit.h"
#include "matches/classify/matches_metrics.h"

/*
 * Compute Probabilities
 */
double archive_score_logit(
    const archive_score_logit_coeff_t* const logit_coeff,
    const matches_predictors_t* const matches_predictors,
    const matches_classification_t* const matches_classification) {
  // Compute lr-factor
  const double lr_factor =
      logit_coeff->coeff_intercept
      +
      logit_coeff->coeff_primary_edit_distance_norm *
      (double)matches_predictors->primary_edit_distance_norm
      +
      logit_coeff->coeff_primary_event_distance_norm *
      (double)matches_predictors->primary_event_distance_norm
      +
      logit_coeff->coeff_primary_swg_score_norm *
      (double)matches_predictors->primary_swg_score_norm
      +
      logit_coeff->coeff_primary_error_quality *
      (double)matches_predictors->primary_error_quality
      +
      logit_coeff->coeff_primary_template_size_sigma_norm *
      (double)matches_predictors->primary_template_size_sigma_norm
      +
      logit_coeff->coeff_mapq_end1 *
      (double)matches_predictors->mapq_end1
      +
      logit_coeff->coeff_mapq_end2 *
      (double)matches_predictors->mapq_end2
      +
      logit_coeff->coeff_subdominant_edit_distance_norm *
      (double)matches_predictors->subdominant_edit_distance_norm
      +
      logit_coeff->coeff_subdominant_event_distance_norm *
      (double)matches_predictors->subdominant_event_distance_norm
      +
      logit_coeff->coeff_subdominant_swg_score_norm *
      (double)matches_predictors->subdominant_swg_score_norm
      +
      logit_coeff->coeff_subdominant_template_size_sigma_norm *
      (double)matches_predictors->subdominant_template_size_sigma_norm
      +
      logit_coeff->coeff_candidates_subdominant *
      (double)matches_predictors->candidates_subdominant
      +
      logit_coeff->coeff_candidates_accepted *
      (double)matches_predictors->candidates_accepted
      +
      logit_coeff->coeff_mcs_end1 *
      (double)matches_predictors->mcs_end1
      +
      logit_coeff->coeff_mcs_end2 *
      (double)matches_predictors->mcs_end2
      +
      logit_coeff->coeff_delta_group *
      (double)matches_classification->delta_group
      +
      logit_coeff->coeff_wdelta_group *
      (double)matches_classification->wdelta_group
      +
      logit_coeff->coeff_sequence_length_norm *
      (double)matches_predictors->sequence_length_norm
      +
      logit_coeff->coeff_avg_region_length_norm *
      (double)matches_predictors->avg_region_length_norm
      +
      logit_coeff->coeff_max_region_length_norm *
      (double)matches_predictors->max_region_length_norm
      +
      logit_coeff->coeff_kmer_frequency *
      (double)matches_predictors->kmer_frequency;
  // Return probability
  return 1.0 / (1.0 + (1.0/exp(lr_factor)));
}
/*
 * Matches Classify Probabilities (using Logistic Regression)
 */
double archive_score_logit_unique(
    const archive_score_logit_model_t* const logit_model,
    const matches_predictors_t* const matches_predictors,
    const matches_classification_t* const matches_classification) {
  return archive_score_logit(&logit_model->unique_logit_coeff,matches_predictors,matches_classification);
}
double archive_score_logit_mmap(
    const archive_score_logit_model_t* const logit_model,
    const matches_predictors_t* const matches_predictors,
    const matches_classification_t* const matches_classification) {
  return archive_score_logit(&logit_model->mmap_logit_coeff,matches_predictors,matches_classification);
}
