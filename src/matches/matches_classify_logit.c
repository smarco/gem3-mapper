/*
 * PROJECT: GEMMapper
 * FILE: matches_classify_logit.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: TODO
 */

#include "matches/matches_classify_logit.h"
#include "matches/matches_metrics.h"

/*
 * Compute Probabilities
 */
double matches_classify_logit(
    const matches_classify_logit_coeff_t* const logit_coeff,
    const matches_predictors_t* const matches_predictors,
    const match_predictors_t* const match_predictors) {
  // Compute lr-factor
  const double lr_factor =
      logit_coeff->coeff_intercept
      +
      logit_coeff->coeff_map_edit_distance_norm *
      (double)match_predictors->map_edit_distance_norm
      +
      logit_coeff->coeff_map_event_distance_norm *
      (double)match_predictors->map_event_distance_norm
      +
      logit_coeff->coeff_map_swg_score_norm *
      (double)match_predictors->map_swg_score_norm
      +
      logit_coeff->coeff_map_template_size_sigma *
      (double)match_predictors->map_template_size_sigma
      +
      logit_coeff->coeff_mapq_end1 *
      (double)match_predictors->mapq_end1
      +
      logit_coeff->coeff_mapq_end2 *
      (double)match_predictors->mapq_end2
      +
      logit_coeff->coeff_best_map_edit_distance_norm *
      (double)matches_predictors->best_map_edit_distance_norm
      +
      logit_coeff->coeff_best_map_event_distance_norm *
      (double)matches_predictors->best_map_event_distance_norm
      +
      logit_coeff->coeff_best_map_swg_score_norm *
      (double)matches_predictors->best_map_swg_score_norm
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
      logit_coeff->coeff_best_map_template_size_sigma *
      (double)matches_predictors->best_map_template_size_sigma
      +
      logit_coeff->coeff_subdominant_template_size_sigma *
      (double)matches_predictors->subdominant_template_size_sigma
      +
      logit_coeff->coeff_matches_accepted *
      (double)matches_predictors->matches_accepted
      +
      logit_coeff->coeff_mcs_end1 *
      (double)matches_predictors->mcs_end1
      +
      logit_coeff->coeff_mcs_end2 *
      (double)matches_predictors->mcs_end2
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
double matches_classify_logit_unique(
    const matches_classify_logit_model_t* const logit_model,
    const matches_predictors_t* const matches_predictors,
    const match_predictors_t* const match_predictors) {
  return matches_classify_logit(&logit_model->unique_logit_coeff,matches_predictors,match_predictors);
}
double matches_classify_logit_mmaps(
    const matches_classify_logit_model_t* const logit_model,
    const matches_predictors_t* const matches_predictors,
    const match_predictors_t* const match_predictors) {
  return matches_classify_logit(&logit_model->mmaps_logit_coeff,matches_predictors,match_predictors);
}
double matches_classify_logit_ties(
    const matches_classify_logit_model_t* const logit_model,
    const matches_predictors_t* const matches_predictors,
    const match_predictors_t* const match_predictors) {
  return matches_classify_logit(&logit_model->ties_logit_coeff,matches_predictors,match_predictors);
}
