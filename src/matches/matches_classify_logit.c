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
    const matches_predictors_t* const predictors,
    const matches_classify_logit_coeff_t* const logit_coeff) {
  // Compute lr-factor
  const double lr_factor = logit_coeff->coeff_intercept
      +
       logit_coeff->coeff_first_map_edit_distance_norm *
      (double)predictors->first_map_edit_distance_norm
      +
       logit_coeff->coeff_first_map_event_distance_norm *
      (double)predictors->first_map_event_distance_norm
      +
       logit_coeff->coeff_first_map_swg_score_norm *
      (double)predictors->first_map_swg_score_norm
      +
       logit_coeff->coeff_subdominant_edit_distance_norm *
      (double)predictors->subdominant_edit_distance_norm
      +
       logit_coeff->coeff_subdominant_event_distance_norm *
      (double)predictors->subdominant_event_distance_norm
      +
       logit_coeff->coeff_subdominant_swg_score_norm *
      (double)predictors->subdominant_swg_score_norm
      +
      logit_coeff->coeff_first_stratum_matches *
     (double)predictors->first_stratum_matches
      +
       logit_coeff->coeff_mcs_end1 *
      (double)predictors->mcs_end1
      +
      logit_coeff->coeff_mcs_end2 *
     (double)predictors->mcs_end2
      +
       logit_coeff->coeff_accepted_candidates_end1 *
      (double)predictors->accepted_candidates_end1
      +
       logit_coeff->coeff_accepted_candidates_end2 *
      (double)predictors->accepted_candidates_end2
      +
       logit_coeff->coeff_max_region_length_norm *
      (double)predictors->max_region_length_norm
      +
       logit_coeff->coeff_mappability_p *
      (double)predictors->mappability_p
      +
       logit_coeff->coeff_mappability_2p *
      (double)predictors->mappability_2p
      +
       logit_coeff->coeff_first_map_template_size_sigma *
      (double)predictors->first_map_template_size_sigma
      +
       logit_coeff->coeff_subdominant_template_size_sigma *
      (double)predictors->subdominant_template_size_sigma
      +
       logit_coeff->coeff_mapq_end1 *
      (double)predictors->mapq_end1
      +
       logit_coeff->coeff_mapq_end2 *
      (double)predictors->mapq_end2;
  // Return probability
  return 1.0 / (1.0 + (1.0/exp(lr_factor)));
}
/*
 * Matches Classify Probabilities (using Logistic Regression)
 */
double matches_classify_logit_unique(
    const matches_predictors_t* const predictors,
    const matches_classify_logit_model_t* const logit_model) {
  // Unique: Probability of the first position-match (primary match) of being a true positive
  return matches_classify_logit(predictors,&logit_model->unique_logit_coeff);
}
double matches_classify_logit_mmaps(
    const matches_predictors_t* const predictors,
    const matches_classify_logit_model_t* const logit_model) {
  // Classify MMaps wrt the probability of the first position-match (primary match) of being a true positive
  return matches_classify_logit(predictors,&logit_model->mmaps_logit_coeff);
}
double matches_classify_logit_ties(
    const matches_predictors_t* const predictors,
    const matches_classify_logit_model_t* const logit_model) {
  // Classify ties wrt the probability of being a true positive
  return matches_classify_logit(predictors,&logit_model->ties_logit_coeff);
}
