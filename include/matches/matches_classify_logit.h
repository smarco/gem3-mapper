/*
 * PROJECT: GEMMapper
 * FILE: matches_classify_logit.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MATCHES_CLASSIFY_LOGIT_H_
#define MATCHES_CLASSIFY_LOGIT_H_

#include "utils/essentials.h"
#include "align/align_swg.h"
#include "matches/matches.h"
#include "matches/matches_predictors.h"
#include "matches/matches_classify.h"
#include "matches/paired_matches.h"

/*
 * Classify logit-coefficients
 */
typedef struct {
  double coeff_intercept;
  // First Map (Primary Match)
  double coeff_first_map_edit_distance_norm;
  double coeff_first_map_event_distance_norm;
  double coeff_first_map_swg_score_norm;
  // Sub-dominant Matches (best in each distance metric)
  double coeff_subdominant_edit_distance_norm;
  double coeff_subdominant_event_distance_norm;
  double coeff_subdominant_swg_score_norm;
  // Search Scope
  double coeff_first_stratum_matches;
  double coeff_mcs_end1;
  double coeff_mcs_end2;
  double coeff_accepted_candidates_end1;
  double coeff_accepted_candidates_end2;
  // Mappability
  double coeff_max_region_length_norm;
  double coeff_mappability_p;
  double coeff_mappability_2p;
  // PE specific
  double coeff_first_map_template_size_sigma;
  double coeff_subdominant_template_size_sigma;
  // MAPQ
  double coeff_mapq_end1;
  double coeff_mapq_end2;
} matches_classify_logit_coeff_t;
typedef struct {
  matches_classify_logit_coeff_t unique_logit_coeff;
  matches_classify_logit_coeff_t mmaps_logit_coeff;
  matches_classify_logit_coeff_t ties_logit_coeff;
} matches_classify_logit_model_t;

/*
 * Compute Probabilities
 */
double matches_classify_logit(
    const matches_predictors_t* const predictors,
    const matches_classify_logit_coeff_t* const logit_coeff);

/*
 * Matches Classify Probabilities (using Logistic Regression)
 */
double matches_classify_logit_unique(
    const matches_predictors_t* const predictors,
    const matches_classify_logit_model_t* const logit_model);
double matches_classify_logit_mmaps(
    const matches_predictors_t* const predictors,
    const matches_classify_logit_model_t* const logit_model);
double matches_classify_logit_ties(
    const matches_predictors_t* const predictors,
    const matches_classify_logit_model_t* const logit_model);

#endif /* MATCHES_CLASSIFY_LOGIT_H_ */
