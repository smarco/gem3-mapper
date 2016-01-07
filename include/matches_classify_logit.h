/*
 * PROJECT: GEMMapper
 * FILE: matches_classify_logit.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MATCHES_CLASSIFY_LOGIT_H_
#define MATCHES_CLASSIFY_LOGIT_H_

#include "essentials.h"
#include "align_swg.h"
#include "matches.h"
#include "paired_matches.h"
#include "matches_classify.h"

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
  double coeff_subdominant_stratum_matches;
  double coeff_mcs;
  double coeff_max_region_length_norm;
  // Subdominant candidates
  double coeff_subdominant_candidates_end1;
  double coeff_subdominant_candidates_end2;
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
