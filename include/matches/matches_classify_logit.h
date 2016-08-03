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
  // Match
  double coeff_map_edit_distance_norm;
  double coeff_map_event_distance_norm;
  double coeff_map_swg_score_norm;
  // Template Size
  double coeff_map_template_size_sigma;
  // MAPQ
  double coeff_mapq_end1;
  double coeff_mapq_end2;
  // Best Match (best in each distance metric)
  double coeff_best_map_edit_distance_norm;
  double coeff_best_map_event_distance_norm;
  double coeff_best_map_swg_score_norm;
  // Sub-dominant Match (sub-dominant in each distance metric)
  double coeff_subdominant_edit_distance_norm;
  double coeff_subdominant_event_distance_norm;
  double coeff_subdominant_swg_score_norm;
  // Template Size
  double coeff_best_map_template_size_sigma;
  double coeff_subdominant_template_size_sigma;
  // Search Scope
  double coeff_matches_accepted;
  double coeff_mcs_end1;
  double coeff_mcs_end2;
  // Mappability
  double coeff_max_region_length_norm;
  double coeff_kmer_frequency;
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
    const matches_classify_logit_coeff_t* const logit_coeff,
    const matches_predictors_t* const matches_predictors,
    const match_predictors_t* const match_predictors);

/*
 * Matches Classify Probabilities (using Logistic Regression)
 */
double matches_classify_logit_unique(
    const matches_classify_logit_model_t* const logit_model,
    const matches_predictors_t* const matches_predictors,
    const match_predictors_t* const match_predictors);
double matches_classify_logit_mmaps(
    const matches_classify_logit_model_t* const logit_model,
    const matches_predictors_t* const matches_predictors,
    const match_predictors_t* const match_predictors);
double matches_classify_logit_ties(
    const matches_classify_logit_model_t* const logit_model,
    const matches_predictors_t* const matches_predictors,
    const match_predictors_t* const match_predictors);

#endif /* MATCHES_CLASSIFY_LOGIT_H_ */
