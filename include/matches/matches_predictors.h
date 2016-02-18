/*
 * PROJECT: GEMMapper
 * FILE: matches_predictors.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MATCHES_PREDICTORS_H_
#define MATCHES_PREDICTORS_H_

#include "utils/essentials.h"
#include "approximate_search/approximate_search_metrics.h"
#include "matches/matches.h"
#include "matches/paired_matches.h"
#include "matches/matches_classify.h"

/*
 * Matches-Predictors (Calculated metrics from a given set of matches & sorting)
 */
typedef struct {
  /* First Map (Primary Match) */
  uint64_t first_map_edit_distance;
  uint64_t first_map_event_distance;
  int32_t first_map_swg_score;
  double first_map_edit_distance_norm;
  double first_map_event_distance_norm;
  double first_map_swg_score_norm;
  /* Sub-dominant Matches (best in each distance metric) */
  uint64_t subdominant_edit_distance;
  uint64_t subdominant_event_distance;
  int32_t subdominant_swg_score;
  double subdominant_edit_distance_norm;
  double subdominant_event_distance_norm;
  double subdominant_swg_score_norm;
  /* Search Scope */
  uint64_t first_stratum_matches;
  uint64_t mcs_end1;
  uint64_t mcs_end2;
  uint64_t accepted_candidates_end1;
  uint64_t accepted_candidates_end2;
  /* Mappability */
  double max_region_length_norm;
  double mappability_p;
  double mappability_2p;
  /* Template Size */
  double first_map_template_size_sigma;
  double subdominant_template_size_sigma;
  /* MAPQ Score */
  uint8_t mapq_end1;
  uint8_t mapq_end2;
} matches_predictors_t;

/*
 * SE Compute Predictors
 */
void matches_predictors_compute(
    matches_t* const matches,
    matches_predictors_t* const predictors,
    approximate_search_metrics_t* const search_metrics,
    const uint64_t mcs);

/*
 * PE Compute Predictors
 */
void paired_matches_predictors_compute(
    paired_matches_t* const paired_matches,
    matches_predictors_t* const predictors,
    approximate_search_metrics_t* const search_metrics_end1,
    approximate_search_metrics_t* const search_metrics_end2,
    const uint64_t mcs_end1,
    const uint64_t mcs_end2);

/*
 * Display
 */
void matches_predictors_se_print(
    FILE* const stream,
    const char* const sequence_tag,
    const matches_class_t matches_class,
    matches_predictors_t* const predictors);
void matches_predictors_pe_print(
    FILE* const stream,
    const char* const sequence_tag,
    const paired_matches_class_t paired_matches_class,
    matches_predictors_t* const predictors);

#endif /* MATCHES_PREDICTORS_H_ */
