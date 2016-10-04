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
#include "matches/matches.h"
#include "matches/paired_matches.h"
#include "matches/classify/matches_classify.h"

/*
 * Predictors
 */
typedef struct {
  /* Primary Match */
  double primary_edit_distance_norm;
  double primary_event_distance_norm;
  double primary_swg_score_norm;
  double primary_template_size_sigma_norm;
  /* Subdominant Match */
  double subdominant_edit_distance_norm;
  double subdominant_event_distance_norm;
  double subdominant_swg_score_norm;
  double subdominant_template_size_sigma_norm;
  /* Search Scope */
  double accepted_candidates;
  double accepted_matches;
  uint64_t mcs_end1;
  uint64_t mcs_end2;
  /* Mappability */
  double max_region_length_norm;
  double kmer_frequency;
  /* MAPQ */
  uint8_t mapq_end1;
  uint8_t mapq_end2;
} matches_predictors_t;

/*
 * Matches Compute Predictors
 */
void matches_predictors_compute_se(
    matches_predictors_t* const predictors,
    matches_t* const matches);
void matches_predictors_compute_pe(
    matches_predictors_t* const predictors,
    paired_matches_t* const paired_matches);

/*
 * Display
 */
void matches_predictors_se_print(
    FILE* const stream,
    const char* const sequence_tag,
    const matches_class_t matches_class,
    matches_predictors_t* const matches_predictors);
void matches_predictors_pe_print(
    FILE* const stream,
    const char* const sequence_tag,
    const paired_matches_class_t paired_matches_class,
    matches_predictors_t* const matches_predictors);

#endif /* MATCHES_PREDICTORS_H_ */
