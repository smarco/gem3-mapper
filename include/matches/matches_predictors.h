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
#include "matches/matches_classify.h"

/*
 * Predictors
 */
typedef struct {
  /* Match */
  double map_edit_distance_norm;
  double map_event_distance_norm;
  double map_swg_score_norm;
  /* Template Size */
  double map_template_size_sigma;
  /* MAPQ Score */
  uint8_t mapq_end1;
  uint8_t mapq_end2;
} match_predictors_t;
typedef struct {
  /* Best Match (Best in each distance metric) */
  double best_map_edit_distance_norm;
  double best_map_event_distance_norm;
  double best_map_swg_score_norm;
  /* Sub-dominant Match (Best in each distance metric) */
  double subdominant_edit_distance_norm;
  double subdominant_event_distance_norm;
  double subdominant_swg_score_norm;
  /* Template Size */
  double best_map_template_size_sigma;
  double subdominant_template_size_sigma;
  /* Search Scope */
  uint64_t matches_accepted;
  uint64_t mcs_end1;
  uint64_t mcs_end2;
  /* Mappability */
  double max_region_length_norm;
  double kmer_frequency;
} matches_predictors_t;

/*
 * Match Compute Predictors
 */
void match_predictors_compute_se(
    match_predictors_t* const match_predictors,
    matches_t* const matches,
    match_trace_t* const match);
void match_predictors_compute_pe(
    match_predictors_t* const match_predictors,
    paired_matches_t* const paired_matches,
    paired_map_t* const paired_map);

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
    matches_predictors_t* const matches_predictors,
    match_predictors_t* const match_predictors);
void matches_predictors_pe_print(
    FILE* const stream,
    const char* const sequence_tag,
    const paired_matches_class_t paired_matches_class,
    matches_predictors_t* const matches_predictors,
    match_predictors_t* const match_predictors);

#endif /* MATCHES_PREDICTORS_H_ */
