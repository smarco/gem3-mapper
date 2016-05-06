/*
 * PROJECT: GEMMapper
 * FILE: matches_classify.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MATCHES_METRICS_H_
#define MATCHES_METRICS_H_

#include "utils/essentials.h"

/*
 * Constants
 */
#define MAX_TEMPLATE_LENGTH_SIGMAS 10.0

/*
 * Matches-Metrics (Current metrics status -always being updated-)
 */
typedef struct {
  /* Aggregated */
  uint64_t total_matches_sampled;    // Total matches sampled
  uint64_t accepted_candidates;      // Number of accepted candidates
  /* Minimums */
  uint64_t min1_counter_value;       // Minimum non-zero counter position (for the distance metric the counters use)
  uint64_t min2_counter_value;       // Second minimum non-zero counter position (for the distance metric the counters use)
  uint64_t min1_edit_distance;       // Minimum edit distance among all found matches
  uint64_t min2_edit_distance;       // Second minimum edit distance among all found matches
  int32_t max1_swg_score;            // Maximum smith-waterman-gotoh score among all found matches
  int32_t max2_swg_score;            // Second maximum smith-waterman-gotoh score among all found matches
  /* Template length */
  double min1_template_length_sigma; // Minimum number of sigma deviation from the mean template-length
  double min2_template_length_sigma; // Second minimum number of sigma deviation from the mean template-length
  /* MAPQ */
  uint8_t mapq;                      // MAPQ score of the primary alignment
} matches_metrics_t;

/*
 * Setup
 */
void matches_metrics_init(matches_metrics_t* const metrics);

/*
 * Accessors
 */
uint64_t matches_metrics_get_min_distance(matches_metrics_t* const metrics);
uint64_t matches_metrics_get_min_edit_distance(matches_metrics_t* const metrics);
int32_t matches_metrics_get_max_swg_score(matches_metrics_t* const metrics);

void matches_metrics_add_accepted_candidates(matches_metrics_t* const metrics,const uint64_t num_candidates);
void matches_metrics_set_mapq(matches_metrics_t* const metrics,const uint8_t mapq);

void matches_metrics_update(
    matches_metrics_t* const matches_metrics,
    const uint64_t distance,
    const uint64_t edit_distance,
    const int32_t swg_score);
void paired_matches_metrics_update(
    matches_metrics_t* const matches_metrics,
    const uint64_t distance,
    const uint64_t edit_distance,
    const int32_t swg_score,
    const double template_length_sigma);

/*
 * Display
 */
void matches_metrics_print(FILE* const stream,matches_metrics_t* const matches_metrics);

#endif /* MATCHES_METRICS_H_ */
