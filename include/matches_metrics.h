/*
 * PROJECT: GEMMapper
 * FILE: matches_classify.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MATCHES_METRICS_H_
#define MATCHES_METRICS_H_

#include "essentials.h"

/*
 * Constants
 */
#define MAX_TEMPLATE_LENGTH_SIGMAS 10.0

/*
 * Matches-Metrics (Current metrics status -always being updated-)
 */
typedef struct {
  uint64_t total_matches_sampled;    // Total matches sampled
  uint64_t min1_counter_value;       // Minimum non-zero counter position (for the distance metric the counters use)
  uint64_t min2_counter_value;       // Second minimum non-zero counter position (for the distance metric the counters use)
  uint64_t max_counter_value;        // Maximum non-zero counter position (for the distance metric the counters use)
  uint64_t min1_edit_distance;       // Minimum edit distance among all found matches
  uint64_t min2_edit_distance;       // Second minimum edit distance among all found matches
  int32_t max1_swg_score;            // Maximum smith-waterman-gotoh score among all found matches
  int32_t max2_swg_score;            // Second maximum smith-waterman-gotoh score among all found matches
  uint64_t subdominant_candidates;   // Number of subdominant candidates (discarded during post-filtering)
  // PE specific
  double min1_template_length_sigma; // Minimum number of sigma deviation from the mean template-length
  double min2_template_length_sigma; // Second minimum number of sigma deviation from the mean template-length
} matches_metrics_t;

/*
 * Matches-Predictors (Calculated metrics from a given set of matches & sorting)
 */
typedef struct {
  // First Map (Primary Match)
  uint64_t first_map_edit_distance;
  uint64_t first_map_event_distance;
  int32_t first_map_swg_score;
  double first_map_edit_distance_norm;
  double first_map_event_distance_norm;
  double first_map_swg_score_norm;
  // Sub-dominant Matches (best in each distance metric)
  uint64_t subdominant_edit_distance;
  uint64_t subdominant_event_distance;
  int32_t subdominant_swg_score;
  double subdominant_edit_distance_norm;
  double subdominant_event_distance_norm;
  double subdominant_swg_score_norm;
  // Search Scope
  uint64_t first_stratum_matches;
  uint64_t subdominant_stratum_matches;
  uint64_t mcs;
  uint64_t max_region_length;
  double max_region_length_norm;
  // Subdominant candidates
  uint64_t subdominant_candidates_end1;
  uint64_t subdominant_candidates_end2;
  // PE specific
  double first_map_template_size_sigma;
  double subdominant_template_size_sigma;
  uint8_t mapq_end1;
  uint8_t mapq_end2;
} matches_predictors_t;

/*
 * Setup
 */
GEM_INLINE void matches_metrics_init(matches_metrics_t* const metrics);

/*
 * Accessors
 */
GEM_INLINE uint64_t matches_metrics_get_min_distance(matches_metrics_t* const metrics);
GEM_INLINE uint64_t matches_metrics_get_max_distance(matches_metrics_t* const metrics);
GEM_INLINE uint64_t matches_metrics_get_min_edit_distance(matches_metrics_t* const metrics);
GEM_INLINE int32_t matches_metrics_get_max_swg_score(matches_metrics_t* const metrics);

/*
 * Update
 */
GEM_INLINE void matches_metrics_update(
    matches_metrics_t* const matches_metrics,const uint64_t distance,
    const uint64_t edit_distance,const int32_t swg_score);
GEM_INLINE void matches_metrics_pe_update(
    matches_metrics_t* const matches_metrics,const uint64_t distance,
    const uint64_t edit_distance,const int32_t swg_score,
    const double template_length_sigma);

GEM_INLINE void matches_metrics_inc_subdominant_candidates(matches_metrics_t* const metrics);
GEM_INLINE void matches_metrics_dec_subdominant_candidates(matches_metrics_t* const metrics);

/*
 * Display
 */
GEM_INLINE void matches_predictors_print(
    matches_predictors_t* const predictors,
    const char* const read_tag,const uint8_t mapq_score);
GEM_INLINE void paired_matches_predictors_print(
    matches_predictors_t* const predictors,
    const char* const read_tag,const uint8_t mapq_score_pair,
    const uint8_t mapq_score_end1,const uint8_t mapq_score_end2);

#endif /* MATCHES_METRICS_H_ */
