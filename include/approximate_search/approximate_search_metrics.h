/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_metrics.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef APPROXIMATE_SEARCH_METRICS_H_
#define APPROXIMATE_SEARCH_METRICS_H_

#include "utils/essentials.h"

typedef struct {
  /* Search Magnitudes */
  double proper_length;
  uint64_t read_length;
  int32_t swg_match_score;
  /* Mappability */
  uint64_t max_region_length;        // Maximum region length
  uint64_t num_zero_regions;         // Number of zero regions
  double mappability_p;              // Approx. Mappability-0 at proper-length
  double mappability_2p;             // Approx. Mappability-0 at 2 times the proper-length
} approximate_search_metrics_t;

/*
 * Setup
 */
void approximate_search_metrics_init(
    approximate_search_metrics_t* const search_metrics,
    const double proper_length,
    const uint64_t read_length,
    const int32_t swg_match_score);

/*
 * Accessors
 */
void approximate_search_metrics_set_max_region_length(
    approximate_search_metrics_t* const search_metrics,
    const uint64_t max_region_length);
void approximate_search_metrics_set_num_zero_regions(
    approximate_search_metrics_t* const search_metrics,
    const uint64_t num_zero_regions);
void approximate_search_metrics_set_mappability(
    approximate_search_metrics_t* const search_metrics,
    const double mappability_p,
    const double mappability_2p);

/*
 * Display
 */
void approximate_search_metrics_print(
    FILE* const stream,
    approximate_search_metrics_t* const search_metrics);

#endif /* APPROXIMATE_SEARCH_METRICS_H_ */
