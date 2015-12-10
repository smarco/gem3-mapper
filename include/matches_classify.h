/*
 * PROJECT: GEMMapper
 * FILE: matches_classify.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MATCHES_CLASSIFY_H_
#define MATCHES_CLASSIFY_H_

#include "essentials.h"
#include "align_swg.h"
#include "matches.h"
#include "paired_matches.h"

/*
 * Constants
 */
#define MATCHES_MIN_CI    0.5
#define MATCHES_UNIQUE_CI 0.998
#define MATCHES_MMAPS_CI  0.995
#define MATCHES_TIES_CI   0.95

/*
 * Matches Classes
 */
typedef enum {
  matches_class_unmapped = 0,
  matches_class_tie_indistinguishable = 1,
  matches_class_tie_swg_score = 2,
  matches_class_tie_edit_distance = 3,
  matches_class_tie_event_distance = 4,
  matches_class_mmap = 5,
  matches_class_unique = 6,
} matches_class_t;
extern const char* matches_class_label[7];

/*
 * Compute predictors
 */
void matches_classify_compute_predictors_unmapped(
    matches_predictors_t* const predictors,matches_metrics_t* const metrics,
    const uint64_t read_length,const uint64_t max_region_length,
    const uint64_t proper_length,const uint64_t max_complete_stratum,
    const uint64_t num_zero_regions);
void matches_classify_compute_predictors_mapped(
    matches_predictors_t* const predictors,matches_metrics_t* const metrics,
    const uint64_t primary_map_distance,const uint64_t primary_map_edit_distance,
    const int32_t primary_map_swg_score,const swg_penalties_t* const swg_penalties,
    const uint64_t read_length,const uint64_t max_region_length,
    const uint64_t proper_length,const uint64_t max_complete_stratum,
    const uint64_t num_zero_regions);

/*
 * SE Classify
 */
matches_class_t matches_classify(matches_t* const matches);
void matches_classify_compute_predictors(
    matches_t* const matches,matches_predictors_t* const predictors,
    const swg_penalties_t* const swg_penalties,const uint64_t read_length,
    const uint64_t max_region_length,uint64_t const proper_length,
    uint64_t const overriding_mcs,const uint64_t num_zero_regions);

double matches_classify_unique(matches_predictors_t* const predictors);
double matches_classify_mmaps(matches_predictors_t* const predictors);
double matches_classify_ties(matches_predictors_t* const predictors);

#endif /* MATCHES_CLASSIFY_H_ */
