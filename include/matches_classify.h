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
  matches_class_unmapped,
  matches_class_tie_indistinguishable,
  matches_class_tie_swg_score,
  matches_class_tie_edit_distance,
  matches_class_tie_event_distance,
  matches_class_mmap,
  matches_class_unique,
} matches_class_t;

/*
 * Compute predictors
 */
GEM_INLINE void matches_classify_compute_predictors_unmapped(
    matches_predictors_t* const predictors,matches_metrics_t* const metrics,
    const uint64_t read_length,const uint64_t max_region_length,
    const uint64_t proper_length,const uint64_t max_complete_stratum,
    const uint64_t num_zero_regions);
GEM_INLINE void matches_classify_compute_predictors_mapped(
    matches_predictors_t* const predictors,matches_metrics_t* const metrics,
    const uint64_t primary_map_distance,const uint64_t primary_map_edit_distance,
    const int32_t primary_map_swg_score,const swg_penalties_t* const swg_penalties,
    const uint64_t read_length,const uint64_t max_region_length,
    const uint64_t proper_length,const uint64_t max_complete_stratum,
    const uint64_t num_zero_regions);

/*
 * SE Classify
 */
GEM_INLINE matches_class_t matches_classify(matches_t* const matches);
GEM_INLINE void matches_classify_compute_predictors(
    matches_t* const matches,matches_predictors_t* const predictors,
    const swg_penalties_t* const swg_penalties,const uint64_t read_length,
    const uint64_t max_region_length,uint64_t const proper_length,
    uint64_t const overriding_mcs,const uint64_t num_zero_regions);

GEM_INLINE double matches_classify_unique(matches_predictors_t* const predictors);
GEM_INLINE double matches_classify_mmaps(matches_predictors_t* const predictors);
GEM_INLINE double matches_classify_ties(matches_predictors_t* const predictors);

#endif /* MATCHES_CLASSIFY_H_ */
