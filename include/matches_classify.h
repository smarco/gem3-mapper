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

typedef enum {
  matches_tie_none = 0,
  matches_tie_event_distance_delta0,
  matches_tie_event_distance_delta1,
  matches_tie_swg_score,
  matches_tie_edit_distance,
} matches_tie_type;

/*
 * Setup
 */
GEM_INLINE void matches_classify_metrics_init(matches_metrics_t* const matches_metrics);

/*
 * Accessors
 */
GEM_INLINE uint64_t matches_classify_metrics_get_min_distance(matches_t* const matches);
GEM_INLINE uint64_t matches_classify_metrics_get_max_distance(matches_t* const matches);
GEM_INLINE uint64_t matches_classify_metrics_get_min_edit_distance(matches_t* const matches);
GEM_INLINE int32_t matches_classify_metrics_get_max_swg_score(matches_t* const matches);

/*
 * Metrics
 */
GEM_INLINE void matches_classify_metrics_update(
    matches_metrics_t* const matches_metrics,const uint64_t distance,
    const uint64_t edit_distance,const int32_t swg_score);
GEM_INLINE void matches_classify_metrics_recompute(matches_t* const matches);

GEM_INLINE void matches_classify_compute_mvalues(
    matches_t* const matches,const swg_penalties_t* const swg_penalties,
    const uint64_t read_length,const uint64_t max_region_length,
    uint64_t const proper_length,uint64_t const overriding_mcs);

/*
 * Classifier
 */
GEM_INLINE matches_tie_type matches_classify_is_tie(matches_t* const matches);
GEM_INLINE bool matches_classify_is_unique(matches_t* const matches);

GEM_INLINE double matches_classify_unique(matches_t* const matches);
GEM_INLINE double matches_classify_mmaps(matches_t* const matches);
GEM_INLINE double matches_classify_delta1(matches_t* const matches);

/*
 * Display
 */
GEM_INLINE void matches_classify_metrics_print(const char* const read_tag,matches_t* const matches);

#endif /* MATCHES_CLASSIFY_H_ */
