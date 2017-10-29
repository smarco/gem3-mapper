/*
 * PROJECT: GEM-Tools library
 * FILE: gt_counters_utils.h
 * DATE: 20/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_COUNTERS_UTILS_H_
#define GT_COUNTERS_UTILS_H_

#include "gt_essentials.h"

/*
 * Accessors
 */
GT_INLINE uint64_t gt_counters_get_num_counters(gt_vector* const counters);
GT_INLINE uint64_t gt_counters_get_counter(gt_vector* const counters,const uint64_t stratum);
GT_INLINE void gt_counters_dynamically_allocate_counter(gt_vector* const counters,const uint64_t stratum);
GT_INLINE void gt_counters_set_counter(gt_vector* const counters,const uint64_t stratum,const uint64_t value);
GT_INLINE void gt_counters_dec_counter(gt_vector* const counters,const uint64_t stratum);
GT_INLINE void gt_counters_inc_counter(gt_vector* const counters,const uint64_t stratum);

/*
 * General Operation with counters
 */
GT_INLINE int64_t gt_counters_get_uniq_degree(gt_vector* const counters);
GT_INLINE bool gt_counters_get_next_matching_strata(
    gt_vector* const counters,const uint64_t begin_strata,
    uint64_t* const next_matching_strata,uint64_t* const num_maps);
GT_INLINE int64_t gt_counters_get_min_matching_strata(gt_vector* const counters);
GT_INLINE void gt_counters_calculate_num_maps(
    gt_vector* const counters,const uint64_t min_decoded_strata,const uint64_t max_decoded_matches,
    uint64_t* num_strata,uint64_t* num_matches);
GT_INLINE uint64_t gt_counters_reduce_sum(gt_vector* const counters);

#endif /* GT_COUNTERS_UTILS_H_ */
