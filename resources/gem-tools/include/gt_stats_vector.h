/*
 * PROJECT: GEM-Tools library
 * FILE: gt_stats_vector.h
 * DATE: 10/12/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_STATS_VECTOR_H_
#define GT_STATS_VECTOR_H_

#include "gt_commons.h"
#include "gt_compact_dna_string.h"
#include "gt_alignment_utils.h"
#include "gt_template_utils.h"

/*
 * Checkers
 */
#define GT_STATS_VECTOR_CHECK(gt_stats_vector) \
  GT_NULL_CHECK(gt_stats_vector); \
  GT_NULL_CHECK((gt_stats_vector)->counters); \
  GT_HASH_CHECK((gt_stats_vector)->out_values); \
  GT_ZERO_CHECK((gt_stats_vector)->num_values); \
  switch ((gt_stats_vector)->type) { \
    case GT_STATS_VECTOR_CUSTOMED_RANGE: \
      GT_NULL_CHECK(customed_range_values); \
      break; \
    case GT_STATS_VECTOR_STEP_RANGE: \
      gt_cond_fatal_error(min_value > max_value,VSTATS_INVALID_MIN_MAX); \
      break; \
    case GT_STATS_VECTOR_RAW: break; \
    default: \
      GT_INVALID_CASE(); \
      break; \
  }

/*
 * Constants
 */
#define GT_STATS_VECTOR_OUT_OF_RANGE UINT64_MAX

/*
 * Stats Vector
 */
typedef enum { GT_STATS_VECTOR_CUSTOMED_RANGE, GT_STATS_VECTOR_STEP_RANGE, GT_STATS_VECTOR_RAW } gt_stats_vector_t;
typedef struct _gt_stats_vector gt_stats_vector; // Forward declaration of gt_stats_vector
struct _gt_stats_vector {
  gt_stats_vector_t type;
  uint64_t* counters;
  /* GT_STATS_VECTOR_STEP_RANGE */
  uint64_t min_value;
  uint64_t max_value;
  uint64_t step;
  /* GT_STATS_VECTOR_CUSTOMED_RANGE */
  uint64_t* customed_range_values;
  uint64_t num_values;
  /* Values out of range */
  uint64_t out_of_range_bucket_size;
  gt_ihash* out_values;
  /* Nested */
  gt_stats_vector* template_vector; // Pattern as to spawn new vectors
  gt_stats_vector* nested_vectors;  // The nested vectors themselves
};
typedef struct {
  gt_stats_vector* stats_vector;
  /* Sliced iteration */
  uint64_t start_index;
  uint64_t end_index;
  /* Iteration state */
  bool eoi; // End-of-iteration
  /* Array State */
  bool is_index_in_array;
  uint64_t array_index;
  /* Hash State */
  gt_ihash_iterator* ihash_iterator;
} gt_stats_vector_iterator;

/*
 * Customed Range Stat-Vector. Eg.
 *   @range[] = { 0, 1, 2, 3 };
 *     @stats_vector[0] => [0,1)
 *     @stats_vector[1] => [1,2)
 *     @stats_vector[2] => [2,3)
 */
GT_INLINE gt_stats_vector* gt_stats_vector_customed_range_new(
    const uint64_t* const customed_range_values,const uint64_t num_values,
    const uint64_t out_of_range_bucket_size);
/*
 * Step Range Stat-Vector. Eg.
 *   @min_value = 0;
 *   @max_value = 6;
 *   @step = 3
 *     @stats_vector[0] => {0,1,2}
 *     @stats_vector[1] => {3,4,5}
 *     @stats_vector[2] => {6,7,8}
 */
GT_INLINE gt_stats_vector* gt_stats_vector_step_range_new(
    const uint64_t min_value,const uint64_t max_value,const uint64_t step,
    const uint64_t out_of_range_bucket_size);
/*
 * Raw Stat-Vector. Eg.
 *   @num_values = 3
 *     @stats_vector[0] => {0}
 *     @stats_vector[1] => {1}
 *     @stats_vector[2] => {2}
 */
GT_INLINE gt_stats_vector* gt_stats_vector_raw_new(
    const uint64_t num_values,const uint64_t out_of_range_bucket_size);

/*
 * Creates a gt_stats_vector from a template (doesn't copy the data, just the type)
 */
GT_INLINE gt_stats_vector* gt_stats_vector_new_from_template(gt_stats_vector* const stats_vector);
/*
 * Add another dimension to the vector (vector[@stats_vector][@nested_stats_vector])
 */
GT_INLINE void gt_stats_vector_add_nested(
    gt_stats_vector* const stats_vector,gt_stats_vector* const nested_stats_vector);

GT_INLINE void gt_stats_vector_clear(gt_stats_vector* const stats_vector);
GT_INLINE void gt_stats_vector_delete(gt_stats_vector* const stats_vector);

/*
 * Increment/Add bucket counter
 */
GT_INLINE void gt_stats_vector_inc(gt_stats_vector* const stats_vector,const uint64_t value);
GT_INLINE void gt_stats_vector_inc_nested_va(gt_stats_vector* const stats_vector,...);
GT_INLINE void gt_stats_vector_add(gt_stats_vector* const stats_vector,const uint64_t quantity,const uint64_t value);
GT_INLINE void gt_stats_vector_add_nested_va(gt_stats_vector* const stats_vector,const uint64_t quantity,...);

/*
 * Inverse. Given the stats_vector index returns the corresponding value/range.
 */
GT_INLINE uint64_t gt_stats_vector_get_value(gt_stats_vector* const stats_vector,const uint64_t index);
GT_INLINE void gt_stats_vector_get_value_range(
    gt_stats_vector* const stats_vector,const uint64_t index,
    uint64_t* const lo_value,uint64_t* const hi_value);

/*
 * Bucket counters getters (Individual buckets)
 */
GT_INLINE uint64_t gt_stats_vector_get_count(
    gt_stats_vector* const stats_vector,const uint64_t value);
GT_INLINE uint64_t gt_stats_vector_get_count_nested_va(
    gt_stats_vector* const stats_vector,const uint64_t value,...);

/*
 * Bucket counters getters (Accumulated ranges)
 */
GT_INLINE uint64_t gt_stats_vector_get_accumulated_count(
    gt_stats_vector* const stats_vector,const uint64_t value_from,const uint64_t value_to);
GT_INLINE uint64_t gt_stats_vector_get_accumulated_count_nested_va(
    gt_stats_vector* const stats_vector,const uint64_t value_from,const uint64_t value_to,...);

/*
 * Merge 2 stats-vector (adding bucket counting)
 */
GT_INLINE void gt_stats_vector_merge(gt_stats_vector* const stats_dst,gt_stats_vector* const stats_src);

/*
 * Display (Printers)
 */
GT_INLINE void gt_stats_vector_print_raw(gt_stats_vector* const stats_vector);
GT_INLINE void gt_stats_vector_print_json(gt_stats_vector* const stats_vector);
/*
 * Iterator
 */
GT_INLINE gt_stats_vector_iterator* gt_stats_vector_iterator_new(gt_stats_vector* const stats_vector);
GT_INLINE gt_stats_vector_iterator* gt_stats_vector_iterator_range_new(
    gt_stats_vector* const stats_vector,const uint64_t value_from,const uint64_t value_to);
GT_INLINE void gt_stats_vector_iterator_delete(gt_stats_vector_iterator* const stats_vector_iterator);

GT_INLINE bool gt_stats_vector_iterator_next(gt_stats_vector_iterator* const stats_vector_iterator);
GT_INLINE uint64_t gt_stats_vector_iterator_get_count(gt_stats_vector_iterator* const stats_vector_iterator);
GT_INLINE void gt_stats_vector_iterator_get_range(
    gt_stats_vector_iterator* const stats_vector_iterator,
    uint64_t* const lo_value,uint64_t* const hi_value);

#endif /* GT_STATS_VECTOR_H_ */
