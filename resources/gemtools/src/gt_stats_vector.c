/*
 * PROJECT: GEM-Tools library
 * FILE: gt_stats_vector.c
 * DATE: 10/12/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODo
 */

#include "gt_stats_vector.h"

/*
 * Constructors
 */
GT_INLINE gt_stats_vector* gt_stats_vector_customed_range_new(
    const uint64_t* const customed_range_values,const uint64_t num_values,
    const uint64_t out_of_range_bucket_size) {
  GT_NULL_CHECK(customed_range_values);
  gt_cond_fatal_error(num_values>=2,INVALID_VALUE,"'num_values'",">2");
  GT_ZERO_CHECK(out_of_range_bucket_size);
  // Allocate handler
  gt_stats_vector* const stats_vector = gt_alloc(gt_stats_vector);
  // Init
  stats_vector->type = GT_STATS_VECTOR_CUSTOMED_RANGE;
  stats_vector->counters = gt_calloc(num_values,uint64_t,true);
  // Customed range vetor
  GT_NULL_CHECK(customed_range_values);
  stats_vector->customed_range_values = customed_range_values;
  stats_vector->num_values = num_values-1;
  stats_vector->out_of_range_bucket_size = out_of_range_bucket_size;
  stats_vector->out_values = gt_ihash_new();
  // Nested
  stats_vector->template_vector = NULL;
  stats_vector->nested_vectors = NULL;
  return stats_vector;
}
GT_INLINE gt_stats_vector* gt_stats_vector_step_range_new(
    const uint64_t min_value,const uint64_t max_value,const uint64_t step,
    const uint64_t out_of_range_bucket_size) {
  GT_ZERO_CHECK(step);
  gt_cond_fatal_error(min_value>max_value,VSTATS_INVALID_MIN_MAX);
  GT_ZERO_CHECK(out_of_range_bucket_size);
  const uint64_t range = max_value-min_value+1;
  const uint64_t num_values = (range+(step-1))/step;
  // Allocate handler
  gt_stats_vector* const stats_vector = gt_alloc(gt_stats_vector);
  // Init
  stats_vector->type = GT_STATS_VECTOR_STEP_RANGE;
  stats_vector->counters = gt_calloc(num_values,uint64_t,true);
  // Step Range
  stats_vector->num_values = num_values;
  stats_vector->min_value = min_value;
  stats_vector->max_value = max_value;
  stats_vector->step = step;
  stats_vector->out_of_range_bucket_size = out_of_range_bucket_size;
  stats_vector->out_values = gt_ihash_new();
  // Nested
  stats_vector->template_vector = NULL;
  stats_vector->nested_vectors = NULL;
  return stats_vector;
}
GT_INLINE gt_stats_vector* gt_stats_vector_raw_new(
    const uint64_t num_values,const uint64_t out_of_range_bucket_size) {
  GT_ZERO_CHECK(out_of_range_bucket_size);
  // Allocate handler
  gt_stats_vector* const stats_vector = gt_alloc(gt_stats_vector);
  // Init
  stats_vector->type = GT_STATS_VECTOR_RAW;
  stats_vector->counters = (num_values) ? gt_calloc(num_values,uint64_t,true) : NULL;
  // Raw
  stats_vector->num_values = num_values;
  stats_vector->out_of_range_bucket_size = out_of_range_bucket_size;
  stats_vector->out_values = gt_ihash_new();
  // Nested
  stats_vector->template_vector = NULL;
  stats_vector->nested_vectors = NULL;
  return stats_vector;
}
GT_INLINE gt_stats_vector* gt_stats_vector_new_from_template(gt_stats_vector* const stats_vector) {
  // Allocate handler
  gt_stats_vector* const stats_vector_copy = gt_alloc(gt_stats_vector);
  // Copy template
  stats_vector_copy->type = stats_vector->type;
  stats_vector_copy->min_value = stats_vector->min_value;
  stats_vector_copy->max_value = stats_vector->max_value;
  stats_vector_copy->step = stats_vector->step;
  stats_vector_copy->customed_range_values = stats_vector->customed_range_values;
  stats_vector_copy->num_values = stats_vector->num_values;
  stats_vector_copy->out_of_range_bucket_size = stats_vector->out_of_range_bucket_size;
  // Init
  stats_vector_copy->counters = gt_calloc(stats_vector->num_values,uint64_t,true);
  stats_vector_copy->out_values = gt_ihash_new();
  // Nested (No copy)
  stats_vector_copy->template_vector = NULL;
  stats_vector_copy->nested_vectors = NULL;
  return stats_vector_copy;
}
GT_INLINE void gt_stats_vector_add_nested(
    gt_stats_vector* const stats_vector,gt_stats_vector* const nested_stats_vector) {
  GT_STATS_VECTOR_CHECK(gt_stats_vector);
  stats_vector->template_vector = nested_stats_vector;
  const uint64_t size_of_uint64_t = sizeof(uint64_t);
  const uint64_t size_of_gt_stats_vector_ptr = sizeof(gt_stats_vector*);
  // Allocate nested vectors container
  if (size_of_uint64_t==size_of_gt_stats_vector_ptr) {
    stats_vector->nested_vectors = stats_vector->counters;
  } else {
    gt_free(stats_vector->counters);
    stats_vector->nested_vectors = gt_calloc(num_values,gt_stats_vector*,true);
    stats_vector->counters = stats_vector->nested_vectors;
  }
}
GT_INLINE void gt_stats_vector_clear(gt_stats_vector* const stats_vector) {
  GT_STATS_VECTOR_CHECK(gt_stats_vector);
  memset(stats_vector->counters,0,stats_vector->num_values);
  gt_ihash_clear(stats_vector->out_values,true);
}
GT_INLINE void gt_stats_vector_delete(gt_stats_vector* const stats_vector) {
  GT_STATS_VECTOR_CHECK(gt_stats_vector);
  gt_free(stats_vector->counters);
  gt_ihash_delete(stats_vector->out_values);
  if (stats_vector->template_vector) gt_stats_vector_delete(stats_vector->template_vector);
  gt_free(stats_vector);
}
/*
 * Calculate index
 */
GT_INLINE uint64_t gt_stats_cvector_get_index(gt_stats_vector* const stats_vector,const uint64_t value) {
  GT_STATS_VECTOR_CHECK(gt_stats_vector);
  uint64_t* const range_values = stats_vector->customed_range_values;
  uint64_t lo = 0;
  uint64_t hi = stats_vector->num_values;
  if (value < range_values[lo] || range_values[hi] < value) {
    return GT_STATS_VECTOR_OUT_OF_RANGE;
  } else {
    while (hi-lo > 1) { // Search for the bucket
      const uint64_t mid = (hi+lo)/2;
      if (value < range_values[mid]) {
        hi = mid;
      } else {
        lo = mid;
      }
    }
    gt_cond_fatal_error(lo+1==hi && range_values[lo] <= value && value < range_values[hi],ALG_INCONSISNTENCY);
    return lo;
  }
}
GT_INLINE uint64_t gt_stats_svector_get_index(gt_stats_vector* const stats_vector,const uint64_t value) {
  GT_STATS_VECTOR_CHECK(gt_stats_vector);
  if (value < stats_vector->min_value || stats_vector->max_value < value) {
    return GT_STATS_VECTOR_OUT_OF_RANGE;
  } else {
    return (value-stats_vector->min_value)/stats_vector->step;
  }
}
GT_INLINE uint64_t gt_stats_rvector_get_index(gt_stats_vector* const stats_vector,const uint64_t value) {
  GT_STATS_VECTOR_CHECK(gt_stats_vector);
  if (value >= stats_vector->num_values) {
    return GT_STATS_VECTOR_OUT_OF_RANGE;
  } else {
    return value;
  }
}
/*
 * Vector's Buckets getters
 */
GT_INLINE uint64_t* gt_stats_hvector_get_counter(gt_stats_vector* const stats_vector,const uint64_t value) {
  GT_STATS_VECTOR_CHECK(gt_stats_vector);
  const uint64_t bucket_index = value/stats_vector->out_of_range_bucket_size;
  // Fetch counter
  uint64_t* counter = gt_ihash_get_element(stats_vector->out_values,bucket_index);
  if (counter!=NULL) return counter;
  // Allocate new counter
  counter = gt_alloc(uint64_t);
  *counter = 0;
  gt_ihash_insert(stats_vector->out_values,bucket_index,counter,uint64_t);
  return counter;
}
GT_INLINE uint64_t* gt_stats_vector_get_counter(gt_stats_vector* const stats_vector,const uint64_t value) {
  GT_STATS_VECTOR_CHECK(gt_stats_vector);
  uint64_t bucket_index;
  switch (stats_vector->type) {
    case GT_STATS_VECTOR_CUSTOMED_RANGE:
      bucket_index = gt_stats_cvector_get_index(stats_vector,value);
      break;
    case GT_STATS_VECTOR_STEP_RANGE:
      bucket_index = gt_stats_svector_get_index(stats_vector,value);
      break;
    case GT_STATS_VECTOR_RAW:
      bucket_index = gt_stats_rvector_get_index(stats_vector,value);
      break;
    default:
      GT_INVALID_CASE();
      break;
  }
  // Return counter
  if (bucket_index==GT_STATS_VECTOR_OUT_OF_RANGE) {
    return gt_stats_hvector_get_counter(stats_vector,value);
  } else {
    return stats_vector->counters+bucket_index;
  }
}
/*
 * Nested vector getters
 */
GT_INLINE gt_stats_vector* gt_stats_hvector_get_nested(gt_stats_vector* const stats_vector,const uint64_t value) {
  GT_STATS_VECTOR_CHECK(gt_stats_vector);
  const uint64_t bucket_value = value/stats_vector->out_of_range_bucket_size;
  // Fetch counter
  gt_stats_vector* nested_vector = gt_ihash_get_element(stats_vector->out_values,bucket_value);
  if (nested_vector!=NULL) return nested_vector;
  // Allocate new counter
  nested_vector = gt_stats_vector_copy_template(stats_vector->template_vector);
  gt_ihash_insert(stats_vector->out_values,bucket_value,counter,uint64_t);
  return nested_vector;
}
GT_INLINE gt_stats_vector* gt_stats_vector_get_nested(gt_stats_vector* const stats_vector,const uint64_t value) {
  GT_STATS_VECTOR_CHECK(gt_stats_vector);
  uint64_t bucket_index;
  switch (stats_vector->type) {
    case GT_STATS_VECTOR_CUSTOMED_RANGE:
      bucket_index = gt_stats_cvector_get_index(stats_vector,value);
      break;
    case GT_STATS_VECTOR_STEP_RANGE:
      bucket_index = gt_stats_svector_get_index(stats_vector,value);
      break;
    case GT_STATS_VECTOR_RAW:
      bucket_index = gt_stats_rvector_get_index(stats_vector,value);
      break;
    default:
      GT_INVALID_CASE();
      break;
  }
  // Return counter
  if (bucket_index==GT_STATS_VECTOR_OUT_OF_RANGE) {
    return gt_stats_hvector_get_nested(stats_vector,value);
  } else {
    return stats_vector->nested_vectors+bucket_index;
  }
}
GT_INLINE void gt_stats_vector_get_nested_v(
    gt_stats_vector* const stats_vector,
    gt_stats_vector** const nested_vector,uint64_t* const nested_vector_value,va_list v_args) {
  GT_STATS_VECTOR_CHECK(gt_stats_vector);
  // Fetch innermost nested vector
  gt_stats_vector* const nested = stats_vector;
  do {
    uint64_t const value = va_arg(v_args,uint64_t);
    if (nested->nested_vectors == NULL) {
      *nested_vector = nested;
      *nested_vector_value = value;
      return;
    } else {
      nested = gt_stats_vector_get_nested(nested,value);
    }
  } while (true);
}
GT_INLINE void gt_stats_vector_get_nested_va(
    gt_stats_vector* const stats_vector,
    gt_stats_vector** const nested_vector,uint64_t* const nested_vector_value,...) {
  // Init va
  va_list v_args;
  va_start(v_args,nested_vector_value);
  // Fetch innermost nested vector
  gt_stats_vector_get_nested_v(stats_vector,nested_vector,nested_vector_value,v_args);
  // End va
  va_end(v_args);
}
/*
 * Increment/Add bucket counter
 */
GT_INLINE void gt_stats_vector_inc_nested_va(gt_stats_vector* const stats_vector,...) {
  GT_STATS_VECTOR_CHECK(gt_stats_vector);
  // Init va
  va_list v_args;
  va_start(v_args,value);
  // Fetch innermost nested vector & the last value from va
  gt_stats_vector* nested_vector = NULL;
  uint64_t nested_vector_value = 0;
  gt_stats_vector_get_nested_va(stats_vector,&nested_vector,&nested_vector_value);
  // Inc
  GT_STATS_VECTOR_CHECK(nested_vector);
  gt_stats_vector_inc(nested_vector,nested_vector_value);
  // End va
  va_end(v_args);
}
GT_INLINE void gt_stats_vector_inc(gt_stats_vector* const stats_vector,const uint64_t value) {
  GT_STATS_VECTOR_CHECK(gt_stats_vector);
  ++(*gt_stats_vector_get_counter(stats_vector,value));
}
GT_INLINE void gt_stats_vector_add_nested_va(gt_stats_vector* const stats_vector,const uint64_t quantity,...) {
  GT_STATS_VECTOR_CHECK(gt_stats_vector);
  // Init va
  va_list v_args;
  va_start(v_args,value);
  // Fetch innermost nested vector & the last value from va
  gt_stats_vector* nested_vector = NULL;
  uint64_t nested_vector_value = 0;
  gt_stats_vector_get_nested_va(stats_vector,&nested_vector,&nested_vector_value);
  // Inc
  GT_STATS_VECTOR_CHECK(nested_vector);
  gt_stats_vector_add(nested_vector,nested_vector_value);
  // End va
  va_end(v_args);
}
GT_INLINE void gt_stats_vector_add(gt_stats_vector* const stats_vector,const uint64_t quantity,const uint64_t value) {
  GT_STATS_VECTOR_CHECK(gt_stats_vector);
  *gt_stats_vector_get_counter(stats_vector,value) += quantity;
}
/*
 * Inverse. Given the stats_vector index returns the corresponding value/range.
 */
GT_INLINE uint64_t gt_stats_vector_get_value(gt_stats_vector* const stats_vector,const uint64_t index) {
  GT_STATS_VECTOR_CHECK(gt_stats_vector);
  if (index < stats_vector->num_values) {
    switch (stats_vector->type) {
      case GT_STATS_VECTOR_CUSTOMED_RANGE:
        return stats_vector->customed_range_values[index];
        break;
      case GT_STATS_VECTOR_STEP_RANGE:
        return stats_vector->min_value + index*stats_vector->step;
        break;
      case GT_STATS_VECTOR_RAW:
        return index;
        break;
      default:
        GT_INVALID_CASE();
        break;
    }
  } else {
    // Hash value
    return index*stats_vector->out_of_range_bucket_size;
  }
}
GT_INLINE void gt_stats_vector_get_value_range(
    gt_stats_vector* const stats_vector,const uint64_t index,
    uint64_t* const lo_value,uint64_t* const hi_value) {
  GT_STATS_VECTOR_CHECK(gt_stats_vector);
  if (index < stats_vector->num_values) {
    switch (stats_vector->type) {
      case GT_STATS_VECTOR_CUSTOMED_RANGE:
        *lo_value = stats_vector->customed_range_values[index];
        *hi_value = stats_vector->customed_range_values[index+1];
        break;
      case GT_STATS_VECTOR_STEP_RANGE:
        *lo_value = stats_vector->min_value + index*stats_vector->step;
        *hi_value = *lo_value + stats_vector->step;
        break;
      case GT_STATS_VECTOR_RAW:
        *lo_value = index;
        *hi_value = index+1;
        break;
      default:
        GT_INVALID_CASE();
        break;
    }
  } else {
    // Hash value
    *lo_value = index*stats_vector->out_of_range_bucket_size;
    *hi_value = *lo_value+stats_vector->out_of_range_bucket_size;
  }
}
/*
 * Bucket counters getters (Individual buckets & Accumulated ranges)
 */
GT_INLINE uint64_t gt_stats_vector_get_count(gt_stats_vector* const stats_vector,const uint64_t value) {
  GT_STATS_VECTOR_CHECK(gt_stats_vector);
  return *gt_stats_vector_get_counter(stats_vector,value);
}
GT_INLINE uint64_t gt_stats_vector_get_count_nested_va(
    gt_stats_vector* const stats_vector,const uint64_t value,...) {
  GT_STATS_VECTOR_CHECK(gt_stats_vector);
  // Init va
  va_list v_args;
  va_start(v_args,value);
  // Fetch innermost nested vector & the last value from va
  gt_stats_vector* nested_vector = NULL;
  uint64_t nested_vector_value = 0;
  gt_stats_vector_get_nested_va(stats_vector,&nested_vector,&nested_vector_value);
  // Inc
  GT_STATS_VECTOR_CHECK(nested_vector);
  uint64_t counter = gt_stats_vector_get_count(nested_vector,nested_vector_value);
  // End va
  va_end(v_args);
  // Return
  return counter;
}
/*
 * Bucket counters getters (Accumulated ranges)
 */
GT_INLINE uint64_t gt_stats_vector_get_accumulated_count(
    gt_stats_vector* const stats_vector,const uint64_t value_from,const uint64_t value_to) {
  GT_STATS_VECTOR_CHECK(gt_stats_vector);
  // TODO
}
GT_INLINE uint64_t gt_stats_vector_get_accumulated_count_nested_va(
    gt_stats_vector* const stats_vector,const uint64_t value_from,const uint64_t value_to,...) {
  GT_STATS_VECTOR_CHECK(gt_stats_vector);
  // TODO
}
/*
 * Merge 2 stats-vector (adding bucket counting)
 */
GT_INLINE void gt_stats_vector_merge(gt_stats_vector* const stats_dst,gt_stats_vector* const stats_src) {
  // TODO
}
/*
 * Display (Printers)
 */
GT_INLINE void gt_stats_vector_print_raw(gt_stats_vector* const stats_vector) {
  // TODO
}
GT_INLINE void gt_stats_vector_print_json(gt_stats_vector* const stats_vector) {
  // TODO
}
/*
 * Iterator
 */
#define gt_string_cmp_wrapper(arg1,arg2) gt_string_cmp((char*)arg1,(char*)arg2)

GT_INLINE gt_stats_vector_iterator* gt_stats_vector_iterator_new(gt_stats_vector* const stats_vector) {
  GT_STATS_VECTOR_CHECK(gt_stats_vector);
  // Allocate
  gt_stats_vector_iterator* const stats_vector_iterator = gt_alloc(gt_stats_vector_iterator);
  // Init
  stats_vector_iterator->stats_vector = stats_vector;
  stats_vector_iterator->start_index = 0;
  stats_vector_iterator->end_index = UINT64_MAX;
  stats_vector_iterator->eoi = (stats_vector->num_values>0 || gt_ihash_get_num_elements(stats_vector->out_values)>0);
  // Array iteration
  stats_vector_iterator->is_index_in_array = (stats_vector->num_values>0);
  stats_vector_iterator->array_index = 0;
  // Hash iteration
  gt_ihash_sort_by_key(stats_vector->out_values);
  stats_vector_iterator->ihash_iterator = gt_ihash_iterator_new(stats_vector->out_values);
}
GT_INLINE gt_stats_vector_iterator* gt_stats_vector_iterator_range_new(
    gt_stats_vector* const stats_vector,const uint64_t value_from,const uint64_t value_to) {
  GT_STATS_VECTOR_CHECK(gt_stats_vector);
  // Allocate
  gt_stats_vector_iterator* const stats_vector_iterator = gt_alloc(gt_stats_vector_iterator);
  // Init
  stats_vector_iterator->stats_vector = stats_vector;
  const uint64_t idx_from = ;
  const uint64_t idx_to = ;
  stats_vector_iterator->start_index = 0;
  stats_vector_iterator->end_index = UINT64_MAX;
  stats_vector_iterator->eoi = (stats_vector->num_values>0 || gt_ihash_get_num_elements(stats_vector->out_values)>0);
  // Array iteration
  stats_vector_iterator->is_index_in_array = (stats_vector->num_values>0);
  stats_vector_iterator->array_index = 0;
  // Hash iteration
  gt_ihash_sort_by_key(stats_vector->out_values);
  stats_vector_iterator->ihash_iterator = gt_ihash_iterator_new(stats_vector->out_values);
}
GT_INLINE void gt_stats_vector_iterator_delete(gt_stats_vector_iterator* const stats_vector_iterator) {
  GT_NULL_CHECK(stats_vector_iterator);
  GT_STATS_VECTOR_CHECK(stats_vector_iterator->stats_vector);

}
GT_INLINE bool gt_stats_vector_iterator_next(gt_stats_vector_iterator* const stats_vector_iterator) {
  GT_NULL_CHECK(stats_vector_iterator);
  GT_STATS_VECTOR_CHECK(stats_vector_iterator->stats_vector);
}
GT_INLINE uint64_t gt_stats_vector_iterator_get_count(gt_stats_vector_iterator* const stats_vector_iterator) {
  GT_NULL_CHECK(stats_vector_iterator);
  GT_STATS_VECTOR_CHECK(stats_vector_iterator->stats_vector);
}
GT_INLINE void gt_stats_vector_iterator_get_range(
    gt_stats_vector_iterator* const stats_vector_iterator,
    uint64_t* const lo_value,uint64_t* const hi_value) {
  GT_NULL_CHECK(stats_vector_iterator);
  GT_STATS_VECTOR_CHECK(stats_vector_iterator->stats_vector);
}


