/*
 * PROJECT: GEMMapper
 * FILE: interval_set.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "interval_set.h"


GEM_INLINE void interval_set_new(interval_set_t* const interval_set) {

}
GEM_INLINE void interval_set_clear(interval_set_t* const interval_set) {

}
GEM_INLINE void interval_set_delete(interval_set_t* const interval_set) {

}

/*
 * Counts the number of candidates in the result set of intervals
 * ** CHECKED[4/8/2011]
 */
#define COUNT_CANDIDATES(candidates, results) \
  candidates=0; \
  INTERVAL_ITERATE(results) { \
    candidates+=interval->hi-interval->lo; \
  } END_INTERVAL_ITERATE
#define COUNT_SET_INTERVALS_CANDIDATES(candidates, result_vector, init_int, end_int, max_misms) { \
  register uint64_t it; \
  register interval_t* result_interval = vector_get_mem(result_vector) + init_int; \
  candidates=0; \
  for (it=init_int; it<end_int; ++it, ++result_interval) { \
    if (result_interval->misms <= max_misms) { \
      candidates+=result_interval->hi-result_interval->lo; \
    } \
  } \
}

GEM_INLINE uint64_t interval_set_count_intervals(interval_set_t* const interval_set) {
  // TODO
  return 0;
}
GEM_INLINE uint64_t interval_set_count_intervals_length(interval_set_t* const interval_set) {
  // TODO
  return 0;
}
GEM_INLINE uint64_t interval_set_count_intervals_length_thresholded(
    interval_set_t* const interval_set,const uint64_t max_error) {

  // TODO
  return 0;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////


//
//
///*
// * Appends to @exclusion_set the intervals contained into @result_set (union set)
// */
//GEM_INLINE void interval_set_union(vector* const exclusion_set, vector* const result_set) {
//  register const uint64_t new_size_excl = vector_get_used(exclusion_set) + vector_get_used(result_set);
//  vector_reserve(exclusion_set, new_size_excl);
//  register interval_t* int_excl = vector_get_mem_next_elm(exclusion_set,interval_t);
//  INTERVAL_ITERATE(result_set) {
//    int_excl->lo = interval->lo;
//    int_excl->hi = interval->hi;
//    int_excl->misms = interval->misms;
//    ++int_excl;
//  } END_INTERVAL_ITERATE;
//  vector_set_used(exclusion_set, new_size_excl);
//}
///*
// * Subtracts to @result_set the intervals contained in @exclusion_set (difference set)
// */
//GEM_INLINE void interval_set_subtract(vector* const result_set, vector* const exclusion_set) {
//  register const uint64_t exclusion_set_size = vector_get_used(exclusion_set);
//  register uint64_t result_set_size = vector_get_used(result_set);
//  register uint64_t i, j;
//  for (i=0; i<result_set_size; ++i) {
//    register interval_t* int_res = (interval_t*)vector_get_mem(result_set) + i;
//    register interval_t* int_excl = (interval_t*)vector_get_mem(exclusion_set);
//    for (j=0; j<exclusion_set_size; ++j, ++int_excl) {
//      register const idx_t hi1 = int_res->hi;
//      register const idx_t lo1 = int_res->lo;
//      register const idx_t hi2 = int_excl->hi;
//      register const idx_t lo2 = int_excl->lo;
//      if (hi1 <= lo2 || hi2 <= lo1) { /* Disjoint intervals */
//        continue;
//      } else {
//        if (lo2 <= lo1 && hi1 <= hi2) { /* Full subtraction */
//          int_res->lo = int_res->hi; // We close the interval
//        } else if (lo1 < lo2 && hi2 < hi1) { /* Exclusion inside result */
//          register idx_t misms1 = int_res->misms;
//          vector_reserve(result_set, result_set_size+1);
//          int_res = (interval_t*)vector_get_mem(result_set) + result_set_size;
//          int_res->lo=hi2; int_res->hi=hi1; int_res->misms=misms1;
//          int_res = (interval_t*)vector_get_mem(result_set) + i;
//          /* int_res->lo=lo1; */ int_res->hi=lo2;
//          ++result_set_size;
//        } else if (lo2 == lo1 && lo1 < hi2) { /* Exclusion overlaps left side of result */
//          int_res->lo=hi2;
//        } else /* if (lo2 < hi1 && hi1 == hi2) */ { /* Exclusion overlaps right side of result */
//          int_res->hi=lo2;
//        }
//      }
//    }
//  }
//  vector_set_used(result_set,result_set_size);
//}
//
//
//
//
//
//
//
//















