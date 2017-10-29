/*
 * PROJECT: GEM-Tools library
 * FILE: gt_vector.h
 * DATE: 01/02/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef _GT_VECTOR_H_GUARD_
#define _GT_VECTOR_H_GUARD_

#include "gt_commons.h"
#include "gt_error.h"

// Codes gt_status
#define GT_VECTOR_OK 0
#define GT_VECTOR_FAIL 1

/*
 * Checkers
 */
#define GT_VECTOR_CHECK(vector) \
  GT_NULL_CHECK(vector); \
  GT_NULL_CHECK(vector->memory); \
  GT_ZERO_CHECK(vector->element_size)
#define GT_VECTOR_RANGE_CHECK(vector,position) \
  GT_VECTOR_CHECK(vector); \
  gt_fatal_check(position>=(vector)->used||position<0,POSITION_OUT_OF_RANGE, \
      (uint64_t)position,(uint64_t)0,((int64_t)(vector)->used)-1);

// Simple linear vector for generic type elements
typedef struct {
  void* memory;
  size_t used;
  size_t element_size;
  size_t elements_allocated;
} gt_vector;

// Get the content of the vector
#define gt_vector_get_mem(vector,type) ((type*)((vector)->memory))
#define gt_vector_get_last_elm(vector,type) (gt_vector_get_mem(vector,type)+(vector)->used-1)
#define gt_vector_get_free_elm(vector,type) (gt_vector_get_mem(vector,type)+(vector)->used)
// Getters/Setter used elements in the vector
#define gt_vector_get_used(vector) ((vector)->used)
#define gt_vector_set_used(vector,total_used) (vector)->used=(total_used)
#define gt_vector_inc_used(vector) (++((vector)->used))
#define gt_vector_dec_used(vector) (--((vector)->used))
#define gt_vector_add_used(vector,additional) gt_vector_set_used(vector,gt_vector_get_used(vector)+additional)
#define gt_vector_clear(vector) (vector)->used=0
#define gt_vector_is_empty(vector) (gt_vector_get_used(vector)==0)
// Initialization and allocation
#define gt_vector_reserve_additional(vector,additional) gt_vector_reserve(vector,gt_vector_get_used(vector)+additional,false)
#define gt_vector_prepare(vector,data_type,num_elements) \
  gt_vector_cast__clear(vector,sizeof(data_type)); \
  gt_vector_reserve(vector,num_elements,false)
// Macro generic iterator
//  GT_VECTOR_ITERATE(vector_of_ints,elm_iterator,elm_counter,int) {
//    ..code..
//  }
#define GT_VECTOR_ITERATE(vector,element,counter,type) \
  uint64_t counter; \
  type* element = gt_vector_get_mem(vector,type); \
  for (counter=0;counter<gt_vector_get_used(vector);++element,++counter)
#define GT_VECTOR_ITERATE_OFFSET(vector,element,counter,offset,type) \
  uint64_t counter; \
  type* element = gt_vector_get_elm(vector,offset,type); \
  for (counter=offset;counter<gt_vector_get_used(vector);++element,++counter)
// Add element to the vector (at the end)
#define gt_vector_insert(vector,element,type) { \
  gt_vector_reserve_additional(vector,1); \
  *(gt_vector_get_free_elm(vector,type))=element; \
  gt_vector_inc_used(vector); \
}

#ifdef GT_NO_CONSISTENCY_CHECKS
#define gt_vector_get_elm(vector,position,type) (gt_vector_get_mem(vector,type)+position)
#else
#define gt_vector_get_elm(vector,position,type) ((type*)gt_vector_get_mem_element(vector,position,sizeof(type)))
#endif

#define gt_vector_set_elm(vector,position,type,elm) (*gt_vector_get_elm(vector,position,type) = elm)

GT_INLINE gt_vector* gt_vector_new(size_t num_initial_elements,size_t element_size); // FIXME: wrt to type MACRO
GT_INLINE gt_status gt_vector_reserve(gt_vector* vector,size_t num_elements,bool zero_mem);
GT_INLINE gt_status gt_vector_resize__clear(gt_vector* vector,size_t num_elements);

GT_INLINE void gt_vector_cast__clear(gt_vector* vector,size_t element_size);
GT_INLINE void gt_vector_delete(gt_vector* vector);
GT_INLINE void gt_vector_copy(gt_vector* vector_to,gt_vector* vector_from);
GT_INLINE gt_vector* gt_vector_dup(gt_vector* vector);

GT_INLINE void* gt_vector_get_mem_element(gt_vector* vector,size_t position,size_t element_size);

#endif /* _GT_VECTOR_H_GUARD_ */
