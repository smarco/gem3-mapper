/*
 * PROJECT: GEM-Tools library
 * FILE: gt_vector.c
 * DATE: 01/02/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_vector.h"
#include "gt_error.h"
#include "gt_mm.h"

#define GT_VECTOR_EXPAND_FACTOR (3.0/2.0)

// FIXME: Keep errors here (we don't want to check this error val all over the code)
// FIXME: This shoud be a type, not a size
GT_INLINE gt_vector* gt_vector_new(size_t num_initial_elements,size_t element_size) {
  GT_ZERO_CHECK(element_size);
  gt_vector* vector=gt_alloc(gt_vector);
  vector->element_size=element_size;
  vector->elements_allocated=num_initial_elements;
  vector->memory=gt_malloc(num_initial_elements*element_size);
  vector->used=0;
  return vector;
}
GT_INLINE gt_status gt_vector_reserve(gt_vector* vector,size_t num_elements,bool zero_mem) {
  GT_VECTOR_CHECK(vector);
  if (vector->elements_allocated < num_elements) {
    size_t proposed=(float)vector->elements_allocated*GT_VECTOR_EXPAND_FACTOR;
    vector->elements_allocated=num_elements>proposed?num_elements:proposed;
    vector->memory=realloc(vector->memory,vector->elements_allocated*vector->element_size);
    if (!vector->memory) return GT_VECTOR_FAIL;
  }
  if (gt_expect_false(zero_mem)) {
    memset(vector->memory+vector->used*vector->element_size,0,
        (vector->elements_allocated-vector->used)*vector->element_size);
  }
  return GT_VECTOR_OK;
}
GT_INLINE gt_status gt_vector_resize__clear(gt_vector* vector,size_t num_elements) {
  GT_VECTOR_CHECK(vector);
  if (vector->elements_allocated < num_elements) {
    size_t proposed=(float)vector->elements_allocated*GT_VECTOR_EXPAND_FACTOR;
    vector->elements_allocated=num_elements>proposed?num_elements:proposed;
    gt_free(vector->memory);
    vector->memory=gt_malloc_nothrow(vector->elements_allocated,vector->element_size,0,0);
    if (!vector->memory) return GT_VECTOR_FAIL;
  }
  vector->used=0;
  return GT_VECTOR_OK;
}

GT_INLINE void gt_vector_cast__clear(gt_vector* vector,size_t element_size) {
  GT_VECTOR_CHECK(vector); GT_ZERO_CHECK(element_size);
  vector->elements_allocated=(vector->elements_allocated*vector->element_size)/element_size;
  vector->element_size=element_size;
  vector->used=0;
}
GT_INLINE void gt_vector_delete(gt_vector* vector) {
  GT_VECTOR_CHECK(vector);
  gt_free(vector->memory);
  gt_free(vector);
}
GT_INLINE void gt_vector_copy(gt_vector* vector_to,gt_vector* vector_from) {
  GT_VECTOR_CHECK(vector_to); GT_VECTOR_CHECK(vector_from);
  gt_vector_cast__clear(vector_to,vector_from->element_size);
  gt_vector_reserve(vector_to,vector_from->used,false);
  gt_vector_set_used(vector_to,gt_vector_get_used(vector_from));
  memcpy(vector_to->memory,vector_from->memory,vector_from->used*vector_from->element_size);
}
GT_INLINE gt_vector* gt_vector_dup(gt_vector* vector) {
  GT_VECTOR_CHECK(vector);
  gt_vector* const vector_cpy = gt_vector_new(vector->used,vector->element_size);
  gt_vector_set_used(vector_cpy,vector->used);
  memcpy(vector_cpy->memory,vector->memory,vector->used*vector->element_size);
  return vector_cpy;
}

GT_INLINE void* gt_vector_get_mem_element(gt_vector* vector,size_t position,size_t element_size) {
  GT_VECTOR_CHECK(vector); GT_ZERO_CHECK(element_size);
  GT_VECTOR_RANGE_CHECK(vector,position);
  return vector->memory+(position*element_size);
}
