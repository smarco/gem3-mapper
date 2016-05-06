/*
 * PROJECT: GEMMapper
 * FILE: vector.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Simple linear vector for generic type elements
 */

#include "utils/vector.h"
#include "system/errors.h"
#include "system/mm.h"

#define VECTOR_EXPAND_FACTOR (3.0/2.0)

/*
 * Setup
 */
vector_t* vector_new_(const uint64_t num_initial_elements,const uint64_t element_size) {
  GEM_CHECK_ZERO(element_size);
  vector_t* const vector_buffer=mm_alloc(vector_t);
  vector_buffer->element_size=element_size;
  vector_buffer->elements_allocated=num_initial_elements;
  vector_buffer->memory=mm_malloc_nothrow(1,num_initial_elements*element_size,false,0);
  gem_cond_fatal_error(!vector_buffer->memory,VECTOR_NEW,num_initial_elements*element_size);
  vector_buffer->used=0;
  return vector_buffer;
}
void vector_reserve(vector_t* const vector,const uint64_t num_elements,const bool zero_mem) {
  VECTOR_CHECK(vector);
  if (vector->elements_allocated < num_elements) {
    const uint64_t proposed=(float)vector->elements_allocated*VECTOR_EXPAND_FACTOR;
    vector->elements_allocated=num_elements>proposed?num_elements:proposed;
    vector->memory=mm_realloc_nothrow(vector->memory,vector->elements_allocated*vector->element_size);
    if (!vector->memory) {
      gem_cond_fatal_error(!vector->memory,VECTOR_RESERVE,vector->elements_allocated*vector->element_size);
    }
  }
  if (gem_expect_false(zero_mem)) {
    memset(vector->memory+vector->used*vector->element_size,0,(vector->elements_allocated-vector->used)*vector->element_size);
  }
}
void vector_resize__clear(vector_t* const vector,const uint64_t num_elements) {
  VECTOR_CHECK(vector);
  if (vector->elements_allocated < num_elements) {
    const uint64_t proposed=(float)vector->elements_allocated*VECTOR_EXPAND_FACTOR;
    vector->elements_allocated=(num_elements>proposed)?num_elements:proposed;
    // Free previous chunk (no need to pay the cost of reallocating memory)
    mm_free(vector->memory);
    // Allocate new block of memory
    vector->memory=mm_malloc_nothrow(vector->elements_allocated,vector->element_size,0,0);
    gem_cond_fatal_error(!vector->memory,VECTOR_RESERVE,vector->elements_allocated*vector->element_size);
  }
  vector->used=0;
}
void vector_cast__clear_(vector_t* const vector,const uint64_t element_size) {
  VECTOR_CHECK(vector);
  GEM_CHECK_ZERO(element_size);
  vector->elements_allocated=(vector->elements_allocated*vector->element_size)/element_size;
  vector->element_size=element_size;
  vector->used=0;
}
void vector_delete(vector_t* const vector) {
  VECTOR_CHECK(vector);
  mm_free(vector->memory);
  mm_free(vector);
}
/*
 * Accessors
 */
#ifdef GEM_DEBUG
void* vector_get_mem_element(vector_t* const vector,const uint64_t position,const uint64_t element_size) {
  VECTOR_CHECK(vector);
  GEM_CHECK_ZERO(element_size);
  if (position >= (vector)->used) {
    VECTOR_RANGE_CHECK(vector,position);
  }
  return vector->memory+(position*element_size);
}
#endif
/*
 * Miscellaneous
 */
void vector_copy(vector_t* const vector_to,vector_t* const vector_from) {
  VECTOR_CHECK(vector_to);
  VECTOR_CHECK(vector_from);
  // Prepare
  vector_cast__clear_(vector_to,vector_from->element_size);
  vector_reserve(vector_to,vector_from->used,false);
  // Copy
  vector_set_used(vector_to,vector_from->used);
  memcpy(vector_to->memory,vector_from->memory,vector_from->used*vector_from->element_size);
}
vector_t* vector_dup(vector_t* const vector_src) {
  VECTOR_CHECK(vector_src);
  vector_t* const vector_cpy = vector_new_(vector_src->used,vector_src->element_size);
  // Copy
  vector_set_used(vector_cpy,vector_src->used);
  memcpy(vector_cpy->memory,vector_src->memory,vector_src->used*vector_src->element_size);
  return vector_cpy;
}

