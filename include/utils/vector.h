/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *  Copyright (c) 2011-2017 by Paolo Ribeca  <paolo.ribeca@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 *            Paolo Ribeca <paolo.ribeca@gmail.com>
 * DESCRIPTION: Simple linear vector for generic type elements
 */

#ifndef VECTOR_H_
#define VECTOR_H_

#include "system/commons.h"

/*
 * Checkers
 */
#define VECTOR_CHECK(vector) \
  GEM_CHECK_NULL(vector); \
  GEM_CHECK_NULL(vector->memory); \
  GEM_CHECK_ZERO(vector->element_size)
#define VECTOR_RANGE_CHECK(vector,position) \
  VECTOR_CHECK(vector); \
  gem_fatal_check(position>=(vector)->used,POSITION_OUT_OF_RANGE, \
      (uint64_t)position,(uint64_t)0,((int64_t)(vector)->used)-1);
/*
 * Data Structures
 */
typedef struct {
  void* memory;
  uint64_t used;
  uint64_t element_size;
  uint64_t elements_allocated;
} vector_t;

/*
 * Vector Setup (Initialization & Allocation)
 */
#define vector_new(num_initial_elements,type) vector_new_(num_initial_elements,sizeof(type))
vector_t* vector_new_(const uint64_t num_initial_elements,const uint64_t element_size);
void vector_reserve(vector_t* const vector,const uint64_t num_elements,const bool zero_mem);
void vector_resize__clear(vector_t* const vector,const uint64_t num_elements);
#define vector_cast__clear(vector,type) vector_cast__clear_s(vector,sizeof(type))
void vector_cast__clear_(vector_t* const vector,const uint64_t element_size);
#define vector_clear(vector) (vector)->used=0   // TODO Implement a REAP mechanism
void vector_delete(vector_t* const vector);
#define vector_is_empty(vector) (vector_get_used(vector)==0)
#define vector_reserve_additional(vector,additional) vector_reserve(vector,vector_get_used(vector)+additional,false)
#define vector_prepare(vector,num_elements,type) \
  vector_cast__clear(vector,sizeof(type)); \
  vector_reserve(vector,num_elements,false);

/*
 * Element Getters/Setters
 */
#define vector_get_mem(vector,type) ((type*)((vector)->memory))
#define vector_get_last_elm(vector,type) (vector_get_mem(vector,type)+(vector)->used-1)
#define vector_get_free_elm(vector,type) (vector_get_mem(vector,type)+(vector)->used)
#ifndef GEM_DEBUG
  #define vector_get_elm(vector,position,type) (vector_get_mem(vector,type)+position)
#else
  void* vector_get_mem_element(vector_t* const vector,const uint64_t position,const uint64_t element_size);
  #define vector_get_elm(vector,position,type) ((type*)vector_get_mem_element(vector,position,sizeof(type)))
#endif
#define vector_set_elm(vector,position,type,elm) *vector_get_elm(vector,position,type) = elm

/*
 * Used elements Getters/Setters
 */
#define vector_get_used(vector) ((vector)->used)
#define vector_set_used(vector,total_used) (vector)->used=(total_used)
#define vector_inc_used(vector) (++((vector)->used))
#define vector_dec_used(vector) (--((vector)->used))
#define vector_add_used(vector,additional) vector_set_used(vector,vector_get_used(vector)+additional)
#define vector_update_used(vector,pointer_to_next_free_element) \
  (vector)->used = (pointer_to_next_free_element) - ((__typeof__(pointer_to_next_free_element))((vector)->memory))


/*
 * Vector Allocate/Insert (Get a new element or Add an element to the end of the vector)
 */
#define vector_alloc_new(vector,type,return_element_pointer) { \
  vector_reserve_additional(vector,1); \
  return_element_pointer = vector_get_free_elm(vector,type); \
  vector_inc_used(vector); \
}
#define vector_insert(vector,element,type) { \
  vector_reserve_additional(vector,1); \
  *(vector_get_free_elm(vector,type))=element; \
  vector_inc_used(vector); \
}

/*
 * Macro generic iterator
 *   VECTOR_ITERATE(vector_of_ints,elm_iterator,elm_counter,int) {
 *     ..code..
 *   }
 */
#define VECTOR_ITERATE(vector,element,counter,type) \
  const uint64_t vector_##element##_used = vector_get_used(vector); \
  type* element = vector_get_mem(vector,type); \
  uint64_t counter; \
  for (counter=0;counter<vector_##element##_used;++element,++counter)
#define VECTOR_ITERATE_OFFSET(vector,element,counter,offset,type) \
  const uint64_t vector_##element##_used = vector_get_used(vector); \
  type* element = vector_get_mem(vector,type)+offset; \
  uint64_t counter; \
  for (counter=offset;counter<vector_##element##_used;++counter,++element)
#define VECTOR_ITERATE_CONST(vector,element,counter,type) \
  const uint64_t vector_##element##_used = vector_get_used(vector); \
  const type* element = vector_get_mem(vector,type); \
  uint64_t counter; \
  for (counter=0;counter<vector_##element##_used;++element,++counter)
#define VECTOR_ITERATE_ELASTIC(vector,element,counter,type) \
  type* element = vector_get_mem(vector,type); \
  uint64_t counter; \
  for (counter=0;counter<vector_get_used(vector);++element,++counter)

/*
 * Miscellaneous
 */
void vector_copy(vector_t* const vector_to,vector_t* const vector_from);
vector_t* vector_dup(vector_t* const vector_src);

/*
 * Error Messages
 */
#define GEM_ERROR_VECTOR_NEW "Could not create new vector (%"PRIu64" bytes requested)"
#define GEM_ERROR_VECTOR_RESERVE "Could not reserve vector (%"PRIu64" bytes requested)"

#endif /* VECTOR_H_ */
