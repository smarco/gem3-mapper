/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
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
 * DESCRIPTION: Simple Hash Implementation (Interger Key, Generic Value)
 */

#ifndef IHASH_H_
#define IHASH_H_

#include "system/commons.h"
#include "system/mm_allocator.h"
#include "resources/uthash/uthash.h"

/*
 * Integer Key Hash
 */
typedef struct {
  int64_t key;
  void* element;
  UT_hash_handle hh;
} ihash_element_t;
typedef struct {
  ihash_element_t* head;
  mm_allocator_t* mm_allocator;
} ihash_t;
typedef struct {
  ihash_element_t* current;
  ihash_element_t* next;
} ihash_iterator_t;

/*
 * Constructor
 */
ihash_t* ihash_new(mm_allocator_t* const mm_allocator);
void ihash_clear(ihash_t* const ihash);
void ihash_delete(ihash_t* const ihash);

/*
 * Basic (Type-unsafe) Accessors
 */
ihash_element_t* ihash_get_ihash_element(ihash_t* const ihash,const int64_t key);
void ihash_insert_element(ihash_t* const ihash,const int64_t key,void* const element);
void ihash_remove_element(ihash_t* const ihash,const int64_t key);
void* ihash_get_element(ihash_t* const ihash,const int64_t key);

/*
 * Type-safe Accessors
 */
#define ihash_insert(ihash,key,element) ihash_insert_element(ihash,key,(void*)(element))
#define ihash_remove(ihash,key) ihash_remove_element(ihash,key)
#define ihash_get(ihash,key,type) ((type*)ihash_get_element(ihash,key))
bool ihash_is_contained(ihash_t* const ihash,const int64_t key);
uint64_t ihash_get_num_elements(ihash_t* const ihash);
uint64_t ihash_get_size(ihash_t* const ihash);

/*
 * Miscellaneous
 */
void ihash_sort_by_key(ihash_t* const ihash);

/*
 * Iterator
 */
#define IHASH_BEGIN_ITERATE(ihash,it_ikey,it_element,type) { \
  ihash_element_t *ihash_##ih_element, *ihash_##tmp; \
  HASH_ITER(hh,ihash->head,ihash_##ih_element,ihash_##tmp) { \
    type* const it_element = (type*)(ihash_##ih_element->element);
#define IHASH_END_ITERATE }}

#define IHASH_BEGIN_ITERATE__KEY(ihash,it_ikey,it_element,type) { \
  ihash_element_t *ihash_##ih_element, *ihash_##tmp; \
  HASH_ITER(hh,ihash->head,ihash_##ih_element,ihash_##tmp) { \
    type* const it_element = (type*)(ihash_##ih_element->element); \
    int64_t const it_ikey = ihash_##ih_element->key;
#define IHASH_END_ITERATE__KEY }}

ihash_iterator_t* ihash_iterator_new(ihash_t* const ihash);
void ihash_iterator_delete(ihash_iterator_t* const ihash_iterator);

bool ihash_iterator_eoi(ihash_iterator_t* const iterator);
bool ihash_iterator_next(ihash_iterator_t* const ihash_iterator);
int64_t ihash_iterator_get_key(ihash_iterator_t* const ihash_iterator);
void* ihash_iterator_get_element(ihash_iterator_t* const ihash_iterator);

#endif /* IHASH_H_ */
