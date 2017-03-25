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

#include "utils/hash.h"
#include "system/mm.h"

/*
 * Constants
 */
#define IHASH_SIZE_PER_ELEMENT 56

/*
 * Allocators
 */
#undef  uthash_malloc
#define uthash_malloc(sz)   (ihash->mm_allocator!=NULL ? mm_allocator_malloc(ihash->mm_allocator,sz) : malloc(sz))
#undef  uthash_free
#define uthash_free(ptr,sz) if (ihash->mm_allocator==NULL) free(ptr);

/*
 * Constructor
 */
ihash_t* ihash_new(mm_allocator_t* const mm_allocator) {
  ihash_t* const ihash = mm_alloc(ihash_t);
  ihash->head = NULL; // uthash initializer
  ihash->mm_allocator = mm_allocator;
  return ihash;
}
void ihash_clear(ihash_t* const ihash) {
  ihash_element_t *ihash_element, *tmp;
  if (ihash->mm_allocator != NULL) {
    ihash->head = NULL; // uthash initializer
  } else {
    HASH_ITER(hh,ihash->head,ihash_element,tmp) {
      HASH_DEL(ihash->head,ihash_element);
      uthash_free(ihash_element,0);
    }
  }
}
void ihash_delete(ihash_t* const ihash) {
  ihash_clear(ihash);
  mm_free(ihash);
}
/*
 * Basic (Type-unsafe) Accessors
 */
ihash_element_t* ihash_get_ihash_element(ihash_t* const ihash,const int64_t key) {
  ihash_element_t *ihash_element;
  HASH_FIND_INT(ihash->head,&key,ihash_element);
  return ihash_element;
}
void ihash_insert_element(ihash_t* ihash,const int64_t key,void* const element) {
  ihash_element_t* ihash_element = ihash_get_ihash_element(ihash,key);
  if (gem_expect_true(ihash_element==NULL)) {
    ihash_element = uthash_malloc(sizeof(ihash_element_t));
    ihash_element->key = key;
    ihash_element->element = element;
    HASH_ADD_INT(ihash->head,key,ihash_element);
  } else {
    // No removal of replaced element
    ihash_element->element = element;
  }
}
void ihash_remove_element(ihash_t* ihash,const int64_t key) {
  ihash_element_t* ihash_element = ihash_get_ihash_element(ihash,key);
  if (ihash_element) {
    HASH_DEL(ihash->head,ihash_element);
  }
}
void* ihash_get_element(ihash_t* const ihash,const int64_t key) {
  ihash_element_t* const ihash_element = ihash_get_ihash_element(ihash,key);
  return gem_expect_true(ihash_element!=NULL) ? ihash_element->element : NULL;
}
/*
 * Type-safe Accessors
 */
bool ihash_is_contained(ihash_t* const ihash,const int64_t key) {
  return (ihash_get_ihash_element(ihash,key)!=NULL);
}
uint64_t ihash_get_num_elements(ihash_t* const ihash) {
  return (uint64_t)HASH_COUNT(ihash->head);
}
uint64_t ihash_get_size(ihash_t* const ihash) {
  /*
   * The hash handle consumes about 32 bytes per item on a 32-bit system, or 56 bytes per item on a 64-bit system.
   * The other overhead costs (the buckets and the table) are negligible in comparison.
   * You can use HASH_OVERHEAD to get the overhead size, in bytes, for a hash table.
   */
  return ihash_get_num_elements(ihash)*IHASH_SIZE_PER_ELEMENT;
}
/*
 * Miscellaneous
 */
int ihash_cmp_keys(int64_t* a,int64_t* b) {
  /*
   * return (int) -1 if (a < b)
   * return (int)  0 if (a == b)
   * return (int)  1 if (a > b)
   */
  return *a-*b;
}
#define ihash_cmp_keys_wrapper(arg1,arg2) ihash_cmp_keys((int64_t*)arg1,(int64_t*)arg2)
void ihash_sort_by_key(ihash_t* ihash) {
  HASH_SORT(ihash->head,ihash_cmp_keys_wrapper); // Sort
}
/*
 * Iterator
 */
ihash_iterator_t* ihash_iterator_new(ihash_t* const ihash) {
  // Allocate
  ihash_iterator_t* const iterator = mm_alloc(ihash_iterator_t);
  // Init
  iterator->current = NULL;
  iterator->next = ihash->head;
  return iterator;
}
void ihash_iterator_delete(ihash_iterator_t* const iterator) {
  mm_free(iterator);
}
bool ihash_iterator_eoi(ihash_iterator_t* const iterator) {
  return (iterator->next!=NULL);
}
bool ihash_iterator_next(ihash_iterator_t* const iterator) {
  if (gem_expect_true(iterator->next!=NULL)) {
    iterator->current = iterator->next;
    iterator->next = iterator->next->hh.next;
    return true;
  } else {
    iterator->current = NULL;
    return false;
  }
}
int64_t ihash_iterator_get_key(ihash_iterator_t* const iterator) {
  return iterator->current->key;
}
void* ihash_iterator_get_element(ihash_iterator_t* const iterator) {
  return iterator->current->element;
}
