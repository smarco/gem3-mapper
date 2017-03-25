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
 * DESCRIPTION: Simple Hash Implementation (String Key, Generic Value)
 */

#include "utils/hash.h"
#include "system/mm.h"

/*
 * Allocators
 */
#undef  uthash_malloc
#define uthash_malloc(sz)   (shash->mm_allocator!=NULL ? mm_allocator_malloc(shash->mm_allocator,sz) : malloc(sz))
#undef  uthash_free
#define uthash_free(ptr,sz) if (shash->mm_allocator==NULL) free(ptr);

/*
 * Constructor
 */
shash_t* shash_new(mm_allocator_t* const mm_allocator) {
  shash_t* const shash = mm_alloc(shash_t);
  shash->head = NULL; // uthash initializer
  shash->mm_allocator = mm_allocator;
  return shash;
}
void shash_clear(shash_t* shash) {
  shash_element_t *shash_element, *tmp;
  if (shash->mm_allocator != NULL) {
    shash->head = NULL; // uthash initializer
  } else {
    HASH_ITER(hh,shash->head,shash_element,tmp) {
      HASH_DEL(shash->head,shash_element);
      mm_free(shash_element);
    }
  }
}
void shash_delete(shash_t* shash) {
  shash_clear(shash);
  mm_free(shash);
}
/*
 * Basic (Type-unsafe) Accessors
 */
shash_element_t* shash_get_shash_element(shash_t* const shash,char* const key) {
  shash_element_t* shash_element;
  HASH_FIND_STR(shash->head,key,shash_element);
  return shash_element;
}
void shash_insert_element(
    shash_t* const shash,
    char* const key,
    const uint64_t key_length,
    void* const element) {
  shash_element_t* shash_element = shash_get_shash_element(shash,key);
  if (gem_expect_true(shash_element==NULL)) {
    shash_element = uthash_malloc(sizeof(shash_element_t));
    shash_element->key = key;
    shash_element->element = element;
    HASH_ADD_KEYPTR(hh,shash->head,shash_element->key,key_length,shash_element);
  } else {
    // No removal of replaced element
    shash_element->element = element;
  }
}
void shash_remove(shash_t* shash,char* const key) {
  shash_element_t* shash_element = shash_get_shash_element(shash,key);
  if (shash_element) {
    HASH_DEL(shash->head,shash_element);
  }
}
void* shash_get_element(shash_t* const shash,char* const key) {
  shash_element_t* shash_element = shash_get_shash_element(shash,key);
  return gem_expect_true(shash_element!=NULL) ? shash_element->element : NULL;
}
/*
 * Type-safe Accessors
 */
bool shash_is_contained(shash_t* const shash,char* const key) {
  return (shash_get_shash_element(shash,key)!=NULL);
}
uint64_t shash_get_num_elements(shash_t* const shash) {
  return (uint64_t)HASH_COUNT(shash->head);
}
/*
 * Iterator
 */
shash_iterator_t* shash_iterator_new(shash_t* const shash) {
  // Allocate
  shash_iterator_t* const iterator = mm_alloc(shash_iterator_t);
  // Init
  iterator->current = NULL;
  iterator->next = shash->head;
  return iterator;
}
void shash_iterator_delete(shash_iterator_t* const iterator) {
  mm_free(iterator);
}
bool shash_iterator_eoi(shash_iterator_t* const iterator) {
  return (iterator->next!=NULL);
}
bool shash_iterator_next(shash_iterator_t* const iterator) {
  if (gem_expect_true(iterator->next!=NULL)) {
    iterator->current = iterator->next;
    iterator->next = iterator->next->hh.next;
    return true;
  } else {
    iterator->current = NULL;
    return false;
  }
}
char* shash_iterator_get_key(shash_iterator_t* const iterator) {
  return iterator->current->key;
}
void* shash_iterator_get_element(shash_iterator_t* const iterator) {
  return iterator->current->element;
}

