/*
 * PROJECT: GEMMapper
 * FILE: shash.c
 * DATE: 10/07/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Simple Hash Implementation (String Key, Generic Value)
 */

#include "hash.h"

#include "mm.h"
#include "shash.h"

/*
 * Constructor
 */
GEM_INLINE shash_t* shash_new(void) {
  shash_t* const shash = mm_alloc(shash_t);
  shash->head = NULL; // uthash initializer
  return shash;
}
GEM_INLINE void shash_clear(shash_t* shash) {
  HASH_CHECK(shash);
  shash_element_t *shash_element, *tmp;
  HASH_ITER(hh,shash->head,shash_element,tmp) {
    HASH_DEL(shash->head,shash_element);
    mm_free(shash_element);
  }
}
GEM_INLINE void shash_delete(shash_t* shash) {
  HASH_CHECK(shash);
  shash_clear(shash);
  mm_free(shash);
}
/*
 * Basic (Type-unsafe) Accessors
 */
GEM_INLINE shash_element_t* shash_get_shash_element(shash_t* const shash,char* const key) {
  HASH_CHECK(shash);
  GEM_CHECK_NULL(key);
  shash_element_t* shash_element;
  HASH_FIND_STR(shash->head,key,shash_element);
  return shash_element;
}
GEM_INLINE void shash_insert_element(shash_t* shash,char* const key,const uint64_t key_length,void* const element) {
  HASH_CHECK(shash);
  GEM_CHECK_NULL(key);
  GEM_CHECK_NULL(element);
  shash_element_t* shash_element = shash_get_shash_element(shash,key);
  if (gem_expect_true(shash_element==NULL)) {
    shash_element = mm_alloc(shash_element_t);
    shash_element->key = key;
    shash_element->element = element;
    HASH_ADD_KEYPTR(hh,shash->head,shash_element->key,key_length,shash_element);
  } else {
    // No removal of replaced element
    shash_element->element = element;
  }
}
GEM_INLINE void shash_remove(shash_t* shash,char* const key) {
  HASH_CHECK(shash);
  GEM_CHECK_NULL(key);
  shash_element_t* shash_element = shash_get_shash_element(shash,key);
  if (shash_element) {
    HASH_DEL(shash->head,shash_element);
  }
}
GEM_INLINE void* shash_get_element(shash_t* const shash,char* const key) {
  HASH_CHECK(shash);
  GEM_CHECK_NULL(key);
  shash_element_t* shash_element = shash_get_shash_element(shash,key);
  return gem_expect_true(shash_element!=NULL) ? shash_element->element : NULL;
}
/*
 * Type-safe Accessors
 */
GEM_INLINE bool shash_is_contained(shash_t* const shash,char* const key) {
  HASH_CHECK(shash);
  GEM_CHECK_NULL(key);
  return (shash_get_shash_element(shash,key)!=NULL);
}
GEM_INLINE uint64_t shash_get_num_elements(shash_t* const shash) {
  HASH_CHECK(shash);
  return (uint64_t)HASH_COUNT(shash->head);
}
/*
 * Iterator
 */
GEM_INLINE shash_iterator_t* shash_iterator_new(shash_t* const shash) {
  HASH_CHECK(shash);
  // Allocate
  shash_iterator_t* const iterator = mm_alloc(shash_iterator_t);
  // Init
  iterator->current = NULL;
  iterator->next = shash->head;
  return iterator;
}
GEM_INLINE void shash_iterator_delete(shash_iterator_t* const iterator) {
  SHASH_ITERATOR_CHECK(iterator);
  mm_free(iterator);
}
GEM_INLINE bool shash_iterator_eoi(shash_iterator_t* const iterator) {
  SHASH_ITERATOR_CHECK(iterator);
  return (iterator->next!=NULL);
}
GEM_INLINE bool shash_iterator_next(shash_iterator_t* const iterator) {
  SHASH_ITERATOR_CHECK(iterator);
  if (gem_expect_true(iterator->next!=NULL)) {
    iterator->current = iterator->next;
    iterator->next = iterator->next->hh.next;
    return true;
  } else {
    iterator->current = NULL;
    return false;
  }
}
GEM_INLINE char* shash_iterator_get_key(shash_iterator_t* const iterator) {
  SHASH_ITERATOR_CHECK(iterator);
  return iterator->current->key;
}
GEM_INLINE void* shash_iterator_get_element(shash_iterator_t* const iterator) {
  SHASH_ITERATOR_CHECK(iterator);
  return iterator->current->element;
}

