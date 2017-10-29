/*
 * PROJECT: GEM-Tools library
 * FILE: gt_shash.c
 * DATE: 10/07/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_hash.h"
#include "gt_commons.h"
#include "gt_error.h"
#include "gt_mm.h"
#include "gt_string.h"

/*
 * Setup
 */
GT_INLINE void gt_shash_free_element(gt_shash_element* const shash_element) {
  GT_NULL_CHECK(shash_element);
  switch (shash_element->element_type) {
    case GT_HASH_TYPE_REGULAR:
      gt_free(shash_element->element);
      break;
    case GT_HASH_TYPE_OBJECT:
      shash_element->element_setup.element_free_fx(shash_element->element);
      break;
    default:
      gt_fatal_error(SELECTION_NOT_VALID);
      break;
  }
}
GT_INLINE void gt_shash_free_shash_element(gt_shash_element* const shash_element,const bool free_element) {
  GT_NULL_CHECK(shash_element);
  // Free key
  gt_free(shash_element->key);
  // Free element
  if (free_element) gt_shash_free_element(shash_element);
  // Free handler
  gt_free(shash_element);
}

/*
 * Constructor
 */
GT_INLINE gt_shash* gt_shash_new(void) {
  gt_shash* shash = gt_alloc(gt_shash);
  shash->shash_head = NULL; // uthash initializer
  return shash;
}
GT_INLINE void gt_shash_clear(gt_shash* const shash,const bool free_element) {
  GT_HASH_CHECK(shash);
  gt_shash_element *shash_element, *tmp;
  HASH_ITER(hh,shash->shash_head,shash_element,tmp) {
    HASH_DEL(shash->shash_head,shash_element);
    gt_shash_free_shash_element(shash_element,free_element);
  }
}
GT_INLINE void gt_shash_delete(gt_shash* const shash,const bool free_element) {
  GT_HASH_CHECK(shash);
  gt_shash_clear(shash,free_element);
  gt_free(shash);
}
GT_INLINE void gt_shash_destroy(gt_shash* const shash) {
  GT_HASH_CHECK(shash);
  gt_shash_clear(shash,true);
  gt_free(shash);
}

/*
 * Basic (Type-unsafe) Accessors
 */
GT_INLINE gt_shash_element* gt_shash_get_shash_element(gt_shash* const shash,char* const key) {
  GT_HASH_CHECK(shash);
  GT_NULL_CHECK(key);
  gt_shash_element *shash_element;
  HASH_FIND_STR(shash->shash_head,key,shash_element);
  return shash_element;
}
GT_INLINE char* gt_shash_insert_primitive(
    gt_shash* const shash,char* const key,void* const element,const int64_t element_size) {
  GT_HASH_CHECK(shash);
  GT_ZERO_CHECK(element_size);
  GT_NULL_CHECK(key); GT_NULL_CHECK(element);
  gt_shash_element *shash_element = gt_shash_get_shash_element(shash,key);
  if (gt_expect_true(shash_element==NULL)) {
    shash_element = gt_alloc(gt_shash_element);
    const uint64_t key_length = strlen(key);
    shash_element->key = gt_strndup(key,key_length);
    shash_element->element = element;
    HASH_ADD_KEYPTR(hh,shash->shash_head,shash_element->key,key_length,shash_element);
  } else {
    gt_shash_free_element(shash_element);
    shash_element->element = element;
  }
  // Set shash element type
  shash_element->element_type = GT_HASH_TYPE_REGULAR;
  shash_element->element_size = element_size;
  // Return key
  return shash_element->key;
}
GT_INLINE char* gt_shash_insert_object(
    gt_shash* const shash,char* const key,
    void* const object,void* (*element_dup_fx)(),void (*element_free_fx)()) {
  GT_HASH_CHECK(shash);
  GT_NULL_CHECK(key); GT_NULL_CHECK(object);
  GT_NULL_CHECK(element_dup_fx); GT_NULL_CHECK(element_free_fx);
  gt_shash_element *shash_element = gt_shash_get_shash_element(shash,key);
  if (gt_expect_true(shash_element==NULL)) {
    shash_element = gt_alloc(gt_shash_element);
    const uint64_t key_length = strlen(key);
    shash_element->key = gt_strndup(key,key_length);
    shash_element->element = object;
    HASH_ADD_KEYPTR(hh,shash->shash_head,shash_element->key,key_length,shash_element);
  } else {
    gt_shash_free_element(shash_element);
    shash_element->element = object;
  }
  // Set shash element type
  shash_element->element_type = GT_HASH_TYPE_OBJECT;
  shash_element->element_setup.element_dup_fx = element_dup_fx;
  shash_element->element_setup.element_free_fx = element_free_fx;
  // Return key
  return shash_element->key;
}
GT_INLINE char* gt_shash_get_key(gt_shash* const shash,char* const key) {
  GT_HASH_CHECK(shash);
  GT_NULL_CHECK(key);
  gt_shash_element *shash_element = gt_shash_get_shash_element(shash,key);
  return gt_expect_true(shash_element!=NULL) ? shash_element->key : NULL;
}
GT_INLINE void* gt_shash_get_element(gt_shash* const shash,char* const key) {
  GT_HASH_CHECK(shash);
  GT_NULL_CHECK(key);
  gt_shash_element *shash_element = gt_shash_get_shash_element(shash,key);
  return gt_expect_true(shash_element!=NULL) ? shash_element->element : NULL;
}
GT_INLINE void gt_shash_remove(gt_shash* const shash,char* const key,const bool free_element) {
  GT_HASH_CHECK(shash);
  GT_NULL_CHECK(key);
  gt_shash_element *shash_element = gt_shash_get_shash_element(shash,key);
  if (shash_element) {
    HASH_DEL(shash->shash_head,shash_element);
    gt_shash_free_shash_element(shash_element,free_element);
  }
}

/*
 * Type-safe Accessors
 */
GT_INLINE bool gt_shash_is_contained(gt_shash* const shash,char* const key) {
  GT_HASH_CHECK(shash);
  GT_NULL_CHECK(key);
  return (gt_shash_get_shash_element(shash,key)!=NULL);
}
GT_INLINE uint64_t gt_shash_get_num_elements(gt_shash* const shash) {
  GT_HASH_CHECK(shash);
  return (uint64_t)HASH_COUNT(shash->shash_head);
}

/*
 * Miscellaneous
 */
GT_INLINE gt_shash* gt_shash_dup(gt_shash* const shash) {
  gt_shash* const shash_cp =  gt_shash_new();
  gt_shash_copy(shash_cp,shash);
  return shash_cp;
}
GT_INLINE void gt_shash_copy(gt_shash* const shash_dst,gt_shash* const shash_src) {
  GT_SHASH_BEGIN_ITERATE(shash_src,skey,shash_element,void) {
    // Insert element into the copy
    switch (shash_sh_element->element_type) {
      case GT_HASH_TYPE_REGULAR: {
        // Copy element
        void* const shash_element_cp = gt_malloc_(1,shash_sh_element->element_size,false,0);
        memcpy(shash_element_cp,shash_element,shash_sh_element->element_size);
        gt_shash_insert_primitive(shash_dst,skey,shash_element_cp,shash_sh_element->element_size);
        break;
      }
      case GT_HASH_TYPE_OBJECT:
        gt_shash_insert_object(shash_dst,skey,shash_sh_element->element_setup.element_dup_fx(shash_element),
            shash_sh_element->element_setup.element_dup_fx,shash_sh_element->element_setup.element_free_fx);
        break;
      default:
        gt_fatal_error(SELECTION_NOT_VALID);
        break;
    }
  } GT_SHASH_END_ITERATE;
}

/*
 * Iterator
 */
GT_INLINE gt_shash_iterator* gt_shash_iterator_new(gt_shash* const shash) {
  GT_HASH_CHECK(shash);
  // Allocate
  gt_shash_iterator* const shash_iterator = gt_alloc(gt_shash_iterator);
  // Init
  shash_iterator->shash = shash;
  shash_iterator->next = shash->shash_head;
  return shash_iterator;
}
GT_INLINE void gt_shash_iterator_delete(gt_shash_iterator* const shash_iterator) {
  GT_HASH_CHECK(shash_iterator->shash);
  gt_free(shash_iterator);
}
GT_INLINE bool gt_shash_iterator_next(gt_shash_iterator* const shash_iterator) {
  GT_HASH_CHECK(shash_iterator->shash);
  if (gt_expect_false(shash_iterator->next==NULL)) return false;
  shash_iterator->next = shash_iterator->next->hh.next;
  return true;
}
GT_INLINE char* gt_shash_iterator_get_key(gt_shash_iterator* const shash_iterator) {
  GT_HASH_CHECK(shash_iterator->shash);
  GT_HASH_CHECK(shash_iterator->next);
  return shash_iterator->next->key;
}
GT_INLINE void* gt_shash_iterator_get_element(gt_shash_iterator* const shash_iterator) {
  GT_HASH_CHECK(shash_iterator->shash);
  GT_HASH_CHECK(shash_iterator->next);
  return shash_iterator->next->element;
}

