/*
 * PROJECT: GEM-Tools library
 * FILE: gt_ihash.c
 * DATE: 2/09/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_hash.h"
#include "gt_commons.h"
#include "gt_error.h"
#include "gt_mm.h"

/*
 * Setup
 */
GT_INLINE void gt_ihash_free_element(gt_ihash_element *ihash_element) {
  GT_NULL_CHECK(ihash_element);
  switch (ihash_element->element_type) {
    case GT_HASH_TYPE_REGULAR:
      gt_free(ihash_element->element);
      break;
    case GT_HASH_TYPE_OBJECT:
      ihash_element->element_setup.element_free_fx(ihash_element->element);
      break;
    default:
      gt_fatal_error(SELECTION_NOT_VALID);
      break;
  }
}
GT_INLINE void gt_ihash_free_ihash_element(gt_ihash_element *ihash_element,const bool free_element) {
  GT_NULL_CHECK(ihash_element);
  // Free element
  if (free_element) gt_ihash_free_element(ihash_element);
  // Free handler
  gt_free(ihash_element);
}

/*
 * Constructor
 */
GT_INLINE gt_ihash* gt_ihash_new(void) {
  gt_ihash* ihash = gt_alloc(gt_ihash);
  ihash->ihash_head = NULL; // uthash initializer
  return ihash;
}
GT_INLINE void gt_ihash_clear(gt_ihash* const ihash,const bool free_element) {
  GT_HASH_CHECK(ihash);
  gt_ihash_element *ihash_element, *tmp;
  HASH_ITER(hh,ihash->ihash_head,ihash_element,tmp) {
    HASH_DEL(ihash->ihash_head,ihash_element);
    gt_ihash_free_ihash_element(ihash_element,free_element);
  }
}
GT_INLINE void gt_ihash_delete(gt_ihash* const ihash,const bool free_element) {
  GT_HASH_CHECK(ihash);
  gt_ihash_clear(ihash,free_element);
  gt_free(ihash);
}
GT_INLINE void gt_ihash_destroy(gt_ihash* const ihash) {
  GT_HASH_CHECK(ihash);
  gt_ihash_clear(ihash,true);
  gt_free(ihash);
}

/*
 * Basic (Type-unsafe) Accessors
 */
GT_INLINE gt_ihash_element* gt_ihash_get_ihash_element(gt_ihash* const ihash,const int64_t key) {
  GT_HASH_CHECK(ihash);
  gt_ihash_element *ihash_element;
  HASH_FIND_INT(ihash->ihash_head,&key,ihash_element);
  return ihash_element;
}
GT_INLINE void gt_ihash_insert_primitive(
    gt_ihash* const ihash,const int64_t key,void* const element,const int64_t element_size) {
  GT_HASH_CHECK(ihash);
  GT_ZERO_CHECK(element_size);
  GT_NULL_CHECK(element);
  gt_ihash_element* ihash_element = gt_ihash_get_ihash_element(ihash,key);
  if (gt_expect_true(ihash_element==NULL)) {
    ihash_element = gt_alloc(gt_ihash_element);
    ihash_element->key = key;
    ihash_element->element = element;
    HASH_ADD_INT(ihash->ihash_head,key,ihash_element);
  } else {
    gt_ihash_free_element(ihash_element);
    ihash_element->element = element;
  }
  // Set shash element type
  ihash_element->element_type = GT_HASH_TYPE_REGULAR;
  ihash_element->element_size = element_size;
}
GT_INLINE void gt_ihash_insert_object(
    gt_ihash* const ihash,const int64_t key,
    void* const object,void* (*element_dup_fx)(),void (*element_free_fx)()) {
  GT_HASH_CHECK(ihash);
  GT_NULL_CHECK(object);
  GT_NULL_CHECK(element_dup_fx); GT_NULL_CHECK(element_free_fx);
  gt_ihash_element* ihash_element = gt_ihash_get_ihash_element(ihash,key);
  if (gt_expect_true(ihash_element==NULL)) {
    ihash_element = gt_alloc(gt_ihash_element);
    ihash_element->key = key;
    ihash_element->element = object;
    HASH_ADD_INT(ihash->ihash_head,key,ihash_element);
  } else {
    gt_ihash_free_element(ihash_element);
    ihash_element->element = object;
  }
  // Set ihash element type
  ihash_element->element_type = GT_HASH_TYPE_OBJECT;
  ihash_element->element_setup.element_dup_fx = element_dup_fx;
  ihash_element->element_setup.element_free_fx = element_free_fx;
}
GT_INLINE void* gt_ihash_get_element(gt_ihash* const ihash,const int64_t key) {
  GT_HASH_CHECK(ihash);
  gt_ihash_element* const ihash_element = gt_ihash_get_ihash_element(ihash,key);
  return gt_expect_true(ihash_element!=NULL) ? ihash_element->element : NULL;
}
GT_INLINE void gt_ihash_remove(gt_ihash* const ihash,const int64_t key,const bool free_element) {
  GT_HASH_CHECK(ihash);
  gt_ihash_element *ihash_element = gt_ihash_get_ihash_element(ihash,key);
  if (ihash_element) {
    HASH_DEL(ihash->ihash_head,ihash_element);
    gt_ihash_free_ihash_element(ihash_element,free_element);
  }
}

/*
 * Type-safe Accessors
 */
GT_INLINE bool gt_ihash_is_contained(gt_ihash* const ihash,const int64_t key) {
  GT_HASH_CHECK(ihash);
  return (gt_ihash_get_ihash_element(ihash,key)!=NULL);
}
GT_INLINE uint64_t gt_ihash_get_num_elements(gt_ihash* const ihash) {
  GT_HASH_CHECK(ihash);
  return (uint64_t)HASH_COUNT(ihash->ihash_head);
}

/*
 * Miscellaneous
 */
GT_INLINE gt_ihash* gt_ihash_dup(gt_ihash* const ihash) {
  gt_ihash* const ihash_cp =  gt_ihash_new();
  gt_ihash_copy(ihash_cp,ihash);
  return ihash_cp;
}
GT_INLINE void gt_ihash_copy(gt_ihash* const ihash_dst,gt_ihash* const ihash_src) {
  GT_IHASH_BEGIN_ITERATE(ihash_src,ikey,ihash_element,void) {
    // Insert element into the copy
    switch (ihash_ih_element->element_type) {
      case GT_HASH_TYPE_REGULAR: {
        // Copy element
        void* const ihash_element_cp = gt_malloc_(1,ihash_ih_element->element_size,false,0);
        memcpy(ihash_element_cp,ihash_element,ihash_ih_element->element_size);
        gt_ihash_insert_primitive(ihash_dst,ikey,ihash_element_cp,ihash_ih_element->element_size);
        break;
      }
      case GT_HASH_TYPE_OBJECT:
        gt_ihash_insert_object(ihash_dst,ikey,ihash_ih_element->element_setup.element_dup_fx(ihash_element),
            ihash_ih_element->element_setup.element_dup_fx,ihash_ih_element->element_setup.element_free_fx);
        break;
      default:
        gt_fatal_error(SELECTION_NOT_VALID);
        break;
    }
  } GT_IHASH_END_ITERATE;
}
int gt_ihash_cmp_keys(int64_t* a,int64_t* b) {
  /*
   * return (int) -1 if (a < b)
   * return (int)  0 if (a == b)
   * return (int)  1 if (a > b)
   */
  return a-b;
}
#define gt_ihash_cmp_keys_wrapper(arg1,arg2) gt_ihash_cmp_keys((int64_t*)arg1,(int64_t*)arg2)
GT_INLINE void gt_ihash_sort_by_key(gt_ihash* const ihash) {
  GT_HASH_CHECK(ihash);
  // Sort
  HASH_SORT(ihash->ihash_head,gt_ihash_cmp_keys_wrapper);
}


/*
 * Iterator
 */
GT_INLINE gt_ihash_iterator* gt_ihash_iterator_new(gt_ihash* const ihash) {
  GT_HASH_CHECK(ihash);
  // Allocate
  gt_ihash_iterator* const ihash_iterator = gt_alloc(gt_ihash_iterator);
  // Init
  ihash_iterator->ihash = ihash;
  ihash_iterator->next = ihash->ihash_head;
  return ihash_iterator;
}
GT_INLINE void gt_ihash_iterator_delete(gt_ihash_iterator* const ihash_iterator) {
  GT_HASH_CHECK(ihash_iterator->ihash);
  gt_free(ihash_iterator);
}
GT_INLINE bool gt_ihash_iterator_next(gt_ihash_iterator* const ihash_iterator) {
  GT_HASH_CHECK(ihash_iterator->ihash);
  if (gt_expect_false(ihash_iterator->next==NULL)) return false;
  ihash_iterator->next = ihash_iterator->next->hh.next;
  return true;
}
GT_INLINE int64_t gt_ihash_iterator_get_key(gt_ihash_iterator* const ihash_iterator) {
  GT_HASH_CHECK(ihash_iterator->ihash);
  GT_HASH_CHECK(ihash_iterator->next);
  return ihash_iterator->next->key;
}
GT_INLINE void* gt_ihash_iterator_get_element(gt_ihash_iterator* const ihash_iterator) {
  GT_HASH_CHECK(ihash_iterator->ihash);
  GT_HASH_CHECK(ihash_iterator->next);
  return ihash_iterator->next->element;
}

