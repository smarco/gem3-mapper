/*
 * PROJECT: GEM-Tools library
 * FILE: gt_ihash.h
 * DATE: 2/09/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_IHASH_H_
#define GT_IHASH_H_

#include "gt_commons.h"
#include "uthash.h"

/*
 * Integer Key Hash
 */
typedef struct {
  int64_t key;
  void* element;
  gt_hash_element_type element_type;
  union {
    size_t element_size;
    gt_hash_element_setup element_setup;
  };
  UT_hash_handle hh;
} gt_ihash_element;
typedef struct {
  gt_ihash_element* ihash_head;
} gt_ihash;
typedef struct {
  gt_ihash* ihash;
  gt_ihash_element* next;
} gt_ihash_iterator;

/*
 * Constructor
 */
GT_INLINE gt_ihash* gt_ihash_new(void);
GT_INLINE void gt_ihash_clear(gt_ihash* const ihash,const bool free_element);
GT_INLINE void gt_ihash_delete(gt_ihash* const ihash,const bool free_element);
GT_INLINE void gt_ihash_destroy(gt_ihash* const ihash);

/*
 * Basic (Type-unsafe) Accessors
 */
GT_INLINE gt_ihash_element* gt_ihash_get_ihash_element(gt_ihash* const ihash,const int64_t key);
GT_INLINE void gt_ihash_insert_primitive(
    gt_ihash* const ihash,const int64_t key,void* const element,const int64_t element_size);
GT_INLINE void gt_ihash_insert_object(
    gt_ihash* const ihash,const int64_t key,
    void* const object,void* (*element_dup_fx)(),void (*element_free_fx)());
GT_INLINE void* gt_ihash_get_element(gt_ihash* const ihash,const int64_t key);
GT_INLINE void gt_ihash_remove(gt_ihash* const ihash,const int64_t key,const bool free_element);

/*
 * Type-safe Accessors
 */
#define gt_ihash_get(ihash,integer_key,type) ((type*)gt_ihash_get_element(ihash,integer_key))
#define gt_ihash_insert(ihash,integer_key,element,type) gt_ihash_insert_primitive(ihash,integer_key,(void*)element,sizeof(type))
#define gt_ihash_insert_string(ihash,integer_key,string) gt_ihash_insert_object(ihash,integer_key,(void*)element,gt_string_dup,gt_string_delete)
GT_INLINE bool gt_ihash_is_contained(gt_ihash* const ihash,const int64_t key);
GT_INLINE uint64_t gt_ihash_get_num_elements(gt_ihash* const ihash);

/*
 * Miscellaneous
 */
GT_INLINE gt_ihash* gt_ihash_dup(gt_ihash* const ihash);
GT_INLINE void gt_ihash_copy(gt_ihash* const ihash_dst,gt_ihash* const ihash_src);
GT_INLINE void gt_ihash_sort_by_key(gt_ihash* const ihash);

/*
 * Iterator
 */
#define GT_IHASH_BEGIN_ITERATE(ihash,it_ikey,it_element,type) { \
  gt_ihash_element *ihash_##ih_element, *ihash_##tmp; \
  HASH_ITER(hh,ihash->ihash_head,ihash_##ih_element,ihash_##tmp) { \
    type* const it_element = (type*)(ihash_##ih_element->element); \
    int64_t const it_ikey = ihash_##ih_element->key;
#define GT_IHASH_END_ITERATE }}

GT_INLINE gt_ihash_iterator* gt_ihash_iterator_new(gt_ihash* const ihash);
GT_INLINE void gt_ihash_iterator_delete(gt_ihash_iterator* const ihash_iterator);

GT_INLINE bool gt_ihash_iterator_next(gt_ihash_iterator* const ihash_iterator);
GT_INLINE int64_t gt_ihash_iterator_get_key(gt_ihash_iterator* const ihash_iterator);
GT_INLINE void* gt_ihash_iterator_get_element(gt_ihash_iterator* const ihash_iterator);

#endif /* GT_IHASH_H_ */
