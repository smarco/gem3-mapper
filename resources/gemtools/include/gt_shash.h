/*
 * PROJECT: GEM-Tools library
 * FILE: gt_shash.h
 * DATE: 10/07/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_SHASH_H_
#define GT_SHASH_H_

#include "gt_hash.h"
#include "uthash.h"

/*
 * String Key Hash
 */
typedef struct {
  char* key;
  void* element;
  gt_hash_element_type element_type;
  //union {
    size_t element_size;
    gt_hash_element_setup element_setup;
  //};
  UT_hash_handle hh;
} gt_shash_element;
typedef struct {
  gt_shash_element* shash_head;
} gt_shash;
typedef struct {
  gt_shash* shash;
  gt_shash_element* next;
} gt_shash_iterator;

/*
 * Constructor
 */
GT_INLINE gt_shash* gt_shash_new(void);
GT_INLINE void gt_shash_clear(gt_shash* const shash,const bool free_element);
GT_INLINE void gt_shash_delete(gt_shash* const shash,const bool free_element);
GT_INLINE void gt_shash_destroy(gt_shash* const shash);

/*
 * Basic (Type-unsafe) Accessors
 */
GT_INLINE gt_shash_element* gt_shash_get_shash_element(gt_shash* const shash,char* const key);
GT_INLINE char* gt_shash_insert_primitive(
    gt_shash* const shash,char* const key,void* const element,const int64_t element_size);
GT_INLINE char* gt_shash_insert_object(
    gt_shash* const shash,char* const key,
    void* const object,void* (*element_dup_fx)(),void (*element_free_fx)());
GT_INLINE char* gt_shash_get_key(gt_shash* const shash,char* const key);
GT_INLINE void* gt_shash_get_element(gt_shash* const shash,char* const key);
GT_INLINE void gt_shash_remove(gt_shash* const shash,char* const key,const bool free_element);

/*
 * Type-safe Accessors
 */
#define gt_shash_get(shash,string_key,type) ((type*)gt_shash_get_element(shash,string_key))
#define gt_shash_insert(shash,string_key,element,type) gt_shash_insert_primitive(shash,string_key,(void*)element,sizeof(type))
#define gt_shash_insert_string(shash,string_key,string) gt_shash_insert_object(shash,string_key,(void*)string,(void*(*)())gt_string_dup,(void(*)())gt_string_delete)
GT_INLINE bool gt_shash_is_contained(gt_shash* const shash,char* const key);
GT_INLINE uint64_t gt_shash_get_num_elements(gt_shash* const shash);

/*
 * Miscellaneous
 */
GT_INLINE gt_shash* gt_shash_dup(gt_shash* const shash);
GT_INLINE void gt_shash_copy(gt_shash* const shash_dst,gt_shash* const shash_src);

/*
 * Iterator
 */
#define GT_SHASH_BEGIN_ITERATE(shash,it_skey,it_element,type) { \
  gt_shash_element *shash_##sh_element, *shash_##tmp; \
  HASH_ITER(hh,shash->shash_head,shash_##sh_element,shash_##tmp) { \
    type* const it_element = (type*)(shash_##sh_element->element); \
    char* const it_skey = shash_##sh_element->key;

#define GT_SHASH_BEGIN_ELEMENT_ITERATE(shash,it_element,type) { \
  gt_shash_element *shash_##sh_element, *shash_##tmp; \
  HASH_ITER(hh,shash->shash_head,shash_##sh_element,shash_##tmp) { \
    type* const it_element = (type*)(shash_##sh_element->element);

#define GT_SHASH_BEGIN_KEY_ITERATE(shash,it_skey) { \
  gt_shash_element *shash_##sh_element, *shash_##tmp; \
  HASH_ITER(hh,shash->shash_head,shash_##sh_element,shash_##tmp) { \
    char* const it_skey = shash_##sh_element->key;

#define GT_SHASH_END_ITERATE }}

GT_INLINE gt_shash_iterator* gt_shash_iterator_new(gt_shash* const shash);
GT_INLINE void gt_shash_iterator_delete(gt_shash_iterator* const shash_iterator);

GT_INLINE bool gt_shash_iterator_next(gt_shash_iterator* const shash_iterator);
GT_INLINE char* gt_shash_iterator_get_key(gt_shash_iterator* const shash_iterator);
GT_INLINE void* gt_shash_iterator_get_element(gt_shash_iterator* const shash_iterator);


#endif /* GT_SHASH_H_ */
