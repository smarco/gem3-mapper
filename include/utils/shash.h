/*
 * PROJECT: GEMMapper
 * FILE: shash.h
 * DATE: 10/07/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Simple Hash Implementation (String Key, Generic Value)
 */

#ifndef SHASH_H_
#define SHASH_H_

#include "system/commons.h"
#include "resources/uthash/uthash.h"

/*
 * String Key Hash
 */
typedef struct {
  char* key;
  void* element;
  UT_hash_handle hh;
} shash_element_t;
typedef struct {
  shash_element_t* head;
} shash_t;
typedef struct {
  shash_element_t* current;
  shash_element_t* next;
} shash_iterator_t;

/*
 * Constructor
 */
shash_t* shash_new(void);
void shash_clear(shash_t* const shash);
void shash_delete(shash_t* const shash);

/*
 * Basic (Type-unsafe) Accessors
 */
shash_element_t* shash_get_shash_element(shash_t* const shash,char* const key);
void* shash_get_element(shash_t* const shash,char* const key);
void shash_insert_element(
    shash_t* const shash,
    char* const key,
    const uint64_t key_length,
    void* const element);
void shash_remove_element(shash_t* const shash,char* const key);

/*
 * Type-safe Accessors
 */
#define shash_insert(shash,key,key_length,element) shash_insert_element(shash,key,key_length,(void*)element)
#define shash_remove_element(shash,key) shash_remove_element(shash,key)
#define shash_get(shash,key,type) ((type*)shash_get_element(shash,key))
bool shash_is_contained(shash_t* const shash,char* const key);
uint64_t shash_get_num_elements(shash_t* const shash);

/*
 * Iterator
 */
#define GT_SHASH_BEGIN_ITERATE(shash,it_skey,it_element,type) { \
  shash_element *shash_##sh_element, *shash_##tmp; \
  HASH_ITER(hh,shash,shash_##sh_element,shash_##tmp) { \
    type* const it_element = (type*)(shash_##sh_element->element); \
    char* const it_skey = shash_##sh_element->key;

#define GT_SHASH_BEGIN_ELEMENT_ITERATE(shash,it_element,type) { \
  shash_element *shash_##sh_element, *shash_##tmp; \
  HASH_ITER(hh,shash,shash_##sh_element,shash_##tmp) { \
    type* const it_element = (type*)(shash_##sh_element->element);

#define GT_SHASH_BEGIN_KEY_ITERATE(shash,it_skey) { \
  shash_element *shash_##sh_element, *shash_##tmp; \
  HASH_ITER(hh,shash,shash_##sh_element,shash_##tmp) { \
    char* const it_skey = shash_##sh_element->key;

#define GT_SHASH_END_ITERATE }}

shash_iterator_t* shash_iterator_new(shash_t* const shash);
void shash_iterator_delete(shash_iterator_t* const iterator);

bool shash_iterator_eoi(shash_iterator_t* const iterator);
bool shash_iterator_next(shash_iterator_t* const iterator);
char* shash_iterator_get_key(shash_iterator_t* const iterator);
void* shash_iterator_get_element(shash_iterator_t* const iterator);


#endif /* SHASH_H_ */
