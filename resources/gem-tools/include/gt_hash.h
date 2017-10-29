/*
 * PROJECT: GEM-Tools library
 * FILE: gt_hash.c
 * DATE: 2/09/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_HASH_H_
#define GT_HASH_H_

#include "gt_commons.h"

/*
 * Checkers
 */
#define GT_HASH_CHECK(hash) GT_NULL_CHECK(hash)

/*
 * Internal Type Setup/Handlers
 */
typedef enum { GT_HASH_TYPE_REGULAR, GT_HASH_TYPE_OBJECT } gt_hash_element_type;

typedef struct {
  void* (*element_dup_fx)();
  void (*element_free_fx)();
} gt_hash_element_setup;

GT_INLINE void gt_hash_free_element(void* const element,const gt_hash_element_type element_type);
GT_INLINE void* gt_hash_copy_element(void* const element,const gt_hash_element_type element_type,const int64_t element_size);

/*
 * Key-specific Hash
 */
#include "gt_ihash.h"
#include "gt_shash.h"

#endif /* GT_HASH_H_ */
