/*
 * PROJECT: GEM-Tools library
 * FILE: gt_misms.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_misms.h"

/*
 * Constructors
 */
GT_INLINE gt_misms* gt_misms_new() {
  gt_misms* misms = gt_alloc(gt_misms);
  misms->position = 0;
  misms->size = 0;
  return misms;
}
GT_INLINE void gt_misms_delete(gt_misms* misms) {
  GT_NULL_CHECK(misms);
  gt_free(misms);
}
GT_INLINE void gt_misms_set_mismatch(gt_misms* const misms,const uint64_t position,const char base) {
  GT_NULL_CHECK(misms);
  misms->misms_type = MISMS;
  misms->position = position;
  misms->base = base;
}
GT_INLINE void gt_misms_set_insertion(gt_misms* const misms,const uint64_t position,const uint64_t size) {
  GT_NULL_CHECK(misms);
  misms->misms_type = INS;
  misms->position = position;
  misms->size = size;
}
GT_INLINE void gt_misms_set_deletion(gt_misms* const misms,const uint64_t position,const uint64_t size) {
  GT_NULL_CHECK(misms);
  misms->misms_type = DEL;
  misms->position = position;
  misms->size = size;
}

/*
 * Accessors
 */
GT_INLINE gt_misms_t gt_misms_get_type(gt_misms* const misms) {
  GT_NULL_CHECK(misms);
  return misms->misms_type;
}
GT_INLINE uint64_t gt_misms_get_position(gt_misms* const misms) {
  GT_NULL_CHECK(misms);
  return misms->position;
}
// Mismatches
GT_INLINE char gt_misms_get_base(gt_misms* const misms) {
  GT_NULL_CHECK(misms);
  gt_check(misms->misms_type!=MISMS,MISMS_TYPE);
  return misms->base;
}
// Insertion/Deletion
GT_INLINE uint64_t gt_misms_get_size(gt_misms* const misms) {
  GT_NULL_CHECK(misms);
  gt_check(misms->misms_type==MISMS,MISMS_TYPE);
  return misms->size;
}
