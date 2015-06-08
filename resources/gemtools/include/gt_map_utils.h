/*
 * PROJECT: GEM-Tools library
 * FILE: gt_map_utils.h
 * DATE: 03/09/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_MAP_UTILS_H_
#define GT_MAP_UTILS_H_

#include "gt_essentials.h"
#include "gt_map.h"

GT_INLINE bool gt_map_block_overlap(
    gt_map* const map_block_1,gt_map* const map_block_2,
    uint64_t* const overlap_start,uint64_t* const overlap_end);

#endif /* GT_MAP_UTILS_H_ */
