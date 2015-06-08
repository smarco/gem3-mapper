/*
 * PROJECT: GEM-Tools library
 * FILE: gt_map_align_bpm64.h
 * DATE: 20/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Myers' Bit-vector algorithm. Straightaway implementation on 64bits vectors
 */

#ifndef GT_MAP_ALIGN_BMP_H_
#define GT_MAP_ALIGN_BMP_H_

#include "gt_commons.h"
#include "gt_map.h"
#include "gt_map_align.h"

// Constants
#define GT_BMP_W64_SIZE 8

/*
 * Bit-compressed (Re)alignment
 *   BMP[BitParalellMayers] - Myers' Fast Bit-Vector algorithm (Levenshtein)
 */
GT_INLINE bool gt_map_block_bpm_get_distance(
    gt_bpm_pattern* const bpm_pattern,char* const sequence,const uint64_t sequence_length,
    uint64_t* const position,uint64_t* const distance,const uint64_t max_distance);
GT_INLINE bool gt_map_block_bpm_get_distance__cutoff(
    gt_bpm_pattern* const bpm_pattern,char* const sequence,const uint64_t sequence_length,
    uint64_t* const position,uint64_t* const distance,const uint64_t max_distance);

#endif /* GT_MAP_ALIGN_BMP_H_ */
