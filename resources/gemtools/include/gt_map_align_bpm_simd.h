/*
 * PROJECT: GEM-Tools library
 * FILE: gt_map_align_bpm_simd.h
 * DATE: 20/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Myers' Bit-vector algorithm.
 *   Implementation using SIMD instructions (vectorial extensions)
 */

#ifndef GT_MAP_ALIGN_BMP_SIMD_H_
#define GT_MAP_ALIGN_BMP_SIMD_H_

#include "gt_commons.h"
#include "gt_map_align.h"

/*
 * Bit-compressed (Re)alignment
 *   BMP[BitParalellMayers] - Myers' Fast Bit-Vector algorithm (Levenshtein)
 */
GT_INLINE bool gt_map_block_bpm_get_distance_simd128_unfold64(
    gt_bpm_pattern* const bpm_pattern,char* const sequence,const uint64_t sequence_length,
    uint64_t* const position,uint64_t* const distance);
GT_INLINE bool gt_map_block_bpm_get_distance_simd128(
    gt_bpm_pattern* const bpm_pattern,char* const sequence,const uint64_t sequence_length,
    uint64_t* const position,uint64_t* const distance);

#endif /* GT_MAP_ALIGN_BMP_SIMD_H_ */
