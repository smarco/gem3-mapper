/*
 * PROJECT: GEM-Tools library
 * FILE: gt_map_align_bmp_generic.h
 * DATE: 20/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Myers' Bit-vector algorithm.
 *   Optimized version from the base implementation.
 *   Generic interface as to use any length word size
 */

#ifndef GT_MAP_ALIGN_BMP_GENERIC_H_
#define GT_MAP_ALIGN_BMP_GENERIC_H_

#include "gt_commons.h"
#include "gt_map.h"
#include "gt_map_align.h"

#define GT_BMP_GENERIC_FUNCTION_PASTE(x,y) x ## _ ## y
#define GT_BMP_GENERIC_FUNCTION_EVALUATOR(x,y)  GT_BMP_GENERIC_FUNCTION_PASTE(x,y)
#define GT_BMP_GENERIC_FUNCTION_NAME(function_name) GT_BMP_GENERIC_FUNCTION_EVALUATOR(function_name,gt_bmp_vector_t)

#endif /* GT_MAP_ALIGN_BMP_GENERIC_H_ */

/*
 * Bit-compressed (Re)alignment
 *   BMP[BitParalellMayers] - Myers' Fast Bit-Vector algorithm (Levenshtein)
 */
#ifndef gt_bmp_vector_t
  #error "GT.Generic.Code:: Type gt_bmp_vector_t should be defined"
#else

GT_INLINE bool GT_BMP_GENERIC_FUNCTION_NAME(gt_map_block_bpm_get_distance)(
    gt_bpm_pattern* const bpm_pattern,char* const sequence,const uint64_t sequence_length,
    uint64_t* const position,uint64_t* const distance);
GT_INLINE bool GT_BMP_GENERIC_FUNCTION_NAME(gt_map_block_bpm_get_distance__cutoff)(
    gt_bpm_pattern* const bpm_pattern,char* const sequence,const uint64_t sequence_length,
    uint64_t* const position,uint64_t* const distance,const uint64_t max_distance);

#endif
