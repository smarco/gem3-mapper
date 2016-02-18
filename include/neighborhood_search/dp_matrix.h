/*
 * PROJECT: GEMMapper
 * FILE: dp_matrix.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef DP_MATRIX_H_
#define DP_MATRIX_H_

#include "utils/essentials.h"

/*
 * Constants
 */
#define NS_INF                              UINT32_MAX
#define NS_SET_PRIORITY(cell_distance,prio) ((uint64_t)(cell_distance) | (uint64_t)(prio))
#define NS_HAS_PRIORITY(cell_distance,prio) ((cell_distance)%2==prio)
#define NS_ENCODE_DISTANCE(distance)        ((distance)<<1)
#define NS_DECODE_DISTANCE(distance)        ((distance)>>1)

/*
 * DP-Matrix
 */
typedef struct {
  // Cells
  uint64_t* cells;
  // Band
  uint64_t band_high_offset;
  uint64_t band_low_offset;
} dp_column_t;

typedef struct {
  dp_column_t* columns;
  uint64_t num_columns;
  uint64_t column_length;
} dp_matrix_t;

/*
 * Display
 */
void dp_column_print_summary(
    FILE* const stream,const dp_matrix_t* const dp_matrix,
    const uint64_t column_position,const uint64_t lo,const uint64_t hi);
void dp_matrix_print(
    FILE* const stream,const dp_matrix_t* const dp_matrix,const bool forward_search,
    const uint8_t* const key,const uint64_t key_begin,const uint64_t key_end,
    const uint8_t* const text,const uint64_t text_begin,const uint64_t text_end);

#endif /* DP_MATRIX_H_ */
