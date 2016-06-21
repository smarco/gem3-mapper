/*
 * PROJECT: GEMMapper
 * FILE: fm_index_query.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef FM_INDEX_QUERY_H_
#define FM_INDEX_QUERY_H_

#include "utils/essentials.h"
#include "data_structures/dna_text.h"
#include "fm_index/fm_index.h"

/*
 * Bi-Directional queries
 */
typedef struct {
  /* Forward Interval (searching forward) */
  uint64_t forward_lo;
  uint64_t forward_hi;
  /* Backward Interval (searching backwards) */
  uint64_t backward_lo;
  uint64_t backward_hi;
} fm_2interval_t;

/*
 * Setup
 */
void fm_index_2query_init(
    const fm_index_t* const fm_index,
    fm_2interval_t* const fm_2interval);

/*
 * FM-Index Bidirectional Operators
 */
void fm_index_2query_forward(
    const fm_index_t* const fm_index,
    const uint8_t char_enc,
    fm_2interval_t* const fm_2interval_in,
    fm_2interval_t* const fm_2interval_out);
void fm_index_2query_backward(
    const fm_index_t* const fm_index,
    const uint8_t char_enc,
    fm_2interval_t* const fm_2interval_in,
    fm_2interval_t* const fm_2interval_out);

#endif /* FM_INDEX_QUERY_H_ */
