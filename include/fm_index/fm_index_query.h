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
 * Precomputed ranks
 */
typedef struct {
  uint64_t pranks[DNA_EXT_RANGE];         // Precomputed eranks
} fm_2erank_elms_t;

/*
 * Setup
 */
void fm_index_2query_init(
    const fm_index_t* const fm_index,
    fm_2interval_t* const fm_2interval);

/*
 * FM-Index Bidirectional Operators
 */
void fm_index_2query_forward_query(
    const fm_index_t* const fm_index,
    fm_2interval_t* const fm_2interval_in,
    fm_2interval_t* const fm_2interval_out,
    const uint8_t char_enc);
void fm_index_2query_backward_query(
    const fm_index_t* const fm_index,
    fm_2interval_t* const fm_2interval_in,
    fm_2interval_t* const fm_2interval_out,
    const uint8_t char_enc);

/*
 * FM-Index Bidirectional 2-step Operators
 */
void fm_index_2query_forward_precompute(
    const fm_index_t* const fm_index,
    fm_2interval_t* const fm_2interval,
    fm_2erank_elms_t* const lo_2erank_elms,
    fm_2erank_elms_t* const hi_2erank_elms);
void fm_index_2query_backward_precompute(
    const fm_index_t* const fm_index,
    fm_2interval_t* const fm_2interval,
    fm_2erank_elms_t* const lo_2erank_elms,
    fm_2erank_elms_t* const hi_2erank_elms);
void fm_index_2query_precomputed_forward_query(
    fm_2erank_elms_t* const lo_2erank_elms,
    fm_2erank_elms_t* const hi_2erank_elms,
    fm_2interval_t* const fm_2interval_in,
    fm_2interval_t* const fm_2interval_out,
    const uint8_t char_enc);
void fm_index_2query_precomputed_backward_query(
    fm_2erank_elms_t* const lo_2erank_elms,
    fm_2erank_elms_t* const hi_2erank_elms,
    fm_2interval_t* const fm_2interval_in,
    fm_2interval_t* const fm_2interval_out,
    const uint8_t char_enc);

#endif /* FM_INDEX_QUERY_H_ */
