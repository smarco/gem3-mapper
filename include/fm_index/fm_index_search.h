/*
 * PROJECT: GEMMapper
 * FILE: fm_index_search.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef FM_INDEX_SEARCH_H_
#define FM_INDEX_SEARCH_H_

#include "utils/essentials.h"
#include "fm_index/fm_index.h"

/*
 * Basic FM-Index search (backwards)
 */
void fm_index_bsearch(
    const fm_index_t* const fm_index,
    const uint8_t* const key,
    uint64_t key_length,
    uint64_t* const hi_out,
    uint64_t* const lo_out);

void fm_index_bsearch_pure(
    const fm_index_t* const fm_index,
    const uint8_t* const key,
    uint64_t key_length,
    uint64_t* const hi_out,
    uint64_t* const lo_out);

uint64_t fm_index_bsearch_continue(
    const fm_index_t* const fm_index,
    const char* const key,
    const uint64_t key_length,
    const bool* const allowed_repl,
    uint64_t last_hi,
    uint64_t last_lo,
    uint64_t begin_pos,
    const uint64_t end_pos,
    uint64_t* const res_hi,
    uint64_t* const res_lo);

/*
 * Basic FM-Index search (forward)
 */
void fm_index_reverse_bsearch_pure(
    const fm_index_t* const fm_index,
    const uint8_t* const key,
    const uint64_t key_length,
    uint64_t* const hi_out,
    uint64_t* const lo_out);

/*
 * Debug
 */
void fm_index_bsearch_debug(
    const fm_index_t* const fm_index,
    const uint8_t* const key,
    uint64_t key_length,
    uint64_t* const hi_out,
    uint64_t* const lo_out,
    uint64_t* const steps_out);

#endif /* FM_INDEX_SEARCH_H_ */
