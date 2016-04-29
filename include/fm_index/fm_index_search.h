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
    const fm_index_t* const restrict fm_index,
    const uint8_t* const restrict key,
    uint64_t key_length,
    uint64_t* const restrict hi_out,
    uint64_t* const restrict lo_out);

void fm_index_bsearch_pure(
    const fm_index_t* const restrict fm_index,
    const uint8_t* const restrict key,
    uint64_t key_length,
    uint64_t* const restrict hi_out,
    uint64_t* const restrict lo_out);

uint64_t fm_index_bsearch_continue(
    const fm_index_t* const restrict fm_index,
    const char* const restrict key,
    const uint64_t key_length,
    const bool* const restrict allowed_repl,
    uint64_t last_hi,
    uint64_t last_lo,
    uint64_t begin_pos,
    const uint64_t end_pos,
    uint64_t* const restrict res_hi,
    uint64_t* const restrict res_lo);

/*
 * Basic FM-Index search (forward)
 */
void fm_index_reverse_bsearch_pure(
    const fm_index_t* const restrict fm_index,
    const uint8_t* const restrict key,
    const uint64_t key_length,
    uint64_t* const restrict hi_out,
    uint64_t* const restrict lo_out);

/*
 * Debug
 */
void fm_index_bsearch_debug(
    const fm_index_t* const restrict fm_index,
    const uint8_t* const restrict key,
    uint64_t key_length,
    uint64_t* const restrict hi_out,
    uint64_t* const restrict lo_out,
    uint64_t* const restrict steps_out);

#endif /* FM_INDEX_SEARCH_H_ */
