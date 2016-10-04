/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   FM-Index module provides high-level FM-index functions
 *   for exact query into the index
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
void fm_index_reverse_bsearch_fb(
    const fm_index_t* const fm_index,
    const uint8_t* const key,
    const uint64_t key_length,
    uint64_t* const hi_out,
    uint64_t* const lo_out);
void fm_index_reverse_bsearch_bf(
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
