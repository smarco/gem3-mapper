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

#include "fm_index/fm_index_search.h"
#include "fm_index/fm_index_query.h"
#include "fm_index/rank_mtable.h"

/*
 * Basic FM-Index search (backwards)
 */
void fm_index_bsearch(
    const fm_index_t* const fm_index,
    const uint8_t* const key,
    uint64_t key_length,
    uint64_t* const hi_out,
    uint64_t* const lo_out) {
  uint64_t lo, hi;
  // Rank queries against the lookup table
  const rank_mtable_t* const rank_mtable = fm_index->rank_table;
  rank_mquery_t query;
  rank_mquery_new(&query);
  while (key_length > 0 && !rank_mquery_is_exhausted(&query)) {
    const uint8_t c = key[key_length-1];
    if (c >= ENC_DNA_CHAR_N) break;
    rank_mquery_add_char(rank_mtable,&query,c); // Rank query (calculate offsets)
    --key_length;
  }
  // Query lookup table
  rank_mtable_fetch(rank_mtable,&query,&lo,&hi);
  // Check zero interval
  if (hi==lo || key_length==0) {
    *hi_out=hi;
    *lo_out=lo;
  } else {
    // Continue with ranks against the FM-Index
    while (key_length > 0 && hi > lo) {
      const uint8_t c = key[--key_length];
      lo = bwt_erank(fm_index->bwt,c,lo);
      hi = bwt_erank(fm_index->bwt,c,hi);
    }
    // Return results
    *hi_out=hi;
    *lo_out=lo;
  }
}
void fm_index_bsearch_pure(
    const fm_index_t* const fm_index,
    const uint8_t* const key,
    uint64_t key_length,
    uint64_t* const hi_out,
    uint64_t* const lo_out) {
  // Query lookup table
  uint64_t lo=0, hi=fm_index_get_length(fm_index);
  // Continue with ranks against the FM-Index
  while (key_length > 0 && hi > lo) {
    const uint8_t c = key[--key_length];
    lo = bwt_erank(fm_index->bwt,c,lo);
    hi = bwt_erank(fm_index->bwt,c,hi);
  }
  // printf("\n");
  // Return results
  *hi_out=hi;
  *lo_out=lo;
}
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
    uint64_t* const res_lo) {
  // Extends the exact search of a key from a given position and interval
  const bwt_t* const bwt = fm_index->bwt;
  // Continue search
  while (begin_pos > end_pos) {
    if (last_lo==last_hi) {
      *res_lo=last_lo;
      *res_hi=last_hi;
      return end_pos;
    }
    // Query step
    const uint8_t enc_char = key[--begin_pos];
    if (!allowed_repl[enc_char]) {
      ++begin_pos; break;
    }
    last_lo=bwt_erank(bwt,enc_char,last_lo);
    last_hi=bwt_erank(bwt,enc_char,last_hi);
  }
  *res_lo=last_lo;
  *res_hi=last_hi;
  return begin_pos;
}
/*
 * Basic FM-Index search (forward)
 */
void fm_index_reverse_bsearch_pure(
    const fm_index_t* const fm_index,
    const uint8_t* const key,
    const uint64_t key_length,
    uint64_t* const hi_out,
    uint64_t* const lo_out) {
  // Init
  fm_2interval_t fm_2interval;
  fm_2interval.backward_lo = 0;
  fm_2interval.backward_hi = fm_index_get_length(fm_index);
  fm_2interval.forward_lo = 0;
  fm_2interval.forward_hi = fm_index_get_length(fm_index);
  // Search using 2interval and eranks against the FM-Index
//  uint64_t i = 0;
//  while (i < key_length && fm_2interval.forward_lo < fm_2interval.forward_hi) {
//    fm_index_2query_forward(fm_index,key[i++],&fm_2interval,&fm_2interval);
//  }

  int64_t i = key_length-1;
  while (i>=0 && fm_2interval.forward_lo < fm_2interval.forward_hi) {
    fm_index_2query_backward_query(fm_index,&fm_2interval,&fm_2interval,key[i--]);
  }

  // Return results
  *lo_out = fm_2interval.backward_lo;
  *hi_out = fm_2interval.backward_hi;
}
void fm_index_reverse_bsearch_fb(
    const fm_index_t* const fm_index,
    const uint8_t* const key,
    const uint64_t key_length,
    uint64_t* const hi_out,
    uint64_t* const lo_out) {
  // Init
  fm_2interval_t fm_2interval;
  fm_2interval.backward_lo = 0;
  fm_2interval.backward_hi = fm_index_get_length(fm_index);
  fm_2interval.forward_lo = 0;
  fm_2interval.forward_hi = fm_index_get_length(fm_index);
  // Search second half (forward)
  int64_t i;
  for (i=10;i<key_length;++i) {
    fm_index_2query_forward_query(fm_index,&fm_2interval,&fm_2interval,key[i]);
    if (fm_2interval.forward_lo >= fm_2interval.forward_hi) {
      *lo_out = 0; *hi_out = 0; return;
    }
  }
  // Search first half (backwards)
  for (i=9;i>=0;--i) {
    fm_index_2query_backward_query(fm_index,&fm_2interval,&fm_2interval,key[i]);
    if (fm_2interval.forward_lo >= fm_2interval.forward_hi) {
      *lo_out = 0; *hi_out = 0; return;
    }
  }
  // Return results
  *lo_out = fm_2interval.backward_lo;
  *hi_out = fm_2interval.backward_hi;
}
void fm_index_reverse_bsearch_bf(
    const fm_index_t* const fm_index,
    const uint8_t* const key,
    const uint64_t key_length,
    uint64_t* const hi_out,
    uint64_t* const lo_out) {
  // Init
  fm_2interval_t fm_2interval;
  fm_2interval.backward_lo = 0;
  fm_2interval.backward_hi = fm_index_get_length(fm_index);
  fm_2interval.forward_lo = 0;
  fm_2interval.forward_hi = fm_index_get_length(fm_index);
  // Search first half (backwards)
  int64_t i;
  for (i=9;i>=0;--i) {
    fm_index_2query_backward_query(fm_index,&fm_2interval,&fm_2interval,key[i]);
    if (fm_2interval.forward_lo >= fm_2interval.forward_hi) {
      *lo_out = 0; *hi_out = 0; return;
    }
  }
  // Search second half (forward)
  for (i=10;i<key_length;++i) {
    fm_index_2query_forward_query(fm_index,&fm_2interval,&fm_2interval,key[i]);
    if (fm_2interval.forward_lo >= fm_2interval.forward_hi) {
      *lo_out = 0; *hi_out = 0; return;
    }
  }
  // Return results
  *lo_out = fm_2interval.backward_lo;
  *hi_out = fm_2interval.backward_hi;
}
/*
 * Debug
 */
void fm_index_bsearch_debug(
    const fm_index_t* const fm_index,
    const uint8_t* const key,
    uint64_t key_length,
    uint64_t* const hi_out,
    uint64_t* const lo_out,
    uint64_t* const steps_out) {
  // Query lookup table
  uint64_t steps=0, lo=0, hi=fm_index_get_length(fm_index);
  // Query ranks against the FM-Index
  while (key_length > 0 && hi > lo) {
    const uint8_t c = key[--key_length];
    lo = bwt_erank(fm_index->bwt,c,lo);
    hi = bwt_erank(fm_index->bwt,c,hi);
    ++steps;
  }
  // Return results
  *hi_out=hi;
  *lo_out=lo;
  *steps_out=steps;
}

//  bwt_block_locator_t bwt_block_locator_lo, bwt_block_locator_hi;
//  bwt_block_elms_t bwt_block_elms_lo, bwt_block_elms_hi;
//    if (bwt_is_same_bucket(lo,hi)) {
//      bwt_precompute_interval(fm_index->bwt,lo,hi,&bwt_block_locator_lo,&bwt_block_elms_lo);
//      bwt_precomputed_erank_interval(fm_index->bwt,c,&lo,&hi,&bwt_block_locator_lo,&bwt_block_elms_lo);
//    } else {
//      bwt_precompute(fm_index->bwt,lo,&bwt_block_locator_lo,&bwt_block_elms_lo);
//      bwt_precompute(fm_index->bwt,hi,&bwt_block_locator_hi,&bwt_block_elms_hi);
//      lo = bwt_precomputed_erank(fm_index->bwt,c,&bwt_block_locator_lo,&bwt_block_elms_lo);
//      hi = bwt_precomputed_erank(fm_index->bwt,c,&bwt_block_locator_hi,&bwt_block_elms_hi);
//    }

