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
 *   FM-Index module provides high-level FM-index query functionality
 *   Includes bidirectional FM-index queries & all-chars expanded query
 */

#include "fm_index/fm_index_query.h"

/*
 * Setup
 */
void fm_index_2query_init(
    const fm_index_t* const fm_index,
    fm_2interval_t* const fm_2interval) {
  const uint64_t fm_index_length = fm_index_get_length(fm_index);
  fm_2interval->backward_lo = 0;
  fm_2interval->backward_hi = fm_index_length;
  fm_2interval->forward_lo = 0;
  fm_2interval->forward_hi = fm_index_length;
}
/*
 * FM-Index Precomputed Bidirectional Operators
 */
void fm_index_precompute_2erank(
    const fm_index_t* const fm_index,
    fm_2erank_elms_t* const fm_2erank_elms,
    const uint64_t position) {
  // Compute all exclusive ranks
  bwt_block_locator_t block_loc;
  bwt_block_elms_t block_elms;
  bwt_t* const bwt = fm_index->bwt;
  bwt_precompute(bwt,position,&block_loc,&block_elms);
  // Store eranks
  fm_2erank_elms->pranks[ENC_DNA_CHAR_A] = bwt_precomputed_erank(bwt,ENC_DNA_CHAR_A,&block_loc,&block_elms);
  fm_2erank_elms->pranks[ENC_DNA_CHAR_C] = bwt_precomputed_erank(bwt,ENC_DNA_CHAR_C,&block_loc,&block_elms);
  fm_2erank_elms->pranks[ENC_DNA_CHAR_G] = bwt_precomputed_erank(bwt,ENC_DNA_CHAR_G,&block_loc,&block_elms);
  fm_2erank_elms->pranks[ENC_DNA_CHAR_T] = bwt_precomputed_erank(bwt,ENC_DNA_CHAR_T,&block_loc,&block_elms);
  fm_2erank_elms->pranks[ENC_DNA_CHAR_N] = bwt_precomputed_erank(bwt,ENC_DNA_CHAR_N,&block_loc,&block_elms);
  fm_2erank_elms->pranks[ENC_DNA_CHAR_SEP] = bwt_precomputed_erank(bwt,ENC_DNA_CHAR_SEP,&block_loc,&block_elms);
}
/*
 * Compute the lo-delta/hi-delta difference as to maintain coherence between forward & backward search interval
 * Eg. char_enc == C
 *    +---+ -\
 *    | A |  | > lo-delta
 *    +---+ -/
 *    | C |
 *    +---+ -\
 *    | G |  |
 *    +---+  |
 *    | T |  | > hi-delta
 *    +---+  |
 *    | N |  |
 *    +---+ -/
 */
uint64_t fm_2erank_elms_upper_acc_erank(
    fm_2erank_elms_t* const lo_2erank_elms,
    fm_2erank_elms_t* const hi_2erank_elms,
    const uint8_t char_enc) {
  // Switch character
  uint64_t acc = 0;
  switch (char_enc) {
    case ENC_DNA_CHAR_SEP:
      acc += hi_2erank_elms->pranks[ENC_DNA_CHAR_N]-lo_2erank_elms->pranks[ENC_DNA_CHAR_N];
      // no break
    case ENC_DNA_CHAR_N:
      acc += hi_2erank_elms->pranks[ENC_DNA_CHAR_A]-lo_2erank_elms->pranks[ENC_DNA_CHAR_A];
      // no break
    case ENC_DNA_CHAR_T:
      acc += hi_2erank_elms->pranks[ENC_DNA_CHAR_C]-lo_2erank_elms->pranks[ENC_DNA_CHAR_C];
      // no break
    case ENC_DNA_CHAR_G:
      acc += hi_2erank_elms->pranks[ENC_DNA_CHAR_G]-lo_2erank_elms->pranks[ENC_DNA_CHAR_G];
      // no break
    case ENC_DNA_CHAR_C:
      acc += hi_2erank_elms->pranks[ENC_DNA_CHAR_T]-lo_2erank_elms->pranks[ENC_DNA_CHAR_T];
      // no break
    case ENC_DNA_CHAR_A:
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  // Return
  return acc;
}
uint64_t fm_2erank_elms_lower_acc_erank(
    fm_2erank_elms_t* const lo_2erank_elms,
    fm_2erank_elms_t* const hi_2erank_elms,
    const uint8_t char_enc) {
  // Switch character
  uint64_t acc = 0;
  switch (char_enc) {
    case ENC_DNA_CHAR_A:
      acc += hi_2erank_elms->pranks[ENC_DNA_CHAR_G]-lo_2erank_elms->pranks[ENC_DNA_CHAR_G];
      // no break
    case ENC_DNA_CHAR_C:
      acc += hi_2erank_elms->pranks[ENC_DNA_CHAR_C]-lo_2erank_elms->pranks[ENC_DNA_CHAR_C];
      // no break
    case ENC_DNA_CHAR_G:
      acc += hi_2erank_elms->pranks[ENC_DNA_CHAR_A]-lo_2erank_elms->pranks[ENC_DNA_CHAR_A];
      // no break
    case ENC_DNA_CHAR_T:
      acc += hi_2erank_elms->pranks[ENC_DNA_CHAR_N]-lo_2erank_elms->pranks[ENC_DNA_CHAR_N];
      // no break
    case ENC_DNA_CHAR_N:
      acc += hi_2erank_elms->pranks[ENC_DNA_CHAR_SEP]-lo_2erank_elms->pranks[ENC_DNA_CHAR_SEP];
      // no break
    case ENC_DNA_CHAR_SEP:
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  // Return
  return acc;
}
void fm_index_precomputed_2query_forward(
    fm_2interval_t* const fm_2interval_in,
    fm_2interval_t* const fm_2interval_out,
    fm_2erank_elms_t* const lo_2erank_elms,
    fm_2erank_elms_t* const hi_2erank_elms,
    const uint8_t char_enc) {
  // Assign forward values (eranks computed forward)
  const uint8_t mirror_char_enc = dna_encoded_complement(char_enc);
  fm_2interval_out->forward_lo = lo_2erank_elms->pranks[mirror_char_enc];
  fm_2interval_out->forward_hi = hi_2erank_elms->pranks[mirror_char_enc];
  // Assign backward values (eranks computed forward)
  fm_2interval_out->backward_lo = fm_2interval_in->backward_lo +
      fm_2erank_elms_upper_acc_erank(lo_2erank_elms,hi_2erank_elms,char_enc);
  fm_2interval_out->backward_hi = fm_2interval_in->backward_hi -
      fm_2erank_elms_lower_acc_erank(lo_2erank_elms,hi_2erank_elms,char_enc);
}
void fm_index_precomputed_2query_backward(
    fm_2interval_t* const fm_2interval_in,
    fm_2interval_t* const fm_2interval_out,
    fm_2erank_elms_t* const lo_2erank_elms,
    fm_2erank_elms_t* const hi_2erank_elms,
    const uint8_t char_enc) {
  // Assign backward values (eranks computed backward)
  fm_2interval_out->backward_lo = lo_2erank_elms->pranks[char_enc];
  fm_2interval_out->backward_hi = hi_2erank_elms->pranks[char_enc];
  // Assign forward values (eranks computed backward)
  const uint8_t mirror_char_enc = dna_encoded_complement(char_enc);
  fm_2interval_out->forward_lo = fm_2interval_in->forward_lo +
      fm_2erank_elms_upper_acc_erank(lo_2erank_elms,hi_2erank_elms,mirror_char_enc);
  fm_2interval_out->forward_hi = fm_2interval_in->forward_hi -
      fm_2erank_elms_lower_acc_erank(lo_2erank_elms,hi_2erank_elms,mirror_char_enc);
}
/*
 * FM-Index Bidirectional Operators
 */
void fm_index_2query_forward_query(
    const fm_index_t* const fm_index,
    fm_2interval_t* const fm_2interval_in,
    fm_2interval_t* const fm_2interval_out,
    const uint8_t char_enc) {
  // Precompute erank
  fm_2erank_elms_t lo_2erank_elms, hi_2erank_elms;
  fm_index_precompute_2erank(fm_index,&lo_2erank_elms,fm_2interval_in->forward_lo);
  fm_index_precompute_2erank(fm_index,&hi_2erank_elms,fm_2interval_in->forward_hi);
  // Compute next 2interval
  fm_index_precomputed_2query_forward(
      fm_2interval_in,fm_2interval_out,&lo_2erank_elms,&hi_2erank_elms,char_enc);
}
void fm_index_2query_backward_query(
    const fm_index_t* const fm_index,
    fm_2interval_t* const fm_2interval_in,
    fm_2interval_t* const fm_2interval_out,
    const uint8_t char_enc) {
  // Precompute erank
  fm_2erank_elms_t lo_2erank_elms, hi_2erank_elms;
  fm_index_precompute_2erank(fm_index,&lo_2erank_elms,fm_2interval_in->backward_lo);
  fm_index_precompute_2erank(fm_index,&hi_2erank_elms,fm_2interval_in->backward_hi);
  // Compute next 2interval
  fm_index_precomputed_2query_backward(
      fm_2interval_in,fm_2interval_out,&lo_2erank_elms,&hi_2erank_elms,char_enc);
}
/*
 * FM-Index Bidirectional 2-step Operators
 */
void fm_index_2query_forward_precompute(
    const fm_index_t* const fm_index,
    fm_2interval_t* const fm_2interval,
    fm_2erank_elms_t* const lo_2erank_elms,
    fm_2erank_elms_t* const hi_2erank_elms) {
  // Precompute erank
  fm_index_precompute_2erank(fm_index,lo_2erank_elms,fm_2interval->forward_lo);
  fm_index_precompute_2erank(fm_index,hi_2erank_elms,fm_2interval->forward_hi);
}
void fm_index_2query_backward_precompute(
    const fm_index_t* const fm_index,
    fm_2interval_t* const fm_2interval,
    fm_2erank_elms_t* const lo_2erank_elms,
    fm_2erank_elms_t* const hi_2erank_elms) {
  // Precompute erank
  fm_index_precompute_2erank(fm_index,lo_2erank_elms,fm_2interval->backward_lo);
  fm_index_precompute_2erank(fm_index,hi_2erank_elms,fm_2interval->backward_hi);
}
void fm_index_2query_precomputed_forward_query(
    fm_2erank_elms_t* const lo_2erank_elms,
    fm_2erank_elms_t* const hi_2erank_elms,
    fm_2interval_t* const fm_2interval_in,
    fm_2interval_t* const fm_2interval_out,
    const uint8_t char_enc) {
  // Compute next 2interval
  fm_index_precomputed_2query_forward(
      fm_2interval_in,fm_2interval_out,
      lo_2erank_elms,hi_2erank_elms,char_enc);
}
void fm_index_2query_precomputed_backward_query(
    fm_2erank_elms_t* const lo_2erank_elms,
    fm_2erank_elms_t* const hi_2erank_elms,
    fm_2interval_t* const fm_2interval_in,
    fm_2interval_t* const fm_2interval_out,
    const uint8_t char_enc) {
  // Compute next 2interval
  fm_index_precomputed_2query_backward(
      fm_2interval_in,fm_2interval_out,
      lo_2erank_elms,hi_2erank_elms,char_enc);
}
