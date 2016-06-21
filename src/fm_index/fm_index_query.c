/*
 * PROJECT: GEMMapper
 * FILE: fm_index_query.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "fm_index/fm_index_query.h"

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
void fm_index_2query_forward(
    const fm_index_t* const fm_index,
    const uint8_t char_enc,
    fm_2interval_t* const fm_2interval_in,
    fm_2interval_t* const fm_2interval_out) {
  // Precompute erank
  fm_2erank_elms_t lo_2erank_elms, hi_2erank_elms;
  fm_index_precompute_2erank(fm_index,&lo_2erank_elms,fm_2interval_in->forward_lo);
  fm_index_precompute_2erank(fm_index,&hi_2erank_elms,fm_2interval_in->forward_hi);
  // Compute next 2interval
  fm_index_precomputed_2query_forward(
      fm_2interval_in,fm_2interval_out,&lo_2erank_elms,&hi_2erank_elms,char_enc);
}
void fm_index_2query_backward(
    const fm_index_t* const fm_index,
    const uint8_t char_enc,
    fm_2interval_t* const fm_2interval_in,
    fm_2interval_t* const fm_2interval_out) {
  // Precompute erank
  fm_2erank_elms_t lo_2erank_elms, hi_2erank_elms;
  fm_index_precompute_2erank(fm_index,&lo_2erank_elms,fm_2interval_in->backward_lo);
  fm_index_precompute_2erank(fm_index,&hi_2erank_elms,fm_2interval_in->backward_hi);
  // Compute next 2interval
  fm_index_precomputed_2query_backward(
      fm_2interval_in,fm_2interval_out,&lo_2erank_elms,&hi_2erank_elms,char_enc);
}
