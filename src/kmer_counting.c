/*
 * PROJECT: GEMMapper
 * FILE: kmer_counting.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *
 * TODO
 *   Trick:
 *     Check for Ns, if it doesn't have N's then make a much simpler check of the candidate read
 *
 *   Title: "Verification tests"
 *     Can be heuristic, like discarding a candidate region if
 *     if doesn't have enough seed supporting it. Can be exact like k-mer counting.
 *     Both can give an upper/lower bound on the matching distance.
 *       Seed-counting => At least 1 error per seed not matching (it needs at least
 *                        MismsRegions to make a match)
 *       K-mer counting => ...
 *
 *   Title: "Sorting candidate regions towards aligning only the best"
 *
 */

#include "kmer_counting.h"
#include "bpm_align.h"

/*
 * Constants
 */
#define KMER_COUNTING_LENGTH 4
#define KMER_COUNTING_MASK 0x00000000000000FFull // 0x00000000000003FFull
#define KMER_COUNTING_NUM_KMERS POW4(KMER_COUNTING_LENGTH)
#define KMER_COUNTING_MASK_INDEX(kmer_idx) ((kmer_idx) & KMER_COUNTING_MASK)
#define KMER_COUNTING_ADD_INDEX(kmer_idx,enc_char)  kmer_idx = (kmer_idx<<2 | (enc_char))
#define KMER_COUNTING_ADD_INDEX__MASK(kmer_idx,enc_char) kmer_idx = KMER_COUNTING_MASK_INDEX(kmer_idx<<2 | (enc_char))

/*
 * Compile Pattern
 */
GEM_INLINE void kmer_counting_compile(
    kmer_counting_t* const kmer_counting,bool* const allowed_enc,
    uint8_t* const pattern,const uint64_t pattern_length,mm_stack_t* const mm_stack) {
  // Check length
  if (pattern_length < KMER_COUNTING_LENGTH) return;
  // Init
  kmer_counting->kmer_count_text = mm_stack_calloc(mm_stack,KMER_COUNTING_NUM_KMERS,uint16_t,true);
  kmer_counting->kmer_count_pattern = mm_stack_calloc(mm_stack,KMER_COUNTING_NUM_KMERS,uint16_t,true);
  kmer_counting->allowed_enc = allowed_enc;
  kmer_counting->pattern_length = pattern_length;
  // Count kmers in pattern
  uint64_t uncalled_impasse = KMER_COUNTING_LENGTH-1;
  uint64_t pos, kmer_idx;
  for (pos=0,kmer_idx=0;pos<pattern_length;++pos) {
    const uint8_t enc_char = pattern[pos]; // Update kmer-index
    KMER_COUNTING_ADD_INDEX__MASK(kmer_idx,enc_char);
    if (!allowed_enc[enc_char]) { // Check allowed
      uncalled_impasse = KMER_COUNTING_LENGTH-1;
    } else {
      if (gem_expect_false(uncalled_impasse > 0)) { // Check impasse
        --uncalled_impasse;
      } else {
        // Increment kmer-count
        ++(kmer_counting->kmer_count_pattern[KMER_COUNTING_MASK_INDEX(kmer_idx)]);
      }
    }
  }
}
/*
 * Filter text region
 */
GEM_INLINE uint64_t kmer_counting_filter(
    const kmer_counting_t* const kmer_counting,
    const uint8_t* const text,const uint64_t text_length,const uint64_t max_error) {
  // Check length
  if (kmer_counting->pattern_length < KMER_COUNTING_LENGTH) return 0; // Don't know
  const uint64_t kmers_error = KMER_COUNTING_LENGTH*max_error;
  const uint64_t kmers_max = kmer_counting->pattern_length - (KMER_COUNTING_LENGTH-1);
  if (kmers_error >= kmers_max) return 0;
  const uint64_t kmers_required = kmer_counting->pattern_length - (KMER_COUNTING_LENGTH-1) - KMER_COUNTING_LENGTH*max_error;
  const bool* const allowed_enc = kmer_counting->allowed_enc;
  memset(kmer_counting->kmer_count_text,0,KMER_COUNTING_NUM_KMERS*UINT16_SIZE);
  uint16_t* const kmer_count_text = kmer_counting->kmer_count_text;
  uint16_t* const kmer_count_pattern = kmer_counting->kmer_count_pattern;
  /*
   * First count (Load)
   */
  const uint64_t init_chunk = MIN(text_length,kmer_counting->pattern_length);
  uint64_t kmers_in_text = 0;
  uint64_t end_pos, kmer_idx_end, uncalled_impasse_end = KMER_COUNTING_LENGTH-1;
  for (end_pos=0,kmer_idx_end=0;end_pos<init_chunk;++end_pos) {
    const uint8_t enc_char_end = text[end_pos];
    KMER_COUNTING_ADD_INDEX__MASK(kmer_idx_end,enc_char_end);
    if (!allowed_enc[enc_char_end]) { // Check allowed
      uncalled_impasse_end = KMER_COUNTING_LENGTH-1;
    } else {
      if (gem_expect_false(uncalled_impasse_end > 0)) { // Check impasse
        --uncalled_impasse_end;
      } else {
        // Increment kmer-count
        const uint16_t count_pattern = kmer_count_pattern[kmer_idx_end];
        if (count_pattern>0 && (kmer_count_text[kmer_idx_end])++ < count_pattern) ++kmers_in_text;
      }
    }
  }
  // Check filter condition
  if (kmers_in_text >= kmers_required) return 0;
  if (init_chunk == text_length) return ALIGN_COLUMN_INF;
  /*
   * Synch begin-index
   */
  uint64_t begin_pos, kmer_idx_begin, uncalled_impasse_begin = 0;
  for (begin_pos=0,kmer_idx_begin=0;begin_pos<KMER_COUNTING_LENGTH-1;++begin_pos) {
    // Update kmer-index
    const uint8_t enc_char_begin = text[begin_pos];
    KMER_COUNTING_ADD_INDEX(kmer_idx_begin,enc_char_begin);
    // Check allowed
    if (!allowed_enc[enc_char_begin]) {
      uncalled_impasse_begin = KMER_COUNTING_LENGTH-1;
    } else {
      // Check impasse
      if (gem_expect_false(uncalled_impasse_begin > 0)) --uncalled_impasse_begin;
    }
  }
  /*
   * Sliding window count
   */
  for (;end_pos<text_length;++end_pos,++begin_pos) {
    // Begin (Decrement)
    const uint8_t enc_char_begin = text[begin_pos];
    KMER_COUNTING_ADD_INDEX__MASK(kmer_idx_begin,enc_char_begin);
    if (!allowed_enc[enc_char_begin]) { // Check allowed
      uncalled_impasse_begin = KMER_COUNTING_LENGTH-1;
    } else {
      if (gem_expect_false(uncalled_impasse_begin > 0)) { // Check impasse
        --uncalled_impasse_begin;
      } else {
        // Decrement kmer-count
        const uint16_t count_pattern = kmer_count_pattern[kmer_idx_begin];
        if (count_pattern > 0 && kmer_count_text[kmer_idx_begin]-- <= count_pattern) --kmers_in_text;
      }
    }
    // End (Increment)
    const uint8_t enc_char_end = text[end_pos];
    KMER_COUNTING_ADD_INDEX__MASK(kmer_idx_end,enc_char_end);
    if (!allowed_enc[enc_char_end]) { // Check allowed
      uncalled_impasse_end = KMER_COUNTING_LENGTH-1;
    } else {
      if (gem_expect_false(uncalled_impasse_end > 0)) { // Check impasse
        --uncalled_impasse_end;
      } else {
        // Increment kmer-count
        const uint16_t count_pattern = kmer_count_pattern[kmer_idx_end];
        if (count_pattern > 0 && kmer_count_text[kmer_idx_end]++ < count_pattern) ++kmers_in_text;
      }
    }
    // Check filter condition
    if (kmers_in_text >= kmers_required) {
//      uint64_t n, total=0, kmers_ot = 0;
//      for (n=0;n<KMER_COUNTING_NUM_KMERS;++n) {
//        kmers_ot += kmer_count_text[n];
//        if (kmer_count_pattern[n]==0) {
//          fprintf(stderr," ");
//        } else if (kmer_count_text[n]>=kmer_count_pattern[n]) {
//          fprintf(stderr,"-"); total+=kmer_count_pattern[n];
//        } else {
//          fprintf(stderr,"*"); total+=kmer_count_text[n];
//        }
//      }
//      fprintf(stderr,"\n");
//      for (n=0;n<KMER_COUNTING_NUM_KMERS;++n) {
//        fprintf(stderr,"(%"PRIu64"/%"PRIu64")[%"PRIu64"] Pattern=%u Text=%u \n",total,kmers_ot,n,kmer_count_pattern[n],kmer_count_text[n]);
//      }
      return begin_pos - (KMER_COUNTING_LENGTH-1);
    }
  }
  // Not passing the filter
  return ALIGN_COLUMN_INF;
}
///*
// * Filter text region
// */
//GEM_INLINE uint64_t kmer_counting_filter(
//    const kmer_counting_t* const kmer_counting,
//    const uint8_t* const text,const uint64_t text_length,const uint64_t max_error) {
//  // Check length
//  if (kmer_counting->pattern_length < KMER_COUNTING_LENGTH) return 0; // Don't know
//  const uint64_t kmers_error = KMER_COUNTING_LENGTH*max_error;
//  const uint64_t kmers_max = kmer_counting->pattern_length - (KMER_COUNTING_LENGTH-1);
//  if (kmers_error >= kmers_max) return 0;
//  const uint64_t kmers_required = kmer_counting->pattern_length - (KMER_COUNTING_LENGTH-1) - KMER_COUNTING_LENGTH*max_error;
//  const bool* const allowed_enc = kmer_counting->allowed_enc;
//  memset(kmer_counting->kmer_count_text,0,KMER_COUNTING_NUM_KMERS*UINT16_SIZE);
//  uint16_t* const kmer_count_text = kmer_counting->kmer_count_text;
//  uint16_t* const kmer_count_pattern = kmer_counting->kmer_count_pattern;
//  /*
//   * First count (Load)
//   */
//  const uint64_t init_chunk = MIN(text_length,kmer_counting->pattern_length);
//  uint64_t kmers_in_text = 0;
//  uint64_t end_pos, kmer_idx_end;
//  for (end_pos=0,kmer_idx_end=0;end_pos<init_chunk;++end_pos) {
//    const uint8_t enc_char_end = text[end_pos];
//    KMER_COUNTING_ADD_INDEX__MASK(kmer_idx_end,enc_char_end);
//    // Increment kmer-count
//    const uint16_t count_pattern = kmer_count_pattern[kmer_idx_end];
//    if (count_pattern>0 && (kmer_count_text[kmer_idx_end])++ < count_pattern) ++kmers_in_text;
//  }
//  // Check filter condition
//  if (kmers_in_text >= kmers_required) return 0;
//  if (init_chunk == text_length) return BPM_ALIGN_DISTANCE_INF;
//  /*
//   * Synch begin-index
//   */
//  uint64_t begin_pos, kmer_idx_begin;
//  for (begin_pos=0,kmer_idx_begin=0;begin_pos<KMER_COUNTING_LENGTH-1;++begin_pos) {
//    // Update kmer-index
//    const uint8_t enc_char_begin = text[begin_pos];
//    KMER_COUNTING_ADD_INDEX(kmer_idx_begin,enc_char_begin);
//  }
//  /*
//   * Sliding window count
//   */
//  for (;end_pos<text_length;++end_pos,++begin_pos) {
//    // Begin (Decrement)
//    const uint8_t enc_char_begin = text[begin_pos];
//    KMER_COUNTING_ADD_INDEX__MASK(kmer_idx_begin,enc_char_begin);
//    // Decrement kmer-count
//    const uint16_t count_pattern_begin = kmer_count_pattern[kmer_idx_begin];
//    if (count_pattern_begin > 0 && kmer_count_text[kmer_idx_begin]-- <= count_pattern_begin) --kmers_in_text;
//    // End (Increment)
//    const uint8_t enc_char_end = text[end_pos];
//    KMER_COUNTING_ADD_INDEX__MASK(kmer_idx_end,enc_char_end);
//    // Increment kmer-count
//    const uint16_t count_pattern_end = kmer_count_pattern[kmer_idx_end];
//    if (count_pattern_end > 0 && kmer_count_text[kmer_idx_end]++ < count_pattern_end) ++kmers_in_text;
//    // Check filter condition
//    if (kmers_in_text >= kmers_required) {
//      return begin_pos - (KMER_COUNTING_LENGTH-1);
//    }
//  }
//  // Not passing the filter
//  return BPM_ALIGN_DISTANCE_INF;
//}





