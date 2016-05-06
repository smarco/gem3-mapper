/*
 * PROJECT: GEMMapper
 * FILE: kmer_counting.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
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

#include "filtering/kmer_counting.h"
#include "data_structures/dna_text.h"
#include "align/align.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PLOW

/*
 * Constants
 */
//#define KMER_COUNTING_LENGTH 3
//#define KMER_COUNTING_MASK   0x000000000000003Full
//#define KMER_COUNTING_LENGTH 4
//#define KMER_COUNTING_MASK   0x00000000000000FFull
//#define KMER_COUNTING_LENGTH 5
//#define KMER_COUNTING_MASK   0x00000000000003FFull
#define KMER_COUNTING_LENGTH 6
#define KMER_COUNTING_MASK   0x0000000000000FFFull
//#define KMER_COUNTING_LENGTH 7
//#define KMER_COUNTING_MASK   0x0000000000003FFFull

#define KMER_COUNTING_NUM_KMERS POW4(KMER_COUNTING_LENGTH)
#define KMER_COUNTING_MASK_INDEX(kmer_idx) ((kmer_idx) & KMER_COUNTING_MASK)
#define KMER_COUNTING_ADD_INDEX(kmer_idx,enc_char) kmer_idx = (kmer_idx<<2 | (enc_char))
#define KMER_COUNTING_ADD_INDEX__MASK(kmer_idx,enc_char) kmer_idx = KMER_COUNTING_MASK_INDEX(kmer_idx<<2 | (enc_char))

#define KMER_COUNTING_EFFECTIVE_THRESHOLD 12

/*
 * Compile Pattern
 */
void kmer_counting_compile(
    kmer_counting_t* const kmer_counting,
    uint8_t* const pattern,
    const uint64_t pattern_length,
    const uint64_t max_error,
    mm_stack_t* const mm_stack) {
  // Check min-length condition
  if (pattern_length < KMER_COUNTING_LENGTH) {
    kmer_counting->enabled = false;
    return;
  }
  // Check efficiency condition
  if (max_error > 0 && pattern_length/max_error < KMER_COUNTING_EFFECTIVE_THRESHOLD) {
    kmer_counting->enabled = false;
    return;
  }
  // Check max-error condition
  const uint64_t kmers_error = KMER_COUNTING_LENGTH * max_error;
  const uint64_t kmers_max = pattern_length - (KMER_COUNTING_LENGTH-1);
  if (kmers_error >= kmers_max) {
    kmer_counting->enabled = false;
    return;
  }
  // Init
  kmer_counting->enabled = true;
  kmer_counting->kmer_count_text = mm_stack_calloc(mm_stack,KMER_COUNTING_NUM_KMERS,uint16_t,true);
  kmer_counting->kmer_count_pattern = mm_stack_calloc(mm_stack,KMER_COUNTING_NUM_KMERS,uint16_t,true);
  kmer_counting->pattern_length = pattern_length;
  kmer_counting->max_error = max_error;
  // Count kmers in pattern
  uint64_t pos=0, kmer_idx=0, acc=0;
  // Compile all k-mers
  for (;pos<pattern_length;++pos) {
    const uint8_t enc_char = pattern[pos];
    if (!is_dna_canonical_encoded(enc_char)) {
      acc=0; kmer_idx=0; // Reset
    } else {
      KMER_COUNTING_ADD_INDEX__MASK(kmer_idx,enc_char);  // Update kmer-index
      if (acc < KMER_COUNTING_LENGTH-1) {
        ++acc; // Inc accumulator
      } else {
        ++(kmer_counting->kmer_count_pattern[kmer_idx]); // Increment kmer-count
      }
    }
  }
}
/*
 * Filter text region
 */
uint64_t kmer_counting_filter(
    const kmer_counting_t* const kmer_counting,
    const uint8_t* const text,
    const uint64_t text_length) {
  PROFILE_START(GP_FC_KMER_COUNTER_FILTER,PROFILE_LEVEL);
  // Check filter enabled
  if (!kmer_counting->enabled) {
    PROFILE_STOP(GP_FC_KMER_COUNTER_FILTER,PROFILE_LEVEL);
    PROF_INC_COUNTER(GP_FC_KMER_COUNTER_FILTER_NA);
    return 0; // Don't filter
  }
  // Prepare filter
  const uint64_t max_error = kmer_counting->max_error;
  const uint64_t kmers_required = kmer_counting->pattern_length -
      (KMER_COUNTING_LENGTH-1) - KMER_COUNTING_LENGTH*max_error;
  memset(kmer_counting->kmer_count_text,0,KMER_COUNTING_NUM_KMERS*UINT16_SIZE);
  uint16_t* const kmer_count_text = kmer_counting->kmer_count_text;
  uint16_t* const kmer_count_pattern = kmer_counting->kmer_count_pattern;
  /*
   * First count (Load)
   */
  const uint64_t total_kmers_text = MAX(text_length,kmer_counting->pattern_length);
  const uint64_t init_chunk = MIN(text_length,kmer_counting->pattern_length);
  uint64_t kmers_left = total_kmers_text, kmers_in_text = 0;
  uint64_t begin_pos, kmer_idx_begin = 0, acc_begin = 0;
  uint64_t end_pos, kmer_idx_end=0, acc_end = 0;
  // Initial window fill
  for (end_pos=0;end_pos<init_chunk;++end_pos,--kmers_left) {
    const uint8_t enc_char_end = text[end_pos];
    if (!is_dna_canonical_encoded(enc_char_end)) {
      acc_end = 0; kmer_idx_end = 0; // Reset
    } else {
      KMER_COUNTING_ADD_INDEX__MASK(kmer_idx_end,enc_char_end); // Update kmer-index
      if (acc_end < KMER_COUNTING_LENGTH-1) {
        ++acc_end;
      } else {
        // Increment kmer-count
        const uint16_t count_pattern = kmer_count_pattern[kmer_idx_end];
        if (count_pattern>0 && (kmer_count_text[kmer_idx_end])++ < count_pattern) ++kmers_in_text;
        // Check filter condition
        if (kmers_in_text >= kmers_required) {
          PROFILE_STOP(GP_FC_KMER_COUNTER_FILTER,PROFILE_LEVEL);
          return 0; // Don't filter
        } else if (kmers_required-kmers_in_text > kmers_left) { // Quick abandon
          PROFILE_STOP(GP_FC_KMER_COUNTER_FILTER,PROFILE_LEVEL);
          return ALIGN_DISTANCE_INF; // Filter out
        }
      }
    }
  }
  /*
   * Sliding window count
   */
  // Initial fill
  for (begin_pos=0;begin_pos<KMER_COUNTING_LENGTH-1;++begin_pos) {
    const uint8_t enc_char_begin = text[begin_pos];
    if (!is_dna_canonical_encoded(enc_char_begin)) {
      acc_begin = 0; kmer_idx_begin = 0; // Reset
    } else {
      KMER_COUNTING_ADD_INDEX__MASK(kmer_idx_begin,enc_char_begin);
      ++acc_begin;
    }
  }
  for (;end_pos<text_length;++end_pos,++begin_pos,--kmers_left) {
    // Begin (Decrement kmer-count)
    const uint8_t enc_char_begin = text[begin_pos];
    if (!is_dna_canonical_encoded(enc_char_begin)) {
      acc_begin = 0; kmer_idx_begin = 0; // Reset
    } else {
      KMER_COUNTING_ADD_INDEX__MASK(kmer_idx_begin,enc_char_begin);
      if (acc_begin < KMER_COUNTING_LENGTH-1) {
        ++acc_begin;
      } else {
        const uint16_t count_pattern_begin = kmer_count_pattern[kmer_idx_begin];
        if (count_pattern_begin > 0 && kmer_count_text[kmer_idx_begin]-- <= count_pattern_begin) --kmers_in_text;
      }
    }
    // End (Increment kmer-count)
    const uint8_t enc_char_end = text[end_pos];
    if (!is_dna_canonical_encoded(enc_char_end)) {
      acc_end = 0; kmer_idx_end = 0; // Reset
    } else {
      KMER_COUNTING_ADD_INDEX__MASK(kmer_idx_end,enc_char_end);
      if (acc_end < KMER_COUNTING_LENGTH-1) {
        ++acc_end;
      } else {
        const uint16_t count_pattern_end = kmer_count_pattern[kmer_idx_end];
        if (count_pattern_end > 0 && kmer_count_text[kmer_idx_end]++ < count_pattern_end) ++kmers_in_text;
      }
    }
    // Check filter condition
    if (kmers_in_text >= kmers_required) {
      PROFILE_STOP(GP_FC_KMER_COUNTER_FILTER,PROFILE_LEVEL);
      return 0; // Don't filter
    } else if (kmers_required-kmers_in_text > kmers_left) { // Quick abandon
      PROFILE_STOP(GP_FC_KMER_COUNTER_FILTER,PROFILE_LEVEL);
      return ALIGN_DISTANCE_INF; // Filter out
    }
  }
  // Not passing the filter
  PROFILE_STOP(GP_FC_KMER_COUNTER_FILTER,PROFILE_LEVEL);
  return ALIGN_DISTANCE_INF; // Filter out
}

