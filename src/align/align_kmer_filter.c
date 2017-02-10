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
 *   Filter based on general k-mer counting as to quickly filter out
 *   candidates that cannot align against its region-text
 *
 *   Verification tests can be heuristic, like discarding a candidate region if
 *   if doesn't have enough seed supporting it. Out it can be exact like k-mer counting.
 *   Both can give an upper/lower bound on the matching distance.
 *     Seed-counting => At least 1 error per seed not matching (it needs at least
 *                      MismsRegions to make a match)
 *     K-mer counting => ...
 */

#include "align/alignment.h"
#include "text/dna_text.h"
#include "align/align_kmer_filter.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PLOW

/*
 * Constants
 */
#define KMER_COUNTING_MASK_3   0x000000000000003Full /* 256 B */
#define KMER_COUNTING_MASK_4   0x00000000000000FFull /*  1 KB */
#define KMER_COUNTING_MASK_5   0x00000000000003FFull /*  4 KB */
#define KMER_COUNTING_MASK_6   0x0000000000000FFFull /* 16 KB */
#define KMER_COUNTING_MASK_7   0x0000000000003FFFull /* 64 KB */

#define KMER_COUNTING_MASK_INDEX(kmer_idx) ((kmer_idx) & kmer_counting->kmer_mask)
#define KMER_COUNTING_ADD_INDEX(kmer_idx,enc_char) kmer_idx = (kmer_idx<<2 | (enc_char))
#define KMER_COUNTING_ADD_INDEX__MASK(kmer_idx,enc_char) kmer_idx = KMER_COUNTING_MASK_INDEX(kmer_idx<<2 | (enc_char))

/*
 * Compile Pattern
 */
void kmer_counting_compile(
    kmer_counting_t* const kmer_counting,
    uint8_t* const key,
    const uint64_t key_length,
    const uint64_t kmer_length,
    mm_stack_t* const mm_stack) {
  // Check kmer length
  kmer_counting->kmer_length = kmer_length;
  switch (kmer_length) {
    case 3: kmer_counting->kmer_mask = KMER_COUNTING_MASK_3; break;
    case 4: kmer_counting->kmer_mask = KMER_COUNTING_MASK_4; break;
    case 5: kmer_counting->kmer_mask = KMER_COUNTING_MASK_5; break;
    case 6: kmer_counting->kmer_mask = KMER_COUNTING_MASK_6; break;
    default:
      gem_warn_msg("K-mer counting. Invalid proposed k-mer length. Set to k=7");
      kmer_counting->kmer_length = 7;
      kmer_counting->kmer_mask = KMER_COUNTING_MASK_7;
      break;
  }
  // Init
  kmer_counting->num_kmers = POW4(kmer_counting->kmer_length);
  kmer_counting->kmer_count_text = mm_stack_calloc(mm_stack,kmer_counting->num_kmers,uint16_t,true);
  kmer_counting->kmer_count_pattern = mm_stack_calloc(mm_stack,kmer_counting->num_kmers,uint16_t,true);
  kmer_counting->pattern_length = key_length;
  // Count kmers in pattern
  uint64_t pos=0, kmer_idx=0, acc=0;
  // Compile all k-mers
  for (;pos<key_length;++pos) {
    const uint8_t enc_char = key[pos];
    if (!is_dna_canonical_encoded(enc_char)) {
      acc=0; kmer_idx=0; // Reset
    } else {
      KMER_COUNTING_ADD_INDEX__MASK(kmer_idx,enc_char);  // Update kmer-index
      if (acc < kmer_counting->kmer_length-1) {
        ++acc; // Inc accumulator
      } else {
        ++(kmer_counting->kmer_count_pattern[kmer_idx]); // Increment kmer-count
      }
    }
  }
}
/*
 * Filter
 */
uint64_t kmer_counting_min_bound(
    kmer_counting_t* const kmer_counting,
    const uint8_t* const text,
    const uint64_t text_length) {
  PROFILE_START(GP_FC_KMER_COUNTER_FILTER,PROFILE_LEVEL);
  // Check min-length condition
  if (kmer_counting->pattern_length < kmer_counting->kmer_length) {
    PROFILE_STOP(GP_FC_KMER_COUNTER_FILTER,PROFILE_LEVEL);
    PROF_INC_COUNTER(GP_FC_KMER_COUNTER_FILTER_NA);
    return 0; // Don't filter
  }
  // Prepare filter
  const uint64_t kmers_max = kmer_counting->pattern_length - (kmer_counting->kmer_length-1);
  memset(kmer_counting->kmer_count_text,0,kmer_counting->num_kmers*UINT16_SIZE);
  uint16_t* const kmer_count_text = kmer_counting->kmer_count_text;
  uint16_t* const kmer_count_pattern = kmer_counting->kmer_count_pattern;
  /*
   * First count (Load)
   */
  const uint64_t total_kmers_text = MAX(text_length,kmer_counting->pattern_length);
  const uint64_t init_chunk = MIN(text_length,kmer_counting->pattern_length);
  uint64_t kmers_left = total_kmers_text, kmers_in_text = 0, max_kmers_in_text = 0;
  uint64_t begin_pos, kmer_idx_begin = 0, acc_begin = 0;
  uint64_t end_pos, kmer_idx_end=0, acc_end = 0;
  // Initial window fill
  for (end_pos=0;end_pos<init_chunk;++end_pos,--kmers_left) {
    const uint8_t enc_char_end = text[end_pos];
    if (!is_dna_canonical_encoded(enc_char_end)) {
      acc_end = 0; kmer_idx_end = 0; // Reset
    } else {
      KMER_COUNTING_ADD_INDEX__MASK(kmer_idx_end,enc_char_end); // Update kmer-index
      if (acc_end < kmer_counting->kmer_length-1) {
        ++acc_end;
      } else {
        // Increment kmer-count
        const uint16_t count_pattern = kmer_count_pattern[kmer_idx_end];
        if (count_pattern>0 && (kmer_count_text[kmer_idx_end])++ < count_pattern) {
          ++kmers_in_text;
          max_kmers_in_text = MAX(max_kmers_in_text,kmers_in_text);
        }
      }
    }
  }
  /*
   * Sliding window count
   */
  // Initial fill
  for (begin_pos=0;begin_pos<kmer_counting->kmer_length-1;++begin_pos) {
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
      if (acc_begin < kmer_counting->kmer_length-1) {
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
      if (acc_end < kmer_counting->kmer_length-1) {
        ++acc_end;
      } else {
        const uint16_t count_pattern_end = kmer_count_pattern[kmer_idx_end];
        if (count_pattern_end > 0 && kmer_count_text[kmer_idx_end]++ < count_pattern_end) {
          ++kmers_in_text;
          max_kmers_in_text = MAX(max_kmers_in_text,kmers_in_text);
        }
      }
    }
  }
  // Compute min-error bound
  PROFILE_STOP(GP_FC_KMER_COUNTER_FILTER,PROFILE_LEVEL);
  const uint64_t min_error = (max_kmers_in_text < kmers_max) ?
      (kmers_max-max_kmers_in_text)/kmer_counting->kmer_length : 0;
  return min_error; // Filter out
}
bool kmer_counting_filter(
    kmer_counting_t* const kmer_counting,
    const uint8_t* const text,
    const uint64_t text_length,
    const uint64_t max_error) {
  PROFILE_START(GP_FC_KMER_COUNTER_FILTER,PROFILE_LEVEL);
  // Check min-length condition
  if (kmer_counting->pattern_length < kmer_counting->kmer_length) {
    PROFILE_STOP(GP_FC_KMER_COUNTER_FILTER,PROFILE_LEVEL);
    PROF_INC_COUNTER(GP_FC_KMER_COUNTER_FILTER_NA);
    return true; // Accept
  }
  // Check max-error condition
  const uint64_t kmers_error = kmer_counting->kmer_length * max_error;
  const uint64_t kmers_max = kmer_counting->pattern_length - (kmer_counting->kmer_length-1);
  if (kmers_error >= kmers_max) {
    PROFILE_STOP(GP_FC_KMER_COUNTER_FILTER,PROFILE_LEVEL);
    PROF_INC_COUNTER(GP_FC_KMER_COUNTER_FILTER_NA);
    return true; // Accept
  }
  // Prepare filter
  const uint64_t kmers_required = kmers_max - kmer_counting->kmer_length*max_error;
  memset(kmer_counting->kmer_count_text,0,kmer_counting->num_kmers*UINT16_SIZE);
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
      if (acc_end < kmer_counting->kmer_length-1) {
        ++acc_end;
      } else {
        // Increment kmer-count
        const uint16_t count_pattern = kmer_count_pattern[kmer_idx_end];
        if (count_pattern>0 && (kmer_count_text[kmer_idx_end])++ < count_pattern) ++kmers_in_text;
        // Check filter condition
        if (kmers_in_text >= kmers_required) {
          PROFILE_STOP(GP_FC_KMER_COUNTER_FILTER,PROFILE_LEVEL);
          return true;  // Accept
        } else if (kmers_required-kmers_in_text > kmers_left) { // Quick abandon
          PROFILE_STOP(GP_FC_KMER_COUNTER_FILTER,PROFILE_LEVEL);
          return false; // Filter out
        }
      }
    }
  }
  /*
   * Sliding window count
   */
  // Initial fill
  for (begin_pos=0;begin_pos<kmer_counting->kmer_length-1;++begin_pos) {
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
      if (acc_begin < kmer_counting->kmer_length-1) {
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
      if (acc_end < kmer_counting->kmer_length-1) {
        ++acc_end;
      } else {
        const uint16_t count_pattern_end = kmer_count_pattern[kmer_idx_end];
        if (count_pattern_end > 0 && kmer_count_text[kmer_idx_end]++ < count_pattern_end) ++kmers_in_text;
      }
    }
    // Check filter condition
    if (kmers_in_text >= kmers_required) {
      PROFILE_STOP(GP_FC_KMER_COUNTER_FILTER,PROFILE_LEVEL);
      return true;  // Accept
    } else if (kmers_required-kmers_in_text > kmers_left) { // Quick abandon
      return false; // Filter out
    }
  }
  // Not passing the filter
  PROFILE_STOP(GP_FC_KMER_COUNTER_FILTER,PROFILE_LEVEL);
  return false; // Filter out
}

