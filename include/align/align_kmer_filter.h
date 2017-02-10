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
 */

#ifndef ALIGN_KMER_FILTER_H_
#define ALIGN_KMER_FILTER_H_

#include "utils/essentials.h"

/*
 * Kmer counting filter
 */
typedef struct {
  // Filter parameters
  uint64_t kmer_length;
  uint64_t kmer_mask;
  uint64_t num_kmers;
  // Pattern parameters
  uint64_t pattern_length;
  // Tables
  uint16_t* kmer_count_text;
  uint16_t* kmer_count_pattern;
} kmer_counting_t;

/*
 * Setup
 */
void kmer_counting_compile(
    kmer_counting_t* const kmer_counting,
    uint8_t* const key,
    const uint64_t key_length,
    const uint64_t kmer_length,
    mm_stack_t* const mm_stack);

/*
 * Filter
 */
uint64_t kmer_counting_min_bound(
    kmer_counting_t* const kmer_counting,
    const uint8_t* const text,
    const uint64_t text_length);
bool kmer_counting_filter(
    kmer_counting_t* const kmer_counting,
    const uint8_t* const text,
    const uint64_t text_length,
    const uint64_t max_error);

#endif /* ALIGN_KMER_FILTER_H_ */
