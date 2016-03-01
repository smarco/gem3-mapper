/*
 * PROJECT: GEMMapper
 * FILE: kmer_counting.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef KMER_COUNTING_H_
#define KMER_COUNTING_H_

#include "utils/essentials.h"

/*
 * Kmer counting filter
 */
typedef struct {
  bool enabled;
  uint16_t* kmer_count_text;
  uint16_t* kmer_count_pattern;
  uint64_t pattern_length;
  uint64_t max_error;
} kmer_counting_t;

/*
 * Compile Pattern
 */
void kmer_counting_compile(
    kmer_counting_t* const kmer_counting,
    uint8_t* const pattern,
    const uint64_t pattern_length,
    const uint64_t max_error,
    mm_stack_t* const mm_stack);

/*
 * Filter text region
 */
uint64_t kmer_counting_filter(
    const kmer_counting_t* const kmer_counting,
    const uint8_t* const text,
    const uint64_t text_length);

#endif /* KMER_COUNTING_H_ */
