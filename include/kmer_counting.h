/*
 * PROJECT: GEMMapper
 * FILE: kmer_counting.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef KMER_COUNTING_H_
#define KMER_COUNTING_H_

#include "essentials.h"

/*
 * Kmer counting filter
 */
typedef struct {
  // Kmer count
  uint16_t* kmer_count_text;
  uint16_t* kmer_count_pattern;
  // Constraints
  uint64_t pattern_length;
  bool* allowed_enc;
} kmer_counting_t;

/*
 * Compile Pattern
 */
void kmer_counting_compile(
    kmer_counting_t* const kmer_counting,bool* const allowed_enc,
    uint8_t* const pattern,const uint64_t pattern_length,mm_stack_t* const mm_stack);

/*
 * Filter text region
 */
uint64_t kmer_counting_filter(
    const kmer_counting_t* const kmer_counting,
    const uint8_t* const text,const uint64_t text_length,const uint64_t max_error);

#endif /* KMER_COUNTING_H_ */
