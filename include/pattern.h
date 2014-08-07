/*
 * PROJECT: GEMMapper
 * FILE: pattern.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef PATTERN_H_
#define PATTERN_H_

#include "essentials.h"
#include "bpm_align.h"

// Approximate Search Pattern
typedef struct {
  /* Processed Search Pattern */
  uint8_t* key;            // Encoded Pattern
  uint8_t* quality_mask;   // Quality Mask
  uint64_t key_length;     // Total Length
  /* Pattern Properties */
  uint64_t num_wildcards;
  uint64_t num_low_quality_bases;
  uint64_t max_effective_filtering_error;
  /* Pattern BitVector-Encoded (Myers-DP) */
  bpm_pattern_t bpm_pattern;
} pattern_t;


#endif /* PATTERN_H_ */
