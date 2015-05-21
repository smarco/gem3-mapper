/*
 * PROJECT: GEMMapper
 * FILE: matches_elements.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MATCH_ELEMENTS_H_
#define MATCH_ELEMENTS_H_

#include "essentials.h"

/*
 * Checkers
 */
#define MATCH_TRACE_CHECK(match) GEM_CHECK_NULL(match)

/*
 * Alignment CIGAR (Mismatches/Indels/...)
 */
typedef enum { cigar_match=0, cigar_mismatch=1, cigar_ins=2, cigar_del=3, cigar_soft_trim=4, cigar_null=5 } cigar_t;
typedef struct {
  int32_t indel_length;   // Indel length
  uint8_t* indel_text;    // Indel-text (Reference text in case of insertion, NULL otherwise)
} cigar_indel_t;
typedef struct {
  cigar_t type;           // Match, Mismatch, insertion or deletion
  union {
    int32_t match_length; // Match length
    uint8_t mismatch;     // Mismatch base
    cigar_indel_t indel;  // Pointer to the insert
  };
} cigar_element_t;
typedef struct {
  uint64_t match_position;
  uint64_t cigar_offset;
  uint64_t cigar_length;
  int64_t effective_length;
  int32_t score;
} match_alignment_t;

/*
 * Region Matching
 *   A chain of regions matching determines the shape of the final
 *   alignment CIGAR. If allowed, the chain of matching regions will
 *   be followed as to rule out regions (key<->text) that match on the high level
 *   This regions can be:
 *     1. Exact-matches (identically mapping)
 *     2. Approximate-matching regions (mapping with a precise indel/mismatch structure => @cigar_buffer_offset)
 *     3. Clipped-regions (excluded from being aligned)
 *        Eg. Adaptors, bad-quality regions, local alignments
 *
 */
typedef enum {
  region_matching_exact,
  region_matching_approximate,
  region_clipped
} region_matching_type;
typedef struct {
  /* Region matching type */
  region_matching_type matching_type;
  /* Error of the region matching */
  uint64_t error;
  uint64_t cigar_buffer_offset;
  uint64_t cigar_length;
  /* Coordinates of the region */
  uint64_t key_begin;
  uint64_t key_end;
  uint64_t text_begin;
  uint64_t text_end;
} region_matching_t;

#endif /* MATCH_ELEMENTS_H_ */
