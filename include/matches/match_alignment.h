/*
 * PROJECT: GEMMapper
 * FILE: match_alignment.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MATCH_ALIGNMENT_H_
#define MATCH_ALIGNMENT_H_

#include "utils/essentials.h"

/*
 * Alignment CIGAR (Mismatches/Indels/...)
 */
typedef enum {
  cigar_null = 0,
  cigar_match = 1,
  cigar_mismatch = 2,
  cigar_ins = 3,
  cigar_del = 4,
} cigar_t;
typedef enum {
  cigar_attr_none = 0,
  cigar_attr_trim = 1,
  cigar_attr_homopolymer = 2,
} cigar_attr_t;
typedef struct {
  cigar_t type;              // Match, Mismatch, insertion or deletion
  cigar_attr_t attributes;   // Attributes
  union {
    int32_t length;          // Match length
    uint8_t mismatch;        // Mismatch base
  };
} cigar_element_t;
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
/*
 * Region Alignment
 */
typedef struct {
  uint64_t match_distance;                  // Distance
  uint64_t text_begin_offset;               // Text begin offset
  uint64_t text_end_offset;                 // Text end offset
} region_alignment_tile_t;
typedef struct {
  uint64_t num_tiles;                       // Total number of tiles
  uint64_t distance_min_bound;              // Distance min-bound (Sum all tile distances)
  region_alignment_tile_t* alignment_tiles; // Alignment of all tiles
} region_alignment_t;
/*
 * Match Alignment
 */
typedef struct {
  uint64_t match_text_offset; // Match text offset (wrt beginning of text-candidate)
  uint64_t match_position;    // Match position
  uint64_t cigar_offset;      // CIGAR offset in buffer
  uint64_t cigar_length;      // CIGAR length
  int64_t effective_length;   // Match effective length
  int32_t score;              // Score assigned by the aligner
} match_alignment_t;
/*
 * Alignment Model
 */
typedef enum {
  alignment_model_none,
  alignment_model_hamming,
  alignment_model_levenshtein,
  alignment_model_gap_affine
} alignment_model_t;

/*
 * Region Matching
 */
uint64_t region_matching_text_coverage(region_matching_t* const region_matching);
uint64_t region_matching_text_distance(
    region_matching_t* const region_matching_a,
    region_matching_t* const region_matching_b);
bool region_matching_text_overlap(
    region_matching_t* const region_matching_a,
    region_matching_t* const region_matching_b);
int region_matching_key_cmp(
    region_matching_t* const region_matching_a,
    region_matching_t* const region_matching_b);

/*
 * Display
 */
void match_alignment_print_pretty(
    FILE* const stream,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector,
    uint8_t* const key,
    const uint64_t key_length,
    uint8_t* const text,
    const uint64_t text_length,
    mm_stack_t* const mm_stack);

#endif /* MATCH_ALIGNMENT_H_ */
