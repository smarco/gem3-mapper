/*
 * PROJECT: GEMMapper
 * FILE: matches_cigar.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Data structure to store alignment matches {sequence,position,strand,...}
 */

#ifndef MATCHES_CIGAR_H_
#define MATCHES_CIGAR_H_

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
 * CIGAR Buffer Handling
 */
void matches_cigar_buffer_add_cigar_element(
    cigar_element_t** const cigar_buffer_sentinel,
    const cigar_t cigar_element_type,
    const uint64_t element_length);
void matches_cigar_buffer_add_cigar_element__attr(
    cigar_element_t** const cigar_buffer_sentinel,
    const cigar_t cigar_element_type,
    const uint64_t element_length,
    const cigar_attr_t attributes);
void matches_cigar_buffer_add_mismatch(
    cigar_element_t** const cigar_buffer_sentinel,
    const uint8_t mismatch);

/*
 * CIGAR Vector Handling
 */
void matches_cigar_vector_append_insertion(
    vector_t* const cigar_vector,
    uint64_t* const current_cigar_length,
    const uint64_t indel_length,
    const cigar_attr_t attributes);
void matches_cigar_vector_append_deletion(
    vector_t* const cigar_vector,
    uint64_t* const current_cigar_length,
    const uint64_t indel_length,
    const cigar_attr_t attributes);
void matches_cigar_vector_append_match(
    vector_t* const cigar_vector,
    uint64_t* const current_cigar_length,
    const uint64_t match_length,
    const cigar_attr_t attributes);
void matches_cigar_vector_append_mismatch(
    vector_t* const cigar_vector,
    uint64_t* const current_cigar_length,
    const uint8_t mismatch,
    const cigar_attr_t attributes);
void matches_cigar_vector_append_cigar_element(
    vector_t* const cigar_vector,
    uint64_t* const cigar_length,
    cigar_element_t* const cigar_element);

/*
 * CIGAR Vector Utils
 */
void matches_cigar_reverse(
    vector_t* const cigar_vector,
    const uint64_t cigar_buffer_offset,
    const uint64_t cigar_length);
void matches_cigar_reverse_colorspace(
    vector_t* const cigar_vector,
    const uint64_t cigar_buffer_offset,
    const uint64_t cigar_length);

uint64_t matches_cigar_compute_event_distance(
    vector_t* const cigar_vector,
    const uint64_t cigar_buffer_offset,
    const uint64_t cigar_length);
uint64_t matches_cigar_compute_edit_distance(
    vector_t* const cigar_vector,
    const uint64_t cigar_buffer_offset,
    const uint64_t cigar_length);
uint64_t matches_cigar_compute_edit_distance__excluding_clipping(
    vector_t* const cigar_vector,
    const uint64_t cigar_buffer_offset,
    const uint64_t cigar_length);
uint64_t matches_cigar_compute_matching_bases(
    vector_t* const cigar_vector,
    const uint64_t cigar_buffer_offset,
    const uint64_t cigar_length);

int64_t matches_cigar_element_effective_length(
    const cigar_element_t* const cigar_element);
int64_t matches_cigar_effective_length(
    vector_t* const cigar_vector,
    const uint64_t cigar_offset,
    const uint64_t cigar_length);

/*
 * CIGAR Vector Compare
 */
int matches_cigar_cmp(
    vector_t* const cigar0_vector,
    const uint64_t cigar0_offset,
    const uint64_t cigar0_length,
    vector_t* const cigar1_vector,
    const uint64_t cigar1_offset,
    const uint64_t cigar1_length);

/*
 * Display
 */
void match_cigar_print(
    FILE* const stream,
    vector_t* const cigar_vector,
    const uint64_t cigar_buffer_offset,
    const uint64_t cigar_length);

#endif /* MATCHES_CIGAR_H_ */
