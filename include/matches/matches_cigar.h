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
void matches_cigar_vector_insert_cigar_element(
    vector_t* const cigar_vector,
    const int64_t position,
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
int64_t matches_cigar_compute_effective_length(
    vector_t* const cigar_vector,
    const uint64_t cigar_offset,
    const uint64_t cigar_length);

float matches_cigar_compute_error_quality(
    vector_t* const cigar_vector,
    const uint64_t cigar_offset,
    const uint64_t cigar_length,
    uint8_t* const quality_mask,
    const uint64_t quality_mask_length);

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
