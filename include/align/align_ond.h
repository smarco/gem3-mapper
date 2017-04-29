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
 *   Alignment module using O(nd)-algorithm to compute levenshtein distance/alignment
 *   (Myers' O(nd)-algorithm to compute levenshtein distance/alignment)
 */

#ifndef ALIGN_OND_H_
#define ALIGN_OND_H_

#include "utils/essentials.h"
#include "matches/align/match_alignment.h"

typedef struct {
  int32_t** contour;
  uint64_t lcs_distance;
  uint64_t match_end_column;
} align_ond_contours_t;

/*
 * O(ND) Compute LCS
 */
void align_ond_compute_lcs_distance(
    const uint8_t* const key,
    const int32_t key_length,
    const uint8_t* const text,
    const int32_t text_length,
    uint64_t* const lcs_distance,
    uint64_t* const match_end_column,
    mm_allocator_t* const mm_allocator);

/*
 * O(ND) Align
 */
void align_ond_compute_contours(
    const uint8_t* const key,
    const int32_t key_length,
    const uint8_t* const text,
    const int32_t text_length,
    const int32_t max_distance,
    align_ond_contours_t* const align_ond_contours,
    mm_allocator_t* const mm_allocator);
void align_ond_backtrace_contours(
    const uint8_t* const key,
    const int32_t key_length,
    const uint8_t* const text,
    const int32_t text_length,
    align_ond_contours_t* const align_ond_contours,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector);
void align_ond_match(
    const uint8_t* const key,
    const int32_t key_length,
    const uint8_t* const text,
    const int32_t text_length,
    const int32_t max_distance,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector,
    mm_allocator_t* const mm_allocator);

/*
 * Display
 */
void align_ond_print_contour(
    FILE* const stream,
    const int32_t* const contour,
    const int32_t begin_contour,
    const int32_t end_contour,
    const int32_t distance);

#endif /* ALIGN_OND_H_ */
