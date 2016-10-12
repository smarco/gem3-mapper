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
 *   Smith-Waterman-Gotoh (SWG) scoring module using custom SWG-penalties
 *   Provides functions to handle scores & score alignment CIGARs
 */

#ifndef ALIGN_SWG_SCORE_H_
#define ALIGN_SWG_SCORE_H_

#include "utils/essentials.h"
#include "text/dna_text.h"
#include "matches/align/match_alignment.h"
#include "matches/matches_cigar.h"

/*
 * Constants
 */
#define SWG_SCORE_MIN (INT16_MIN)

/*
 * SWG Penalties
 */
typedef int32_t swg_matching_score_t[DNA__N_RANGE][DNA__N_RANGE];
typedef struct {
  int32_t gap_open_score;
  int32_t gap_extension_score;
  int32_t generic_match_score;
  int32_t generic_mismatch_score;
  swg_matching_score_t matching_score;
} swg_penalties_t;

/*
 * SWG Score
 */
int32_t align_swg_score_deletion(const swg_penalties_t* const swg_penalties,const int32_t length);
int32_t align_swg_score_insertion(const swg_penalties_t* const swg_penalties,const int32_t length);
int32_t align_swg_score_mismatch(const swg_penalties_t* const swg_penalties);
int32_t align_swg_score_match(const swg_penalties_t* const swg_penalties,const int32_t match_length);

int32_t align_swg_score_cigar_element(
    const swg_penalties_t* const swg_penalties,
    const cigar_element_t* const cigar_element);
int32_t align_swg_score_cigar(
    const swg_penalties_t* const swg_penalties,
    vector_t* const cigar_vector,
    const uint64_t cigar_offset,
    const uint64_t cigar_length);
int32_t align_swg_score_cigar_excluding_deletions(
    const swg_penalties_t* const swg_penalties,
    vector_t* const cigar_vector,
    const uint64_t cigar_offset,
    const uint64_t cigar_length);
int32_t align_swg_score_cigar_excluding_clipping(
    const swg_penalties_t* const swg_penalties,
    vector_t* const cigar_vector,
    const uint64_t cigar_offset,
    const uint64_t cigar_length);

/*
 * Bounding scores
 */
int32_t align_swg_score_compute_min_score_bound(
    const swg_penalties_t* const swg_penalties,
    const uint64_t edit_distance,
    const uint64_t key_length);
int32_t align_swg_score_compute_max_score_bound(
    const swg_penalties_t* const swg_penalties,
    const uint64_t edit_distance,
    const uint64_t key_length);

/*
 * Bounding edit distance
 */
int32_t align_swg_score_compute_max_edit_bound(
    const swg_penalties_t* const swg_penalties,
    const uint64_t swg_score,
    const uint64_t key_length);

#endif /* ALIGN_SWG_SCORE_H_ */
