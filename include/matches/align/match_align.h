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

#ifndef MATCH_ALIGN_H_
#define MATCH_ALIGN_H_

#include "utils/essentials.h"
#include "matches/matches.h"
#include "matches/align/match_alignment.h"
#include "matches/scaffold/match_scaffold.h"

/*
 * Match Clipping
 */
void match_aling_add_clipping(
    match_trace_t* const match_trace,
    vector_t* const cigar_vector,
    const uint64_t sequence_clip_left,
    const uint64_t sequence_clip_right);

/*
 * Exact Alignment
 */
void match_align_exact(
    matches_t* const matches,
    match_trace_t* const match_trace,
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters);

/*
 * Pseudo-Alignment
 */
void match_align_pseudoalignment(
    matches_t* const matches,
    match_trace_t* const match_trace,
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters);

/*
 * Hamming Alignment (Only mismatches)
 */
void match_align_hamming(
    matches_t* const matches,
    match_trace_t* const match_trace,
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters);

/*
 * Levenshtein Alignment (edit distance)
 */
void match_align_levenshtein(
    matches_t* const matches,
    match_trace_t* const match_trace,
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    mm_stack_t* const mm_stack);

/*
 * SWG Alignment (Gap-affine)
 */
void match_align_smith_waterman_gotoh(
    matches_t* const matches,
    match_trace_t* const match_trace,
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    match_scaffold_t* const match_scaffold,
    mm_stack_t* const mm_stack);

#endif /* MATCH_ALIGN_H_ */
