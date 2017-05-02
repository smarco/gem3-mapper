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

#ifndef MATCH_SCAFFOLD_LEVENSHTEIN_H_
#define MATCH_SCAFFOLD_LEVENSHTEIN_H_

#include "utils/essentials.h"
#include "matches/scaffold/match_scaffold.h"
#include "matches/matches.h"

/*
 * Setup
 */
void match_scaffold_levenshtein_allocate(
    match_scaffold_t* const match_scaffold,
    const uint64_t key_length,
    const uint64_t matching_min_length,
    mm_allocator_t* const mm_allocator);

/*
 * Compose the scaffolding
 */
void match_scaffold_levenshtein_compose_alignment(
    match_scaffold_t* const match_scaffold,
    const match_alignment_t* const match_alignment,
    uint64_t key_offset,
    const uint64_t matching_min_length,
    matches_t* const matches);

/*
 * Levenshtein Scaffold Tiled
 */
void match_scaffold_levenshtein_tile(
    match_scaffold_t* const match_scaffold,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    alignment_tile_t* const alignment_tile,
    pattern_tile_t* const pattern_tile,
    const uint64_t matching_min_length,
    matches_t* const matches,
    mm_allocator_t* const mm_allocator);

/*
 * Levenshtein Scaffold
 */
bool match_scaffold_levenshtein(
    match_scaffold_t* const match_scaffold,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    alignment_t* const alignment,
    const uint64_t matching_min_length,
    matches_t* const matches,
    mm_allocator_t* const mm_allocator);

#endif /* MATCH_SCAFFOLD_LEVENSHTEIN_H_ */
