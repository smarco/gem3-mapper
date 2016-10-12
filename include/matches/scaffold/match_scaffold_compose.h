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

#ifndef MATCH_SCAFFOLD_COMPOSE_H_
#define MATCH_SCAFFOLD_COMPOSE_H_

#include "utils/essentials.h"
#include "matches/scaffold/match_scaffold.h"
#include "matches/align/match_alignment.h"

/*
 * Add alignment-region
 */
void match_scaffold_compose_add_approximate_match(
    match_scaffold_t* const match_scaffold,
    const uint64_t key_begin,
    const uint64_t key_end,
    const uint64_t text_begin,
    const uint64_t text_end,
    const uint64_t cigar_offset,
    const uint64_t cigar_length,
    const uint64_t score);
match_alignment_region_t* match_scaffold_compose_add_exact_match(
    match_scaffold_t* const match_scaffold,
    uint64_t* const key_offset,
    uint64_t* const text_offset,
    const uint64_t cigar_offset,
    const uint64_t match_length);
void match_scaffold_compose_add_mismatch(
    match_scaffold_t* const match_scaffold,
    match_alignment_region_t* const last_alignment_region,
    uint64_t* const key_offset,
    uint64_t* const text_offset,
    const uint64_t match_length);

#endif /* MATCH_SCAFFOLD_COMPOSE_H_ */
