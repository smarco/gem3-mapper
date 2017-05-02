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

#ifndef MATCH_SCAFFOLD_CHAIN_H_
#define MATCH_SCAFFOLD_CHAIN_H_

#include "utils/essentials.h"
#include "matches/scaffold/match_scaffold.h"
#include "matches/matches.h"

/*
 * Alignment-region chain scaffolding
 */
void match_scaffold_chain(
    match_scaffold_t* const match_scaffold,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    const bool exact_extend,
    mm_allocator_t* const mm_allocator);

#endif /* MATCH_SCAFFOLD_CHAIN_H_ */
