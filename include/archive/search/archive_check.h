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
 *   Archive checker module providing basic functions to verify correctness
 *   and completeness of archive-searches (SE/PE)
 */

#ifndef ARCHIVE_CHECK_H_
#define ARCHIVE_CHECK_H_

#include "archive/archive.h"
#include "archive/search/archive_search_se_parameters.h"
#include "text/sequence.h"
#include "matches/align/match_align.h"
#include "matches/paired_matches.h"

/*
 * Check Matches
 */
void archive_check_se_matches(
    archive_t* const archive,
    const match_alignment_model_t match_alignment_model,
    swg_penalties_t* swg_penalties,
    sequence_t* const sequence,
    matches_t* const matches,
    const archive_check_type check_type,
    mm_allocator_t* const mm_allocator);
void archive_check_pe_matches(
    archive_t* const archive,
    const match_alignment_model_t match_alignment_model,
    swg_penalties_t* swg_penalties,
    sequence_t* const sequence_end1,
    sequence_t* const sequence_end2,
    paired_matches_t* const paired_matches,
    const archive_check_type check_type,
    mm_allocator_t* const mm_allocator);

#endif /* ARCHIVE_CHECK_H_ */
