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
 *   Archive-Search Paired-End module encapsulating basic
 *   PE-search stages
 */

#ifndef ARCHIVE_SEARCH_PE_STAGES_H_
#define ARCHIVE_SEARCH_PE_STAGES_H_

#include "utils/essentials.h"
#include "archive/search/archive_search.h"
#include "matches/paired_matches.h"

/*
 * Archive Search PE Stages
 */
void archive_search_pe_begin(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2);

void archive_search_pe_search_end1(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches);
void archive_search_pe_search_end2(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches);

bool archive_search_pe_use_shortcut_extension(
    archive_search_t* const archive_search_extended,
    archive_search_t* const archive_search_candidate,
    matches_t* const matches);

void archive_search_pe_find_pairs(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches);

void archive_search_pe_end(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches);

#endif /* ARCHIVE_SEARCH_PE_STAGES_H_ */
