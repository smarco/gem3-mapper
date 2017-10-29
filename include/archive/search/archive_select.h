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
 *   Archive select module provides functions to select and discard
 *   matches according with the user specified parameters
 */

#ifndef ARCHIVE_SELECT_H_
#define ARCHIVE_SELECT_H_

#include "archive/search/archive_search.h"
#include "archive/search/archive_select_parameters.h"
#include "matches/matches.h"
#include "matches/paired_matches.h"

/*
 * Setup
 */
void archive_select_configure_se(archive_search_t* const archive_search);
void archive_select_configure_pe(archive_search_t* const archive_search);

/*
 * Select Paired-Matches
 */
void archive_select_se_matches(
    select_parameters_t* const select_parameters,
    matches_t* const matches);
void archive_select_pe_matches(
    select_parameters_t* const select_parameters,
    mapper_stats_t* const mapper_stats,
    paired_matches_t* const paired_matches);

#endif /* ARCHIVE_SELECT_H_ */
