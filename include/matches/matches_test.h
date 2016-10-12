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

#ifndef MATCHES_TEST_H_
#define MATCHES_TEST_H_

#include "utils/essentials.h"
#include "archive/search/archive_search_se_parameters.h"
#include "matches/matches.h"

/*
 * Matches Condition Tests
 */
bool matches_test_max_matches_reached(
    matches_t* const matches,
    const uint64_t mcs,
    const uint64_t key_length,
    search_parameters_t* const search_parameters);
bool matches_test_accuracy_reached(
    matches_t* const matches,
    const uint64_t mcs,
    const uint64_t key_length,
    search_parameters_t* const search_parameters,
    const uint64_t max_complete_error,
    uint64_t* const max_complete_error_required);

#endif /* MATCHES_TEST_H_ */
