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

#ifndef NSEARCH_HAMMING_H_
#define NSEARCH_HAMMING_H_

#include "utils/essentials.h"
#include "fm_index/fm_index.h"
#include "approximate_search/approximate_search.h"
#include "filtering/region_profile/region_profile.h"
#include "neighborhood_search/nsearch_schedule.h"

/*
 * Hamming Brute Force
 */
void nsearch_hamming_brute_force(
    approximate_search_t* const search,
    matches_t* const matches);

/*
 * Perform scheduled search
 */
void nsearch_hamming_scheduled_search(
    nsearch_schedule_t* const nsearch_schedule);

/*
 * Hamming Neighborhood Search
 */
void nsearch_hamming(
    approximate_search_t* const search,
    matches_t* const matches);
void nsearch_hamming_preconditioned(
    approximate_search_t* const search,
    matches_t* const matches);

#endif /* NSEARCH_HAMMING_H_ */
