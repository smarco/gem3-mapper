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
 *   Region-Profile module provides functions to generate candidate
 *   bwt-positions from a region-profile.
 */

#ifndef REGION_PROFILE_SCHEDULE_H_
#define REGION_PROFILE_SCHEDULE_H_

#include "utils/essentials.h"
#include "filtering/region_profile/region_profile.h"

/*
 * Exact Scheduling
 */
void region_profile_schedule_exact(
    region_profile_t* const region_profile,
    const uint64_t candidates_threshold);
void region_profile_schedule_exact_best(
    region_profile_t* const region_profile,
    const uint64_t num_regions,
    const uint64_t candidates_threshold);

/*
 * Approximate Scheduling
 */
void region_profile_schedule_approximate(
    region_profile_t* const region_profile,
    const uint64_t region_error,
    const uint64_t candidates_threshold);

#endif /* REGION_PROFILE_SCHEDULE_H_ */
