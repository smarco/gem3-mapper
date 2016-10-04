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

#include "filtering/region_profile/region_profile_schedule.h"

/*
 * Region Profile Scheduling
 */
void region_profile_schedule_filtering_exact(region_profile_t* const region_profile) {
  // Parameters
  const uint64_t num_regions = region_profile->num_filtering_regions;
  region_search_t* const filtering_region = region_profile->filtering_region;
  // Assign degree-zero to all regions
  uint64_t i;
  for (i=0;i<num_regions;++i) {
    filtering_region[i].degree = REGION_FILTER_DEGREE_ZERO;
  }
}
