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
 *   Region-Profile module provides functions to generate an key partition
 *   based of the entropy of the regions. An adaptive profile determines a
 *   key partition into regions that have few matches in the index.
 *     Note that if the algorithm cannot find any region
 *     could be due to the following reasons
 *       - There are wildcards which prevents region generation
 *       - There are too many exact matches (preventing unique regions)
 */

#ifndef REGION_PROFILE_ADAPTIVE_H_
#define REGION_PROFILE_ADAPTIVE_H_

#include "utils/essentials.h"
#include "filtering/region_profile/region_profile.h"

/*
 * Region Profile Generator & Query
 */
typedef struct {
  // Region Profile
  region_profile_t* region_profile;
  // Region State
  uint64_t last_cut;
  uint64_t lo_cut;
  uint64_t hi_cut;
  uint64_t expected_count;
  uint64_t max_steps;
  bool allow_zero_regions;
  // Query
  fm_index_t* fm_index;
  const uint8_t* key;
  uint64_t key_length;
  // Query state
  uint64_t key_position;
  rank_mquery_t rank_mquery;
  uint64_t lo;
  uint64_t hi;
} region_profile_generator_t;

/*
 * Region Profile Adaptive Iterator
 */
void region_profile_generator_init(
    region_profile_generator_t* const generator,
    region_profile_t* const region_profile,
    fm_index_t* const fm_index,
    const uint8_t* const key,
    const uint64_t key_length,
    const bool allow_zero_regions);
bool region_profile_generator_next_region(
    region_profile_t* const region_profile,
    region_profile_generator_t* const generator,
    const region_profile_model_t* const profile_model);

/*
 * Region Profile Generation
 */
void region_profile_generate_adaptive(
    region_profile_t* const region_profile,
    fm_index_t* const fm_index,
    const uint8_t* const key,
    const uint64_t key_length,
    const region_profile_model_t* const profile_model,
    const uint64_t max_regions,
    const bool allow_zero_regions);
void region_profile_generate_adaptive_limited(
    region_profile_t* const region_profile,
    fm_index_t* const fm_index,
    const uint8_t* const key,
    const uint64_t key_length,
    const region_profile_model_t* const profile_model);

#endif /* REGION_PROFILE_ADAPTIVE_H_ */
