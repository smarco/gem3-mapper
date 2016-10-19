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
 *   Mapper-profile module provides functions to display the
 *   GEM-profiling counters involved in the mapping search
 */

#ifndef MAPPER_PROFILE_SEARCH_H_
#define MAPPER_PROFILE_SEARCH_H_

#include "utils/essentials.h"
#include "mapper/mapper_profile_counters.h"

/*
 * Region Profile
 */
void mapper_profile_print_region_profile(FILE* const stream);

/*
 * Candidates Generation
 */
void mapper_profile_print_candidate_generation(FILE* const stream);

/*
 * Candidate Verification
 */
void mapper_profile_print_candidate_verification(FILE* const stream);

/*
 * Candidate realign
 */
void mapper_profile_print_candidate_realign(FILE* const stream);
void mapper_profile_print_candidate_realign_local(FILE* const stream);

/*
 * Candidate Scaffold
 */
void mapper_profile_print_candidate_scaffold(FILE* const stream);

/*
 * Candidate misc
 */
void mapper_profile_print_candidate_misc(FILE* const stream);

/*
 * Neighborhood Search
 */
void mapper_profile_print_neighborhood_search(FILE* const stream);

/*
 * Approximate Search
 */
void mapper_profile_print_approximate_search(FILE* const stream);

/*
 * Approximate Search Profile Summary
 */
void mapper_profile_print_approximate_search_summary(FILE* const stream);

#endif /* MAPPER_PROFILE_SEARCH_H_ */
