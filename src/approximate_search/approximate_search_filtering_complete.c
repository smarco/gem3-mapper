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
 *   Approximate-String-Matching (ASM) using modified adaptive-filtering techniques
 *   as to produce complete results (complete AF)
 */

#include "approximate_search/approximate_search_filtering_complete.h"
#include "approximate_search/approximate_search_neighborhood.h"
#include "approximate_search/approximate_search_region_profile.h"
#include "approximate_search/approximate_search_generate_candidates.h"
#include "approximate_search/approximate_search_stages.h"
#include "filtering/region_profile/region_profile_schedule.h"
#include "filtering/candidates/filtering_candidates_process.h"
#include "filtering/candidates/filtering_candidates_verify.h"
#include "filtering/candidates/filtering_candidates_align.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Approximate Complete-Search based on filtering
 */
void approximate_search_filtering_complete(
    approximate_search_t* const search,
    matches_t* const matches) {
  // Parameters
  pattern_t* const pattern = &search->pattern;
  region_profile_t* const region_profile = &search->region_profile;
  // Exact Search
  if (search->current_max_complete_error==0) {
    approximate_search_neighborhood_exact_search(search,matches);
    return;
  }
  // Compute the region profile
  approximate_search_region_profile_adaptive(search,region_profile_adaptive_limited);
  if (search->processing_state==asearch_processing_state_no_regions) return;
  // Generate exact-candidates
  region_profile_schedule_filtering_exact(region_profile);
  approximate_search_generate_candidates_exact(search);
  // Verify candidates
  filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
  filtering_candidates_process_candidates(filtering_candidates,pattern,true);
  filtering_candidates_verify_candidates(filtering_candidates,pattern);
  // Align candidates
  filtering_candidates_align_candidates(filtering_candidates,pattern,false,false,matches);
  // Finish search
  approximate_search_end(search,matches);
}
