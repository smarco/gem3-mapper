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
 *   Filtering module provides functions to produce local-alignments
 *   from the discarded filtering-regions and its match alignment-regions
 */

#include "align/alignment.h"
#include "filtering/candidates/filtering_candidates_align_local.h"
#include "filtering/candidates/filtering_candidates_align.h"
#include "filtering/candidates/filtering_candidates_verify.h"
#include "filtering/region/filtering_region_verify.h"

/*
 * Debug
 */
#define DEBUG_FILTERING_CANDIDATES  GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Filtering Candidates (Re)Alignment
 */
void filtering_candidates_align_local_discarded(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    matches_t* const matches) {
  // Parameters
  search_parameters_t* const search_parameters = filtering_candidates->search_parameters;
  select_parameters_t* const select_parameters = &search_parameters->select_parameters;
  const uint64_t max_searched_matches = select_parameters->max_searched_matches;
  // Exact extend all discarded regions
  filtering_candidates_extend_discarded_candidates(filtering_candidates,pattern,matches);
  // Local-align the most promising regions
  filtering_region_t** const regions_discarded = filtering_candidates_get_discarded_regions(filtering_candidates);
  const uint64_t num_regions_discarded = filtering_candidates_get_num_discarded_regions(filtering_candidates);
  uint64_t i;
  for (i=0;i<num_regions_discarded;++i) {
    // Cut-off max-reported matches
    filtering_region_t* const filtering_region = regions_discarded[i];
    if (matches_get_num_match_traces(matches) >= max_searched_matches) break;
    // Align Region
    PROF_INC_COUNTER(GP_CANDIDATE_REGION_LOCAL_ALIGNED);
    filtering_candidates_align_local_region(
        filtering_candidates,filtering_region,pattern,matches);
  }
  // Clear discarded-candidates
  filtering_candidates_clear_discarded_regions(filtering_candidates,true);
  PROF_ADD_COUNTER(GP_CANDIDATE_REGION_LOCAL,num_regions_discarded);
}
void filtering_candidates_align_local(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    matches_t* const matches) {
  PROFILE_START(GP_FC_REALIGN_LOCAL_CANDIDATE_REGIONS,PROFILE_LEVEL);
  // Add pending local matches (found so far)
  locator_t* const locator = filtering_candidates->archive->locator;
  matches_local_pending_add_to_regular_matches(matches,locator);
  // Check total alignments found
  search_parameters_t* const search_parameters = filtering_candidates->search_parameters;
  select_parameters_t* const select_parameters = &search_parameters->select_parameters;
  const uint64_t max_searched_matches = select_parameters->max_searched_matches;
  uint64_t total_matches = matches_get_num_match_traces(matches);
  if (total_matches >= max_searched_matches) {
    PROFILE_STOP(GP_FC_REALIGN_LOCAL_CANDIDATE_REGIONS,PROFILE_LEVEL);
    return;
  }
  // Local align all discarded regions
  filtering_candidates_align_local_discarded(filtering_candidates,pattern,matches);
  // DEBUG
  gem_cond_debug_block(DEBUG_FILTERING_CANDIDATES) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Filtering.Candidates (local_align)\n");
    tab_global_inc();
    filtering_candidates_print_regions(gem_log_get_stream(),filtering_candidates,false);
    tab_global_dec();
  }
  PROFILE_STOP(GP_FC_REALIGN_LOCAL_CANDIDATE_REGIONS,PROFILE_LEVEL);
}
