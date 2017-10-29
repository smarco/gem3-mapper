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

#include "neighborhood_search/nsearch_filtering.h"
#include "neighborhood_search/nsearch_levenshtein_control.h"
#include "filtering/candidates/filtering_candidates_process.h"
#include "filtering/candidates/filtering_candidates_verify.h"
#include "filtering/candidates/filtering_candidates_align.h"
#include "matches/classify/matches_classify.h"

/*
 * NS Filter candidates
 */
void nsearch_filtering(nsearch_schedule_t* const nsearch_schedule) {
#ifndef NSEARCH_ENUMERATE
  // Parameters
  filtering_candidates_t* const filtering_candidates = nsearch_schedule->filtering_candidates;
  pattern_t* const pattern = nsearch_schedule->pattern;
  matches_t* const matches = nsearch_schedule->matches;
  // Quick check to avoid further processing
  if (filtering_candidates_get_num_positions(filtering_candidates) > 0) {
    // Process+Verify candidates
    PROF_START(GP_NS_VERIFICATION);
    filtering_candidates_process_candidates(filtering_candidates,pattern,false);
    filtering_candidates_verify_candidates(filtering_candidates,pattern,matches);
    PROF_STOP(GP_NS_VERIFICATION);
    // Align
    PROF_START(GP_NS_ALIGN);
    filtering_candidates_align_candidates(filtering_candidates,pattern,matches);
    PROF_STOP(GP_NS_ALIGN);
    // Check search-depth
    search_parameters_t* const search_parameters = nsearch_schedule->search_parameters;
    nsearch_parameters_t* const nsearch_parameters = &search_parameters->nsearch_parameters;
    if (nsearch_parameters->matches_accuracy_cutoff) {
      matches_classify(matches);
      nsearch_schedule->quick_abandon = matches_classify_min_depth_searched(&matches->classification);
    }
    // Check max-matches
    if (!nsearch_schedule->quick_abandon && nsearch_parameters->matches_max_searched_cutoff) {
      nsearch_schedule->quick_abandon =
          matches_classify_max_matches_searched(matches,search_parameters,pattern);
    }
  }
#endif
}
