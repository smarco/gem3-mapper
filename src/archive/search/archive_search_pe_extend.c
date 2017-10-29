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
 *   Archive-Search Paired-End module encapsulating basic
 *   PE-search stages
 */

#include "archive/search/archive_search_pe_extend.h"
#include "archive/search/archive_select.h"
#include "archive/score/archive_score_se.h"
#include "approximate_search/approximate_search_verify_candidates.h"
#include "matches/classify/matches_classify.h"
#include "matches/classify/paired_matches_classify.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PHIGH

/*
 * Perform extensions
 */
void archive_search_pe_extend_matches_ends(
    archive_search_t* const candidate_archive_search,
    const sequence_end_t candidate_end,
    search_parameters_t* const search_parameters,
    mapper_stats_t* const mapper_stats,
    matches_t* const extended_matches,
    matches_t* const candidate_matches,
    paired_matches_t* const paired_matches) {
  PROF_INC_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_CANDIDATES_TOTAL);
  PROFILE_START(GP_ARCHIVE_SEARCH_PE_EXTEND_CANDIDATES,PROFILE_LEVEL);
  // Parameters
  search_paired_parameters_t* const search_paired_parameters = &search_parameters->search_paired_parameters;
  approximate_search_t* const approximate_search = &candidate_archive_search->approximate_search;
  filtering_candidates_t* const filtering_candidates = approximate_search->filtering_candidates;
  pattern_t* const candidate_pattern = &candidate_archive_search->approximate_search.pattern;
  // Check orientation
  pair_relation_t* const pair_orientation = search_paired_parameters->pair_orientation;
  if (pair_orientation[pair_orientation_FR] != pair_relation_concordant) return;
  // Iterate over all matches of the extended end
  const uint64_t num_extended_match_traces = matches_get_num_match_traces(extended_matches);
  match_trace_t** const extended_match_traces = matches_get_match_traces(extended_matches);
  uint64_t i;
  for (i=0;i<num_extended_match_traces;++i) {
    // Fetch match
    match_trace_t* const extended_match = extended_match_traces[i];
    // Skip extended or subdominant ends
    if (extended_match->type == match_type_extended) continue;
    if (paired_matches_classify_subdominant_end(paired_matches,candidate_matches,extended_match)) continue;
    // Extend (filter nearby region)
    PROF_INC_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_NUM_MATCHES);
    vector_clear(candidate_matches->match_traces_extended);
    approximate_search_verify_extend_candidate(
        filtering_candidates,candidate_pattern,
        extended_match,mapper_stats,
        paired_matches,candidate_end);
    // (Re)Score Matches
    const uint64_t num_matches_result = matches_get_num_match_traces_extended(candidate_matches);
    if (num_matches_result > 0 || candidate_matches->match_replaced) {
      archive_score_matches_se(candidate_archive_search,candidate_matches);
      // Recompute metrics in case updated SE match was paired
      if (candidate_matches->match_replaced) {
        paired_matches_recompute_metrics(paired_matches);
      }
      candidate_matches->match_replaced = false;
    }
    // Cross-Pair extended matches
    match_trace_t** const matches_result = vector_get_mem(candidate_matches->match_traces_extended,match_trace_t*);
    uint64_t j;
    for (j=0;j<num_matches_result;++j) {
      if (candidate_end == paired_end2) {
        paired_matches_cross_pair(
            paired_matches,search_parameters,mapper_stats,
            extended_match,matches_result[j]);
      } else {
        paired_matches_cross_pair(
            paired_matches,search_parameters,mapper_stats,
            matches_result[j],extended_match);
      }
    }
    // Check Search Accomplished (Quick abandon condition)
    if (paired_matches_classify_search_accomplished(paired_matches)) return;
  }
  // Update matches extended
  extended_matches->matches_extended = true;
  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_EXTEND_CANDIDATES,PROFILE_LEVEL);
}
/*
 * Extend matches
 */
void archive_search_pe_extend_matches(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches,
    const sequence_end_t candidate_end) {
  // Parameters
  search_parameters_t* const search_parameters = &archive_search_end1->search_parameters;
  matches_t* const matches_end1 = paired_matches->matches_end1;
  matches_t* const matches_end2 = paired_matches->matches_end2;
  mapper_stats_t* const mapper_stats = archive_search_end1->mapper_stats;
  // Configure extension
  archive_search_t* candidate_archive_search;
  matches_t* extended_matches, *candidate_matches;
  if (candidate_end==paired_end2) {
    // Extend towards end/2
    candidate_archive_search = archive_search_end2;
    extended_matches = matches_end1;
    candidate_matches = matches_end2;
  } else {
    // Extend towards end/1
    candidate_archive_search = archive_search_end1;
    extended_matches = matches_end2;
    candidate_matches = matches_end1;
  }
  // Perform extension
  archive_search_pe_extend_matches_ends(
      candidate_archive_search,candidate_end,
      search_parameters,mapper_stats,
      extended_matches,candidate_matches,
      paired_matches);
}
