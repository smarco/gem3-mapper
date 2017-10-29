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
 *   Approximate-String-Matching (ASM) using neighborhood-search (NS)
 */

#include "approximate_search/approximate_search_neighborhood.h"
#include "approximate_search/approximate_search_stages.h"
#include "fm_index/fm_index_search.h"
#include "filtering/candidates/filtering_candidates.h"
#include "filtering/candidates/filtering_candidates_accessors.h"
#include "filtering/candidates/filtering_candidates_process.h"
#include "filtering/candidates/filtering_candidates_verify.h"
#include "filtering/candidates/filtering_candidates_align.h"
#include "align/pattern/pattern.h"
#include "neighborhood_search/nsearch_hamming.h"
#include "neighborhood_search/nsearch_levenshtein.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Exact Search
 */
void approximate_search_neighborhood_exact_search(
    approximate_search_t* const search,
    matches_t* const matches) {
  pattern_t* const pattern = &search->pattern;
  uint64_t hi, lo;
  // FM-Index basic exact search
  fm_index_bsearch(search->archive->fm_index,pattern->key,pattern->key_length,&hi,&lo);
  // fm_index_bsearch_pure(search->archive->fm_index,pattern->key,pattern->key_length,&hi,&lo);
  // fm_index_reverse_bsearch_pure(search->archive->fm_index,pattern->key,pattern->key_length,&hi,&lo);
  // fm_index_reverse_bsearch_fb(search->archive->fm_index,pattern->key,pattern->key_length,&hi,&lo);
  // fm_index_reverse_bsearch_bf(search->archive->fm_index,pattern->key,pattern->key_length,&hi,&lo);
  // Add interval
  filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
  filtering_candidates_add_positions_from_interval(
      filtering_candidates,search->search_parameters,
      pattern,lo,hi,0,pattern->key_length,0);
  // Process+Verify candidates
  PROF_START(GP_NS_VERIFICATION);
  filtering_candidates_process_candidates(search->filtering_candidates,&search->pattern,false);
  filtering_candidates_verify_candidates(filtering_candidates,pattern,matches);
  PROF_STOP(GP_NS_VERIFICATION);
  // Align
  PROF_START(GP_NS_ALIGN);
  filtering_candidates_align_candidates(filtering_candidates,pattern,matches);
  PROF_STOP(GP_NS_ALIGN);
  // Finish Search
  matches_update_mcs(matches,1);
  approximate_search_end(search,matches);
}
/*
 * Neighborhood Search Brute Force
 */
void approximate_search_neighborhood_search_brute_force(
    approximate_search_t* const search,
    matches_t* const matches) {
  PROFILE_START(GP_AS_NEIGHBORHOOD_SEARCH,PROFILE_LEVEL);
  // Parameters
  search_parameters_t* const search_parameters = search->search_parameters;
  filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
  pattern_t* const pattern = &search->pattern;
  // Generate Candidates (Select Alignment Model)
  if (search_parameters->match_alignment_model == match_alignment_model_hamming) {
    nsearch_hamming_brute_force(search,matches);
  } else {
    nsearch_levenshtein_brute_force(search,true,matches);
  }
  // Process+Verify candidates
  PROF_START(GP_NS_VERIFICATION);
  filtering_candidates_process_candidates(search->filtering_candidates,&search->pattern,false);
  filtering_candidates_verify_candidates(filtering_candidates,pattern,matches);
  PROF_STOP(GP_NS_VERIFICATION);
  // Align
  PROF_START(GP_NS_ALIGN);
  filtering_candidates_align_candidates(filtering_candidates,pattern,matches);
  PROF_STOP(GP_NS_ALIGN);
  // Finish Search
  matches_update_mcs(matches,search->max_search_error+1);
  approximate_search_end(search,matches);
  PROFILE_STOP(GP_AS_NEIGHBORHOOD_SEARCH,PROFILE_LEVEL);
}
/*
 * Neighborhood Search (Using partitions)
 */
void approximate_search_neighborhood_search_partition(
    approximate_search_t* const search,
    matches_t* const matches) {
  PROFILE_START(GP_AS_NEIGHBORHOOD_SEARCH,PROFILE_LEVEL);
  // Parameters
  search_parameters_t* const search_parameters = search->search_parameters;
  region_profile_clear(&search->region_profile); // Clear for scoring purposes
  // Generate Candidates (Select Alignment Model)
  if (search_parameters->match_alignment_model == match_alignment_model_hamming) {
    nsearch_hamming(search,matches);
  } else {
    nsearch_levenshtein(search,matches);
  }
  PROFILE_STOP(GP_AS_NEIGHBORHOOD_SEARCH,PROFILE_LEVEL);
}
/*
 * Neighborhood Search (Using partitions + region-profile preconditioned)
 */
void approximate_search_neighborhood_search_partition_preconditioned(
    approximate_search_t* const search,
    matches_t* const matches) {
  PROFILE_START(GP_AS_NEIGHBORHOOD_SEARCH,PROFILE_LEVEL);
  // Parameters
  search_parameters_t* const search_parameters = search->search_parameters;
  region_profile_t* const region_profile = &search->region_profile;
  const uint64_t max_search_error = search->max_search_error;
  const uint64_t mcs = matches->max_complete_stratum;
  // Compute error limits
  region_profile_compute_error_limits(region_profile,mcs,max_search_error);
  // Generate Candidates (Select Alignment Model)
  if (search_parameters->match_alignment_model == match_alignment_model_hamming) {
    nsearch_hamming_preconditioned(search,matches);
  } else {
    nsearch_levenshtein_preconditioned(search,matches);
  }
  PROFILE_STOP(GP_AS_NEIGHBORHOOD_SEARCH,PROFILE_LEVEL);
}
