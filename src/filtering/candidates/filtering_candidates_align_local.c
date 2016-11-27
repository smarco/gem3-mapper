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
#include "filtering/region/filtering_region_align_local.h"
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
void filtering_candidates_align_local(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    matches_t* const matches) {
  PROFILE_START(GP_FC_REALIGN_LOCAL_CANDIDATE_REGIONS,PROFILE_LEVEL);
  // Parameters
  archive_text_t* const archive_text = filtering_candidates->archive->text;
  locator_t* const locator = filtering_candidates->archive->locator;
  text_collection_t* const text_collection = &filtering_candidates->text_collection;
  mm_stack_t* const mm_stack = filtering_candidates->mm->mm_general;
  /*
   * Add pending local matches (found so far)
   */
  matches_add_pending_local_matches(matches,locator);
  // Check total alignments found
  search_parameters_t* const search_parameters = filtering_candidates->search_parameters;
  select_parameters_t* const select_parameters = &search_parameters->select_parameters;
  const uint64_t max_searched_matches = select_parameters->max_searched_matches;
  uint64_t total_matches = matches_get_num_match_traces(matches);
  if (total_matches >= max_searched_matches) {
    PROFILE_STOP(GP_FC_REALIGN_LOCAL_CANDIDATE_REGIONS,PROFILE_LEVEL);
    return;
  }
  /*
   * Exact extend all discarded regions
   */
  const uint64_t num_regions_discarded = filtering_candidates_get_num_discarded_regions(filtering_candidates);
  PROF_ADD_COUNTER(GP_CANDIDATE_REGION_LOCAL,num_regions_discarded);
  filtering_region_t** const regions_discarded = filtering_candidates_get_discarded_regions(filtering_candidates);
  match_align_input_t align_input;
  match_align_parameters_t align_parameters = { .allowed_enc = search_parameters->allowed_enc };
  uint64_t i;
  for (i=0;i<num_regions_discarded;++i) {
    filtering_region_t* const filtering_region = regions_discarded[i];
    // Retrieve Text
    filtering_region_retrieve_text(filtering_region,pattern,archive_text,text_collection);
    text_trace_t* const text_trace = text_collection_get_trace(text_collection,filtering_region->text_trace_offset);
    // Configure
    align_input.key = pattern->key;
    align_input.key_length = pattern->key_length;
    align_input.text = text_trace->text;
    align_input.text_length = text_trace->text_length;
    // Exact extend
    match_scaffold_exact(&filtering_region->match_scaffold,
        &align_input,&align_parameters,matches,mm_stack);
  }
  /*
   * Local-align the most promising regions
   */
  // Sort by scaffold-coverage
  filtering_region_cache_clear(&filtering_candidates->filtering_region_cache); // Clear cache
  filtering_candidates_sort_discarded_by_scaffold_coverage(filtering_candidates);
  const uint64_t max_regions_considered = MIN(num_regions_discarded,search_parameters->alignment_local_max_candidates);
  for (i=0;i<max_regions_considered;++i) {
    filtering_region_t* const filtering_region = regions_discarded[i];
    // Cut-off max-reported matches
    if (total_matches >= max_searched_matches) break;
    // Align Region
    PROF_INC_COUNTER(GP_CANDIDATE_REGION_LOCAL_ALIGNED);
    const bool match_added = filtering_candidates_align_region(
        filtering_candidates,filtering_region,pattern,true,false,matches);
    if (match_added) ++total_matches;
  }
  // Clear discarded-candidates
  filtering_candidates_clear_discarded_regions(filtering_candidates);
  // DEBUG
  gem_cond_debug_block(DEBUG_FILTERING_CANDIDATES) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Filtering.Candidates (local_align)\n");
    tab_global_inc();
    filtering_candidates_print_regions(gem_log_get_stream(),filtering_candidates,false);
    tab_global_dec();
  }
  PROFILE_STOP(GP_FC_REALIGN_LOCAL_CANDIDATE_REGIONS,PROFILE_LEVEL);
}
