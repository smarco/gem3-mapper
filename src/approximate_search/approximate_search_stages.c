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
 *   Approximate-String-Matching (ASM) module encapsulating
 *   the basic search-stages that many ASM approaches use.
 */

#include "approximate_search/approximate_search_stages.h"
#include "approximate_search/approximate_search_region_profile.h"
#include "approximate_search/approximate_search_generate_candidates.h"
#include "approximate_search/approximate_search_verify_candidates.h"
#include "approximate_search/approximate_search_neighborhood.h"
#include "filtering/region_profile/region_profile.h"
#include "filtering/region_profile/region_profile_adaptive.h"
#include "filtering/region_profile/region_profile_schedule.h"
#include "filtering/candidates/filtering_candidates_accessors.h"
#include "filtering/candidates/filtering_candidates_process.h"
#include "filtering/candidates/filtering_candidates_verify.h"
#include "filtering/candidates/filtering_candidates_align.h"
#include "filtering/candidates/filtering_candidates_align_local.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Basic Cases Selector
 */
void approximate_search_begin(approximate_search_t* const search) {
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,"[GEM]>ASM::Basic Cases\n");
    tab_global_inc();
  }
  // Parameters
  pattern_t* const pattern = &search->pattern;
  const uint64_t key_length = pattern->key_length;
  const uint64_t num_wildcards = pattern->num_wildcards;
  // All characters are wildcards
  if (key_length==num_wildcards || key_length==0) {
    search->search_stage = asearch_stage_end;
    return;
  }
  // Otherwise, go to standard exact filtering
  search->search_stage = asearch_stage_filtering_adaptive;
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); }
  return;
}
/*
 * Exact Filtering Adaptive
 */
void approximate_search_exact_filtering_adaptive(
    approximate_search_t* const search,
    matches_t* const matches) {
  PROFILE_START(GP_AS_FILTERING_ADATIVE,PROFILE_LEVEL);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,"[GEM]>ASM::Adaptive Filtering (Exact)\n");
    tab_global_inc();
    tab_global_inc();
  }
  // Region Profile (a.k.a. seeding stage)
  approximate_search_region_profile(search);
  if (search->processing_state==asearch_processing_state_no_regions) return; // Corner case
  // Verify Candidates (if needed)
  const bool verify_candidates = (matches != NULL);
  if (verify_candidates) {
    // Generate candidates
    approximate_search_generate_candidates(search);
    // Process candidates
    filtering_candidates_process_candidates(search->filtering_candidates,&search->pattern,true);
    // Verify candidates
    approximate_search_verify_candidates(search,matches);
    search->processing_state = asearch_processing_state_candidates_verified;
  } else {
    search->processing_state = asearch_processing_state_candidates_processed;
  }
  // Update MCS
  matches_update_mcs(matches,search->region_profile.num_filtered_regions); // (+ pattern->num_wildcards)
  matches_update_limited_exact_matches(matches,search->num_limited_exact_matches);
  // DEBUG
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); tab_global_dec(); }
  PROFILE_STOP(GP_AS_FILTERING_ADATIVE,PROFILE_LEVEL);
}
void approximate_search_exact_filtering_adaptive_cutoff(
    approximate_search_t* const search,
    matches_t* const matches) {
  PROFILE_START(GP_AS_FILTERING_ADATIVE,PROFILE_LEVEL);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,"[GEM]>ASM::Adaptive Filtering (Fast Cut-off Exact)\n");
    tab_global_inc();
  }
  // Parameters
  search_parameters_t* const search_parameters = search->search_parameters;
  archive_t* const archive = search->archive;
  fm_index_t* const fm_index = archive->fm_index;
  pattern_t* const pattern = &search->pattern;
  region_profile_t* const region_profile = &search->region_profile;
  filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
  // Select Key (Regular/RL)
  uint8_t* key;
  uint64_t key_length;
  if (pattern->run_length) {
    key = pattern->rl_key;
    key_length = pattern->rl_key_length;
  } else {
    key = pattern->key;
    key_length = pattern->key_length;
  }
  // Iterate process of region-candidates-verification
  const region_profile_model_t* const profile_model = &search_parameters->region_profile_model;
  region_profile_generator_t generator;
  region_profile_generator_init(&generator,region_profile,fm_index,key,key_length,false);
  while (region_profile_generator_next_region(region_profile,&generator,profile_model)) {
    // Cut-off
    if (region_profile->num_filtering_regions >= region_profile->max_expected_regions) break;
    // Generate candidates for the last region found
    PROFILE_START(GP_AS_GENERATE_CANDIDATES,PROFILE_LEVEL);
    region_search_t* const last_region = region_profile->filtering_region + (region_profile->num_filtering_regions-1);
    filtering_candidates_add_positions_from_interval(
        filtering_candidates,search_parameters,pattern,
        last_region->lo,last_region->hi,last_region->begin,last_region->end,0);
    PROFILE_STOP(GP_AS_GENERATE_CANDIDATES,PROFILE_LEVEL);
    // Verify candidates
    PROFILE_START(GP_AS_GENERATE_CANDIDATES_DYNAMIC_FILTERING,PROFILE_LEVEL);
    filtering_candidates_process_candidates(filtering_candidates,pattern,true);
    filtering_candidates_verify_candidates(filtering_candidates,pattern,matches);
    filtering_candidates_align_candidates(filtering_candidates,pattern,matches);
    PROFILE_STOP(GP_AS_GENERATE_CANDIDATES_DYNAMIC_FILTERING,PROFILE_LEVEL);
    // Cut-off condition
    ++(region_profile->num_filtered_regions);
    // TODO if (asearch_control_fulfilled(search,matches)) break;
  }
  // Check region-profile status
  if (region_profile->num_filtering_regions==0) {
    search->processing_state = asearch_processing_state_no_regions;
  } else {
    search->processing_state = asearch_processing_state_candidates_verified;
  }
  PROFILE_STOP(GP_AS_FILTERING_ADATIVE,PROFILE_LEVEL);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); }
}
/*
 * Filtering Verification (+ realign)
 */
void approximate_search_verify(
    approximate_search_t* const search,
    matches_t* const matches) {
  // Verify
  const uint64_t num_accepted_regions =
      filtering_candidates_verify_candidates(
          search->filtering_candidates,&search->pattern,matches);
  // Realign
  if (num_accepted_regions > 0) {
    filtering_candidates_align_candidates(
        search->filtering_candidates,&search->pattern,matches);
  }
  // Update state
  search->processing_state = asearch_processing_state_candidates_verified;
}
/*
 * Neighborhood Search
 */
void approximate_search_neighborhood(
    approximate_search_t* const search,
    matches_t* const matches) {
  PROFILE_START(GP_AS_NEIGHBORHOOD_SEARCH,PROFILE_LEVEL);
  // Parameters
  pattern_t* const pattern = &search->pattern;
  region_profile_t* const region_profile = &search->region_profile;
  // Prepare region-profile (fill gaps)
  region_profile_fill_gaps(region_profile,
      pattern->key,pattern->key_length,pattern->num_wildcards);
  region_profile_merge_small_regions(region_profile,search->archive->fm_index->proper_length);
  // NS
  approximate_search_neighborhood_search_partition_preconditioned(search,matches);
  PROFILE_STOP(GP_AS_NEIGHBORHOOD_SEARCH,PROFILE_LEVEL);
}
/*
 * Unbound Filtering (+ realign)
 */
void approximate_search_align_local(
    approximate_search_t* const search,
    matches_t* const matches) {
  // DEBUG
  PROFILE_START(GP_AS_LOCAL_ALIGN,PROFILE_LEVEL);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,"[GEM]>ASM::Local-Alignment Filtering (local align)\n");
    tab_global_inc();
  }
  // Align-local discarded candidates
  filtering_candidates_align_local(search->filtering_candidates,&search->pattern,matches);
  // DEBUG
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); }
  PROFILE_STOP(GP_AS_LOCAL_ALIGN,PROFILE_LEVEL);
}
/*
 * End of the search
 */
void approximate_search_end(
    approximate_search_t* const search,
    matches_t* const matches) {
  // DEBUG
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,"[GEM]>ASM::Search END\n");
  }
  // Update metrics
  const double proper_length = fm_index_get_proper_length(search->archive->fm_index);
  const search_parameters_t* const search_parameters = search->search_parameters;
  matches_metrics_set_proper_length(&matches->metrics,proper_length);
  matches_metrics_set_read_length(&matches->metrics,search->pattern.key_length);
  matches_metrics_set_swg_match_score(&matches->metrics,search_parameters->swg_penalties.generic_match_score);
  region_profile_t* const region_profile = &search->region_profile;
  matches_metrics_set_region_profile_metrics(
      &matches->metrics,
      region_profile->avg_region_length,
      region_profile->max_region_length,
      region_profile->kmer_frequency);
  // PROF
  PROF_ADD_COUNTER(GP_AS_MCS,matches->max_complete_stratum);
}
