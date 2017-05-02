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
 *   Filtering module provides functions to produce full-alignments
 *   from the accepted filtering-regions
 */

#include "filtering/candidates/filtering_candidates_align.h"
#include "filtering/candidates/filtering_candidates_process.h"
#include "filtering/candidates/filtering_candidates_classify.h"
#include "filtering/region/filtering_region_align.h"

/*
 * Debug
 */
#define DEBUG_FILTERING_CANDIDATES  GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Filtering Candidates Cache
 */
bool filtering_candidates_align_search_cache(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const region,
    const uint64_t run_length,
    match_trace_t* const match_trace) {
  // Search the cache
  match_trace_t* const match_trace_cache =
      filtering_region_transient_cache_search(
          &filtering_candidates->filtering_region_cache,region);
  if (match_trace_cache==NULL) return false;
  // Clone the match-trace found in the cache
  filtering_region_align_clone(match_trace_cache,match_trace,region);
  return true;
}
/*
 * Filtering Candidates Region Align
 */
void filtering_candidates_align_region(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const region,
    pattern_t* const pattern,
    const bool local_alignment,
    matches_t* const matches,
    bool* const region_accepted,
    bool* const match_accepted) {
  // Parameters
  archive_t* const archive = filtering_candidates->archive;
  locator_t* const locator = archive->locator;
  archive_text_t* const archive_text = archive->text;
  // Retrieve Candidate (if needed)
  filtering_region_retrieve_text(region,
      pattern,archive_text,filtering_candidates->mm_allocator);
  // Search Cache (Before jumping into aligning the region)
  match_trace_t match_trace;
  bool match_trace_aligned =
      filtering_candidates_align_search_cache(
          filtering_candidates,region,pattern->run_length,&match_trace);
  // Align the region
  if (!match_trace_aligned) {
    match_trace_aligned = filtering_region_align(filtering_candidates,
        region,pattern,local_alignment,matches,&match_trace);
    if (!match_trace_aligned) { // Not aligned or subdominant
      *region_accepted = false;
      *match_accepted = false;
      return;
    }
  }
  // Add to matches
  if (match_trace.type == match_type_local && !local_alignment) {
    // Add Local Alignment (Pending)
    matches_add_local_match_pending(matches,&match_trace);
    *region_accepted = true;
    *match_accepted = false;
    return;
  } else {
    // Add Global Alignment
    bool match_replaced;
    match_trace_t* const match_trace_added =
        matches_add_match_trace(matches,locator,&match_trace,&match_replaced);
    if (match_trace_added==NULL) {
      *region_accepted = true;
      *match_accepted = false;
      return;
    }
    filtering_region_transient_cache_add(
        &filtering_candidates->filtering_region_cache,region,match_trace_added);
    *region_accepted = true;
    *match_accepted = !match_replaced; // Repeated?
    return;
  }
}
void filtering_candidates_align_extended_region(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const region,
    pattern_t* const pattern,
    matches_t* const matches,
    bool* const region_accepted,
    bool* const match_accepted) {
  // Parameters
  archive_t* const archive = filtering_candidates->archive;
  locator_t* const locator = archive->locator;
  archive_text_t* const archive_text = archive->text;
  // Retrieve Candidate (if needed)
  filtering_region_retrieve_text(region,
      pattern,archive_text,filtering_candidates->mm_allocator);
  // Align the region
  match_trace_t match_trace;
  bool match_trace_aligned =
      filtering_region_align(
          filtering_candidates,region,pattern,
          false,matches,&match_trace);
  if (!match_trace_aligned) { // Not aligned
    *region_accepted = false;
    *match_accepted = false;
    return;
  }
  // Add to matches (Global/Local Alignment)
  bool match_replaced;
  match_trace_t* const match_trace_added =
      matches_add_match_trace(matches,locator,&match_trace,&match_replaced);
  if (match_trace_added==NULL) {
    *region_accepted = true;
    *match_accepted = false;
    return;
  }
  // Add to extended-matches (Global/Local Alignment)
  match_trace_added->type = match_type_extended;
  if (!match_replaced) {
    vector_insert(matches->match_traces_extended,match_trace_added,match_trace_t*);
  }
  *region_accepted = true;
  *match_accepted = !match_replaced; // Repeated?
  return;
}
/*
 * Filtering Candidates (Re)Alignment
 */
uint64_t filtering_candidates_align_candidates(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    const bool local_alignment,
    matches_t* const matches) {
  // DEBUG
  gem_cond_debug_block(DEBUG_FILTERING_CANDIDATES) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Filtering.Candidates (align_acepted_regions)\n");
    tab_global_inc();
  }
  // Hint to matches
  const uint64_t num_filtering_regions = filtering_candidates_get_num_regions(filtering_candidates);
  if (num_filtering_regions==0) return 0;
  PROFILE_START(GP_FC_REALIGN_CANDIDATE_REGIONS,PROFILE_LEVEL);
  // Prepare Candidate Vectors
  filtering_candidates_sort_regions_by_align_distance(filtering_candidates); // Sort wrt align_distance
  filtering_region_t** const regions_in = filtering_candidates_get_regions(filtering_candidates);
  filtering_region_t** const regions_discarded =
      filtering_candidates_reserve_discarded_regions(filtering_candidates,num_filtering_regions);
  // Clear cache
  filtering_region_cache_clear(&filtering_candidates->filtering_region_cache);
  // Traverse all accepted candidates (text-space)
  uint64_t n, num_accepted_matches = 0;
  bool region_accepted, match_accepted;
  for (n=0;n<num_filtering_regions;++n) {
    // Fetch
    filtering_region_t* const filtering_region = regions_in[n];
    // Check if candidate is subdominant (check distance bounds)
    bool candidate_subdominant =
        filtering_candidates_classify_subdominant_match(
            filtering_candidates,filtering_region,pattern,matches);
    if (candidate_subdominant) {
      PROF_INC_COUNTER(GP_FC_SELECT_PRUNE_HIT);
      matches_metrics_set_limited_candidates(&matches->metrics,true);
      filtering_region->status = filtering_region_accepted_subdominant;
      continue;
    }
    // Align Region
    filtering_candidates_align_region(
        filtering_candidates,filtering_region,pattern,local_alignment,
        matches,&region_accepted,&match_accepted);
    if (!region_accepted) {
      filtering_region->status = filtering_region_accepted_subdominant;
    } else {
      filtering_region->status = filtering_region_accepted;
      if (match_accepted) ++num_accepted_matches;
    }
  }
  // Clean
  uint64_t num_regions_out = 0, num_regions_discarded = 0;
  for (n=0;n<num_filtering_regions;++n) {
    filtering_region_t* const filtering_region = regions_in[n];
    if (filtering_region->status==filtering_region_accepted) {
      filtering_candidates_free_region(filtering_candidates,filtering_region); // Free
    } else {
      regions_discarded[num_regions_discarded++] = filtering_region;
    }
  }
  // Update used
  matches_metrics_add_accepted_candidates(&matches->metrics,num_filtering_regions);
  filtering_candidates_set_num_regions(filtering_candidates,num_regions_out);
  filtering_candidates_add_num_discarded_regions(filtering_candidates,num_regions_discarded);
  // DEBUG
  gem_cond_debug_block(DEBUG_FILTERING_CANDIDATES) {
    tab_global_dec();
    tab_fprintf(gem_log_get_stream(),"[GEM]>Filtering.Candidates (align_acepted_regions - aftermath)\n");
    tab_global_inc();
    filtering_candidates_print_regions(gem_log_get_stream(),filtering_candidates,false);
    tab_global_dec();
  }
  // Return total accepted regions
  PROFILE_STOP(GP_FC_REALIGN_CANDIDATE_REGIONS,PROFILE_LEVEL);
  return num_accepted_matches;
}
uint64_t filtering_candidates_align_extended_candidates(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    matches_t* const matches) {
  // DEBUG
  gem_cond_debug_block(DEBUG_FILTERING_CANDIDATES) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Filtering.Candidates (align_acepted_regions)\n");
    tab_global_inc();
  }
  // Check number of filtering regions
  const uint64_t num_filtering_regions = filtering_candidates_get_num_regions(filtering_candidates);
  if (num_filtering_regions==0) return 0;
  PROFILE_START(GP_FC_REALIGN_CANDIDATE_REGIONS,PROFILE_LEVEL);
  // Sort
  filtering_candidates_sort_regions_by_align_distance(filtering_candidates); // Sort wrt align_distance
  // Traverse all accepted candidates
  filtering_region_t** const regions_in = filtering_candidates_get_regions(filtering_candidates);
  uint64_t n, num_accepted_matches = 0;
  for (n=0;n<num_filtering_regions;++n) {
    // Fetch Region
    filtering_region_t* const filtering_region = regions_in[n];
    // Align Region
    bool region_accepted, match_accepted;
    filtering_candidates_align_extended_region(
        filtering_candidates,filtering_region,pattern,
        matches,&region_accepted,&match_accepted);
    if (!region_accepted) {
      filtering_region->status = filtering_region_accepted_subdominant;
    } else {
      filtering_region->status = filtering_region_accepted;
      if (match_accepted) ++num_accepted_matches;
    }
    // Free
    filtering_candidates_free_region(filtering_candidates,filtering_region);
  }
  // Clean
  filtering_candidates_set_num_regions(filtering_candidates,0);

  matches_metrics_add_accepted_candidates(&matches->metrics,num_filtering_regions);

  // Return total accepted regions
  PROFILE_STOP(GP_FC_REALIGN_CANDIDATE_REGIONS,PROFILE_LEVEL);
  return num_accepted_matches;
}
