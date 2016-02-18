/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates_align.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering/filtering_candidates_align_local.h"
#include "filtering/filtering_candidates_align.h"

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
void filtering_candidates_align_local_accepted_subdominant(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    const bool emulated_rc_search,
    matches_t* const matches) {
  // Parameters
  const uint64_t num_regions = vector_get_used(filtering_candidates->discarded_regions);
  filtering_region_t* regions_in = vector_get_mem(filtering_candidates->discarded_regions,filtering_region_t);
  filtering_region_t* regions_out = regions_in;
  // Traverse all regions
  uint64_t i;
  for (i=0;i<num_regions;++i,++regions_in) {
    // Skip regions
    if (regions_in->status != filtering_region_accepted_subdominant) {
      *regions_out = *regions_in;
      ++regions_out;
      continue;
    }
    // Align Region
    filtering_candidates_align_region(filtering_candidates,
        regions_in,pattern,emulated_rc_search,true,false,matches);
  }
  // Update used
  vector_update_used(filtering_candidates->discarded_regions,regions_out);
}
void filtering_candidates_align_local(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    const bool emulated_rc_search,
    matches_t* const matches) {
  PROFILE_START(GP_FC_REALIGN_LOCAL,PROFILE_LEVEL);
  // Align-Local accepted-subdominant
  filtering_candidates_align_local_accepted_subdominant(
      filtering_candidates,pattern,emulated_rc_search,matches);
  // Add pending local matches (found so far)
  locator_t* const locator = filtering_candidates->archive->locator;
  mm_stack_t* const mm_stack = filtering_candidates->mm_stack;
  matches_add_pending_local_matches(matches,locator,mm_stack);
//  // Check total alignments found
//  if (!matches_is_mapped(matches)) {
//    // Resort to align verify-discarded
//     TODO
//  }
  // DEBUG
  gem_cond_debug_block(DEBUG_FILTERING_CANDIDATES) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Filtering.Candidates (local_align)\n");
    tab_global_inc();
    filtering_candidates_print_regions(gem_log_get_stream(),filtering_candidates,false,false);
    tab_global_dec();
  }
  PROFILE_STOP(GP_FC_REALIGN_LOCAL,PROFILE_LEVEL);
}
