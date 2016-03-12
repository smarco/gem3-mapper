/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates_verify.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering/filtering_candidates_verify.h"
#include "filtering/filtering_candidates_process.h"
#include "filtering/filtering_region_verify.h"

/*
 * Debug
 */
#define DEBUG_FILTERING_CANDIDATES  GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Candidate Verification
 */
uint64_t filtering_candidates_verify_filtering_regions(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern) {
  PROFILE_START(GP_FC_VERIFY_CANDIDATES,PROFILE_LEVEL);
  // Traverse all regions (text-space)
  const uint64_t num_filtering_regions = vector_get_used(filtering_candidates->filtering_regions);
  vector_reserve_additional(filtering_candidates->verified_regions,num_filtering_regions);
  vector_reserve_additional(filtering_candidates->discarded_regions,num_filtering_regions);
  verified_region_t* regions_verified = vector_get_free_elm(filtering_candidates->verified_regions,verified_region_t);
  filtering_region_t* regions_discarded = vector_get_free_elm(filtering_candidates->discarded_regions,filtering_region_t);
  filtering_region_t* regions_in = vector_get_mem(filtering_candidates->filtering_regions,filtering_region_t);
  filtering_region_t* regions_out = regions_in;
  uint64_t n, num_regions_accepted = 0;
  for (n=0;n<num_filtering_regions;++n,++regions_in) {
    // Check region status (Skip other than unverified)
    if (regions_in->status!=filtering_region_unverified) {
      *regions_out = *regions_in;
      ++regions_out;
    } else {
      // Verify region
      if (filtering_region_verify(filtering_candidates,regions_in,pattern,true)) {
        *regions_out = *regions_in;
        ++regions_out;
        ++num_regions_accepted;
      } else {
        *regions_discarded = *regions_in;
        ++regions_discarded;
      }
      // Add to verify regions
      regions_verified->begin_position = regions_in->text_begin_position;
      regions_verified->end_position = regions_in->text_end_position;
      ++regions_verified;
    }
  }
  // Update Used
  vector_update_used(filtering_candidates->filtering_regions,regions_out);
  vector_update_used(filtering_candidates->discarded_regions,regions_discarded);
  vector_update_used(filtering_candidates->verified_regions,regions_verified);
  // DEBUG
  gem_cond_debug_block(DEBUG_FILTERING_CANDIDATES) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Filtering.Candidates (verify_regions)\n");
    tab_global_inc();
    filtering_candidates_print_regions(gem_log_get_stream(),filtering_candidates,false);
    tab_global_dec();
  }
  PROFILE_STOP(GP_FC_VERIFY_CANDIDATES,PROFILE_LEVEL);
  return num_regions_accepted;
}
uint64_t filtering_candidates_verify_candidates(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern) {
  // Check number of filtering regions
  uint64_t pending_candidates = vector_get_used(filtering_candidates->filtering_regions);
  if (pending_candidates==0) return 0;
  PROFILE_START(GP_FC_VERIFICATION,PROFILE_LEVEL);
  // Verify candidates
  pending_candidates = filtering_candidates_verify_filtering_regions(filtering_candidates,pattern);
  PROFILE_STOP(GP_FC_VERIFICATION,PROFILE_LEVEL);
  return pending_candidates;
}
