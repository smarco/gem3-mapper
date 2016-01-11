/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates_verify.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering_candidates_verify.h"
#include "filtering_candidates_process.h"
#include "filtering_region_verify.h"

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
    filtering_candidates_t* const filtering_candidates,text_collection_t* const text_collection,
    const pattern_t* const pattern,const as_parameters_t* const as_parameters,
    matches_t* const matches) {
  PROFILE_START(GP_FC_VERIFY_CANDIDATE_REGIONS,PROFILE_LEVEL);
  // Parameters
  search_parameters_t* const search_parameters = as_parameters->search_parameters;
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
      if (filtering_region_verify(regions_in,text_collection,search_parameters,pattern)) {
        *regions_out = *regions_in;
        ++regions_out;
        ++num_regions_accepted;
      } else {
        *regions_discarded = *regions_in;
        ++regions_discarded;
      }
      // Add to verify regions
      regions_verified->begin_position = regions_in->begin_position;
      regions_verified->end_position = regions_in->end_position;
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
    filtering_candidates_print_regions(gem_log_get_stream(),filtering_candidates,text_collection,false,false);
    tab_global_dec();
  }
  PROFILE_STOP(GP_FC_VERIFY_CANDIDATE_REGIONS,PROFILE_LEVEL);
  return num_regions_accepted;
}
uint64_t filtering_candidates_verify_filtering_regions_multiple_hits(
    filtering_candidates_t* const filtering_candidates,text_collection_t* const text_collection,
    const pattern_t* const pattern,const as_parameters_t* const as_parameters,
    matches_t* const matches) {
  // Parameters
  search_parameters_t* const search_parameters = as_parameters->search_parameters;
  // Traverse all regions (text-space)
  const uint64_t num_filtering_regions = vector_get_used(filtering_candidates->filtering_regions);
  vector_reserve_additional(filtering_candidates->verified_regions,num_filtering_regions);
  vector_reserve_additional(filtering_candidates->discarded_regions,num_filtering_regions);
  verified_region_t* regions_verified = vector_get_free_elm(filtering_candidates->verified_regions,verified_region_t);
  filtering_region_t* regions_discarded = vector_get_free_elm(filtering_candidates->discarded_regions,filtering_region_t);
  filtering_region_t* regions_in = vector_get_mem(filtering_candidates->filtering_regions,filtering_region_t);
  uint64_t n, total_regions_accepted = 0;
  for (n=0;n<num_filtering_regions;++n,++regions_in) {
    // Check region status (Skip other than unverified)
    if (regions_in->status!=filtering_region_unverified) {
      filtering_region_t* regions_out;
      vector_alloc_new(filtering_candidates->filtering_regions,filtering_region_t,regions_out);
      *regions_out = *regions_in;
    } else {
      // Verify region
      uint64_t num_regions_accepted = filtering_region_verify_multiple_hits(
          filtering_candidates->filtering_regions,regions_in,text_collection,search_parameters,pattern);
      if (num_regions_accepted > 0) {
        total_regions_accepted += num_regions_accepted;
      } else {
        *regions_discarded = *regions_in;
        ++regions_discarded;
      }
      // Add to verify regions
      regions_verified->begin_position = regions_in->begin_position;
      regions_verified->end_position = regions_in->end_position;
      ++regions_verified;
    }
  }
  // Update used
  vector_update_used(filtering_candidates->discarded_regions,regions_discarded);
  vector_update_used(filtering_candidates->verified_regions,regions_verified);
  // DEBUG
  gem_cond_debug_block(DEBUG_FILTERING_CANDIDATES) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Filtering.Candidates (verify_regions_multiple_hits)\n");
    tab_global_inc();
    filtering_candidates_print_regions(gem_log_get_stream(),filtering_candidates,text_collection,false,false);
    tab_global_dec();
  }
  return total_regions_accepted;
}
uint64_t filtering_candidates_verify_candidates(
    filtering_candidates_t* const filtering_candidates,archive_t* const archive,
    text_collection_t* const text_collection,const pattern_t* const pattern,
    const as_parameters_t* const as_parameters,
    matches_t* const matches,mm_stack_t* const mm_stack) {
  // Check number of filtering regions
  uint64_t pending_candidates = vector_get_used(filtering_candidates->filtering_regions);
  if (pending_candidates==0) return 0;
  PROFILE_START(GP_FC_VERIFICATION,PROFILE_LEVEL);
  // Retrieve text-candidates
  filtering_candidates_retrieve_filtering_regions(filtering_candidates,archive->text,text_collection,mm_stack);
  // Verify candidates
  pending_candidates = filtering_candidates_verify_filtering_regions(
      filtering_candidates,text_collection,pattern,as_parameters,matches);
  PROFILE_STOP(GP_FC_VERIFICATION,PROFILE_LEVEL);
  return pending_candidates;
}
