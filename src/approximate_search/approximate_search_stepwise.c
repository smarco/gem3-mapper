/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_filtering_stages_stepwise.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "approximate_search/approximate_search_stepwise.h"
#include "approximate_search/approximate_search_region_profile.h"
#include "approximate_search/approximate_search_verify_candidates.h"
#include "approximate_search/approximate_search_generate_candidates.h"
#include "approximate_search/approximate_search_filtering_adaptive.h"
#include "approximate_search/approximate_search_filtering_stages.h"
#include "approximate_search/approximate_search_control.h"
#include "approximate_search/approximate_search_neighborhood.h"
#include "filtering/filtering_candidates_align.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Control
 */
void asearch_control_next_state_filtering_adaptive(
    approximate_search_t* const search,
    matches_t* const matches) {
  // Select state
  switch (search->processing_state) {
    case asearch_processing_state_no_regions:
      search->search_stage = asearch_stage_end;
      break;
    case asearch_processing_state_candidates_verified:
      PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MAPPED,matches_is_mapped(matches)?1:0);
      PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MCS,search->max_complete_stratum);
      // Local alignment
      if (search->search_parameters->local_alignment==local_alignment_never) {
        search->search_stage = asearch_stage_end;
      } else { // local_alignment_if_unmapped
        if (matches_is_mapped(matches)) {
          search->search_stage = asearch_stage_end;
        } else {
          search->search_stage = asearch_stage_local_alignment;
        }
      }
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
/*
 * AM Stepwise :: Region Profile
 */
void approximate_search_stepwise_region_profile_generate_static(approximate_search_t* const search) {
  while (true) {
    switch (search->search_stage) {
      case asearch_stage_begin: // Search Start. Check basic cases
        approximate_search_filtering_adaptive_basic_cases(search);
        break;
      case asearch_stage_read_recovery:
      case asearch_stage_filtering_adaptive:
        approximate_search_region_partition_fixed(search);
        return;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
}
void approximate_search_stepwise_region_profile_generate_adaptive(approximate_search_t* const search) {
  while (true) {
    switch (search->search_stage) {
      case asearch_stage_begin: // Search Start. Check basic cases
        approximate_search_filtering_adaptive_basic_cases(search);
        break;
      case asearch_stage_read_recovery:
      case asearch_stage_filtering_adaptive:
        // Adaptive region-profile
        approximate_search_exact_filtering_adaptive_lightweight(search,NULL);
        // Check no-regions
        if (search->processing_state == asearch_processing_state_no_regions) return;
        // Set as region-profiled
        search->processing_state = asearch_processing_state_region_profiled;
        // Check exact matches (limit the number of matches)
        region_profile_t* const region_profile = &search->region_profile;
        if (region_profile_has_exact_matches(region_profile)) {
          search_parameters_t* const search_parameters = search->search_parameters;
          select_parameters_t* const select_parameters = &search_parameters->select_parameters_align;
          region_search_t* const filtering_region = region_profile->filtering_region;
          const uint64_t total_candidates = filtering_region->hi - filtering_region->lo;
          if (select_parameters->min_reported_strata_nominal==0 &&
              total_candidates > select_parameters->max_reported_matches) {
            filtering_region->hi = filtering_region->lo + select_parameters->max_reported_matches;
            region_profile->total_candidates = select_parameters->max_reported_matches;
          }
        }
        return;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
}
void approximate_search_stepwise_region_profile_copy(
    approximate_search_t* const search,
    gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search) {
  if (search->processing_state == asearch_processing_state_region_partitioned) {
    approximate_search_region_profile_buffered_copy(search,gpu_buffer_fmi_search);
  }
}
void approximate_search_stepwise_region_profile_retrieve(
    approximate_search_t* const search,
    gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search) {
  if (search->processing_state == asearch_processing_state_region_partitioned) {
    approximate_search_region_profile_buffered_retrieve(search,gpu_buffer_fmi_search);
  }
  // Check unsuccessful cases
  if (search->processing_state == asearch_processing_state_no_regions) {
    PROF_START(GP_ASSW_REGION_PROFILE_UNSUCCESSFUL);
    // Re-Compute Region Profile Adaptively
    approximate_search_region_profile_buffered_recompute(search);
    PROF_STOP(GP_ASSW_REGION_PROFILE_UNSUCCESSFUL);
    if (search->processing_state==asearch_processing_state_no_regions) return; // Corner cases
    search->processing_state = asearch_processing_state_region_profiled;
  }
  // Check exact matches (limit the number of matches)
  region_profile_t* const region_profile = &search->region_profile;
  if (region_profile_has_exact_matches(region_profile)) {
    search_parameters_t* const search_parameters = search->search_parameters;
    select_parameters_t* const select_parameters = &search_parameters->select_parameters_align;
    region_search_t* const filtering_region = region_profile->filtering_region;
    const uint64_t total_candidates = filtering_region->hi - filtering_region->lo;
    if (select_parameters->min_reported_strata_nominal==0 &&
        total_candidates > select_parameters->max_reported_matches) {
      filtering_region->hi = filtering_region->lo + select_parameters->max_reported_matches;
      region_profile->total_candidates = select_parameters->max_reported_matches;
    }
  }
}
/*
 * AM Stepwise :: Decode Candidates
 */
void approximate_search_stepwise_decode_candidates_copy(
    approximate_search_t* const search,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  if (search->processing_state == asearch_processing_state_region_profiled) {
    approximate_search_generate_candidates_buffered_copy(search,gpu_buffer_fmi_decode);
  }
}
void approximate_search_stepwise_decode_candidates_retrieve(
    approximate_search_t* const search,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  if (search->processing_state == asearch_processing_state_region_profiled) {
    approximate_search_generate_candidates_buffered_retrieve(search,gpu_buffer_fmi_decode);
  }
}
/*
 * AM Stepwise :: Verify Candidates
 */
void approximate_search_stepwise_verify_candidates_copy(
    approximate_search_t* const search,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  if (search->processing_state == asearch_processing_state_candidates_processed) {
    approximate_search_verify_candidates_buffered_copy(search,gpu_buffer_align_bpm);
  }
}
void approximate_search_stepwise_verify_candidates_retrieve(
    approximate_search_t* const search,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    matches_t* const matches) {
  if (search->processing_state == asearch_processing_state_candidates_processed) {
    approximate_search_verify_candidates_buffered_retrieve(search,gpu_buffer_align_bpm,matches);
  }
}
/*
 * AM Stepwise :: Finish Search
 */
void approximate_search_stepwise_finish(
    approximate_search_t* const search,
    matches_t* const matches) {
  if (search->search_stage == asearch_stage_read_recovery) {
    search->search_stage = asearch_stage_end;
  } else  if (search->search_stage == asearch_stage_filtering_adaptive) {
    asearch_control_next_state_filtering_adaptive(search,matches); // Next State
  }
  // Finish search using regular workflow
  approximate_search_filtering_adaptive(search,matches);
}
