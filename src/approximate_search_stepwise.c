/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_filtering_stages_stepwise.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "approximate_search_stepwise.h"
#include "approximate_search_region_profile.h"
#include "approximate_search_verify_candidates.h"
#include "approximate_search_generate_candidates.h"
#include "approximate_search_filtering_adaptive.h"
#include "approximate_search_filtering_stages.h"
#include "approximate_search_control.h"
#include "approximate_search_neighborhood.h"
#include "filtering_candidates_align.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Control
 */
//GEM_INLINE void asearch_control_next_state_generate_candidates_recovery(approximate_search_t* const search) {
//  // Stats
//  PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MCS,search->max_complete_stratum);
//  // Select state
//  switch (search->search_state) {
//    case asearch_no_regions:
//      search->search_state = asearch_end;
//      break;
//    case asearch_exact_matches:
//    case asearch_verify_candidates:
//      break;
//    default:
//      GEM_INVALID_CASE();
//      break;
//  }
//}
//GEM_INLINE void asearch_control_next_state_generate_candidates(approximate_search_t* const search) {
//  // Stats
//  PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MCS,search->max_complete_stratum);
//  // Select state
//  switch (search->search_state) {
//    case asearch_no_regions:
//      search->search_state = asearch_exact_filtering_boost;
//      break;
//    case asearch_exact_matches:
//    case asearch_verify_candidates:
//      break;
//    default:
//      GEM_INVALID_CASE();
//      break;
//  }
//}
//GEM_INLINE void asearch_control_next_state_candidates_verified(approximate_search_t* const search,matches_t* const matches) {
//  // Select state
//  switch (search->search_state) {
//    case asearch_exact_matches:
//      PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MAPPED,1);
//      PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MCS,1);
//      return;
//    case asearch_candidates_verified:
//      PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MAPPED,matches_is_mapped(matches)?1:0);
//      PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MCS,search->max_complete_stratum);
//      search->search_state = asearch_control_trigger_boost(search,matches) ?
//          asearch_exact_filtering_boost : asearch_end;
//      break;
//    case asearch_exact_filtering_boost:
//    case asearch_end:
//      break;
//    default:
//      GEM_INVALID_CASE();
//      break;
//  }
//}
/*
 * AM Stepwise :: Region Profile
 */
GEM_INLINE void approximate_search_stepwise_region_profile_generate(approximate_search_t* const search) {
  while (true) {
    switch (search->search_stage) {
      case asearch_stage_begin: // Search Start. Check basic cases
        approximate_search_filtering_adaptive_basic_cases(search);
        break;
      case asearch_stage_filtering_adaptive:
        approximate_search_region_partition_fixed(search);
        return;
      case asearch_stage_read_recovery: return;
      case asearch_stage_neighborhood: return;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
}
GEM_INLINE void approximate_search_stepwise_region_profile_copy(
    approximate_search_t* const search,gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search) {
  if (search->processing_state == asearch_processing_state_region_partitioned) {
    approximate_search_region_profile_buffered_copy(search,gpu_buffer_fmi_search);
  }
}
GEM_INLINE void approximate_search_stepwise_region_profile_retrieve(
    approximate_search_t* const search,gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search) {
  if (search->processing_state == asearch_processing_state_region_partitioned) {
    approximate_search_region_profile_buffered_retrieve(search,gpu_buffer_fmi_search);
    // Check results
    if (search->processing_state == asearch_processing_state_no_regions) {
      approximate_search_exact_filtering_adaptive_lightweight(search,NULL);
      GEM_NOT_IMPLEMENTED(); // TODO Split until processing
    }
  }
}
/*
 * AM Stepwise :: Decode Candidates
 */
GEM_INLINE void approximate_search_stepwise_decode_candidates_generate(approximate_search_t* const search) {
  // ???
}
GEM_INLINE void approximate_search_stepwise_decode_candidates_copy(
    approximate_search_t* const search,gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  if (search->processing_state == asearch_processing_state_region_profiled) {
    approximate_search_generate_exact_candidates_buffered_copy(search,gpu_buffer_fmi_decode);
  }
}
GEM_INLINE void approximate_search_stepwise_decode_candidates_retrieve(
    approximate_search_t* const search,gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  if (search->processing_state == asearch_processing_state_region_profiled) {
    approximate_search_generate_exact_candidates_buffered_retrieve(search,gpu_buffer_fmi_decode);
  }
}
/*
 * AM Stepwise :: Verify Candidates
 */
GEM_INLINE void approximate_search_stepwise_verify_candidates_generate(approximate_search_t* const search) {
  // ???
}
GEM_INLINE void approximate_search_stepwise_verify_candidates_copy(
    approximate_search_t* const search,gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  if (search->processing_state == asearch_processing_state_candidates_processed) {
    approximate_search_verify_candidates_buffered_copy(search,gpu_buffer_align_bpm);
  }
}
GEM_INLINE void approximate_search_stepwise_verify_candidates_retrieve(
    approximate_search_t* const search,gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    matches_t* const matches) {
  if (search->processing_state == asearch_processing_state_candidates_processed) {
    approximate_search_verify_candidates_buffered_retrieve(search,gpu_buffer_align_bpm,matches);
  }
}
/*
 * AM Stepwise :: Finish Search
 */
GEM_INLINE void approximate_search_stepwise_finish(approximate_search_t* const search,matches_t* const matches) {
  // Next State
  // FIXME asearch_control_next_state_candidates_verified(search,matches);
  // Finish search using regular workflow
  approximate_search_filtering_adaptive(search,matches);
}
