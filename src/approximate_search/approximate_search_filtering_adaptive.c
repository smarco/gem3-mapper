/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_filtering_adaptive.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include <approximate_search/approximate_search_stages.h>
#include "approximate_search/approximate_search_filtering_adaptive.h"
#include "approximate_search/approximate_search_control.h"
#include "approximate_search/approximate_search_neighborhood.h"

/*
 * Control
 */
void as_filtering_control_begin(approximate_search_t* const search) {
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
  PROF_INC_COUNTER(GP_AS_FILTERING_ADATIVE_CALL);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); }
  return;
}
void as_filtering_control_filtering_adaptive_next_state(
    approximate_search_t* const search,
    matches_t* const matches) {
  PROF_ADD_COUNTER(GP_AS_FILTERING_ADATIVE_MCS,search->region_profile.num_filtered_regions);
  // Select state
  switch (search->processing_state) {
    case asearch_processing_state_no_regions:
      search->search_stage = asearch_stage_end;
      break;
    case asearch_processing_state_candidates_verified: {
      // Local alignment
      search_parameters_t* const search_parameters = search->search_parameters;
      if (search_parameters->local_alignment==local_alignment_never || matches_is_mapped(matches)) {
        search->search_stage = asearch_stage_end;
      } else {  // local_alignment_if_unmapped
        PROF_INC_COUNTER(GP_AS_LOCAL_ALIGN_CALL);
        search->search_stage = asearch_stage_local_alignment;
      }
      break;
    }
    default:
      GEM_INVALID_CASE();
      break;
  }
}
/*
 * Adaptive mapping [GEM-workflow 4.0]
 */
void approximate_search_filtering_adaptive(
    approximate_search_t* const search,
    matches_t* const matches) {
  // Process proper search-stage
  while (true) {
    switch (search->search_stage) {
      case asearch_stage_begin: // Search begin
        as_filtering_control_begin(search);
        break;
      case asearch_stage_filtering_adaptive: // Exact-Filtering (Adaptive)
        // approximate_search_exact_filtering_adaptive_cutoff(search,matches);
        approximate_search_exact_filtering_adaptive(search,matches);
        as_filtering_control_filtering_adaptive_next_state(search,matches); // Next State
        break;
      case asearch_stage_filtering_adaptive_finished:
        as_filtering_control_filtering_adaptive_next_state(search,matches); // Next State
        break;
      case asearch_stage_local_alignment: // Local alignments
        approximate_search_align_local(search,matches);
        search->search_stage = asearch_stage_end; // Next State
        break;
      case asearch_stage_end:
        approximate_search_end(search,matches);
        return;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
}
