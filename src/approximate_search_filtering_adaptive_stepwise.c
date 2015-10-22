/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_filtering_stages_stepwise.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "approximate_search_filtering_adaptive_stepwise.h"
#include "approximate_search_filtering_adaptive.h"
#include "approximate_search_filtering_base.h"
#include "approximate_search_filtering_stages.h"
#include "approximate_search_filtering_control.h"
#include "approximate_search_neighborhood.h"

/*
 * Control
 */
GEM_INLINE void asearch_control_next_state_generate_candidates_recovery(approximate_search_t* const search) {
  // Stats
  PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MCS,search->max_complete_stratum);
  // Select state
  switch (search->search_state) {
    case asearch_no_regions:
      search->search_state = asearch_end;
      break;
    case asearch_exact_matches:
    case asearch_verify_candidates:
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
GEM_INLINE void asearch_control_next_state_generate_candidates(approximate_search_t* const search) {
  // Stats
  PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MCS,search->max_complete_stratum);
  // Select state
  switch (search->search_state) {
    case asearch_no_regions:
      search->search_state = asearch_exact_filtering_boost;
      break;
    case asearch_exact_matches:
    case asearch_verify_candidates:
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
GEM_INLINE void asearch_control_next_state_candidates_verified(approximate_search_t* const search,matches_t* const matches) {
  // Select state
  switch (search->search_state) {
    case asearch_exact_matches:
      PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MAPPED,1);
      PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MCS,1);
      return;
    case asearch_candidates_verified:
      PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MAPPED,matches_is_mapped(matches)?1:0);
      PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MCS,search->max_complete_stratum);
      search->search_state = asearch_control_trigger_boost(search,matches) ?
          asearch_exact_filtering_boost : asearch_end;
      break;
    case asearch_exact_filtering_boost:
    case asearch_end:
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
/*
 * Approximate Search based on Adaptive filtering (Stepwise)
 */
GEM_INLINE void approximate_search_filtering_adaptive_stepwise_generate_regions(approximate_search_t* const search) {
  approximate_search_exact_filtering_fixed(search);
}
GEM_INLINE void approximate_search_filtering_adaptive_stepwise_generate_candidates(approximate_search_t* const search) {
  // Process proper search-stage
  while (search->search_state != asearch_end) {
    switch (search->search_state) {
      case asearch_begin: // Search Start. Check basic cases
        approximate_search_filtering_adaptive_basic_cases(search);
        break;
      case asearch_read_recovery: // Read recovery
        approximate_search_exact_filtering_adaptive_recovery(search,NULL);
        asearch_control_next_state_generate_candidates_recovery(search); // Next State
        return; // Return (candidates generated)
      case asearch_exact_filtering_adaptive: // Exact-Filtering (Adaptive)
        approximate_search_exact_filtering_adaptive_lightweight(search,NULL);
        asearch_control_next_state_generate_candidates(search); // Next State
        return; // Return (candidates generated)
      case asearch_neighborhood:
        return; // Return (no candidates generated)
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
}
GEM_INLINE void approximate_search_filtering_adaptive_stepwise_finish(approximate_search_t* const search,matches_t* const matches) {
  // Next State
  asearch_control_next_state_candidates_verified(search,matches);
  // Finish search using regular workflow
  approximate_search_filtering_adaptive(search,matches);
}
