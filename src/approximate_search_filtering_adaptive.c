/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_filtering_stages.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "approximate_search_filtering_adaptive.h"
#include "approximate_search_filtering_base.h"
#include "approximate_search_filtering_stages.h"
#include "approximate_search_filtering_control.h"
#include "approximate_search_neighborhood.h"

/*
 * Approximate Search based on Adaptive filtering
 */
GEM_INLINE void approximate_search_filtering_adaptive(approximate_search_t* const search,matches_t* const matches) {
  // Parameters
  const as_parameters_t* const actual_parameters = search->as_parameters;
  const search_parameters_t* const search_parameters = actual_parameters->search_parameters;
  pattern_t* const pattern = &search->pattern;
  // Process proper search-stage
  while (search->search_state != asearch_end) {
    switch (search->search_state) {
      case asearch_begin: // Search Start. Check basic cases
        approximate_search_adaptive_mapping_basic_cases(search);
        break;
      case asearch_read_recovery: // Read recovery
        approximate_search_exact_filtering_adaptive_recovery(search,matches);
        break;
      case asearch_exact_filtering_adaptive: // Exact-Filtering (Adaptive)
        if (search_parameters->mapping_mode==mapping_adaptive_filtering_fast) {
          approximate_search_exact_filtering_adaptive_cutoff(search,matches);
        } else {
          approximate_search_exact_filtering_adaptive_lightweight(search,matches);
        }
        break;
      case asearch_verify_candidates: // Verify Candidates
        approximate_search_verify_candidates(search,matches);
        break;
      case asearch_candidates_verified:
        asearch_control_next_state(search,asearch_exact_filtering_adaptive,matches);
        break;
      case asearch_exact_filtering_boost:
        approximate_search_exact_filtering_boost(search,matches);
        break;
      case asearch_inexact_filtering: // Inexact-Filtering
        approximate_search_inexact_filtering(search,matches);
        break;
      case asearch_neighborhood:
        approximate_search_neighborhood_search(search,matches);
        break;
      case asearch_exact_matches: // Exact Matches
        gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_fprintf(stderr,"[GEM]>ASM::Exact-Matches\n"); }
        matches_add_interval_match(matches,search->lo_exact_matches,
            search->hi_exact_matches,pattern->key_length,0,search->emulated_rc_search); // Add interval
        approximate_search_update_mcs(search,1);
        search->search_state = asearch_end;
        break;
      case asearch_unbounded_alignment: // Unbounded alignments
        approximate_search_unbounded_align(search,matches);
        break;
      case asearch_no_regions: // No regions found
        gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_fprintf(stderr,"[GEM]>ASM::No-Regions\n"); }
        approximate_search_update_mcs(search,pattern->num_wildcards);
        search->search_state = asearch_end;
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
  // Done!!
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    if (search->search_state == asearch_end) {
      tab_fprintf(stderr,"[GEM]>ASM::END ASM +++\n");
    }
  }
  PROF_ADD_COUNTER(GP_AS_ADAPTIVE_MCS,search->max_complete_stratum);
}
GEM_INLINE void approximate_search_filtering_adaptive_generate_regions(approximate_search_t* const search) {
  approximate_search_exact_filtering_fixed(search);
}
GEM_INLINE void approximate_search_filtering_adaptive_generate_candidates(approximate_search_t* const search) {
  // Process proper search-stage
  while (search->search_state != asearch_end) {
    switch (search->search_state) {
      case asearch_begin: // Search Start. Check basic cases
        approximate_search_adaptive_mapping_basic_cases(search);
        break;
      case asearch_read_recovery: // Read recovery
        approximate_search_exact_filtering_adaptive_recovery(search,NULL);
        return; // Return (candidates generated)
      case asearch_exact_filtering_adaptive: // Exact-Filtering (Adaptive)
        approximate_search_exact_filtering_adaptive_lightweight(search,NULL);
        return; // Return (candidates generated)
      case asearch_no_regions:
      case asearch_exact_matches:
      case asearch_neighborhood:
        return; // Return (no candidates generated)
      case asearch_verify_candidates:
      case asearch_candidates_verified:
      case asearch_exact_filtering_boost:
      case asearch_inexact_filtering:
      case asearch_unbounded_alignment:
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
  // Done!!
  PROF_ADD_COUNTER(GP_AS_ADAPTIVE_MCS,search->max_complete_stratum);
}
/*
 * Test
 */
GEM_INLINE void approximate_search_test(approximate_search_t* const search,matches_t* const matches) {
  pattern_t* const pattern = &search->pattern;
  // FM-Index basic exact search
//  fm_index_bsearch(search->archive->fm_index,pattern->key,
//    pattern->key_length,&search->hi_exact_matches,&search->lo_exact_matches);
//  fm_index_reverse_bsearch(search->archive->fm_index,pattern->key,
//      pattern->key_length,&search->hi_exact_matches,&search->lo_exact_matches);
//  fm_index_bsearch_pure(search->archive->fm_index,pattern->key,
//      pattern->key_length,&search->hi_exact_matches,&search->lo_exact_matches);
//  fm_index_reverse_bsearch_pure(search->archive->fm_index,pattern->key,
//      pattern->key_length,&search->hi_exact_matches,&search->lo_exact_matches);
  // Add to matches
  matches_add_interval_match(matches,search->lo_exact_matches,
      search->hi_exact_matches,pattern->key_length,0,search->emulated_rc_search);
  // Update MCS
  approximate_search_update_mcs(search,1);
  // Update next state
  search->search_state = asearch_end;
}
