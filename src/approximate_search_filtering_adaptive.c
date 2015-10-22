/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_filtering_adaptive.h
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
 * Control
 */
GEM_INLINE void asearch_control_next_state_read_recovery(approximate_search_t* const search,matches_t* const matches) {
  switch (search->search_state) {
    case asearch_no_regions:
      search->search_state = asearch_end;
      break;
    case asearch_exact_matches:
      PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MAPPED,1);
      PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MCS,1);
      break;
    case asearch_candidates_verified:
      PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MAPPED,matches_is_mapped(matches)?1:0);
      PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MCS,search->max_complete_stratum);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
GEM_INLINE void asearch_control_next_state_exact_filtering_adaptive(approximate_search_t* const search,matches_t* const matches) {
  // Stats
  // Select state
  switch (search->search_state) {
    case asearch_no_regions:
      search->search_state = asearch_exact_filtering_boost;
      break;
    case asearch_exact_matches:
      PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MAPPED,1);
      PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MCS,1);
      break;
    case asearch_candidates_verified:
      PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MAPPED,matches_is_mapped(matches)?1:0);
      PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MCS,search->max_complete_stratum);
      search->search_state = asearch_control_trigger_boost(search,matches) ?
          asearch_exact_filtering_boost : asearch_end;
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
GEM_INLINE void asearch_control_next_state_exact_filtering_boost(approximate_search_t* const search,matches_t* const matches) {
  search_parameters_t* const search_parameters = search->as_parameters->search_parameters;
  if (matches_is_mapped(matches) || search_parameters->unbounded_alignment==unbounded_alignment_never) {
    search->search_state = asearch_end;
  } else {
    search->search_state = asearch_unbounded_alignment;
  }
}
/*
 * Adaptive mapping Initial Basic Cases Selector
 */
GEM_INLINE void approximate_search_filtering_adaptive_basic_cases(approximate_search_t* const search) {
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,"[GEM]>ASM::Basic Cases\n");
    tab_global_inc();
  }
  // Parameters
  search_parameters_t* const search_parameters = search->as_parameters->search_parameters;
  pattern_t* const pattern = &search->pattern;
  const uint64_t key_length = pattern->key_length;
  const uint64_t num_wildcards = pattern->num_wildcards;
  // All characters are wildcards
  if (key_length==num_wildcards) {
    search->hi_exact_matches = 0;
    search->lo_exact_matches = 0;
    approximate_search_update_mcs(search,key_length);
    search->search_state = asearch_end;
    return;
  }
  /*
   * Recovery
   *   If the number of wildcards (or errors required) is too high we try to recover as
   *   many matches as possible. We extract feasible regions from the read and filter
   *   them trying to recover anything out of bad quality reads.
   */
  if (num_wildcards > 0) {
    switch (search_parameters->mapping_mode) {
      case mapping_adaptive_filtering_fast:
      case mapping_adaptive_filtering_thorough:
        if (num_wildcards >= search->as_parameters->alignment_max_error_nominal) {
          search->search_state = asearch_read_recovery;
          return;
        }
        break;
      default:
        if (num_wildcards >= search->max_complete_error) {
          search->search_state = asearch_read_recovery;
          return;
        }
        break;
    }
  }
  // Exact search
  if (search->max_complete_error==0) {
    search->search_state = asearch_neighborhood;
    return;
  }
  // Very short reads (Neighborhood search)
  if (key_length <= RANK_MTABLE_SEARCH_DEPTH || key_length < search->archive->fm_index->proper_length) {
    search->search_state = asearch_neighborhood;
    return;
  }
  // Otherwise, go to standard exact filtering
  search->search_state = asearch_exact_filtering_adaptive;
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); }
  return;
}
/*
 * Adaptive mapping [GEM-workflow 4.0]
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
        approximate_search_filtering_adaptive_basic_cases(search);
        break;
      case asearch_read_recovery: // Read recovery
        approximate_search_exact_filtering_adaptive_recovery(search,matches);
        asearch_control_next_state_read_recovery(search,matches); // Next State
        break;
      case asearch_exact_filtering_adaptive: // Exact-Filtering (Adaptive)
        if (search_parameters->mapping_mode==mapping_adaptive_filtering_fast) {
          approximate_search_exact_filtering_adaptive_cutoff(search,matches);
        } else {
          approximate_search_exact_filtering_adaptive_lightweight(search,matches);
        }
        asearch_control_next_state_exact_filtering_adaptive(search,matches); // Next State
        break;
      case asearch_exact_filtering_boost:
        approximate_search_exact_filtering_boost(search,matches);
        asearch_control_next_state_exact_filtering_boost(search,matches);
        break;
//      case asearch_inexact_filtering: // Inexact-Filtering
//        approximate_search_inexact_filtering(search,matches);
//        break;
      case asearch_neighborhood:
        approximate_search_neighborhood_search(search,matches);
        search->search_state = asearch_end; // Next State
        break;
      case asearch_exact_matches: // Exact Matches
        gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_fprintf(stderr,"[GEM]>ASM::Exact-Matches\n"); }
        matches_add_interval_match(matches,search->lo_exact_matches,
            search->hi_exact_matches,pattern->key_length,0,search->emulated_rc_search); // Add interval
        approximate_search_update_mcs(search,1);
        search->search_state = asearch_end; // Next State
        break;
      case asearch_unbounded_alignment: // Unbounded alignments
        approximate_search_unbounded_align(search,matches);
        search->search_state = asearch_end; // Next State
        break;
      case asearch_no_regions: // No regions found
        gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_fprintf(stderr,"[GEM]>ASM::No-Regions\n"); }
        approximate_search_update_mcs(search,pattern->num_wildcards);
        search->search_state = asearch_end; // Next State
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
  // Set MCS
  matches->max_complete_stratum = MIN(matches->max_complete_stratum,search->max_complete_stratum);
  // Done!!
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    if (search->search_state == asearch_end) {
      tab_fprintf(stderr,"[GEM]>ASM::END ASM +++\n");
    }
  }
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
