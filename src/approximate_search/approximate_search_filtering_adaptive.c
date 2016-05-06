/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_filtering_adaptive.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "approximate_search/approximate_search_filtering_adaptive.h"
#include "approximate_search/approximate_search_filtering_stages.h"
#include "approximate_search/approximate_search_control.h"
#include "approximate_search/approximate_search_neighborhood.h"

/*
 * Control
 */
void asearch_control_next_state_read_recovery(
    approximate_search_t* const search,
    matches_t* const matches) {
  switch (search->processing_state) {
    case asearch_processing_state_no_regions:
      break;
    case asearch_processing_state_candidates_verified:
      PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MAPPED,matches_is_mapped(matches)?1:0);
      PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MCS,search->max_complete_stratum);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  // Finish
  search->search_stage = asearch_stage_end;
}
void asearch_control_next_state_exact_filtering_adaptive(
    approximate_search_t* const search,
    matches_t* const matches) {
  // Select state
  search_parameters_t* const search_parameters = search->search_parameters;
  switch (search->processing_state) {
    case asearch_processing_state_no_regions:
      search->search_stage = asearch_stage_end;
      break;
    case asearch_processing_state_candidates_verified:
      PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MAPPED,matches_is_mapped(matches)?1:0);
      PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MCS,search->max_complete_stratum);
      // Local alignment
      if (search_parameters->local_alignment==local_alignment_never) {
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
 * Adaptive mapping Initial Basic Cases Selector
 */
void approximate_search_filtering_adaptive_basic_cases(approximate_search_t* const search) {
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,"[GEM]>ASM::Basic Cases\n");
    tab_global_inc();
  }
  // Parameters
  pattern_t* const pattern = &search->pattern;
  const uint64_t key_length = pattern->key_length;
  const uint64_t num_wildcards = pattern->num_wildcards;
  // All characters are wildcards
  if (key_length==num_wildcards) {
    approximate_search_update_mcs(search,key_length);
    search->search_stage = asearch_stage_end;
    return;
  }
  /*
   * Recovery
   *   If the number of wildcards (or errors required) is too high we try to recover as
   *   many matches as possible. We extract feasible regions from the read and filter
   *   them trying to recover anything out of bad quality reads.
   */
  if (num_wildcards > 0 && num_wildcards >= search->search_parameters->alignment_max_error_nominal) {
    search->search_stage = asearch_stage_read_recovery;
    return;
  }
//  // Exact search
//  if (search->max_complete_error==0) {
//    search->search_stage = asearch_stage_neighborhood;
//    return;
//  }
//  // Very short reads (Neighborhood search)
//  if (key_length <= RANK_MTABLE_SEARCH_DEPTH || key_length < search->archive->fm_index->proper_length) {
//    search->search_stage = asearch_stage_neighborhood;
//    return;
//  }
  // Otherwise, go to standard exact filtering
  search->search_stage = asearch_stage_filtering_adaptive;
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); }
  return;
}
/*
 * Adaptive mapping [GEM-workflow 4.0]
 */
void approximate_search_filtering_adaptive(
    approximate_search_t* const search,
    matches_t* const matches) {
  // Parameters
  const search_parameters_t* const search_parameters = search->search_parameters;
  // Process proper search-stage
  while (true) {
    switch (search->search_stage) {
      case asearch_stage_begin: // Search Start. Check basic cases
        approximate_search_filtering_adaptive_basic_cases(search);
        break;
      case asearch_stage_read_recovery: // Read recovery
        approximate_search_exact_filtering_adaptive_heavyweight(search,matches);
        asearch_control_next_state_read_recovery(search,matches); // Next State
        break;
      case asearch_stage_filtering_adaptive: // Exact-Filtering (Adaptive)
        if (search_parameters->mapping_mode==mapping_adaptive_filtering_fast) {
          approximate_search_exact_filtering_adaptive_cutoff(search,matches);
        } else {
          approximate_search_exact_filtering_adaptive_lightweight(search,matches);
        }
        asearch_control_next_state_exact_filtering_adaptive(search,matches); // Next State
        break;
//      case asearch_stage_neighborhood:
//        approximate_search_neighborhood_search(search,matches);
//        search->search_stage = asearch_stage_end; // Next State
//        break;
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
/*
 * Test // FIXME Delete ME
 */
void approximate_search_filtering_adaptive_test(
    approximate_search_t* const search,
    matches_t* const matches) {
//  pattern_t* const pattern = &search->pattern;
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
//  // TODO
//  matches_add_interval_match(
//      matches,search->lo_exact_matches,search->hi_exact_matches,
//      pattern->key_length,0,search->emulated_rc_search);
  // Update MCS
  approximate_search_update_mcs(search,1);
  // Update next state
  search->search_stage = asearch_stage_end;
}
