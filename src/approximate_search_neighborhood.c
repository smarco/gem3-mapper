/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_filtering.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "approximate_search_neighborhood.h"
#include "nsearch.h"
#include "fm_index_search.h"

/*
 * Neighborhood Generation (Exact Search)
 */
GEM_INLINE void approximate_search_neighborhood_exact_search(approximate_search_t* const search,matches_t* const matches) {
  PROF_START(GP_AS_EXACT_SEARCH);
  pattern_t* const pattern = &search->pattern;
  // FM-Index basic exact search
//  fm_index_bsearch(search->archive->fm_index,pattern->key,
//    pattern->key_length,&search->hi_exact_matches,&search->lo_exact_matches);
//  fm_index_bsearch_pure(search->archive->fm_index,pattern->key,
//      pattern->key_length,&search->hi_exact_matches,&search->lo_exact_matches);
//  fm_index_reverse_bsearch_pure(search->archive->fm_index,pattern->key,
//      pattern->key_length,&search->hi_exact_matches,&search->lo_exact_matches);
  fm_index_bsearch(search->archive->fm_index,pattern->key,
      pattern->key_length,&search->hi_exact_matches,&search->lo_exact_matches);
  // Add to matches
  matches_add_interval_match(matches,search->lo_exact_matches,
      search->hi_exact_matches,pattern->key_length,0,search->emulated_rc_search);
  // Update MCS
  approximate_search_update_mcs(search,1);
  // Update next state
  search->search_state = asearch_end;
  PROF_STOP(GP_AS_EXACT_SEARCH);
}
/*
 * Neighborhood Generation (Inexact Search)
 */
GEM_INLINE void approximate_search_neighborhood_inexact_search(approximate_search_t* const search,matches_t* const matches) {
  // Parameters, pattern & interval-set
  pattern_t* const pattern = &search->pattern;
  interval_set_t* const intervals_result = search->interval_set;
  // Basic search (Brute force mitigated by mrank_table)
  interval_set_clear(search->interval_set); // Clear
  neighborhood_search(search->archive->fm_index,pattern->key,pattern->key_length,
      search->max_complete_error,intervals_result,search->mm_stack);
  // Add results
  matches_add_interval_set(matches,intervals_result,pattern->key_length,search->emulated_rc_search);
  // Update MCS
  approximate_search_update_mcs(search,search->max_complete_error+1);
  // Update next state
  search->search_state = asearch_end;
}
/*
 * Neighborhood Search
 */
GEM_INLINE void approximate_search_neighborhood_search(approximate_search_t* const search,matches_t* const matches) {
  // Check max-differences allowed
  if (search->max_complete_error==0) {
    approximate_search_neighborhood_exact_search(search,matches);   // Exact Search
  } else {
    approximate_search_neighborhood_inexact_search(search,matches); // Basic brute force search
  }
  search->search_state = asearch_end; // Update search state
}
