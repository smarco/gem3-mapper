/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_filtering.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "approximate_search/approximate_search_neighborhood.h"
#include "neighborhood_search/nsearch.h"
#include "fm_index/fm_index_search.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED


/*
 * Neighborhood Generation (Exact Search)
 */
void approximate_search_neighborhood_exact_search(
    approximate_search_t* const search,
    matches_t* const matches) {
  PROFILE_START(GP_AS_EXACT_SEARCH,PROFILE_LEVEL);
//  pattern_t* const pattern = &search->pattern;
  // FM-Index basic exact search
//  fm_index_bsearch(search->archive->fm_index,pattern->key,
//    pattern->key_length,&search->hi_exact_matches,&search->lo_exact_matches);
//  fm_index_bsearch_pure(search->archive->fm_index,pattern->key,
//      pattern->key_length,&search->hi_exact_matches,&search->lo_exact_matches);
//  fm_index_reverse_bsearch_pure(search->archive->fm_index,pattern->key,
//      pattern->key_length,&search->hi_exact_matches,&search->lo_exact_matches);
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
//  fm_index_bsearch(search->archive->fm_index,pattern->key,
//      pattern->key_length,&search->hi_exact_matches,&search->lo_exact_matches);
//  // Add to matches
//  matches_add_interval_match(matches,search->lo_exact_matches,
//      search->hi_exact_matches,pattern->key_length,0,search->emulated_rc_search);
  // Update MCS
  approximate_search_update_mcs(search,1);
  PROFILE_STOP(GP_AS_EXACT_SEARCH,PROFILE_LEVEL);
}
/*
 * Neighborhood Generation (Inexact Search)
 */
void approximate_search_neighborhood_inexact_search(
    approximate_search_t* const search,
    matches_t* const matches) {
  PROFILE_START(GP_AS_NEIGHBORHOOD_SEARCH,PROFILE_LEVEL);
  // Parameters, pattern & interval-set
  pattern_t* const pattern = &search->pattern;
  interval_set_t* const intervals_result = search->interval_set;
  // Basic search (Brute force mitigated by mrank_table)
  interval_set_clear(search->interval_set); // Clear
  neighborhood_search(search->archive->fm_index,pattern->key,pattern->key_length,
      search->max_complete_error,intervals_result,search->mm_stack);
  // Add results
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
  // matches_add_interval_set(matches,intervals_result,pattern->key_length,search->emulated_rc_search);
  // Update MCS
  approximate_search_update_mcs(search,search->max_complete_error+1);
  PROFILE_STOP(GP_AS_NEIGHBORHOOD_SEARCH,PROFILE_LEVEL);
}
/*
 * Neighborhood Search
 */
void approximate_search_neighborhood_search(
    approximate_search_t* const search,
    matches_t* const matches) {
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,"[GEM]>ASM::Neighborhood Search\n");
    tab_global_inc();
  }
  // Check max-differences allowed
  if (search->max_complete_error==0) {
    approximate_search_neighborhood_exact_search(search,matches);   // Exact Search
  } else {
    approximate_search_neighborhood_inexact_search(search,matches); // Basic brute force search
  }
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); }
}
