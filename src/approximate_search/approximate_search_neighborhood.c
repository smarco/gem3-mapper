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
#include "filtering/filtering_candidates_process.h"
#include "filtering/filtering_candidates_verify.h"
#include "filtering/filtering_candidates_align.h"

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
  pattern_t* const pattern = &search->pattern;
  uint64_t hi, lo;
  // FM-Index basic exact search
  //fm_index_bsearch(search->archive->fm_index,pattern->key,pattern->key_length,&hi,&lo);
  //fm_index_bsearch_pure(search->archive->fm_index,pattern->key,pattern->key_length,&hi,&lo);
  fm_index_reverse_bsearch_pure(search->archive->fm_index,pattern->key,pattern->key_length,&hi,&lo);
  //fm_index_reverse_bsearch_fb(search->archive->fm_index,pattern->key,pattern->key_length,&hi,&lo);
  //fm_index_reverse_bsearch_bf(search->archive->fm_index,pattern->key,pattern->key_length,&hi,&lo);

  printf("%lu\t%lu\n",lo,hi);

  // Add interval
  filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
  filtering_candidates_add_read_interval(filtering_candidates,
      search->search_parameters,lo,hi,pattern->key_length,0);
  // Process candidates
  filtering_candidates_process_candidates(search->filtering_candidates,&search->pattern);
  // Verify Candidates
  filtering_candidates_verify_candidates(filtering_candidates,pattern);
  // Align
  filtering_candidates_align_candidates(filtering_candidates,
      pattern,search->emulated_rc_search,false,false,matches);
  // Update MCS
  approximate_search_update_mcs(search,1);
  // Set MCS
  matches->max_complete_stratum = MIN(matches->max_complete_stratum,search->max_complete_stratum);
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
