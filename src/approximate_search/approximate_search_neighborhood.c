/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_filtering.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "approximate_search/approximate_search_neighborhood.h"
#include "fm_index/fm_index_search.h"
#include "filtering/filtering_candidates_process.h"
#include "filtering/filtering_candidates_verify.h"
#include "filtering/filtering_candidates_align.h"
#include "neighborhood_search/nsearch_hamming.h"
#include "neighborhood_search/nsearch_levenshtein.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Neighborhood Search Brute Force
 */
void approximate_search_neighborhood_search_brute_force(
    approximate_search_t* const search,
    matches_t* const matches) {
  PROFILE_START(GP_AS_NEIGHBORHOOD_SEARCH,PROFILE_LEVEL);
  // Parameters
  search_parameters_t* const search_parameters = search->search_parameters;
  filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
  pattern_t* const pattern = &search->pattern;
  // Generate Candidates (Select Alignment Model)
  if (search_parameters->alignment_model == alignment_model_hamming) {
    nsearch_hamming_brute_force(search,matches);
  } else {
    nsearch_levenshtein_brute_force(search,true,matches);
  }
  // Process candidates
  filtering_candidates_process_candidates(search->filtering_candidates,&search->pattern);
  // Verify Candidates
  filtering_candidates_verify_candidates(filtering_candidates,pattern);
  // Align
  filtering_candidates_align_candidates(filtering_candidates,pattern,false,false,matches);
  // Set MCS
  approximate_search_update_mcs(search,search->current_max_complete_error+1);
  matches->max_complete_stratum = MIN(matches->max_complete_stratum,search->current_max_complete_stratum);
  PROFILE_STOP(GP_AS_NEIGHBORHOOD_SEARCH,PROFILE_LEVEL);
}
/*
 * Neighborhood Search (Using partitions)
 */
void approximate_search_neighborhood_search_partition(
    approximate_search_t* const search,
    matches_t* const matches) {
  PROFILE_START(GP_AS_NEIGHBORHOOD_SEARCH,PROFILE_LEVEL);
  // Parameters
  search_parameters_t* const search_parameters = search->search_parameters;
  filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
  pattern_t* const pattern = &search->pattern;
  // Generate Candidates (Select Alignment Model)
  if (search_parameters->alignment_model == alignment_model_hamming) {
    nsearch_hamming(search,matches);
  } else {
    nsearch_levenshtein(search,matches);
  }
  // Process candidates
  filtering_candidates_process_candidates(search->filtering_candidates,&search->pattern);
  // Verify Candidates
  filtering_candidates_verify_candidates(filtering_candidates,pattern);
  // Align
  filtering_candidates_align_candidates(filtering_candidates,pattern,false,false,matches);
  // Set MCS
  approximate_search_update_mcs(search,search->current_max_complete_error+1);
  matches->max_complete_stratum = MIN(matches->max_complete_stratum,search->current_max_complete_stratum);
  PROFILE_STOP(GP_AS_NEIGHBORHOOD_SEARCH,PROFILE_LEVEL);
}
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
  // fm_index_bsearch(search->archive->fm_index,pattern->key,pattern->key_length,&hi,&lo);
  // fm_index_bsearch_pure(search->archive->fm_index,pattern->key,pattern->key_length,&hi,&lo);
  // fm_index_reverse_bsearch_pure(search->archive->fm_index,pattern->key,pattern->key_length,&hi,&lo);
  // fm_index_reverse_bsearch_fb(search->archive->fm_index,pattern->key,pattern->key_length,&hi,&lo);
  fm_index_reverse_bsearch_bf(search->archive->fm_index,pattern->key,pattern->key_length,&hi,&lo);
  // Add interval
  filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
  filtering_candidates_add_region_interval(filtering_candidates,
      search->search_parameters,pattern,lo,hi,0,pattern->key_length,0);
  // Process candidates
  filtering_candidates_process_candidates(search->filtering_candidates,&search->pattern);
  // Verify Candidates
  filtering_candidates_verify_candidates(filtering_candidates,pattern);
  // Align
  filtering_candidates_align_candidates(filtering_candidates,pattern,false,false,matches);
  // Set MCS
  approximate_search_update_mcs(search,1);
  matches->max_complete_stratum = MIN(matches->max_complete_stratum,search->current_max_complete_stratum);
  PROFILE_STOP(GP_AS_EXACT_SEARCH,PROFILE_LEVEL);
}
/*
 * Neighborhood Generation (Inexact Search)
 */
void approximate_search_neighborhood_inexact_search(
    approximate_search_t* const search,
    matches_t* const matches) {
  PROFILE_START(GP_AS_NEIGHBORHOOD_SEARCH,PROFILE_LEVEL);
  // TODO
  // Set MCS
  approximate_search_update_mcs(search,search->current_max_complete_error+1);
  matches->max_complete_stratum = MIN(matches->max_complete_stratum,search->current_max_complete_stratum);
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
  if (search->current_max_complete_error==0) {
    approximate_search_neighborhood_exact_search(search,matches); // Exact Search
  } else {
    approximate_search_neighborhood_inexact_search(search,matches); // Basic brute force search
  }
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); }
}
