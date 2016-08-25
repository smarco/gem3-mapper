/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_filtering.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "approximate_search/approximate_search_neighborhood.h"
#include "approximate_search/approximate_search_stages.h"
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
 * Exact Search
 */
void approximate_search_neighborhood_exact_search(
    approximate_search_t* const search,
    matches_t* const matches) {
  pattern_t* const pattern = &search->pattern;
  uint64_t hi, lo;
  // FM-Index basic exact search
  fm_index_bsearch(search->archive->fm_index,pattern->key,pattern->key_length,&hi,&lo);
  // fm_index_bsearch_pure(search->archive->fm_index,pattern->key,pattern->key_length,&hi,&lo);
  // fm_index_reverse_bsearch_pure(search->archive->fm_index,pattern->key,pattern->key_length,&hi,&lo);
  // fm_index_reverse_bsearch_fb(search->archive->fm_index,pattern->key,pattern->key_length,&hi,&lo);
  // fm_index_reverse_bsearch_bf(search->archive->fm_index,pattern->key,pattern->key_length,&hi,&lo);
  // Add interval
  bool limited;
  filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
  filtering_candidates_add_region_interval(
      filtering_candidates,search->search_parameters,pattern,
      lo,hi,0,pattern->key_length,0,&limited);
  // Process+Verify candidates
  PROF_START(GP_NS_VERIFICATION);
  filtering_candidates_process_candidates(search->filtering_candidates,&search->pattern,false);
  filtering_candidates_verify_candidates(filtering_candidates,pattern);
  PROF_STOP(GP_NS_VERIFICATION);
  // Align
  PROF_START(GP_NS_ALIGN);
  filtering_candidates_align_candidates(filtering_candidates,pattern,false,false,matches);
  PROF_STOP(GP_NS_ALIGN);
  // Finish Search
  approximate_search_update_mcs(search,1);
  approximate_search_end(search,matches);
}
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
  // Process+Verify candidates
  PROF_START(GP_NS_VERIFICATION);
  filtering_candidates_process_candidates(search->filtering_candidates,&search->pattern,false);
  filtering_candidates_verify_candidates(filtering_candidates,pattern);
  PROF_STOP(GP_NS_VERIFICATION);
  // Align
  PROF_START(GP_NS_ALIGN);
  filtering_candidates_align_candidates(filtering_candidates,pattern,false,false,matches);
  PROF_STOP(GP_NS_ALIGN);
  // Finish Search
  approximate_search_update_mcs(search,search->current_max_complete_error+1);
  approximate_search_end(search,matches);
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
  // Process+Verify candidates
  PROF_START(GP_NS_VERIFICATION);
  filtering_candidates_process_candidates(search->filtering_candidates,&search->pattern,false);
  filtering_candidates_verify_candidates(filtering_candidates,pattern);
  PROF_STOP(GP_NS_VERIFICATION);
  // Align
  PROF_START(GP_NS_ALIGN);
  filtering_candidates_align_candidates(filtering_candidates,pattern,false,false,matches);
  PROF_STOP(GP_NS_ALIGN);
  // Finish Search
  approximate_search_update_mcs(search,search->current_max_complete_error+1);
  approximate_search_end(search,matches);
  PROFILE_STOP(GP_AS_NEIGHBORHOOD_SEARCH,PROFILE_LEVEL);
}
/*
 * Neighborhood Search (Using partitions + region-profile preconditioned)
 */
void approximate_search_neighborhood_search_partition_preconditioned(
    approximate_search_t* const search,
    matches_t* const matches) {
  // Parameters
  search_parameters_t* const search_parameters = search->search_parameters;
  pattern_t* const pattern = &search->pattern;
  region_profile_t* const region_profile = &search->region_profile;
  filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
  const uint64_t max_complete_error = search->current_max_complete_error;
  const uint64_t mcs = search->current_max_complete_stratum;
  // Compute error limits
  region_profile_compute_error_limits(region_profile,mcs,max_complete_error);
  // Generate Candidates (Select Alignment Model)
  if (search_parameters->alignment_model == alignment_model_hamming) {
    nsearch_hamming_preconditioned(search,matches);
  } else {
    nsearch_levenshtein_preconditioned(search,matches);
  }
  // Process+Verify candidates
  PROF_START(GP_NS_VERIFICATION);
  filtering_candidates_process_candidates(search->filtering_candidates,&search->pattern,false);
  filtering_candidates_verify_candidates(filtering_candidates,pattern);
  PROF_STOP(GP_NS_VERIFICATION);
  // Align
  PROF_START(GP_NS_ALIGN);
  filtering_candidates_align_candidates(filtering_candidates,pattern,false,false,matches);
  PROF_STOP(GP_NS_ALIGN);
  // Update MCS
  approximate_search_update_mcs(search,max_complete_error+1);
}
