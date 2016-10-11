/*
 * PROJECT: GEMMapper
 * FILE: nsearch_levenshtein.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "neighborhood_search/nsearch_levenshtein_control.h"
#include "fm_index/fm_index_query.h"
#include "filtering/candidates/filtering_candidates.h"
#include "filtering/candidates/filtering_candidates_process.h"
#include "filtering/candidates/filtering_candidates_verify.h"
#include "filtering/candidates/filtering_candidates_align.h"

/*
 * Search candidates cut-off
 */
bool nsearch_levenshtein_candidates_cutoff(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t num_candidates,
    nsearch_query_t* const nsearch_query,
    nsearch_query_t* const next_nsearch_query) {
  // Parameters
  search_parameters_t* const search_parameters = nsearch_schedule->search_parameters;
  const uint64_t ns_filtering_threshold = search_parameters->region_profile_model.ns_filtering_threshold;
  // Check number of candidates
  if (num_candidates <= ns_filtering_threshold) {
    if (nsearch_query->max_steps==0) {
      next_nsearch_query->max_steps = search_parameters->region_profile_model.max_steps;
      return false;
    } else {
      if (nsearch_query->max_steps == 1) return true;
      next_nsearch_query->max_steps = nsearch_query->max_steps - 1;
      return false;
    }
  } else {
    next_nsearch_query->max_steps = 0;
    return false;
  }
}
/*
 * Standard search terminate search-branch
 */
uint64_t nsearch_levenshtein_terminate(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t text_position,
    uint64_t lo,
    uint64_t hi,
    const uint64_t align_distance) {
  // PROFILE
  PROF_ADD_COUNTER(GP_NS_SEARCH_DEPTH,text_position);
  PROF_ADD_COUNTER(GP_NS_BRANCH_CANDIDATES_GENERATED,(hi-lo));
#ifdef NSEARCH_ENUMERATE
  const uint8_t* const text = nsearch_schedule->pending_searches->text;
  dna_buffer_print(stdout,text,text_position+1,true);
  fprintf(stdout,"\n");
  return 1;
#else
  // Parameters
  filtering_candidates_t* const filtering_candidates = nsearch_schedule->filtering_candidates;
  search_parameters_t* const search_parameters = nsearch_schedule->search_parameters;
  select_parameters_t* const select_parameters = &search_parameters->select_parameters_align;
  pattern_t* const pattern = nsearch_schedule->pattern;
  // FIXME: Depending on the mapping-strategy regulate this (avoid on complete-search)
  // Limit the number of candidates (cases than can exponentially explode)
  if (select_parameters->min_reported_strata_nominal==0) {
    const uint64_t num_candidates = hi - lo;
    if (num_candidates > select_parameters->max_reported_matches) {
      hi = lo + select_parameters->max_reported_matches;
    }
  }
  // Add candidates to filtering
  bool limited;
  filtering_candidates_add_positions_from_interval(
      filtering_candidates,search_parameters,pattern,
      lo,hi,0,pattern->key_length,align_distance,&limited);
  return hi-lo;
#endif
}
/*
 * Scheduled search terminate search-branch
 */
uint64_t nsearch_levenshtein_scheduled_terminate(
    nsearch_schedule_t* const nsearch_schedule,
    nsearch_operation_t* const nsearch_operation,
    const uint64_t text_length,
    nsearch_query_t* const nsearch_query,
    const uint64_t align_distance) {
  // Parameters
  fm_2interval_t* const fm_2interval = &nsearch_query->fm_2interval;
#ifdef NSEARCH_ENUMERATE
  // PROFILE
  PROF_ADD_COUNTER(GP_NS_SEARCH_DEPTH,text_length);
  PROF_ADD_COUNTER(GP_NS_CANDIDATES_GENERATED,1);
  // Search Text
  nsearch_operation->text_position = text_length;
  nsearch_operation_state_print_global_text(stdout,nsearch_operation); // Print
  return 1;
#else
  // Parameters
  const bool forward_search = (nsearch_operation->search_direction==direction_forward);
  filtering_candidates_t* const filtering_candidates = nsearch_schedule->filtering_candidates;
  search_parameters_t* const search_parameters = nsearch_schedule->search_parameters;
  pattern_t* const pattern = nsearch_schedule->pattern;
  nsearch_operation->text_position = text_length;
  // Compute region boundaries
  const uint64_t max_error = nsearch_schedule->max_error;
  const uint64_t global_key_begin = nsearch_operation->global_key_begin;
  const uint64_t global_key_end = nsearch_operation->global_key_end;
  const uint64_t max_region_scope = text_length + max_error;
  uint64_t region_begin, region_end;
  if (forward_search) {
    region_begin = BOUNDED_SUBTRACTION(global_key_begin,max_error,0);
    region_end = global_key_begin + max_region_scope;
  } else {
    region_begin = BOUNDED_SUBTRACTION(global_key_end,max_region_scope,0);
    region_end = global_key_end + max_error;
  }
  // PROFILE
  PROF_ADD_COUNTER(GP_NS_SEARCH_DEPTH,text_length);
  PROF_ADD_COUNTER(GP_NS_BRANCH_CANDIDATES_GENERATED,
      (fm_2interval->backward_hi-fm_2interval->backward_lo));
  // Add interval
  bool limited;
  filtering_candidates_add_positions_from_interval(
      filtering_candidates,search_parameters,pattern,fm_2interval->backward_lo,
      fm_2interval->backward_hi,region_begin,region_end,align_distance,&limited);
  // Dynamic filtering
  if (nsearch_schedule->dynamic_filtering) {
    // Process+Verify candidates
    PROF_START(GP_NS_VERIFICATION);
    filtering_candidates_process_candidates(filtering_candidates,nsearch_schedule->pattern,false);
    filtering_candidates_verify_candidates(filtering_candidates,pattern);
    PROF_STOP(GP_NS_VERIFICATION);
    // Align
    PROF_START(GP_NS_ALIGN);
    filtering_candidates_align_candidates(filtering_candidates,
        pattern,false,false,nsearch_schedule->matches);
    PROF_STOP(GP_NS_ALIGN);
    // Check max-matches quit condition
    // TODO
  }
  // Return
  return fm_2interval->backward_hi-fm_2interval->backward_lo;
#endif
}
