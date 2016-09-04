/*
 * PROJECT: GEMMapper
 * FILE: nsearch_levenshtein.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "neighborhood_search/nsearch_levenshtein.h"
#include "neighborhood_search/nsearch_levenshtein_state.h"
#include "neighborhood_search/nsearch_partition.h"

/*
 * Query
 */
void nsearch_levenshtein_query(
    nsearch_schedule_t* const nsearch_schedule,
    nsearch_operation_t* const nsearch_operation,
    const uint64_t current_position,
    const uint8_t char_enc,
    const uint64_t lo_in,
    const uint64_t hi_in,
    uint64_t* const lo_out,
    uint64_t* const hi_out) {
#ifdef NSEARCH_ENUMERATE
  nsearch_schedule->pending_searches->text[current_position] = char_enc;
  *lo_out = 0; *hi_out = 1;
#else
  fm_index_t* const fm_index = nsearch_schedule->search->archive->fm_index;
  *lo_out = bwt_erank(fm_index->bwt,char_enc,lo_in);
  *hi_out = bwt_erank(fm_index->bwt,char_enc,hi_in);
#endif
}
uint64_t nsearch_levenshtein_terminate(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t text_position,
    uint64_t lo,
    uint64_t hi,
    const uint64_t align_distance) {
  // PROFILE
  PROF_ADD_COUNTER(GP_NS_SEARCH_DEPTH,text_position);
  PROF_ADD_COUNTER(GP_NS_CANDIDATES_GENERATED,(hi-lo));
#ifdef NSEARCH_ENUMERATE
  const uint8_t* const text = nsearch_schedule->pending_searches->text;
  dna_buffer_print(stdout,text,text_position+1,true);
  fprintf(stdout,"\n");
  return 1;
#else
  // Parameters
  filtering_candidates_t* const filtering_candidates = nsearch_schedule->search->filtering_candidates;
  search_parameters_t* const search_parameters = nsearch_schedule->search->search_parameters;
  select_parameters_t* const select_parameters = &search_parameters->select_parameters_align;
  pattern_t* const pattern = &nsearch_schedule->search->pattern;
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
  filtering_candidates_add_region_interval(
      filtering_candidates,search_parameters,pattern,
      lo,hi,0,pattern->key_length,align_distance,&limited);
  return hi-lo;
#endif
}
/*
 * Levenshtein Brute Force
 */
uint64_t nsearch_levenshtein_brute_force_step(
    nsearch_schedule_t* const nsearch_schedule,
    nsearch_operation_t* const nsearch_operation,
    const bool supercondensed,
    const uint64_t lo,
    const uint64_t hi) {
  // Parameters
  const uint64_t max_error = nsearch_schedule->max_error;
  const uint8_t* const key = nsearch_schedule->key;
  const uint64_t key_length = nsearch_schedule->key_length;
  const uint64_t text_position = nsearch_operation->text_position;
  const uint64_t max_text_length = key_length + max_error;
  uint64_t total_matches_found = 0;
  // Expand node for all characters
  uint64_t next_lo, next_hi;
  uint8_t char_enc;
  for (char_enc=0;char_enc<DNA_RANGE;++char_enc) {
    // Compute DP-next
    uint64_t min_val, align_distance;
    nsearch_levenshtein_state_compute_chararacter(
        &nsearch_operation->nsearch_state,false,key,key_length,
        text_position,char_enc,max_error,&min_val,&align_distance);
    if (min_val > max_error) continue;
    // Query
    nsearch_levenshtein_query(
        nsearch_schedule,nsearch_operation,
        text_position,char_enc,lo,hi,&next_lo,&next_hi);
    if (next_lo >= next_hi) continue;
    nsearch_operation->text_position = text_position + 1;
    // Keep searching
    if (supercondensed) {
      // Supercondensed Neighbourhood
      if (align_distance <= max_error) {
        total_matches_found += nsearch_levenshtein_terminate(
            nsearch_schedule,text_position,next_lo,next_hi,align_distance);
      } else if (text_position < max_text_length) {
        total_matches_found += nsearch_levenshtein_brute_force_step(
            nsearch_schedule,nsearch_operation,supercondensed,next_lo,next_hi);
      }
    } else {
      // Full Neighbourhood
      if (min_val <= align_distance && text_position < max_text_length) {
        total_matches_found += nsearch_levenshtein_brute_force_step(
            nsearch_schedule,nsearch_operation,supercondensed,next_lo,next_hi);
      }
      if (align_distance <= max_error) {
        total_matches_found += nsearch_levenshtein_terminate(
            nsearch_schedule,text_position,next_lo,next_hi,align_distance);
      }
    }
  }
  // Return matches found
  NSEARCH_PROF_NODE(nsearch_schedule,total_matches_found); // PROFILE
  return total_matches_found;
}
void nsearch_levenshtein_brute_force(
    approximate_search_t* const search,
    const bool supercondensed,
    matches_t* const matches) {
  PROF_START(GP_NS_GENERATION);
  // Init
  nsearch_schedule_t nsearch_schedule;
  nsearch_schedule_init(&nsearch_schedule,nsearch_model_levenshtein,search,matches);
  nsearch_operation_t* const nsearch_operation = nsearch_schedule.pending_searches;
  nsearch_levenshtein_state_prepare(&nsearch_operation->nsearch_state,supercondensed);
  nsearch_operation->text_position = 0;
#ifdef NSEARCH_ENUMERATE
  const uint64_t init_lo = 0;
  const uint64_t init_hi = 1;
#else
  const uint64_t init_lo = 0;
  const uint64_t init_hi = fm_index_get_length(nsearch_schedule.search->archive->fm_index);
#endif
  // Search
  nsearch_levenshtein_brute_force_step(&nsearch_schedule,nsearch_operation,supercondensed,init_lo,init_hi);
  // PROFILE
  // nsearch_schedule_print_profile(stderr,&nsearch_schedule);
  PROF_ADD_COUNTER(GP_NS_NODES,nsearch_schedule.profile.ns_nodes);
  PROF_ADD_COUNTER(GP_NS_NODES_SUCCESS,nsearch_schedule.profile.ns_nodes_success);
  PROF_ADD_COUNTER(GP_NS_NODES_FAIL,nsearch_schedule.profile.ns_nodes_fail);
  PROF_STOP(GP_NS_GENERATION);
}
/*
 * Levenshtein Neighborhood Search
 */
void nsearch_levenshtein(
    approximate_search_t* const search,
    matches_t* const matches) {
  // Search
  nsearch_schedule_t nsearch_schedule;
  nsearch_schedule_init(&nsearch_schedule,nsearch_model_levenshtein,search,matches);
  nsearch_schedule_search(&nsearch_schedule);
  // nsearch_schedule_print_profile(stderr,&nsearch_schedule); // PROFILE
}
/*
 * Neighborhood Search (Preconditioned by region profile)
 */
void nsearch_levenshtein_preconditioned(
    approximate_search_t* const search,
    matches_t* const matches) {
  // Search
  nsearch_schedule_t nsearch_schedule;
  nsearch_schedule_init(&nsearch_schedule,nsearch_model_levenshtein,search,matches);
  nsearch_schedule_search_preconditioned(&nsearch_schedule);
  // nsearch_schedule_print_profile(stderr,&nsearch_schedule); // PROFILE
}
/*
 * Display
 */
void nsearch_levenshtein_print_status(
    FILE* const stream,
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t pending_searches) {
  nsearch_operation_t* nsearch_operation = NULL;
  const int64_t first_idx_op = pending_searches;
  const int64_t last_idx_op = nsearch_schedule->num_pending_searches-1;
  int64_t i;
  // Print searched string
  for (i=last_idx_op;i>=first_idx_op;--i) {
    nsearch_operation = nsearch_schedule->pending_searches + i;
    fprintf(stream,"[%lu] key[%lu,%lu) Local-text=",last_idx_op-i,
        nsearch_operation->global_key_begin,nsearch_operation->global_key_end);
    nsearch_operation_state_print_local_text(stream,nsearch_operation);
    nsearch_operation_state_print(stream,nsearch_operation,nsearch_schedule->key);
  }
}
