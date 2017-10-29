/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "neighborhood_search/nsearch_levenshtein_control.h"
#include "neighborhood_search/nsearch_filtering.h"
#include "fm_index/fm_index_query.h"
#include "filtering/candidates/filtering_candidates.h"
#include "filtering/candidates/filtering_candidates_process.h"
#include "filtering/candidates/filtering_candidates_verify.h"
#include "filtering/candidates/filtering_candidates_align.h"

/*
 * Search cut-off(s)
 */
bool nsearch_levenshtein_candidates_cutoff(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t num_candidates,
    nsearch_query_t* const nsearch_query,
    nsearch_query_t* const next_nsearch_query) {
  // Parameters
  search_parameters_t* const search_parameters = nsearch_schedule->search_parameters;
  nsearch_parameters_t* const nsearch_parameters = &search_parameters->nsearch_parameters;
  // Enabled
  if (!nsearch_parameters->filtering_cutoff) return false;
  // Direct filtering step
  const uint64_t filtering_quick_th = nsearch_parameters->filtering_quick_th;
  if (num_candidates <= filtering_quick_th) return true; // Cut-off
  // Optimization steps (if number of candidates below threshold)
  const uint64_t filtering_region_opt_th = nsearch_parameters->filtering_region_opt_th;
  if (num_candidates <= filtering_region_opt_th) {
    if (nsearch_query->num_optimization_steps==0) {
      next_nsearch_query->num_optimization_steps = nsearch_parameters->filtering_region_opt_steps;
    } else {
      if (nsearch_query->num_optimization_steps == 1) return true; // Cut-off
      next_nsearch_query->num_optimization_steps = nsearch_query->num_optimization_steps - 1;
    }
  } else {
    next_nsearch_query->num_optimization_steps = 0;
  }
  // Default
  return false;
}
/*
 * Candidate interval limit
 */
void nsearch_levenshtein_control_limit_interval(
    select_parameters_t* const select_parameters,
    uint64_t* const lo,
    uint64_t* const hi) {
  // Check
  if (select_parameters->min_reported_strata_nominal==0) {
    const uint64_t num_candidates = *hi - *lo;
    if (num_candidates > select_parameters->max_searched_matches) {
      *hi = *lo + select_parameters->max_searched_matches;
    }
  }
}
/*
 * Search terminate branch
 */
uint64_t nsearch_levenshtein_terminate(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t text_position,
    uint64_t lo,
    uint64_t hi,
    const uint64_t align_distance,
    const bool full_alignment) {
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
  pattern_t* const pattern = nsearch_schedule->pattern;
  // Limit interval (if full alignment by means of NS)
  if (full_alignment) {
    nsearch_levenshtein_control_limit_interval(&search_parameters->select_parameters,&lo,&hi);
  }
  // Add candidates to filtering
  filtering_candidates_add_positions_from_interval(
      filtering_candidates,search_parameters,pattern,
      lo,hi,0,pattern->key_length,align_distance);
  return hi-lo;
#endif
}
uint64_t nsearch_levenshtein_scheduled_terminate(
    nsearch_schedule_t* const nsearch_schedule,
    nsearch_operation_t* const nsearch_operation,
    const uint64_t text_length,
    nsearch_query_t* const nsearch_query,
    const uint64_t align_distance,
    const bool full_alignment) {
  // Parameters
  fm_2interval_t* const fm_2interval = &nsearch_query->fm_2interval;
#ifdef NSEARCH_ENUMERATE
  // PROFILE
  PROF_ADD_COUNTER(GP_NS_SEARCH_DEPTH,text_length);
  PROF_ADD_COUNTER(GP_NS_BRANCH_CANDIDATES_GENERATED,1);
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
  // Limit interval (if full alignment by means of NS)
  if (full_alignment) {
    nsearch_levenshtein_control_limit_interval(
        &search_parameters->select_parameters,
        &fm_2interval->backward_lo,&fm_2interval->backward_hi);
  }
  // Add interval
  filtering_candidates_add_positions_from_interval(
      filtering_candidates,search_parameters,pattern,fm_2interval->backward_lo,
      fm_2interval->backward_hi,region_begin,region_end,align_distance);
  // Filtering
  nsearch_parameters_t* const nsearch_parameters = &nsearch_schedule->search_parameters->nsearch_parameters;
  const uint64_t max_searched_matches = search_parameters->select_parameters.max_searched_matches;
  const uint64_t max_candidates_acc = nsearch_parameters->filtering_max_candidates_acc;
  const uint64_t num_candidates =
      filtering_candidates_get_num_positions(nsearch_schedule->filtering_candidates);
  if (num_candidates >= max_searched_matches || num_candidates >= max_candidates_acc) {
    nsearch_filtering(nsearch_schedule);
  }
  // PROFILE
  PROF_ADD_COUNTER(GP_NS_SEARCH_DEPTH,text_length);
  PROF_ADD_COUNTER(GP_NS_BRANCH_CANDIDATES_GENERATED,num_candidates);
  // Return
  return fm_2interval->backward_hi-fm_2interval->backward_lo;
#endif
}
