/*
 * PROJECT: GEMMapper
 * FILE: archive_select.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive/archive_select.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Setup
 */
void archive_select_configure_se(archive_search_t* const archive_search) {
  search_parameters_t* const search_parameters = &archive_search->search_parameters;
  select_parameters_t* const select_parameters_report = &search_parameters->select_parameters_report;
  select_parameters_t* const select_parameters_align = &search_parameters->select_parameters_align;
  select_parameters_align->min_reported_strata = select_parameters_report->min_reported_strata;
  select_parameters_align->min_reported_matches = select_parameters_report->min_reported_matches;
  select_parameters_align->max_reported_matches = select_parameters_report->max_reported_matches;
}
void archive_select_configure_pe(archive_search_t* const archive_search) {
  search_parameters_t* const search_parameters = &archive_search->search_parameters;
  select_parameters_t* const select_parameters_align = &search_parameters->select_parameters_align;
  select_parameters_align->min_reported_strata = 0;
  select_parameters_align->min_reported_matches = 5;
  select_parameters_align->max_reported_matches = 5;
}
/*
 * Decoding Matches (Retrieving & Processing matches)
 */
void archive_select_process_trace_matches(
    archive_search_t* const archive_search,
    matches_t* const matches,
    const uint64_t reported_strata,
    uint64_t* const last_stratum_reported_matches) {
  const uint64_t num_initial_matches = vector_get_used(matches->position_matches);
  const uint64_t last_stratum_distance = reported_strata-1;
  // Count already decoded matches & discard unwanted matches
  uint64_t num_matches_last_stratum = 0;
  match_trace_t* position_matches_it = matches_get_match_trace_buffer(matches);
  VECTOR_ITERATE(matches->position_matches,match_trace,match_trace_num,match_trace_t) {
    if (match_trace->distance <= last_stratum_distance) {
      // Count matches last stratum
      if (match_trace->distance == last_stratum_distance) {
        if (num_matches_last_stratum >= *last_stratum_reported_matches) {
          continue; // Too many matches decoded in the last stratum
        }
        ++num_matches_last_stratum;
      }
      // Add the match (Store it in the vector, removing unwanted ones)
      if (position_matches_it != match_trace) *position_matches_it = *match_trace;
      ++position_matches_it;
    }
  }
  vector_update_used(matches->position_matches,position_matches_it); // Update used
  if (num_initial_matches != vector_get_used(matches->position_matches)) {
    // Because the matches had been reallocated, the indexed-positions are no longer valid
    matches_index_rebuild(matches,archive_search->mm_stack);
    matches_recompute_metrics(matches);
  }
  // Return matches to decode in the last stratum
  *last_stratum_reported_matches -= num_matches_last_stratum;
}
/*
 * Select Paired-Matches
 */
void archive_select_se_matches(
    archive_search_t* const archive_search,
    select_parameters_t* const select_parameters,
    matches_t* const matches) {
  PROFILE_START(GP_ARCHIVE_SELECT_SE_MATCHES,PROFILE_LEVEL);
  // Calculate the number of matches to decode wrt input parameters
  uint64_t reported_strata = 0, last_stratum_reported_matches = 0;
  matches_counters_compute_matches_to_decode(
      matches->counters,select_parameters->min_reported_strata_nominal,select_parameters->min_reported_matches,
      select_parameters->max_reported_matches,&reported_strata,&last_stratum_reported_matches);
  if (reported_strata==0) {
    // Remove all matches
    matches_get_clear_match_traces(matches);
  } else {
    // Decode matches (discard unwanted matches)
    archive_select_process_trace_matches(archive_search,
        matches,reported_strata,&last_stratum_reported_matches);
  }
  PROFILE_STOP(GP_ARCHIVE_SELECT_SE_MATCHES,PROFILE_LEVEL);
}
void archive_select_pe_matches(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    select_parameters_t* const select_parameters,
    paired_matches_t* const paired_matches) {
  // Update stats (Check number of paired-matches)
  const uint64_t num_matches = paired_matches_get_num_maps(paired_matches);
  if (num_matches==0) return;
  PROFILE_START(GP_ARCHIVE_SELECT_PE_MATCHES,PROFILE_LEVEL);
  // Sample unique
  if (num_matches==1) {
    const paired_map_t* const paired_map = paired_matches_get_maps(paired_matches);
    if (paired_map->pair_relation==pair_relation_concordant) {
      mapper_stats_template_length_sample(archive_search_end1->mapper_stats,paired_map->template_length);
    }
  } else {
    // Sort by distance (whichever it's selected)
    paired_matches_sort_by_distance(paired_matches);
  }
  // Discard surplus
  const uint64_t num_paired_matches = vector_get_used(paired_matches->paired_maps);
  if (num_paired_matches > select_parameters->max_reported_matches) {
    vector_set_used(paired_matches->paired_maps,select_parameters->max_reported_matches);
  }
  PROFILE_STOP(GP_ARCHIVE_SELECT_PE_MATCHES,PROFILE_LEVEL);
}
