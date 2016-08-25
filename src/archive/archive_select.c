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
  select_parameters_align->max_reported_matches = select_parameters_report->max_reported_matches;
}
void archive_select_configure_pe(archive_search_t* const archive_search) {
  search_parameters_t* const search_parameters = &archive_search->search_parameters;
  select_parameters_t* const select_parameters_align = &search_parameters->select_parameters_align;
  select_parameters_align->min_reported_strata = 0;
  select_parameters_align->max_reported_matches = 100;
}
/*
 * Select Matches
 */
void archive_select_se_matches_discard(
    archive_search_t* const archive_search,
    matches_t* const matches,
    const uint64_t matches_to_report,
    const uint64_t strata_to_report) {
  // Check limits on the maximum stratum & discard unwanted matches
  const uint64_t max_edit_distance = strata_to_report+1;
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  match_trace_t** match_trace_in = matches_get_match_traces(matches);
  match_trace_t** match_trace_out = match_trace_in;
  uint64_t i, num_matches_accepted = 0;
  for (i=0;i<num_matches;++i) {
    // Remove unwanted ones
    if (match_trace_in[i]->edit_distance > max_edit_distance) continue;
    *match_trace_out = *match_trace_in;
    ++match_trace_out; ++num_matches_accepted;
  }
  if (num_matches_accepted != num_matches) {
    matches_metrics_set_limited_candidates(&matches->metrics,true);
    vector_set_used(matches->match_traces,num_matches_accepted); // Update used
  }
}
void archive_select_se_matches(
    archive_search_t* const archive_search,
    select_parameters_t* const select_parameters,
    matches_t* const matches,
    const bool paired_end_search) {
  PROFILE_START(GP_ARCHIVE_SELECT_SE_MATCHES,PROFILE_LEVEL);
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  // Check min-reported-strata constrain
  if (select_parameters->min_reported_strata_nominal > 0) {
    // Calculate the number of matches to decode wrt input parameters
    uint64_t strata_to_report = 0, matches_to_report = 0;
    matches_counters_compute_matches_to_report(matches->counters,
        select_parameters->min_reported_strata_nominal,
        select_parameters->max_reported_matches,
        &matches_to_report,&strata_to_report);
    // Discard unwanted matches
    if (matches_to_report > 0) {
      archive_select_se_matches_discard(
          archive_search,matches,matches_to_report,strata_to_report);
    }
  } else if (num_matches > select_parameters->max_reported_matches) {
    // Discard unwanted
    vector_set_used(matches->match_traces,select_parameters->max_reported_matches);
    matches_metrics_set_limited_candidates(&matches->metrics,true);
  }
  PROFILE_STOP(GP_ARCHIVE_SELECT_SE_MATCHES,PROFILE_LEVEL);
}
/*
 * Select Paired-Matches
 */
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
    paired_matches_sort_by_swg_score(paired_matches);
  }
  // Discard unwanted
  const uint64_t num_paired_matches = vector_get_used(paired_matches->paired_maps);
  if (num_paired_matches > select_parameters->max_reported_matches) {
    vector_set_used(paired_matches->paired_maps,select_parameters->max_reported_matches);
  }
  PROFILE_STOP(GP_ARCHIVE_SELECT_PE_MATCHES,PROFILE_LEVEL);
}
