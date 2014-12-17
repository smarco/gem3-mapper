/*
 * PROJECT: GEMMapper
 * FILE: archive_select_parameters.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive_select_parameters.h"

/*
 * Select Parameters Setup
 */
GEM_INLINE void archive_select_parameters_init(select_parameters_t* const select_parameters) {
  // Reporting
  select_parameters->min_decoded_strata = 0;
  select_parameters->max_decoded_matches = 20;
  select_parameters->min_reported_matches = 1;
  select_parameters->max_reported_matches = 100;
  // Check
  select_parameters->check_matches_mask = check_none;
}
GEM_INLINE void archive_select_configure_reporting(
    select_parameters_t* const select_parameters,
    float min_decoded_strata,uint64_t max_decoded_matches,
    uint64_t min_reported_matches,uint64_t max_reported_matches) {
  // Reporting
  select_parameters->min_decoded_strata = min_decoded_strata;
  select_parameters->max_decoded_matches = max_decoded_matches;
  select_parameters->min_reported_matches = min_reported_matches;
  select_parameters->max_reported_matches = max_reported_matches;
}
GEM_INLINE void archive_select_configure_check(
    select_parameters_t* const select_parameters,const uint64_t check_matches_mask) {
  select_parameters->check_matches_mask = check_matches_mask;
}
GEM_INLINE void archive_select_instantiate_values(
    select_parameters_t* const select_parameters,const uint64_t sequence_length) {
  select_parameters->min_decoded_strata_nominal = integer_proportion(select_parameters->min_decoded_strata,sequence_length);
}
