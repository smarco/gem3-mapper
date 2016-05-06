/*
 * PROJECT: GEMMapper
 * FILE: archive_select_parameters.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive/archive_select_parameters.h"

/*
 * Archive Select Parameters Setup
 */
void select_parameters_init(select_parameters_t* const select_parameters) {
  // Reporting
  select_parameters->min_reported_strata = 0;
  select_parameters->min_reported_matches = 5;
  select_parameters->max_reported_matches = 5;
}
void select_parameters_configure_reporting(
    select_parameters_t* const select_parameters,
    const float min_reported_strata,
    const uint64_t min_reported_matches,
    const uint64_t max_reported_matches) {
  // Reporting
  select_parameters->min_reported_strata = min_reported_strata;
  select_parameters->min_reported_matches = min_reported_matches;
  select_parameters->max_reported_matches = max_reported_matches;
}
void select_parameters_instantiate_values(
    select_parameters_t* const select_parameters,
    const uint64_t sequence_length) {
  select_parameters->min_reported_strata_nominal =
      integer_proportion(select_parameters->min_reported_strata,sequence_length);
}
