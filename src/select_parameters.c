/*
 * PROJECT: GEMMapper
 * FILE: select_parameters.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "select_parameters.h"

/*
 * Select Parameters Setup
 */
void select_parameters_init(select_parameters_t* const select_parameters) {
  // MAPQ Score
  select_parameters->mapq_model = mapq_model_gem;
  select_parameters->mapq_threshold = 0;
  // Reporting
  select_parameters->min_reported_strata = 0;
  select_parameters->min_reported_matches = 10;
  select_parameters->max_reported_matches = 100;
  // Check
  select_parameters->check_type = archive_check_nothing;
}
void select_configure_reporting(
    select_parameters_t* const select_parameters,const float min_reported_strata,
    const uint64_t min_reported_matches,const uint64_t max_reported_matches) {
  // Reporting
  select_parameters->min_reported_strata = min_reported_strata;
  select_parameters->min_reported_matches = min_reported_matches;
  select_parameters->max_reported_matches = max_reported_matches;
}
void select_instantiate_values(
    select_parameters_t* const select_parameters,const uint64_t sequence_length) {
  select_parameters->min_reported_strata_nominal = integer_proportion(select_parameters->min_reported_strata,sequence_length);
}
