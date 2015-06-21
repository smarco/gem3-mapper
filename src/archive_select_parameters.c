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
  // MAPQ Score
  select_parameters->mapq_model = mapq_model_gem;
  select_parameters->mapq_threshold = 0;
  // Reporting
  select_parameters->min_decoded_strata = 0;
  select_parameters->min_reported_matches = 2;
  select_parameters->max_reported_matches = 20;
  // Check
  select_parameters->check_correct = false;
  select_parameters->check_optimum = false;
  select_parameters->check_complete = false;
}
GEM_INLINE void archive_select_configure_reporting(
    select_parameters_t* const select_parameters,const float min_decoded_strata,
    const uint64_t min_reported_matches,const uint64_t max_reported_matches) {
  // Reporting
  select_parameters->min_decoded_strata = min_decoded_strata;
  select_parameters->min_reported_matches = min_reported_matches;
  select_parameters->max_reported_matches = max_reported_matches;
}
GEM_INLINE void archive_select_instantiate_values(
    select_parameters_t* const select_parameters,const uint64_t sequence_length) {
  select_parameters->min_decoded_strata_nominal = integer_proportion(select_parameters->min_decoded_strata,sequence_length);
}
