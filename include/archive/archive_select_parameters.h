/*
 * PROJECT: GEMMapper
 * FILE: archive_select_parameters.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef ARCHIVE_SELECT_PARAMETERS_H_
#define ARCHIVE_SELECT_PARAMETERS_H_

#include "utils/essentials.h"

/*
 * Select Parameters
 */
typedef struct {
  double min_reported_strata;
  uint64_t min_reported_strata_nominal;
  uint64_t min_reported_matches;
  uint64_t max_reported_matches;
} select_parameters_t;

/*
 * Select Parameters Setup
 */
void select_parameters_init(select_parameters_t* const select_parameters);

void select_parameters_configure_reporting(
    select_parameters_t* const select_parameters,
    const float min_decoded_strata,
    const uint64_t min_reported_matches,
    const uint64_t max_reported_matches);

void select_parameters_instantiate_values(
    select_parameters_t* const select_parameters,
    const uint64_t sequence_length);

#endif /* ARCHIVE_SELECT_PARAMETERS_H_ */
