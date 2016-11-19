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
 * DESCRIPTION:
 *   Archive select parameters module encapsulates the parameters to
 *   select/discard matches once the search is accomplish
 */

#ifndef ARCHIVE_SELECT_PARAMETERS_H_
#define ARCHIVE_SELECT_PARAMETERS_H_

#include "utils/essentials.h"

/*
 * Select Parameters
 */
typedef struct {
  double min_reported_strata;            // Minimum number of complete strata to be reported (prop. read-length)
  uint64_t min_reported_strata_nominal;  // Minimum number of complete strata to be reported (nominal)
  uint64_t max_reported_matches;         // Maximum number of reported matches
  uint64_t max_searched_matches;         // Maximum number of searched se-matches (internal)
  uint64_t max_searched_paired_matches;  // Maximum number of searched pe-matches (internal)
} select_parameters_t;

/*
 * Select Parameters Setup
 */
void select_parameters_init(select_parameters_t* const select_parameters);

void select_parameters_configure_reporting(
    select_parameters_t* const select_parameters,
    const float min_reported_strata,
    const uint64_t max_reported_matches);

void select_parameters_instantiate_values(
    select_parameters_t* const select_parameters,
    const uint64_t sequence_length);

#endif /* ARCHIVE_SELECT_PARAMETERS_H_ */
