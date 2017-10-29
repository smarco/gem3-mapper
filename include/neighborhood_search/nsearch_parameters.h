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

#ifndef NSEARCH_PARAMETERS_H_
#define NSEARCH_PARAMETERS_H_

#include "utils/essentials.h"

typedef struct {
  // Filtering cut-offs
  bool matches_accuracy_cutoff;
  bool matches_max_searched_cutoff;
  uint64_t filtering_max_candidates_acc;
  // Filtering search-drop limits
  bool filtering_cutoff;
  uint64_t filtering_quick_th;
  uint64_t filtering_region_opt_th;
  uint64_t filtering_region_opt_steps;
} nsearch_parameters_t;

void nsearch_parameters_init(nsearch_parameters_t* const nsearch_parameters);

#endif /* NSEARCH_PARAMETERS_H_ */
