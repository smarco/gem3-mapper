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

#include "neighborhood_search/nsearch_parameters.h"

void nsearch_parameters_init(nsearch_parameters_t* const nsearch_parameters) {
  // Filtering cut-offs
  nsearch_parameters->matches_accuracy_cutoff = true;
  nsearch_parameters->matches_max_searched_cutoff = true;
  nsearch_parameters->filtering_cutoff = true;
  nsearch_parameters->filtering_max_candidates_acc = 1000;
  nsearch_parameters->filtering_quick_th = 0;
  nsearch_parameters->filtering_region_opt_th = 4;
  nsearch_parameters->filtering_region_opt_steps = 10;
}
