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

#include "utils/region_set.h"

/*
 * Constants
 */
#define REGION_SET_NUN_INITIAL_INTERVALS 20

/*
 * Setup
 */
void region_set_init(region_set_t* const region_set) {
  region_set->region_intervals = vector_new(REGION_SET_NUN_INITIAL_INTERVALS,region_interval_t);
}
void region_set_clear(region_set_t* const region_set) {
  vector_clear(region_set->region_intervals);
}
void region_set_destroy(region_set_t* const region_set) {
  vector_delete(region_set->region_intervals);
}
/*
 * Adding
 */
void region_set_add(
    region_set_t* const region_set,
    const uint64_t begin_position,
    const uint64_t end_position) {
  // Allocate
  region_interval_t* region_interval;
  vector_alloc_new(region_set->region_intervals,region_interval_t,region_interval);
  // Add
  region_interval->begin_position = begin_position;
  region_interval->end_position = end_position;
}
