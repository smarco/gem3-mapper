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

#ifndef REGION_SET_H_
#define REGION_SET_H_

#include "utils/essentials.h"

/*
 * Region Set
 */
typedef struct {
  uint64_t begin_position;
  uint64_t end_position;
} region_interval_t;
typedef struct {
  vector_t* region_intervals; // (region_interval_t*)
} region_set_t;

/*
 * Setup
 */
void region_set_init(region_set_t* const region_set);
void region_set_clear(region_set_t* const region_set);
void region_set_destroy(region_set_t* const region_set);

/*
 * Adding
 */
void region_set_add(
    region_set_t* const region_set,
    const uint64_t begin_position,
    const uint64_t end_position);

/*
 * Sort
 */
// TODO void region_set_sort(region_set_t* const region_set);

/*
 * Set operators
 */
// TODO bool region_set_is_contained(
//  region_set_t* const region_set,
//  const uint64_t begin_position,
//  const uint64_t end_position);

/*
 * Macro iterator
 */
#define REGION_SET_ITERATE(region_set,element) \
  VECTOR_ITERATE(region_set->region_intervals,element,region_set_##position,region_interval_t)


#endif /* REGION_SET_H_ */
