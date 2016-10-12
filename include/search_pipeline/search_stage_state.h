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

#ifndef SEARCH_STAGE_STATE_H_
#define SEARCH_STAGE_STATE_H_

#include "utils/essentials.h"
#include "archive/search/archive_search.h"
#include "archive/search/archive_search_se_stepwise.h"

/*
 * Buffer State
 */
typedef enum {
  search_group_buffer_phase_sending,
  search_group_buffer_phase_retrieving
} search_stage_mode_t;
typedef struct {
  uint64_t num_buffers;           // Total buffers
  uint64_t current_buffer_idx;    // Current buffer index
  uint64_t num_searches;          // Total archive-searches (in current buffer)
  uint64_t current_search_idx;    // Current search index (in current buffer)
} search_stage_iterator_t;

#endif /* SEARCH_STAGE_STATE_H_ */
