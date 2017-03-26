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

#ifndef SEARCH_STATE_DECODE_H_
#define SEARCH_STATE_DECODE_H_


#include "utils/essentials.h"
#include "archive/search/archive_search.h"
#include "search_pipeline/search_stage_state.h"
#include "search_pipeline/search_pipeline_handlers.h"

/*
 * Search-group Decode Candidates
 */
typedef struct {
  // Configuration
  search_stage_mode_t search_stage_mode;                    // Stage Mode (Sending/Receiving)
  // Decode Candidates Buffers
  vector_t* buffers;                                        // Verify Candidates Buffers (search_stage_decode_buffer_t*)
  search_stage_iterator_t iterator;                         // Buffers State
  /* Support Data Structures */
  search_pipeline_handlers_t* search_pipeline_handlers;
} search_stage_decode_t;

/*
 * Setup
 */
search_stage_decode_t* search_stage_decode_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffers_offset,
    const uint64_t num_buffers,
    const uint32_t sampling_rate,
    const bool decode_sa_enabled,
    const bool decode_text_enabled,
    search_pipeline_handlers_t* const search_pipeline_handlers);
void search_stage_decode_clear(
    search_stage_decode_t* const search_stage);
void search_stage_decode_delete(
    search_stage_decode_t* const search_stage);

/*
 * Send Searches (buffered)
 */
bool search_stage_decode_send_se_search(
    search_stage_decode_t* const search_stage,
    archive_search_t* const archive_search);
bool search_stage_decode_send_pe_search(
    search_stage_decode_t* const search_stage,
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2);

/*
 * Retrieve Searches (buffered)
 */
bool search_stage_decode_retrieve_finished(
    search_stage_decode_t* const search_stage);
bool search_stage_decode_retrieve_se_search(
    search_stage_decode_t* const search_stage,
    archive_search_t** const archive_search);
bool search_stage_decode_retrieve_pe_search(
    search_stage_decode_t* const search_stage,
    archive_search_t** const archive_search_end1,
    archive_search_t** const archive_search_end2);


#endif /* SEARCH_STATE_DECODE_H_ */
