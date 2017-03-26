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

#ifndef SEARCH_STAGE_DECODE_BUFFER_H_
#define SEARCH_STAGE_DECODE_BUFFER_H_

#include "utils/essentials.h"
#include "archive/search/archive_search.h"
#include "archive/search/archive_search_cache.h"
#include "gpu/gpu_buffer_fmi_decode.h"

/*
 * Search-Stage Verify Candidates Buffer
 */
typedef struct {
  gpu_buffer_fmi_decode_t* gpu_buffer_fmi_decode; // GPU-BPM-Buffer
  vector_t* archive_searches;                     // Vector of archive-searches (archive_search_t*)
} search_stage_decode_buffer_t;

/*
 * Setup
 */
search_stage_decode_buffer_t* search_stage_decode_buffer_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_no,
    const uint32_t sampling_rate,
    const bool decode_sa_enabled,
    const bool decode_text_enabled);
void search_stage_decode_buffer_clear(
    search_stage_decode_buffer_t* const decode_buffer);
void search_stage_decode_buffer_delete(
    search_stage_decode_buffer_t* const decode_buffer);

/*
 * Occupancy
 */
bool search_stage_decode_buffer_fits(
    search_stage_decode_buffer_t* const decode_buffer,
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2);

/*
 * Send/Receive
 */
void search_stage_decode_buffer_send(
    search_stage_decode_buffer_t* const decode_buffer);
void search_stage_decode_buffer_receive(
    search_stage_decode_buffer_t* const decode_buffer);

/*
 * Accessors
 */
void search_stage_decode_buffer_add(
    search_stage_decode_buffer_t* const decode_buffer,
    archive_search_t* const archive_search);
void search_stage_decode_buffer_retrieve(
    search_stage_decode_buffer_t* const decode_buffer,
    const uint64_t search_idx,
    archive_search_t** const archive_search);

#endif /* SEARCH_STAGE_DECODE_BUFFER_H_ */
