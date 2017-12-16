/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *  Copyright (c) 2013-2017 by Alejandro Chacon <alejandro.chacond@gmail.com>
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
 *            Alejandro Chacon <alejandro.chacond@gmail.com>
 */

#ifndef GPU_BUFFER_BPM_DISTANCE_H_
#define GPU_BUFFER_BPM_DISTANCE_H_

#include "utils/essentials.h"
#include "profiler/profiler_timer.h"
#include "align/pattern/pattern.h"
#include "gpu/gpu_buffer_collection.h"

/*
 * GPU BMP Buffer
 */
typedef struct {
  /* GPU Buffer */
  void* buffer;                       // GPU Generic Buffer
  /* Buffer state */
  bool bpm_distance_enabled;          // Enabled GPU-computing BPM-Distance
  uint32_t current_query_offset;      // Current query offset (Once a pattern is added)
  uint32_t current_candidates_added;  // Current number of candidates added
  /* Buffer Queries */
  uint32_t num_queries;               // Total queries
  uint32_t num_query_entries;         // Total query-entries (BPM encoded chunks)
  /* Buffer Candidates */
  uint32_t num_candidates;            // Total candidates
  /* Stats */
  gem_counter_t query_length;         // Tracks queries' length
  gem_counter_t candidates_per_tile;  // Tracks number of candidates per tile
  uint32_t query_same_length;         // Tracks same-read-length buffers
  gem_timer_t timer;
} gpu_buffer_bpm_distance_t;

/*
 * Setup
 */
gpu_buffer_bpm_distance_t* gpu_buffer_bpm_distance_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_no,
    const bool bpm_distance_enabled);
void gpu_buffer_bpm_distance_clear(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance);
void gpu_buffer_bpm_distance_delete(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance);

/*
 * Occupancy & Limits
 */
uint64_t gpu_buffer_bpm_distance_get_entry_length(void);
uint64_t gpu_buffer_bpm_distance_get_max_candidates(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance);
uint64_t gpu_buffer_bpm_distance_get_max_queries(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance);
uint64_t gpu_buffer_bpm_distance_get_num_candidates(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance);

/*
 * Dimensions
 */
void gpu_buffer_bpm_distance_compute_dimensions(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance,
    pattern_t* const pattern,
    const uint64_t num_candidates,
    uint64_t* const total_queries,
    uint64_t* const total_query_entries,
    uint64_t* const total_candidates);
bool gpu_buffer_bpm_distance_fits_in_buffer(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance,
    const uint64_t total_queries,
    const uint64_t total_query_entries,
    const uint64_t total_candidates);

/*
 * Pattern accessor
 */
void gpu_buffer_bpm_distance_add_pattern(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance,
    pattern_t* const pattern);

/*
 * Candidate accessor
 */
void gpu_buffer_bpm_distance_add_candidate(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance,
    const uint64_t tile_offset,
    const uint64_t candidate_text_position,
    const uint64_t candidate_length);
void gpu_buffer_bpm_distance_get_candidate(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance,
    const uint64_t candidate_offset,
    uint64_t* const candidate_text_position,
    uint32_t* const candidate_length);

/*
 * Alignment distance accessor
 */
void gpu_buffer_bpm_distance_get_distance(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance,
    const uint64_t candidate_offset,
    uint32_t* const levenshtein_distance,
    uint32_t* const levenshtein_match_pos);

/*
 * Send/Receive
 */
void gpu_buffer_bpm_distance_send(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance);
void gpu_buffer_bpm_distance_receive(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance);

#endif /* GPU_BUFFER_BPM_DISTANCE_H_ */
