/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *  Copyright (c) 2011-2017 by Alejandro Chacon <alejandro.chacon@uab.es>
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
 */

#ifndef GPU_BUFFER_BPM_ALIGN_H_
#define GPU_BUFFER_BPM_ALIGN_H_

#include "utils/essentials.h"
#include "profiler/profiler_timer.h"
#include "align/pattern/pattern.h"
#include "gpu/gpu_buffer_collection.h"

/*
 * GPU BMP Buffer
 */
typedef struct {
  /* GPU Generic Buffer*/
  void* buffer;                       // GPU Generic Buffer
  /* Dimensions Hints */
  gem_counter_t query_length;         // Tracks queries' length
  gem_counter_t candidates_per_tile;  // Tracks number of candidates per tile
  uint32_t query_same_length;         // Tracks same-read-length buffers
  /* Buffer state */
  bool bpm_align_enabled;             // Enabled GPU-computing BPM-Distance
  uint32_t current_query_offset;
  uint32_t num_entries;
  uint32_t num_queries;
  uint32_t num_candidates;
  /* Profile */
  gem_timer_t timer;
} gpu_buffer_bpm_align_t;

/*
 * Setup
 */
gpu_buffer_bpm_align_t* gpu_buffer_bpm_align_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_no,
    const bool bpm_align_enabled);
void gpu_buffer_bpm_align_clear(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align);
void gpu_buffer_bpm_align_delete(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align);

/*
 * Occupancy & Limits
 */
uint64_t gpu_buffer_bpm_align_get_entry_length(void);
uint64_t gpu_buffer_bpm_align_get_max_candidates(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align);
uint64_t gpu_buffer_bpm_align_get_max_queries(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align);

uint64_t gpu_buffer_bpm_align_get_num_candidates(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align);
uint64_t gpu_buffer_bpm_align_get_num_queries(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align);

void gpu_buffer_bpm_align_compute_dimensions(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    pattern_t* const pattern,
    const uint64_t num_candidates,
    uint64_t* const total_entries,
    uint64_t* const total_queries,
    uint64_t* const total_candidates);
bool gpu_buffer_bpm_align_fits_in_buffer(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    const uint64_t total_entries,
    const uint64_t total_queries,
    const uint64_t total_candidates);

/*
 * Accessors
 */
void gpu_buffer_bpm_align_add_pattern(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    pattern_t* const pattern);
void gpu_buffer_bpm_align_add_candidate(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    const uint64_t tile_offset,
    const uint64_t candidate_text_position,
    const uint64_t candidate_length);

void gpu_buffer_bpm_align_get_candidate(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    const uint64_t candidate_offset,
    uint64_t* const candidate_text_position,
    uint32_t* const candidate_length);
void gpu_buffer_bpm_align_get_result(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    const uint64_t candidate_offset,
    uint32_t* const levenshtein_distance,
    uint32_t* const levenshtein_match_pos);

/*
 * Hints
 */
void gpu_buffer_bpm_align_record_query_length(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    const uint64_t query_length);
void gpu_buffer_bpm_align_record_candidates_per_tile(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    const uint64_t num_candidates);
uint64_t gpu_buffer_bpm_align_get_mean_query_length(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align);
uint64_t gpu_buffer_bpm_align_get_mean_candidates_per_tile(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align);

/*
 * Send/Receive
 */
void gpu_buffer_bpm_align_send(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align);
void gpu_buffer_bpm_align_receive(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align);

#endif /* GPU_BUFFER_BPM_ALIGN_H_ */
