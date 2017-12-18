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

#ifndef GPU_BUFFER_KMER_FILTER_H_
#define GPU_BUFFER_KMER_FILTER_H_

#include "utils/essentials.h"
#include "profiler/profiler_timer.h"
#include "align/pattern/pattern.h"
#include "align/align_kmer_filter_nway.h"
#include "gpu/gpu_buffer_collection.h"

/*
 * GPU BMP Buffer
 */
typedef struct {
  /* GPU Generic Buffer*/
  void* buffer;                       // GPU Generic Buffer
  /* Buffer state */
  bool kmer_filter_enabled;           // Enabled GPU-computing BPM-Distance
  uint32_t current_query_offset;      // Current query offset (Once a pattern is added)
  uint32_t current_candidates_added;  // Current number of candidates added
  /* Buffer Queries & Candidates */
  uint32_t num_queries;
  uint32_t query_buffer_offset;
  uint32_t num_candidates;
  /* Stats */
  gem_counter_t query_length;         // Tracks queries' length
  gem_counter_t candidates_per_query; // Tracks number of candidates per query
  gem_timer_t timer;
} gpu_buffer_kmer_filter_t;

/*
 * Setup
 */
gpu_buffer_kmer_filter_t* gpu_buffer_kmer_filter_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_no,
    const bool kmer_filter_enabled);
void gpu_buffer_kmer_filter_clear(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter);
void gpu_buffer_kmer_filter_delete(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter);

/*
 * Occupancy & Limits
 */
uint64_t gpu_buffer_kmer_filter_get_max_queries(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter);
uint64_t gpu_buffer_kmer_filter_get_query_buffer_size(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter);
uint64_t gpu_buffer_kmer_filter_get_max_candidates(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter);
uint64_t gpu_buffer_kmer_filter_get_num_queries(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter);
uint64_t gpu_buffer_kmer_filter_get_num_candidates(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter);

/*
 * Dimensions
 */
void gpu_buffer_kmer_filter_compute_dimensions(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter,
    kmer_counting_nway_t* const kmer_counting,
    const uint64_t num_candidates,
    uint64_t* const total_queries,
    uint64_t* const total_queries_length,
    uint64_t* const total_candidates);
bool gpu_buffer_kmer_filter_fits_in_buffer(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter,
    const uint64_t total_queries,
    const uint64_t total_queries_length,
    const uint64_t total_candidates);

/*
 * Accessors
 */
void gpu_buffer_kmer_filter_add_pattern(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter,
    kmer_counting_nway_t* const kmer_counting);
void gpu_buffer_kmer_filter_add_candidate(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter,
    kmer_counting_nway_t* const kmer_counting,
    const uint64_t candidate_text_position);

uint32_t gpu_buffer_kmer_filter_get_min_edit_bound(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter,
    const uint64_t candidate_base_offset,
    const uint64_t num_tiles);

/*
 * Send/Receive
 */
void gpu_buffer_kmer_filter_send(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter);
void gpu_buffer_kmer_filter_receive(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter);

#endif /* GPU_BUFFER_KMER_FILTER_H_ */
