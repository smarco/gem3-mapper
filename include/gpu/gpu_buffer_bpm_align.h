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

#ifndef GPU_BUFFER_BPM_ALIGN_H_
#define GPU_BUFFER_BPM_ALIGN_H_

#include "utils/essentials.h"
#include "profiler/profiler_timer.h"
#include "align/pattern/pattern.h"
#include "gpu/gpu_buffer_collection.h"
#include "gpu/gpu_buffer_bpm_pattern.h"

/*
 * GPU BMP Buffer
 */
typedef struct {
  /* GPU Buffer */
  void* buffer;                         // GPU Generic Buffer
  /* Buffer state */
  bool bpm_align_enabled;               // Enabled GPU-computing BPM-Distance
  uint32_t current_query_idx;           // Current query index (Once a pattern is added)
  uint32_t current_query_buffer_offset; // Current query-buffer offset (Once a pattern is added)
  uint32_t current_candidates_added;    // Current number of candidates added
  /* Buffer Queries */
  uint32_t num_queries;                 // Total queries
  uint32_t num_query_entries;           // Total query-entries (BPM encoded chunks)
  uint32_t query_buffer_offset;         // Query-buffer offset (Plain text)
  /* Buffer Candidates */
  uint32_t num_candidates;              // Total candidates
  uint32_t candidate_buffer_offset;     // Candidate-buffer offset (Plain text)
  /* Buffer Cigars */
  uint32_t num_cigar_entries;           // Total cigar-entries
  /* Stats */
  gem_counter_t query_length;                  // Tracks queries' length
  gem_counter_t candidate_length;              // Tracks candidates' length
  gem_counter_t candidates_per_tile;           // Tracks number of candidates per tile
  gem_counter_t canonical_candidates_per_tile; // Tracks number of canonical candidates per tile
  uint32_t query_same_length;                  // Tracks same-read-length buffers
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
uint64_t gpu_buffer_bpm_align_get_max_queries(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align);
uint64_t gpu_buffer_bpm_align_get_max_query_entries(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align);
uint64_t gpu_buffer_bpm_align_get_max_query_buffer_size(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align);
uint64_t gpu_buffer_bpm_align_get_max_candidates(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align);
uint64_t gpu_buffer_bpm_align_get_max_candidate_length(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align);
uint64_t gpu_buffer_bpm_align_get_max_cigar_entries(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align);
uint64_t gpu_buffer_bpm_align_get_num_candidates(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align);

/*
 * Dimensions
 */
void gpu_buffer_bpm_align_compute_dimensions(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    pattern_t* const pattern,
    const uint64_t num_candidates,
    const uint64_t candidates_length,
    uint64_t* const total_queries,
    uint64_t* const total_query_entries,
    uint64_t* const total_query_length,
    uint64_t* const total_candidates);
bool gpu_buffer_bpm_align_fits_in_buffer(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    const uint64_t total_queries,
    const uint64_t total_query_entries,
    const uint64_t total_query_length,
    const uint64_t total_candidates);

/*
 * Pattern accessor
 */
void gpu_buffer_bpm_align_add_pattern(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    pattern_t* const pattern);

/*
 * Candidate accessor
 */
void gpu_buffer_bpm_align_add_candidate(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    const uint64_t candidate_text_position,
    const uint64_t tile_idx,
    const uint64_t tile_offset,
    const uint64_t candidate_length,
    const bool left_gap_align);

/*
 * Alignment accessor (CIGAR)
 */
void gpu_buffer_bpm_align_retrieve_alignment(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    uint8_t* const candidate_buffer,
    const uint64_t candidate_offset,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector);

/*
 * Send/Receive
 */
void gpu_buffer_bpm_align_send(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align);
void gpu_buffer_bpm_align_receive(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align);

#endif /* GPU_BUFFER_BPM_ALIGN_H_ */
