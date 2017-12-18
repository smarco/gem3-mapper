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
 * DESCRIPTION:
 *   GPU-adaptor module provides support data-structures and functions
 *   to offload the decoding of bwt-positions to text-positions
 *   from a batch of searches to the GPU
 */

#ifndef GPU_BUFFER_FMI_DECODE_H_
#define GPU_BUFFER_FMI_DECODE_H_

#include "utils/essentials.h"
#include "profiler/profiler_timer.h"
#include "gpu/gpu_buffer_collection.h"

/*
 * GPU FM-Index Backward Search Buffer
 */
typedef struct {
  /* GPU Generic Buffer*/
  void* buffer;             // GPU Generic Buffer
  /* Buffer state */
  bool decode_sa_enabled;   // Enabled GPU Computing
  bool decode_text_enabled; // Enabled GPU Computing
  uint32_t num_queries;     // Buffer state
  uint32_t sampling_rate;   // Decode sampling-rate
  /* Profile */
  gem_timer_t timer;        // Internal timer
} gpu_buffer_fmi_decode_t;

/*
 * Setup
 */
gpu_buffer_fmi_decode_t* gpu_buffer_fmi_decode_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_pos,
    const uint32_t sampling_rate,
    const bool decode_sa_enabled,
    const bool decode_text_enabled);
void gpu_buffer_fmi_decode_clear(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode);
void gpu_buffer_fmi_decode_delete(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode);

/*
 * Occupancy & Limits
 */
uint64_t gpu_buffer_fmi_decode_get_max_queries(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode);
uint64_t gpu_buffer_fmi_decode_get_num_queries(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode);
bool gpu_buffer_fmi_decode_fits_in_buffer(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,
    const uint64_t total_queries);

/*
 * Accessors
 */
void gpu_buffer_fmi_decode_add_query(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,
    const uint64_t bwt_position);
void gpu_buffer_fmi_decode_get_query(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,
    const uint64_t buffer_pos,
    uint64_t* const position);

void gpu_buffer_fmi_decode_get_position_sa(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,
    const uint64_t buffer_pos,
    uint64_t* const bwt_sampled_position,
    uint64_t* const lf_steps);
void gpu_buffer_fmi_decode_get_position_text(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,
    const uint64_t buffer_pos,
    uint64_t* const text_position);

/*
 * Send/Receive
 */
void gpu_buffer_fmi_decode_send(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode);
void gpu_buffer_fmi_decode_receive(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode);

#endif /* GPU_BUFFER_FMI_DECODE_H_ */
