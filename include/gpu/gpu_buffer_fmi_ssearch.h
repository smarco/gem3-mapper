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
 *            Alejandro Chacon <alejandro.chacon@uab.es>
 * DESCRIPTION:
 *   GPU-adaptor module provides support data-structures and functions
 *   to offload the static search of key partitions (fixed region-profile)
 *   from a batch of searches to the GPU
 */

#ifndef GPU_BUFFER_FMI_SSEARCH_H_
#define GPU_BUFFER_FMI_SSEARCH_H_

#include "utils/essentials.h"
#include "profiler/profiler_timer.h"
#include "align/pattern/pattern.h"
#include "gpu/gpu_buffer_collection.h"

/*
 * GPU FM-Index Backward Search Buffer
 */
typedef struct {
  /* GPU Generic Buffer*/
  void* buffer;             // GPU Generic Buffer
  /* Buffer state */
  bool fmi_search_enabled;  // Enabled GPU fmi-search
  uint32_t num_queries;     // Buffer state
  /* Profile */
  gem_timer_t timer;
} gpu_buffer_fmi_ssearch_t;

/*
 * Setup
 */
gpu_buffer_fmi_ssearch_t* gpu_buffer_fmi_ssearch_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_pos,
    const bool fmi_search_enabled);
void gpu_buffer_fmi_ssearch_clear(gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch);
void gpu_buffer_fmi_ssearch_delete(gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch);

/*
 * Occupancy & Limits
 */
uint64_t gpu_buffer_fmi_ssearch_get_max_queries(gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch);
uint64_t gpu_buffer_fmi_ssearch_get_num_queries(gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch);
bool gpu_buffer_fmi_ssearch_fits_in_buffer(
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch,
    const uint64_t total_queries);

/*
 * Accessors
 */
void gpu_buffer_fmi_ssearch_add_query(
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch,
    pattern_t* const pattern,
    const uint64_t begin,
    const uint64_t end);
void gpu_buffer_fmi_ssearch_get_result(
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch,
    const uint64_t buffer_pos,
    uint64_t* const hi,
    uint64_t* const lo);

/*
 * Send/Receive
 */
void gpu_buffer_fmi_ssearch_send(gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch);
void gpu_buffer_fmi_ssearch_receive(gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch);

#endif /* GPU_BUFFER_FMI_SSEARCH_H_ */
