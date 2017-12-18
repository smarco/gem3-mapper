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
 *   to offload the adaptive search of key partitions (adaptive region-profile)
 *   from a batch of searches to the GPU
 */

#ifndef GPU_BUFFER_FMI_SEARCH_H_
#define GPU_BUFFER_FMI_SEARCH_H_

#include "utils/essentials.h"
#include "profiler/profiler_timer.h"
#include "align/pattern/pattern.h"
#include "gpu/gpu_buffer_collection.h"

/*
 * GPU FM-Index Backward Search Buffer
 */
typedef struct {
  /* GPU Generic Buffer*/
  void* buffer;                  // GPU Generic Buffer
  bool fmi_search_enabled;       // Enabled GPU fmi-search
  /* Region-Profile Parameters */
  uint32_t occ_min_threshold;
  uint32_t extra_search_steps;
  uint32_t alphabet_size;
  /* Buffer state */
  uint32_t num_queries;          // Current number of queries
  uint32_t num_bases;            // Current number of bases
  uint32_t num_regions;          // Current number of regions
  /* Profile */
  gem_counter_t query_length;    // Tracks queries' length
  gem_timer_t timer;
} gpu_buffer_fmi_asearch_t;

/*
 * Setup
 */
gpu_buffer_fmi_asearch_t* gpu_buffer_fmi_asearch_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_pos,
    const bool fmi_search_enabled,
    const uint32_t occ_min_threshold,
    const uint32_t extra_search_steps,
    const uint32_t alphabet_size);
void gpu_buffer_fmi_asearch_clear(gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch);
void gpu_buffer_fmi_asearch_delete(gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch);

/*
 * Occupancy & Limits
 */
uint64_t gpu_buffer_fmi_asearch_get_max_queries(gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch);
uint64_t gpu_buffer_fmi_asearch_get_max_bases(gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch);
uint64_t gpu_buffer_fmi_asearch_get_max_regions(gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch);

uint64_t gpu_buffer_fmi_asearch_get_num_queries(gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch);
uint64_t gpu_buffer_fmi_asearch_get_num_bases(gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch);
uint64_t gpu_buffer_fmi_asearch_get_num_regions(gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch);

bool gpu_buffer_fmi_asearch_fits_in_buffer(
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch,
    const uint64_t total_queries,
    const uint64_t total_bases,
    const uint64_t total_regions);

uint32_t gpu_buffer_fmi_asearch_get_mean_query_length(
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch);

/*
 * Accessors
 */
uint64_t gpu_buffer_fmi_asearch_add_query(
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch,
    pattern_t* const pattern,
    const uint64_t max_regions);
void gpu_buffer_fmi_asearch_get_result_total_regions(
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch,
    const uint64_t query_offset,
    uint64_t* const regions_offset,
    uint64_t* const num_regions);
void gpu_buffer_fmi_asearch_get_result_region(
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch,
    const uint64_t region_offset,
    uint64_t* const region_begin,
    uint64_t* const region_end,
    uint64_t* const region_hi,
    uint64_t* const region_lo);

/*
 * Send/Receive
 */
void gpu_buffer_fmi_asearch_send(gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch);
void gpu_buffer_fmi_asearch_receive(gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch);

#endif /* GPU_BUFFER_FMI_SEARCH_H_ */
