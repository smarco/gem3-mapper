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

#include "gpu/gpu_buffer_kmer_filter.h"

/*
 * Errors
 */
#define GEM_ERROR_GPU_KMER_FILTER_QUERY_BUFFER_SIZE \
  "GPU.Kmer.Filter. Query pattern (%"PRIu64" bases) exceeds maximum buffer capacity (%"PRIu64" bases)"
#define GEM_ERROR_GPU_KMER_FILTER_MAX_CANDIDATES \
  "GPU.Kmer.Filter. Number of candidates (%"PRIu64") exceeds maximum buffer capacity (%"PRIu64" candidates)"
#define GEM_ERROR_GPU_KMER_FILTER_MAX_QUERIES \
  "GPU.Kmer.Filter. Number of queries (%"PRIu64") exceeds maximum buffer capacity (%"PRIu64" queries)"

/*
 * Constants :: Buffer Hints
 */
#define GPU_KMER_FILTER_MIN_NUM_SAMPLES             1
#define GPU_KMER_FILTER_AVERAGE_QUERY_LENGTH      150
#define GPU_KMER_FILTER_CANDIDATES_PER_QUERY       20
#define GPU_KMER_FILTER_QUERY_PADDING               8

/*
 * CUDA Support
 */
#ifdef HAVE_CUDA
/*
 * Stats
 */
void gpu_buffer_kmer_filter_record_query_length(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter,
    const uint64_t query_length) {
  COUNTER_ADD(&gpu_buffer_kmer_filter->query_length,query_length);
}
void gpu_buffer_kmer_filter_record_candidates_per_query(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter,
    const uint64_t candidates_per_query) {
  // Update candidates per query
  PROF_ADD_COUNTER(GP_GPU_BUFFER_KMER_FILTER_CANDIDATES_PER_QUERY,candidates_per_query);
  COUNTER_ADD(&gpu_buffer_kmer_filter->candidates_per_query,candidates_per_query);
}
uint64_t gpu_buffer_kmer_filter_get_mean_query_length(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter) {
  if (COUNTER_GET_NUM_SAMPLES(&gpu_buffer_kmer_filter->query_length) >= GPU_KMER_FILTER_MIN_NUM_SAMPLES) {
    return (uint64_t)ceil(COUNTER_GET_MEAN(&gpu_buffer_kmer_filter->query_length));
  } else {
    return GPU_KMER_FILTER_AVERAGE_QUERY_LENGTH;
  }
}
uint64_t gpu_buffer_kmer_filter_get_mean_candidates_per_query(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter) {
  if (COUNTER_GET_NUM_SAMPLES(&gpu_buffer_kmer_filter->candidates_per_query) >= GPU_KMER_FILTER_MIN_NUM_SAMPLES) {
    return (uint64_t)ceil(COUNTER_GET_MEAN(&gpu_buffer_kmer_filter->candidates_per_query));
  } else {
    return GPU_KMER_FILTER_CANDIDATES_PER_QUERY;
  }
}
/*
 * Setup
 */
gpu_buffer_kmer_filter_t* gpu_buffer_kmer_filter_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_no,
    const bool kmer_filter_enabled) {
  PROF_START(GP_GPU_BUFFER_KMER_FILTER_ALLOC);
  // Allocate handler
  gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter = mm_alloc(gpu_buffer_kmer_filter_t);
  // Stats
  COUNTER_RESET(&gpu_buffer_kmer_filter->query_length);
  COUNTER_RESET(&gpu_buffer_kmer_filter->candidates_per_query);
  TIMER_RESET(&gpu_buffer_kmer_filter->timer);
  // Buffer state
  gpu_buffer_kmer_filter->kmer_filter_enabled = kmer_filter_enabled;
  gpu_buffer_kmer_filter->buffer = gpu_buffer_collection_get_buffer(gpu_buffer_collection,buffer_no);
  gpu_kmer_filter_init_buffer_(gpu_buffer_kmer_filter->buffer,
      gpu_buffer_kmer_filter_get_mean_query_length(gpu_buffer_kmer_filter),
      gpu_buffer_kmer_filter_get_mean_candidates_per_query(gpu_buffer_kmer_filter));
  gpu_buffer_kmer_filter->current_query_offset = 0;
  gpu_buffer_kmer_filter->current_candidates_added = 0;
  // Buffer Queries & Candidates
  gpu_buffer_kmer_filter->num_queries = 0;
  gpu_buffer_kmer_filter->num_candidates = 0;
  gpu_buffer_kmer_filter->query_buffer_offset = 0;
  // Allocate GPU-buffer
  const int64_t thread_id = gtid(); // Between [1,num_threads] (zero is master)
  gpu_alloc_buffer_(gpu_buffer_kmer_filter->buffer,thread_id);
  // Return
  PROF_STOP(GP_GPU_BUFFER_KMER_FILTER_ALLOC);
  return gpu_buffer_kmer_filter;
}
void gpu_buffer_kmer_filter_clear(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter) {
  // Init buffer
  gpu_kmer_filter_init_buffer_(gpu_buffer_kmer_filter->buffer,
      gpu_buffer_kmer_filter_get_mean_query_length(gpu_buffer_kmer_filter),
      gpu_buffer_kmer_filter_get_mean_candidates_per_query(gpu_buffer_kmer_filter));
  gpu_buffer_kmer_filter->current_query_offset = 0;
  gpu_buffer_kmer_filter->current_candidates_added = 0;
  gpu_buffer_kmer_filter->num_queries = 0;
  gpu_buffer_kmer_filter->query_buffer_offset = 0;
  gpu_buffer_kmer_filter->num_candidates = 0;
}
void gpu_buffer_kmer_filter_delete(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter) {
  mm_free(gpu_buffer_kmer_filter);
}
/*
 * Occupancy & Limits
 */
uint64_t gpu_buffer_kmer_filter_get_max_queries(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter) {
  return gpu_kmer_filter_buffer_get_max_queries_(gpu_buffer_kmer_filter->buffer);
}
uint64_t gpu_buffer_kmer_filter_get_query_buffer_size(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter) {
  return gpu_kmer_filter_buffer_get_max_qry_bases_(gpu_buffer_kmer_filter->buffer);
}
uint64_t gpu_buffer_kmer_filter_get_max_candidates(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter) {
  return gpu_kmer_filter_buffer_get_max_candidates_(gpu_buffer_kmer_filter->buffer);
}
uint64_t gpu_buffer_kmer_filter_get_num_queries(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter) {
  return gpu_buffer_kmer_filter->num_queries;
}
uint64_t gpu_buffer_kmer_filter_get_num_candidates(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter) {
  return gpu_buffer_kmer_filter->num_candidates;
}
/*
 * Dimensions
 */
void gpu_buffer_kmer_filter_compute_dimensions(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter,
    kmer_counting_nway_t* const kmer_counting,
    const uint64_t num_candidates,
    uint64_t* const total_queries,
    uint64_t* const total_queries_length,
    uint64_t* const total_candidates) {
  // Calculate dimensions
  *total_queries += kmer_counting->num_tiles;
  *total_queries_length += kmer_counting->key_length;
  *total_candidates += kmer_counting->num_tiles*num_candidates;
}
bool gpu_buffer_kmer_filter_fits_in_buffer(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter,
    const uint64_t total_queries,
    const uint64_t total_queries_length,
    const uint64_t total_candidates) {
  // Get Limits
  uint64_t max_queries = gpu_buffer_kmer_filter_get_max_queries(gpu_buffer_kmer_filter);
  uint64_t max_candidates = gpu_buffer_kmer_filter_get_max_candidates(gpu_buffer_kmer_filter);
  uint64_t query_buffer_size = gpu_buffer_kmer_filter_get_query_buffer_size(gpu_buffer_kmer_filter);
  // Compute total bases + padding (worst case)
  const uint64_t total_bases_padded = total_queries_length + 2*GPU_KMER_FILTER_QUERY_PADDING;
  // Check available space in buffer for the pattern
  if (gpu_buffer_kmer_filter->num_queries+total_queries > max_queries ||
      gpu_buffer_kmer_filter->num_candidates+total_candidates > max_candidates ||
      gpu_buffer_kmer_filter->query_buffer_offset+total_bases_padded > query_buffer_size) {
    // Check buffer occupancy
    if (gpu_buffer_kmer_filter->num_queries > 0) {
      return false; // Leave it to the next fresh buffer
    }
    // Reallocate buffer
    gpu_kmer_filter_init_and_realloc_buffer_(
        gpu_buffer_kmer_filter->buffer,
        total_bases_padded,total_candidates,total_queries);
    // Check reallocated buffer dimensions (error otherwise)
    max_queries = gpu_buffer_kmer_filter_get_max_queries(gpu_buffer_kmer_filter);
    gem_cond_fatal_error(total_queries > max_queries,
        GPU_KMER_FILTER_MAX_QUERIES,total_queries,max_queries);
    max_candidates = gpu_buffer_kmer_filter_get_max_candidates(gpu_buffer_kmer_filter);
    gem_cond_fatal_error(total_candidates > max_candidates,
        GPU_KMER_FILTER_MAX_CANDIDATES,total_candidates,max_candidates);
    query_buffer_size = gpu_buffer_kmer_filter_get_query_buffer_size(gpu_buffer_kmer_filter);
    gem_cond_fatal_error(total_bases_padded > query_buffer_size,
        GPU_KMER_FILTER_QUERY_BUFFER_SIZE,total_bases_padded,query_buffer_size);
    // Return OK (after reallocation)
    return true;
  }
  // Return OK (fits in buffer)
  return true;
}
/*
 * Accessors
 */
void gpu_buffer_kmer_filter_add_pattern(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter,
    kmer_counting_nway_t* const kmer_counting) {
  // Stats
  if (gpu_buffer_kmer_filter->current_candidates_added != 0) {
    gpu_buffer_kmer_filter_record_candidates_per_query(
        gpu_buffer_kmer_filter,gpu_buffer_kmer_filter->current_candidates_added);
    gpu_buffer_kmer_filter->current_candidates_added = 0;
  }
  // Add query (all tiles)
  void* const gpu_buffer = gpu_buffer_kmer_filter->buffer;
  gpu_kmer_filter_qry_entry_t* const gpu_kmer_qry_entry = gpu_kmer_filter_buffer_get_queries_(gpu_buffer);
  gpu_kmer_filter_qry_info_t* const gpu_kmer_qry_info =
      gpu_kmer_filter_buffer_get_qry_info_(gpu_buffer) + gpu_buffer_kmer_filter->num_queries;
  const uint64_t num_tiles = kmer_counting->num_tiles;
  uint64_t i, query_buffer_offset;
  query_buffer_offset = gpu_buffer_kmer_filter->query_buffer_offset;
  for (i=0;i<num_tiles;++i) {
    // Copy query
    const uint64_t query_size = kmer_counting->key_tiles[i].end-kmer_counting->key_tiles[i].begin;
    memcpy(gpu_kmer_qry_entry+query_buffer_offset,
        kmer_counting->key+kmer_counting->key_tiles[i].begin,query_size);
    gpu_buffer_kmer_filter_record_query_length(gpu_buffer_kmer_filter,query_size); // Stats
    // Set query info
    gpu_kmer_qry_info[i].init_offset = query_buffer_offset;
    gpu_kmer_qry_info[i].query_size = query_size;
    // Account for padding
    query_buffer_offset += query_size;
    const uint64_t num_bases_mod = query_size % GPU_KMER_FILTER_QUERY_PADDING;
    if (num_bases_mod > 0) {
      query_buffer_offset += GPU_KMER_FILTER_QUERY_PADDING - num_bases_mod;
    }
  }
  gpu_buffer_kmer_filter->query_buffer_offset = query_buffer_offset;
  gpu_buffer_kmer_filter->current_query_offset = gpu_buffer_kmer_filter->num_queries;
  gpu_buffer_kmer_filter->num_queries += num_tiles;
}
void gpu_buffer_kmer_filter_add_candidate(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter,
    kmer_counting_nway_t* const kmer_counting,
    const uint64_t candidate_text_position) {
  // Add candidates
  gpu_kmer_filter_cand_info_t* const gpu_kmer_cand_info =
      gpu_kmer_filter_buffer_get_candidates_(gpu_buffer_kmer_filter->buffer) + gpu_buffer_kmer_filter->num_candidates;
  const uint64_t num_tiles = kmer_counting->num_tiles;
  uint64_t i;
  for (i=0;i<num_tiles;++i) {
    gpu_kmer_cand_info[i].position = candidate_text_position + kmer_counting->text_tiles[i].text_begin;
    gpu_kmer_cand_info[i].query = gpu_buffer_kmer_filter->current_query_offset + i;
    gpu_kmer_cand_info[i].size = kmer_counting->text_tiles[i].text_end - kmer_counting->text_tiles[i].text_begin;
  }
  gpu_buffer_kmer_filter->num_candidates += num_tiles;
  gpu_buffer_kmer_filter->current_candidates_added += num_tiles;
  // PROFILE
#ifdef GEM_PROFILE
  PROF_ADD_COUNTER(GP_GPU_BUFFER_KMER_FILTER_NUM_QUERIES,num_tiles);
  for (i=0;i<num_tiles;++i) {
    PROF_ADD_COUNTER(GP_GPU_BUFFER_KMER_FILTER_CANDIDATE_LENGTH,
        kmer_counting->text_tiles[i].text_end - kmer_counting->text_tiles[i].text_begin);
  }
#endif
}
uint32_t gpu_buffer_kmer_filter_get_min_edit_bound(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter,
    const uint64_t candidate_base_offset,
    const uint64_t num_tiles) {
  // Get alignment result(s)
  gpu_kmer_filter_alg_entry_t* const gpu_kmer_alg_entry =
      gpu_kmer_filter_buffer_get_alignments_(gpu_buffer_kmer_filter->buffer) + candidate_base_offset;
  // Add the bound of all tiles
  uint64_t i, min_edit_bound = 0;
  for (i=0;i<num_tiles;++i) {
    min_edit_bound += gpu_kmer_alg_entry[i];
  }
  // Return
  return min_edit_bound;
}
/*
 * Send/Receive
 */
void gpu_buffer_kmer_filter_send(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter) {
  PROF_START(GP_GPU_BUFFER_KMER_FILTER_SEND);
#ifdef GEM_PROFILE
  const uint64_t max_candidates = gpu_buffer_kmer_filter_get_max_candidates(gpu_buffer_kmer_filter);
  const uint64_t max_queries = gpu_buffer_kmer_filter_get_max_queries(gpu_buffer_kmer_filter);
  const uint64_t query_buffer_size = gpu_buffer_kmer_filter_get_query_buffer_size(gpu_buffer_kmer_filter);
  const uint64_t used_candidates = gpu_buffer_kmer_filter->num_candidates;
  const uint64_t used_queries = gpu_buffer_kmer_filter->num_queries;
  const uint64_t query_buffer_offset = gpu_buffer_kmer_filter->query_buffer_offset;
  PROF_ADD_COUNTER(GP_GPU_BUFFER_KMER_FILTER_USAGE_CANDIDATES,(100*used_candidates)/max_candidates);
  PROF_ADD_COUNTER(GP_GPU_BUFFER_KMER_FILTER_USAGE_QUERIES,(100*used_queries)/max_queries);
  PROF_ADD_COUNTER(GP_GPU_BUFFER_KMER_FILTER_USAGE_QUERY_BUFFER,(100*query_buffer_offset)/query_buffer_size);
  TIMER_START(&gpu_buffer_kmer_filter->timer);
#endif
  // Select computing device
  if (gpu_buffer_kmer_filter->kmer_filter_enabled) {
    if (gpu_buffer_kmer_filter->num_candidates > 0) {
      gpu_kmer_filter_send_buffer_(
          gpu_buffer_kmer_filter->buffer,
          gpu_buffer_kmer_filter->query_buffer_offset,
          gpu_buffer_kmer_filter->num_queries,
          gpu_buffer_kmer_filter->num_candidates,
          100); // FIXME
    }
  }
	PROF_STOP(GP_GPU_BUFFER_KMER_FILTER_SEND);
}
void gpu_buffer_kmer_filter_receive(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter) {
  PROF_START(GP_GPU_BUFFER_KMER_FILTER_RECEIVE);
  if (gpu_buffer_kmer_filter->kmer_filter_enabled) {
    if (gpu_buffer_kmer_filter->num_candidates > 0) {
      gpu_kmer_filter_receive_buffer_(gpu_buffer_kmer_filter->buffer);
    }
  }
  PROF_STOP(GP_GPU_BUFFER_KMER_FILTER_RECEIVE);
#ifdef GEM_PROFILE
  TIMER_STOP(&gpu_buffer_kmer_filter->timer);
  COUNTER_ADD(&PROF_GET_TIMER(GP_GPU_BUFFER_KMER_FILTER_DUTY_CYCLE)->time_ns,
      gpu_buffer_kmer_filter->timer.accumulated);
#endif
}
/*
 * CUDA NOT-Supported
 */
#else
/*
 * Setup
 */
gpu_buffer_kmer_filter_t* gpu_buffer_kmer_filter_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_no,
    const bool kmer_filter_enabled) { GEM_CUDA_NOT_SUPPORTED(); return NULL; }
void gpu_buffer_kmer_filter_clear(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter) { GEM_CUDA_NOT_SUPPORTED(); }
void gpu_buffer_kmer_filter_delete(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter) { GEM_CUDA_NOT_SUPPORTED(); }
/*
 * Occupancy & Limits
 */
uint64_t gpu_buffer_kmer_filter_get_max_queries(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
uint64_t gpu_buffer_kmer_filter_get_query_buffer_size(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
uint64_t gpu_buffer_kmer_filter_get_max_candidates(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
uint64_t gpu_buffer_kmer_filter_get_num_queries(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
uint64_t gpu_buffer_kmer_filter_get_num_candidates(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
/*
 * Dimensions
 */
void gpu_buffer_kmer_filter_compute_dimensions(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter,
    kmer_counting_nway_t* const kmer_counting,
    const uint64_t num_candidates,
    uint64_t* const total_queries,
    uint64_t* const total_queries_length,
    uint64_t* const total_candidates) { GEM_CUDA_NOT_SUPPORTED(); }
bool gpu_buffer_kmer_filter_fits_in_buffer(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter,
    const uint64_t total_queries,
    const uint64_t total_queries_length,
    const uint64_t total_candidates) { GEM_CUDA_NOT_SUPPORTED(); return false; }
/*
 * Accessors
 */
void gpu_buffer_kmer_filter_add_pattern(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter,
    kmer_counting_nway_t* const kmer_counting) { GEM_CUDA_NOT_SUPPORTED(); }
void gpu_buffer_kmer_filter_add_candidate(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter,
    kmer_counting_nway_t* const kmer_counting,
    const uint64_t candidate_text_position) { GEM_CUDA_NOT_SUPPORTED(); }
uint32_t gpu_buffer_kmer_filter_get_min_edit_bound(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter,
    const uint64_t candidate_base_offset,
    const uint64_t num_tiles) { GEM_CUDA_NOT_SUPPORTED(); return 0;}
/*
 * Send/Receive
 */
void gpu_buffer_kmer_filter_send(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter) { GEM_CUDA_NOT_SUPPORTED(); }
void gpu_buffer_kmer_filter_receive(
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter) { GEM_CUDA_NOT_SUPPORTED(); }

#endif /* HAVE_CUDA */

