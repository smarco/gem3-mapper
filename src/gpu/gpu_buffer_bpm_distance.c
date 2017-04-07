/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *  Copyright (c) 2011-2017 by Alejandro Chacon <alejandro.chacond@gmail.com>
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
 */

#include "gpu/gpu_buffer_bpm_distance.h"
#include "resources/gpu_modules/gpu_interface.h"
#include "align/align_bpm_distance.h"

/*
 * Errors
 */
#define GEM_ERROR_GPU_BPM_DISTANCE_MAX_PATTERN_LENGTH \
  "GPU.BPM.Distance. Query pattern (%"PRIu64" entries) exceeds maximum buffer capacity (%"PRIu64" entries)"
#define GEM_ERROR_GPU_BPM_DISTANCE_MAX_CANDIDATES \
  "GPU.BPM.Distance. Number of candidates (%"PRIu64") exceeds maximum buffer capacity (%"PRIu64" candidates)"
#define GEM_ERROR_GPU_BPM_DISTANCE_MAX_QUERIES \
  "GPU.BPM.Distance. Number of queries (%"PRIu64") exceeds maximum buffer capacity (%"PRIu64" queries)"

/*
 * Constants :: GPU BMP-Distance
 */
#define GPU_BPM_DISTANCE_NUM_SUB_ENTRIES  GPU_BPM_FILTER_PEQ_SUBENTRIES
#define GPU_BPM_DISTANCE_ENTRY_LENGTH     GPU_BPM_FILTER_PEQ_ENTRY_LENGTH
#define GPU_BPM_DISTANCE_SUBENTRY_LENGTH  GPU_BPM_FILTER_PEQ_SUBENTRY_LENGTH
#define GPU_BPM_DISTANCE_ENTRY_SIZE       (GPU_BPM_DISTANCE_ENTRY_LENGTH/UINT8_LENGTH)
#define GPU_BPM_DISTANCE_ALPHABET_LENGTH  GPU_BPM_FILTER_PEQ_ALPHABET_SIZE

/*
 * Constants :: Buffer Hints
 */
#define GPU_BPM_DISTANCE_MIN_NUM_SAMPLES           1
#define GPU_BPM_DISTANCE_AVERAGE_QUERY_LENGTH      150
#define GPU_BPM_DISTANCE_CANDIDATES_PER_TILE       20

/*
 * CUDA Support
 */
#ifdef HAVE_CUDA
/*
 * Setup
 */
gpu_buffer_bpm_distance_t* gpu_buffer_bpm_distance_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_no,
    const bool bpm_distance_enabled) {
  PROF_START(GP_GPU_BUFFER_BPM_DISTANCE_ALLOC);
  // Alloc
  gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance = mm_alloc(gpu_buffer_bpm_distance_t);
  // Dimensions Hints
  COUNTER_RESET(&gpu_buffer_bpm_distance->query_length);
  COUNTER_RESET(&gpu_buffer_bpm_distance->candidates_per_tile);
  gpu_buffer_bpm_distance->query_same_length = 0;
  // Buffer state
  gpu_buffer_bpm_distance->bpm_distance_enabled = bpm_distance_enabled;
  gpu_buffer_bpm_distance->buffer = gpu_buffer_collection_get_buffer(gpu_buffer_collection,buffer_no);
  gpu_buffer_bpm_distance->num_entries = 0;
  gpu_buffer_bpm_distance->num_queries = 0;
  gpu_buffer_bpm_distance->num_candidates = 0;
  gpu_buffer_bpm_distance->current_query_offset = 0;
  TIMER_RESET(&gpu_buffer_bpm_distance->timer);
  // Init buffer
  const int64_t thread_id = gtid(); // Between [1,num_threads] (zero is master)
  gpu_alloc_buffer_(gpu_buffer_bpm_distance->buffer, thread_id);
  gpu_bpm_filter_init_buffer_(gpu_buffer_bpm_distance->buffer,
      gpu_buffer_bpm_distance_get_mean_query_length(gpu_buffer_bpm_distance),
      gpu_buffer_bpm_distance_get_mean_candidates_per_tile(gpu_buffer_bpm_distance));
  PROF_STOP(GP_GPU_BUFFER_BPM_DISTANCE_ALLOC);
  // Return
  return gpu_buffer_bpm_distance;
}
void gpu_buffer_bpm_distance_clear(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance) {
  // Init buffer
  gpu_bpm_filter_init_buffer_(gpu_buffer_bpm_distance->buffer,
      gpu_buffer_bpm_distance_get_mean_query_length(gpu_buffer_bpm_distance),
      gpu_buffer_bpm_distance_get_mean_candidates_per_tile(gpu_buffer_bpm_distance));
  // Dimensions Hints
  gpu_buffer_bpm_distance->query_same_length = 0;
  gpu_buffer_bpm_distance->num_entries = 0;
  gpu_buffer_bpm_distance->num_queries = 0;
  gpu_buffer_bpm_distance->num_candidates = 0;
  gpu_buffer_bpm_distance->current_query_offset = 0;
}
void gpu_buffer_bpm_distance_delete(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance) {
  mm_free(gpu_buffer_bpm_distance);
}
/*
 * Occupancy & Limits
 */
uint64_t gpu_buffer_bpm_distance_get_entry_length(void) {
  return GPU_BPM_DISTANCE_ENTRY_LENGTH;
}
uint64_t gpu_buffer_bpm_distance_get_max_candidates(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance) {
  return gpu_bpm_filter_buffer_get_max_candidates_(gpu_buffer_bpm_distance->buffer);
}
uint64_t gpu_buffer_bpm_distance_get_max_queries(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance) {
  return gpu_bpm_filter_buffer_get_max_queries_(gpu_buffer_bpm_distance->buffer);
}
uint64_t gpu_buffer_bpm_distance_get_max_entries(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance) {
  return gpu_bpm_filter_buffer_get_max_peq_entries_(gpu_buffer_bpm_distance->buffer);
}
uint64_t gpu_buffer_bpm_distance_get_num_candidates(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance) {
  return gpu_buffer_bpm_distance->num_candidates;
}
uint64_t gpu_buffer_bpm_distance_get_num_queries(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance) {
  return gpu_buffer_bpm_distance->num_queries;
}
void gpu_buffer_bpm_distance_compute_dimensions(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance,
    pattern_t* const pattern,
    const uint64_t num_candidates,
    uint64_t* const total_entries,
    uint64_t* const total_queries,
    uint64_t* const total_candidates) {
  // Calculate dimensions
  const uint64_t num_entries = DIV_CEIL(pattern->key_length,GPU_BPM_DISTANCE_ENTRY_LENGTH);
  *total_entries += num_entries;
  *total_queries += pattern->pattern_tiled.num_tiles;
  *total_candidates += pattern->pattern_tiled.num_tiles*num_candidates;
}
bool gpu_buffer_bpm_distance_fits_in_buffer(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance,
    const uint64_t total_entries,
    const uint64_t total_queries,
    const uint64_t total_candidates) {
  // Get Limits
  uint64_t max_PEQ_entries =  gpu_buffer_bpm_distance_get_max_entries(gpu_buffer_bpm_distance);
  uint64_t max_queries = gpu_buffer_bpm_distance_get_max_queries(gpu_buffer_bpm_distance);
  uint64_t max_candidates = gpu_buffer_bpm_distance_get_max_candidates(gpu_buffer_bpm_distance);
  // Check available space in buffer for the pattern
  if (gpu_buffer_bpm_distance->num_queries+total_queries > max_queries ||
      gpu_buffer_bpm_distance->num_entries+total_entries > max_PEQ_entries ||
      gpu_buffer_bpm_distance->num_candidates+total_candidates > max_candidates) {
    // Check buffer occupancy
    if (gpu_buffer_bpm_distance->num_queries > 0) {
      return false; // Leave it to the next fresh buffer
    }
    // Reallocate buffer
    gpu_bpm_filter_init_and_realloc_buffer_(
        gpu_buffer_bpm_distance->buffer,
        total_entries,total_candidates,total_queries);
    // Check reallocated buffer dimensions (error otherwise)
    max_queries = gpu_buffer_bpm_distance_get_max_queries(gpu_buffer_bpm_distance);
    gem_cond_fatal_error(total_queries > max_queries,
        GPU_BPM_DISTANCE_MAX_QUERIES,total_queries,max_queries);
    max_PEQ_entries = gpu_buffer_bpm_distance_get_max_entries(gpu_buffer_bpm_distance);
    gem_cond_fatal_error(total_entries > max_PEQ_entries,
        GPU_BPM_DISTANCE_MAX_PATTERN_LENGTH,total_entries,max_PEQ_entries);
    max_candidates = gpu_buffer_bpm_distance_get_max_candidates(gpu_buffer_bpm_distance);
    gem_cond_fatal_error(total_candidates > max_candidates,
        GPU_BPM_DISTANCE_MAX_CANDIDATES,total_candidates,max_candidates);
    // Return OK (after reallocation)
    return true;
  }
  // Return OK (fits in buffer)
  return true;
}
/*
 * Pattern Compile/Decompile
 */
void gpu_buffer_bpm_distance_pattern_decompile(
    gpu_bpm_filter_qry_entry_t* const pattern_entry,
    const uint32_t pattern_length,
    bpm_pattern_t* const bpm_pattern,
    mm_allocator_t* const mm_allocator) {
  // Calculate dimensions
  const uint64_t pattern_num_words = DIV_CEIL(pattern_length,UINT64_LENGTH);
  const uint64_t pattern_mod = pattern_length%UINT64_LENGTH;
  // Init fields
  bpm_pattern->pattern_length = pattern_length;
  bpm_pattern->pattern_num_words64 = pattern_num_words;
  bpm_pattern->pattern_mod = pattern_mod;
  // Allocate memory
  const uint64_t gpu_pattern_num_words = DIV_CEIL(pattern_length,GPU_BPM_DISTANCE_ENTRY_LENGTH);
  const uint64_t aux_vector_size = gpu_pattern_num_words*GPU_BPM_DISTANCE_ENTRY_SIZE;
  const uint64_t PEQ_size = DNA__N_RANGE*aux_vector_size;
  const uint64_t score_size = pattern_num_words*UINT64_SIZE;
  uint64_t* PEQ = mm_allocator_malloc(mm_allocator,PEQ_size);
  bpm_pattern->PEQ = PEQ;
  bpm_pattern->P = mm_allocator_malloc(mm_allocator,aux_vector_size);
  bpm_pattern->M = mm_allocator_malloc(mm_allocator,aux_vector_size);
  bpm_pattern->level_mask = mm_allocator_malloc(mm_allocator,aux_vector_size);
  bpm_pattern->score = mm_allocator_malloc(mm_allocator,score_size);
  bpm_pattern->init_score = mm_allocator_malloc(mm_allocator,score_size);
  bpm_pattern->pattern_left = mm_allocator_malloc(mm_allocator,(pattern_num_words+1)*UINT64_SIZE);
  // Init PEQ
  uint8_t enc_char;
  uint64_t i, entry=0, subentry=0;
  for (i=0;i<gpu_pattern_num_words;++i) {
    subentry = 0;
    for (enc_char=0;enc_char<DNA__N_RANGE;++enc_char) {
      *PEQ = ((uint64_t)pattern_entry[entry].bitmap[enc_char][subentry]) |
             ((uint64_t)pattern_entry[entry].bitmap[enc_char][subentry+1] << 32);
      ++PEQ;
    }
    subentry += 2;
    for (enc_char=0;enc_char<DNA__N_RANGE;++enc_char) {
      *PEQ = ((uint64_t)pattern_entry[entry].bitmap[enc_char][subentry]) |
             ((uint64_t)pattern_entry[entry].bitmap[enc_char][subentry+1] << 32);
      ++PEQ;
    }
    ++entry;
  }
  // Init auxiliary data
  uint64_t pattern_left = pattern_length;
  const uint64_t top = pattern_num_words-1;
  memset(bpm_pattern->level_mask,0,aux_vector_size);
  for (i=0;i<top;++i) {
    bpm_pattern->level_mask[i] = BMP_W64_MASK;
    bpm_pattern->init_score[i] = UINT64_LENGTH;
    bpm_pattern->pattern_left[i] = pattern_left;
    pattern_left = (pattern_left > UINT64_LENGTH) ? pattern_left-UINT64_LENGTH : 0;
  }
  for (;i<=pattern_num_words;++i) {
    bpm_pattern->pattern_left[i] = pattern_left;
    pattern_left = (pattern_left > UINT64_LENGTH) ? pattern_left-UINT64_LENGTH : 0;
  }
  if (pattern_mod>0) {
    const uint64_t mask_shift = pattern_mod-1;
    bpm_pattern->level_mask[top] = 1ull<<(mask_shift);
    bpm_pattern->init_score[top] = pattern_mod;
  } else {
    bpm_pattern->level_mask[top] = BMP_W64_MASK;
    bpm_pattern->init_score[top] = UINT64_LENGTH;
  }
}
void gpu_buffer_bpm_distance_pattern_compile(
    gpu_bpm_filter_qry_entry_t* const pattern_entry,
    const bpm_pattern_t* const bpm_pattern) {
  // Copy PEQ pattern
  const uint64_t bpm_pattern_num_words = bpm_pattern->pattern_num_words64;
  const uint64_t pattern_length = bpm_pattern->pattern_length;
  const uint64_t gpu_pattern_num_words = DIV_CEIL(pattern_length,GPU_BPM_DISTANCE_ENTRY_LENGTH)*2;
  const uint32_t* PEQ = (uint32_t*) bpm_pattern->PEQ;
  uint64_t i, entry=0, subentry=0;
  for (i=0;i<bpm_pattern_num_words;++i) {
    // Update location
    if (subentry==GPU_BPM_DISTANCE_NUM_SUB_ENTRIES) {
      subentry = 0; ++entry;
    }
    // Copy pattern
    uint8_t enc_char;
    for (enc_char=0;enc_char<DNA__N_RANGE;++enc_char) {
      pattern_entry[entry].bitmap[enc_char][subentry] = *PEQ; ++PEQ;
      pattern_entry[entry].bitmap[enc_char][subentry+1] = *PEQ; ++PEQ;
    }
    subentry += 2;
  }
  for (;i<gpu_pattern_num_words;++i) {
    // Update location
    if (subentry==GPU_BPM_DISTANCE_NUM_SUB_ENTRIES) {
      subentry = 0; ++entry;
    }
    // Copy pattern
    uint8_t enc_char;
    for (enc_char=0;enc_char<DNA__N_RANGE;++enc_char) {
      pattern_entry[entry].bitmap[enc_char][subentry] = UINT32_ONES;
      pattern_entry[entry].bitmap[enc_char][subentry+1] = UINT32_ONES;
    }
    subentry += 2;
  }
}
/*
 * Accessors
 */
void gpu_buffer_bpm_distance_add_pattern(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance,
    pattern_t* const pattern) {
  // Parameters (Buffer & Dimensions)
  const uint64_t key_length = pattern->key_length;
  pattern_tiled_t* const pattern_tiled = &pattern->pattern_tiled;
  const uint64_t buffer_entries_offset = gpu_buffer_bpm_distance->num_entries;
  void* const gpu_buffer = gpu_buffer_bpm_distance->buffer;
  const uint64_t num_tiles = pattern_tiled->num_tiles;
  const uint64_t num_entries = DIV_CEIL(key_length,GPU_BPM_DISTANCE_ENTRY_LENGTH);
  // Add query (Metadata)
  const uint32_t buffer_pattern_query_offset = gpu_buffer_bpm_distance->num_queries;
  (gpu_buffer_bpm_distance->num_queries) += num_tiles;
  gpu_buffer_bpm_distance->current_query_offset = buffer_pattern_query_offset;
  // Add query (for all tiles)
  gpu_bpm_filter_qry_info_t* buffer_pattern_query =
      gpu_bpm_filter_buffer_get_peq_info_(gpu_buffer) + buffer_pattern_query_offset;
  uint64_t i, buffer_entries_added = 0;
  for (i=0;i<num_tiles;++i,++buffer_pattern_query) {
    // Compute tile dimensions
    const uint64_t tile_length = pattern_tiled->tiles[i].tile_length;
    const uint64_t tile_entries = DIV_CEIL(tile_length,GPU_BPM_DISTANCE_ENTRY_LENGTH);
    // Set entry & size
    buffer_pattern_query->posEntry = buffer_entries_offset + buffer_entries_added;
    buffer_pattern_query->size = tile_length;
    buffer_entries_added += tile_entries;
    // Check tile length
    if (gpu_buffer_bpm_distance->query_same_length != UINT32_MAX) {
      if (gpu_buffer_bpm_distance->query_same_length == 0) {
        gpu_buffer_bpm_distance->query_same_length = tile_entries;
      } else if (gpu_buffer_bpm_distance->query_same_length != tile_entries) {
        gpu_buffer_bpm_distance->query_same_length = UINT32_MAX; // Not the same length
      }
    }
  }
  // [DTO] Compile PEQ pattern (Add pattern entries)
  bpm_pattern_t* const bpm_pattern = &pattern_tiled->bpm_pattern;
  gpu_bpm_filter_qry_entry_t* const buffer_pattern_entry =
      gpu_bpm_filter_buffer_get_peq_entries_(gpu_buffer) + buffer_entries_offset;
  gpu_buffer_bpm_distance->num_entries += num_entries;
  gpu_buffer_bpm_distance_pattern_compile(buffer_pattern_entry,bpm_pattern);
}
void gpu_buffer_bpm_distance_add_candidate(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance,
    const uint64_t tile_offset,
    const uint64_t candidate_text_position,
    const uint64_t candidate_length) {
  // Insert candidate
#ifdef GEM_PROFILE
  gpu_bpm_filter_qry_info_t* const pattern_query =
      gpu_bpm_filter_buffer_get_peq_info_(gpu_buffer_bpm_distance->buffer) +
      (gpu_buffer_bpm_distance->current_query_offset + tile_offset);
  PROF_INC_COUNTER(GP_GPU_BUFFER_BPM_DISTANCE_NUM_QUERIES);
  PROF_ADD_COUNTER(GP_GPU_BUFFER_BPM_DISTANCE_CANDIDATES_LENGTH,candidate_length);
  PROF_ADD_COUNTER(GP_GPU_BUFFER_BPM_DISTANCE_CELLS,candidate_length*pattern_query->size);
#endif
  const uint64_t candidate_offset = gpu_buffer_bpm_distance->num_candidates;
  gpu_bpm_filter_cand_info_t* const candidate =
      gpu_bpm_filter_buffer_get_candidates_(gpu_buffer_bpm_distance->buffer) + candidate_offset;
  candidate->query = gpu_buffer_bpm_distance->current_query_offset + tile_offset;
  candidate->position = candidate_text_position;
  candidate->size = candidate_length;
  ++(gpu_buffer_bpm_distance->num_candidates);
}
void gpu_buffer_bpm_distance_get_candidate(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance,
    const uint64_t candidate_offset,
    uint64_t* const candidate_text_position,
    uint32_t* const candidate_length) {
  // Get candidate
  gpu_bpm_filter_cand_info_t* const candidate =
      gpu_bpm_filter_buffer_get_candidates_(gpu_buffer_bpm_distance->buffer) + candidate_offset;
  *candidate_text_position = candidate->position;
  *candidate_length = candidate->size;
}
void gpu_buffer_bpm_distance_get_result(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance,
    const uint64_t candidate_offset,
    uint32_t* const levenshtein_distance,
    uint32_t* const levenshtein_match_pos) {
  // Get candidate results
  gpu_bpm_filter_alg_entry_t* const result =
      gpu_bpm_filter_buffer_get_alignments_(gpu_buffer_bpm_distance->buffer) + candidate_offset;
  *levenshtein_distance = result->score;
  *levenshtein_match_pos = result->column;
}
void gpu_buffer_bpm_distance_retrieve_pattern(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance,
    const uint64_t candidate_offset,
    bpm_pattern_t* const bpm_pattern,
    mm_allocator_t* const mm_allocator) {
  // Parameters
  void* const gpu_buffer = gpu_buffer_bpm_distance->buffer;
  // Get candidate
  gpu_bpm_filter_cand_info_t* const candidate = gpu_bpm_filter_buffer_get_candidates_(gpu_buffer) + candidate_offset;
  // Get Query & Entry
  gpu_bpm_filter_qry_info_t* const pattern_query =
      gpu_bpm_filter_buffer_get_peq_info_(gpu_buffer) + candidate->query;
  gpu_bpm_filter_qry_entry_t* const pattern_entry =
      gpu_bpm_filter_buffer_get_peq_entries_(gpu_buffer) + pattern_query->posEntry;
  // Decompile
  gpu_buffer_bpm_distance_pattern_decompile(pattern_entry,pattern_query->size,bpm_pattern,mm_allocator);
}
/*
 * Hints
 */
void gpu_buffer_bpm_distance_record_query_length(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance,
    const uint64_t query_length) {
  COUNTER_ADD(&gpu_buffer_bpm_distance->query_length,query_length);
}
void gpu_buffer_bpm_distance_record_candidates_per_tile(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance,
    const uint64_t candidates_per_tile) {
  PROF_ADD_COUNTER(GP_GPU_BUFFER_BPM_DISTANCE_CANDIDATES_PER_TILE,candidates_per_tile);
  COUNTER_ADD(&gpu_buffer_bpm_distance->candidates_per_tile,candidates_per_tile);
}
uint64_t gpu_buffer_bpm_distance_get_mean_query_length(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance) {
  if (COUNTER_GET_NUM_SAMPLES(&gpu_buffer_bpm_distance->query_length) >= GPU_BPM_DISTANCE_MIN_NUM_SAMPLES) {
    return (uint64_t)ceil(COUNTER_GET_MEAN(&gpu_buffer_bpm_distance->query_length));
  } else {
    return GPU_BPM_DISTANCE_AVERAGE_QUERY_LENGTH;
  }
}
uint64_t gpu_buffer_bpm_distance_get_mean_candidates_per_tile(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance) {
  if (COUNTER_GET_NUM_SAMPLES(&gpu_buffer_bpm_distance->candidates_per_tile) >= GPU_BPM_DISTANCE_MIN_NUM_SAMPLES) {
    return (uint64_t)ceil(COUNTER_GET_MEAN(&gpu_buffer_bpm_distance->candidates_per_tile));
  } else {
    return GPU_BPM_DISTANCE_CANDIDATES_PER_TILE;
  }
}
/*
 * Send/Receive
 */
void gpu_buffer_bpm_distance_send(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance) {
  PROF_START(GP_GPU_BUFFER_BPM_DISTANCE_SEND);
#ifdef GEM_PROFILE
  const uint64_t max_candidates = gpu_buffer_bpm_distance_get_max_candidates(gpu_buffer_bpm_distance);
  const uint64_t max_queries = gpu_buffer_bpm_distance_get_max_queries(gpu_buffer_bpm_distance);
  const uint64_t max_entries = gpu_buffer_bpm_distance_get_max_entries(gpu_buffer_bpm_distance);
  const uint64_t used_candidates = gpu_buffer_bpm_distance->num_candidates;
  const uint64_t used_queries = gpu_buffer_bpm_distance->num_queries;
  const uint64_t used_entries = gpu_buffer_bpm_distance->num_entries;
  PROF_ADD_COUNTER(GP_GPU_BUFFER_BPM_DISTANCE_USAGE_CANDIDATES,(100*used_candidates)/max_candidates);
  PROF_ADD_COUNTER(GP_GPU_BUFFER_BPM_DISTANCE_USAGE_QUERIES,(100*used_queries)/max_queries);
  PROF_ADD_COUNTER(GP_GPU_BUFFER_BPM_DISTANCE_USAGE_PEQ_ENTRIES,(100*used_entries)/max_entries);
  TIMER_START(&gpu_buffer_bpm_distance->timer);
#endif
  // Select computing device
  if (gpu_buffer_bpm_distance->bpm_distance_enabled) {
    if (gpu_buffer_bpm_distance->num_candidates > 0) {
      gpu_bpm_filter_send_buffer_(gpu_buffer_bpm_distance->buffer,gpu_buffer_bpm_distance->num_entries,
          gpu_buffer_bpm_distance->num_queries,gpu_buffer_bpm_distance->num_candidates,
          (gpu_buffer_bpm_distance->query_same_length==UINT32_MAX) ?
              0 : gpu_buffer_bpm_distance->query_same_length);
    }
  }
	PROF_STOP(GP_GPU_BUFFER_BPM_DISTANCE_SEND);
}
void gpu_buffer_bpm_distance_receive(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance) {
  PROF_START(GP_GPU_BUFFER_BPM_DISTANCE_RECEIVE);
  if (gpu_buffer_bpm_distance->bpm_distance_enabled) {
    if (gpu_buffer_bpm_distance->num_candidates > 0) {
      gpu_bpm_filter_receive_buffer_(gpu_buffer_bpm_distance->buffer);
    }
  }
  PROF_STOP(GP_GPU_BUFFER_BPM_DISTANCE_RECEIVE);
#ifdef GEM_PROFILE
  TIMER_STOP(&gpu_buffer_bpm_distance->timer);
  COUNTER_ADD(&PROF_GET_TIMER(GP_GPU_BUFFER_BPM_DISTANCE_DUTY_CYCLE)->time_ns,
      gpu_buffer_bpm_distance->timer.accumulated);
#endif
}
/*
 * CUDA NOT-Supported
 */
#else
/*
 * Setup
 */
gpu_buffer_bpm_distance_t* gpu_buffer_bpm_distance_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_no,
    const bool bpm_distance_enabled) { GEM_CUDA_NOT_SUPPORTED(); return NULL; }
void gpu_buffer_bpm_distance_clear(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance) { GEM_CUDA_NOT_SUPPORTED(); }
void gpu_buffer_bpm_distance_delete(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance) { GEM_CUDA_NOT_SUPPORTED(); }
/*
 * Occupancy & Limits
 */
uint64_t gpu_buffer_bpm_distance_get_entry_length() { GEM_CUDA_NOT_SUPPORTED(); return 0; }
uint64_t gpu_buffer_bpm_distance_get_max_candidates(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
uint64_t gpu_buffer_bpm_distance_get_max_queries(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
uint64_t gpu_buffer_bpm_distance_get_num_candidates(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
uint64_t gpu_buffer_bpm_distance_get_num_queries(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
void gpu_buffer_bpm_distance_compute_dimensions(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance,
    pattern_t* const pattern,
    const uint64_t num_candidates,
    uint64_t* const total_entries,
    uint64_t* const total_queries,
    uint64_t* const total_candidates) { GEM_CUDA_NOT_SUPPORTED(); }
bool gpu_buffer_bpm_distance_fits_in_buffer(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance,
    const uint64_t total_entries,
    const uint64_t total_queries,
    const uint64_t total_candidates) { GEM_CUDA_NOT_SUPPORTED(); return false; }
/*
 * Accessors
 */
void gpu_buffer_bpm_distance_add_pattern(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance,
    pattern_t* const pattern) { GEM_CUDA_NOT_SUPPORTED(); }
void gpu_buffer_bpm_distance_add_candidate(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance,
    const uint64_t tile_offset,
    const uint64_t candidate_text_position,
    const uint64_t candidate_length) { GEM_CUDA_NOT_SUPPORTED(); }
void gpu_buffer_bpm_distance_get_candidate(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance,
    const uint64_t candidate_offset,
    uint64_t* const candidate_text_position,
    uint32_t* const candidate_length) { GEM_CUDA_NOT_SUPPORTED(); }
void gpu_buffer_bpm_distance_get_result(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance,
    const uint64_t candidate_offset,
    uint32_t* const levenshtein_distance,
    uint32_t* const levenshtein_match_pos) { GEM_CUDA_NOT_SUPPORTED(); }
void gpu_buffer_bpm_distance_retrieve_pattern(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance,
    const uint64_t candidate_offset,
    bpm_pattern_t* const bpm_pattern,
    mm_allocator_t* const mm_allocator) { GEM_CUDA_NOT_SUPPORTED(); }
/*
 * Hints
 */
void gpu_buffer_bpm_distance_record_query_length(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance,
    const uint64_t query_length) { GEM_CUDA_NOT_SUPPORTED(); }
void gpu_buffer_bpm_distance_record_candidates_per_tile(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance,
    const uint64_t num_candidates) { GEM_CUDA_NOT_SUPPORTED(); }
uint64_t gpu_buffer_bpm_distance_get_mean_query_length(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
uint64_t gpu_buffer_bpm_distance_get_mean_candidates_per_tile(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
/*
 * Send/Receive
 */
void gpu_buffer_bpm_distance_send(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance) { GEM_CUDA_NOT_SUPPORTED(); }
void gpu_buffer_bpm_distance_receive(
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance) { GEM_CUDA_NOT_SUPPORTED(); }
#endif /* HAVE_CUDA */

