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

#include "gpu/gpu_buffer_bpm_align.h"
#include "gpu/gpu_buffer_bpm_pattern.h"
#include "align/align_bpm_distance.h"

/*
 * Errors
 */
#define GEM_ERROR_GPU_BPM_ALIGN_MAX_QUERIES \
  "GPU.BPM.Align. Number of queries (%"PRIu64") exceeds maximum buffer capacity (%"PRIu64" queries)"
#define GEM_ERROR_GPU_BPM_ALIGN_MAX_PATTERN_ENTRIES \
  "GPU.BPM.Align. Query pattern (%"PRIu64" entries) exceeds maximum buffer capacity (%"PRIu64" entries)"
#define GEM_ERROR_GPU_BPM_ALIGN_MAX_PATTERN_LENGTH \
  "GPU.BPM.Align. Query pattern (%"PRIu64" bases) exceeds maximum buffer capacity (%"PRIu64" bases)"
#define GEM_ERROR_GPU_BPM_ALIGN_MAX_CANDIDATES \
  "GPU.BPM.Align. Number of candidates (%"PRIu64" candidates) exceeds maximum buffer capacity (%"PRIu64" candidates)"
#define GEM_ERROR_GPU_BPM_ALIGN_MAX_CANDIDATE_LENGTH \
  "GPU.BPM.Align. Candidate length (%"PRIu64" bases) exceeds maximum buffer capacity (%"PRIu64" bases)"

/*
 * Constants :: GPU BMP-Distance
 */
#define GPU_BPM_ALIGN_NUM_SUB_ENTRIES  GPU_BPM_ALIGN_PEQ_SUBENTRIES
#define GPU_BPM_ALIGN_ENTRY_LENGTH     GPU_BPM_ALIGN_PEQ_ENTRY_LENGTH
#define GPU_BPM_ALIGN_SUBENTRY_LENGTH  GPU_BPM_ALIGN_PEQ_SUBENTRY_LENGTH
#define GPU_BPM_ALIGN_ENTRY_SIZE       (GPU_BPM_ALIGN_ENTRY_LENGTH/UINT8_LENGTH)
#define GPU_BPM_ALIGN_ALPHABET_LENGTH  GPU_BPM_ALIGN_PEQ_ALPHABET_SIZE

/*
 * Constants :: Buffer Hints
 */
#define GPU_BPM_ALIGN_MIN_NUM_SAMPLES           1
#define GPU_BPM_ALIGN_AVERAGE_QUERY_LENGTH      150
#define GPU_BPM_ALIGN_CANDIDATES_PER_TILE       20

/*
 * CUDA Support
 */
#ifdef HAVE_CUDA
/*
 * Stats accessors
 */
void gpu_buffer_bpm_align_record_query_length(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    const uint64_t query_length) {
  COUNTER_ADD(&gpu_buffer_bpm_align->query_length,query_length);
}
void gpu_buffer_bpm_align_record_candidate_length(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    const uint64_t candidate_length) {
  COUNTER_ADD(&gpu_buffer_bpm_align->candidate_length,candidate_length);
}
void gpu_buffer_bpm_align_record_candidates_per_tile(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    const uint64_t candidates_per_tile) {
  PROF_ADD_COUNTER(GP_GPU_BUFFER_BPM_ALIGN_CANDIDATE_PER_TILE,candidates_per_tile);
  COUNTER_ADD(&gpu_buffer_bpm_align->candidates_per_tile,candidates_per_tile);
}
uint64_t gpu_buffer_bpm_align_get_mean_query_length(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align) {
  if (COUNTER_GET_NUM_SAMPLES(&gpu_buffer_bpm_align->query_length) >= GPU_BPM_ALIGN_MIN_NUM_SAMPLES) {
    return (uint64_t)ceil(COUNTER_GET_MEAN(&gpu_buffer_bpm_align->query_length));
  } else {
    return GPU_BPM_ALIGN_AVERAGE_QUERY_LENGTH;
  }
}
uint64_t gpu_buffer_bpm_align_get_mean_candidate_length(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align) {
  if (COUNTER_GET_NUM_SAMPLES(&gpu_buffer_bpm_align->candidate_length) >= GPU_BPM_ALIGN_MIN_NUM_SAMPLES) {
    return (uint64_t)ceil(COUNTER_GET_MEAN(&gpu_buffer_bpm_align->candidate_length));
  } else {
    return GPU_BPM_ALIGN_CANDIDATES_PER_TILE;
  }
}
uint64_t gpu_buffer_bpm_align_get_mean_candidates_per_tile(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align) {
  if (COUNTER_GET_NUM_SAMPLES(&gpu_buffer_bpm_align->candidates_per_tile) >= GPU_BPM_ALIGN_MIN_NUM_SAMPLES) {
    return (uint64_t)ceil(COUNTER_GET_MEAN(&gpu_buffer_bpm_align->candidates_per_tile));
  } else {
    return GPU_BPM_ALIGN_CANDIDATES_PER_TILE;
  }
}
/*
 * Setup
 */
gpu_buffer_bpm_align_t* gpu_buffer_bpm_align_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_no,
    const bool bpm_align_enabled) {
  PROF_START(GP_GPU_BUFFER_BPM_ALIGN_ALLOC);
  // Allocate handler
  gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align = mm_alloc(gpu_buffer_bpm_align_t);
  // Initialize buffer state
  gpu_buffer_bpm_align->buffer = gpu_buffer_collection_get_buffer(gpu_buffer_collection,buffer_no);
  gpu_buffer_bpm_align->bpm_align_enabled = bpm_align_enabled;
  gpu_buffer_bpm_align->current_query_offset = 0;
  gpu_buffer_bpm_align->current_candidates_added = 0;
  // Buffer Queries
  gpu_buffer_bpm_align->num_queries = 0;
  gpu_buffer_bpm_align->num_query_entries = 0;
  gpu_buffer_bpm_align->query_buffer_offset = 0;
  // Buffer Candidates
  gpu_buffer_bpm_align->num_candidates = 0;
  gpu_buffer_bpm_align->candidate_buffer_offset = 0;
  // Init buffer
  const int64_t thread_id = gtid(); // Between [1,num_threads] (zero is master)
  gpu_alloc_buffer_(gpu_buffer_bpm_align->buffer,thread_id);
  gpu_bpm_align_init_buffer_(gpu_buffer_bpm_align->buffer,
      gpu_buffer_bpm_align_get_mean_query_length(gpu_buffer_bpm_align),
      gpu_buffer_bpm_align_get_mean_candidate_length(gpu_buffer_bpm_align),
      gpu_buffer_bpm_align_get_mean_candidates_per_tile(gpu_buffer_bpm_align));
  // Stats
  COUNTER_RESET(&gpu_buffer_bpm_align->query_length);
  COUNTER_RESET(&gpu_buffer_bpm_align->candidate_length);
  COUNTER_RESET(&gpu_buffer_bpm_align->candidates_per_tile);
  gpu_buffer_bpm_align->query_same_length = 0;
  TIMER_RESET(&gpu_buffer_bpm_align->timer);
  // Return
  PROF_STOP(GP_GPU_BUFFER_BPM_ALIGN_ALLOC);
  return gpu_buffer_bpm_align;
}
void gpu_buffer_bpm_align_clear(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align) {
  // Init buffer
  gpu_bpm_align_init_buffer_(gpu_buffer_bpm_align->buffer,
      gpu_buffer_bpm_align_get_mean_query_length(gpu_buffer_bpm_align),
      gpu_buffer_bpm_align_get_mean_candidate_length(gpu_buffer_bpm_align),
      gpu_buffer_bpm_align_get_mean_candidates_per_tile(gpu_buffer_bpm_align));
  // Clear buffer
  gpu_buffer_bpm_align->current_query_offset = 0;
  gpu_buffer_bpm_align->current_candidates_added = 0;
  gpu_buffer_bpm_align->num_queries = 0;
  gpu_buffer_bpm_align->num_query_entries = 0;
  gpu_buffer_bpm_align->query_buffer_offset = 0;
  gpu_buffer_bpm_align->num_candidates = 0;
  gpu_buffer_bpm_align->candidate_buffer_offset = 0;
  gpu_buffer_bpm_align->query_same_length = 0;
}
void gpu_buffer_bpm_align_delete(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align) {
  mm_free(gpu_buffer_bpm_align);
}
/*
 * Limits
 */
uint64_t gpu_buffer_bpm_align_get_entry_length(void) {
  return GPU_BPM_ALIGN_ENTRY_LENGTH;
}
uint64_t gpu_buffer_bpm_align_get_max_queries(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align) {
  return gpu_bpm_align_buffer_get_max_queries_(gpu_buffer_bpm_align->buffer);
}
uint64_t gpu_buffer_bpm_align_get_max_query_entries(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align) {
  return gpu_bpm_align_buffer_get_max_peq_entries_(gpu_buffer_bpm_align->buffer);
}
uint64_t gpu_buffer_bpm_align_get_max_query_buffer_size(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align) {
  return gpu_bpm_align_buffer_get_max_query_bases_(gpu_buffer_bpm_align->buffer);
}
uint64_t gpu_buffer_bpm_align_get_max_candidates(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align) {
  return gpu_bpm_align_buffer_get_max_candidates_(gpu_buffer_bpm_align->buffer);
}
uint64_t gpu_buffer_bpm_align_get_max_candidate_buffer_size(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align) {
  return gpu_bpm_align_buffer_get_max_candidate_size_(gpu_buffer_bpm_align->buffer);
}
uint64_t gpu_buffer_bpm_align_get_num_candidates(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align) {
  return gpu_buffer_bpm_align->num_candidates;
}
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
    uint64_t* const total_candidates,
    uint64_t* const total_candidates_length) {
  // Calculate queries' dimensions
  const uint64_t num_entries = DIV_CEIL(pattern->key_length,GPU_BPM_ALIGN_ENTRY_LENGTH);
  *total_queries += pattern->pattern_tiled.num_tiles;
  *total_query_entries += num_entries;
  *total_query_length += pattern->key_length;
  // Calculate candidates' dimensions
  *total_candidates += pattern->pattern_tiled.num_tiles*num_candidates;
  *total_candidates_length += candidates_length;
}
bool gpu_buffer_bpm_align_fits_in_buffer(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    const uint64_t total_queries,
    const uint64_t total_query_entries,
    const uint64_t total_query_length,
    const uint64_t total_candidates,
    const uint64_t total_candidates_length) {
  // Get Limits
  uint64_t max_queries = gpu_buffer_bpm_align_get_max_queries(gpu_buffer_bpm_align);
  uint64_t max_query_entries =  gpu_buffer_bpm_align_get_max_query_entries(gpu_buffer_bpm_align);
  uint64_t max_query_buffer_size = gpu_buffer_bpm_align_get_max_query_buffer_size(gpu_buffer_bpm_align);
  uint64_t max_candidates = gpu_buffer_bpm_align_get_max_candidates(gpu_buffer_bpm_align);
  uint64_t max_candidate_buffer_size = gpu_buffer_bpm_align_get_max_candidate_buffer_size(gpu_buffer_bpm_align);
  // Check available space in buffer for the pattern
  if (gpu_buffer_bpm_align->num_queries+total_queries > max_queries ||
      gpu_buffer_bpm_align->num_query_entries+total_query_entries > max_query_entries ||
      gpu_buffer_bpm_align->query_buffer_offset+total_query_length > max_query_buffer_size ||
      gpu_buffer_bpm_align->num_candidates+total_candidates > max_candidates ||
      gpu_buffer_bpm_align->candidate_buffer_offset+total_candidates_length > max_candidate_buffer_size) {
    // Check buffer occupancy
    if (gpu_buffer_bpm_align->num_queries > 0) {
      return false; // Leave it to the next fresh buffer
    }
    // Reallocate buffer
    gpu_bpm_align_init_and_realloc_buffer_(
        gpu_buffer_bpm_align->buffer,total_query_entries,total_query_length,
        total_candidates_length,total_queries,total_candidates);
    // Check reallocated buffer dimensions (error otherwise)
    max_queries = gpu_buffer_bpm_align_get_max_queries(gpu_buffer_bpm_align);
    gem_cond_fatal_error(total_queries > max_queries,
        GPU_BPM_ALIGN_MAX_QUERIES,total_queries,max_queries);
    max_query_entries =  gpu_buffer_bpm_align_get_max_query_entries(gpu_buffer_bpm_align);
    gem_cond_fatal_error(total_query_entries > max_query_entries,
        GPU_BPM_ALIGN_MAX_PATTERN_ENTRIES,total_query_entries,max_query_entries);
    max_query_buffer_size = gpu_buffer_bpm_align_get_max_query_buffer_size(gpu_buffer_bpm_align);
    gem_cond_fatal_error(total_query_length > max_query_buffer_size,
        GPU_BPM_ALIGN_MAX_PATTERN_LENGTH,total_query_length,max_query_buffer_size);
    max_candidates = gpu_buffer_bpm_align_get_max_candidates(gpu_buffer_bpm_align);
    gem_cond_fatal_error(total_candidates > max_candidates,
        GPU_BPM_ALIGN_MAX_CANDIDATES,total_candidates,max_candidates);
    max_candidate_buffer_size = gpu_buffer_bpm_align_get_max_candidate_buffer_size(gpu_buffer_bpm_align);
    gem_cond_fatal_error(total_candidates_length > max_candidate_buffer_size,
        GPU_BPM_ALIGN_MAX_CANDIDATE_LENGTH,total_candidates_length,max_candidate_buffer_size);
    // Return OK (after reallocation)
    return true;
  }
  // Return OK (fits in buffer)
  return true;
}
/*
 * Pattern accessor
 */
void gpu_buffer_bpm_align_add_query_info(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    pattern_t* const pattern) {
  // Parameters GPU-Buffer
  void* const gpu_buffer = gpu_buffer_bpm_align->buffer;
  const uint64_t num_query_entries = gpu_buffer_bpm_align->num_query_entries;
  // Parameters Pattern
  pattern_tiled_t* const pattern_tiled = &pattern->pattern_tiled;
  const uint64_t num_tiles = pattern_tiled->num_tiles;
  // Add Query (Metadata)
  const uint32_t num_queries = gpu_buffer_bpm_align->num_queries;
  (gpu_buffer_bpm_align->num_queries) += num_tiles;
  gpu_buffer_bpm_align->current_query_offset = num_queries;
  // Add Query Tiles (Metadata)
  gpu_bpm_align_qry_info_t* gpu_queries_info =
      gpu_bpm_align_buffer_get_queries_info_(gpu_buffer) + num_queries;
  uint64_t i, buffer_entries_added = 0;
  for (i=0;i<num_tiles;++i,++gpu_queries_info) {
    // Compute tile dimensions
    const uint64_t tile_length = pattern_tiled->tiles[i].tile_length;
    const uint64_t tile_entries = DIV_CEIL(tile_length,GPU_BPM_ALIGN_ENTRY_LENGTH);
    // Set entry & size
    gpu_queries_info->posEntryBase = num_query_entries;
    gpu_queries_info->posEntryPEQ = num_query_entries + buffer_entries_added;
    gpu_queries_info->size = tile_length;
    buffer_entries_added += tile_entries;
    // Check tile length
    if (gpu_buffer_bpm_align->query_same_length != UINT32_MAX) {
      if (gpu_buffer_bpm_align->query_same_length == 0) {
        gpu_buffer_bpm_align->query_same_length = tile_entries;
      } else if (gpu_buffer_bpm_align->query_same_length != tile_entries) {
        gpu_buffer_bpm_align->query_same_length = UINT32_MAX; // Not the same length
      }
    }
    // Stats
    gpu_buffer_bpm_align_record_query_length(gpu_buffer_bpm_align,tile_length);
  }
}
void gpu_buffer_bpm_align_add_query_entries(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    pattern_t* const pattern) {
  // Parameters
  const uint64_t key_length = pattern->key_length;
  pattern_tiled_t* const pattern_tiled = &pattern->pattern_tiled;
  const uint64_t num_entries = DIV_CEIL(key_length,GPU_BPM_ALIGN_ENTRY_LENGTH);
  const uint64_t num_query_entries = gpu_buffer_bpm_align->num_query_entries;
  bpm_pattern_t* const bpm_pattern = &pattern_tiled->bpm_pattern;
  // Compile pattern entries (PEQ vector)
  gpu_bpm_align_peq_entry_t* const gpu_queries_entry =
      gpu_bpm_align_buffer_get_peq_entries_(gpu_buffer_bpm_align->buffer) + num_query_entries;
  gpu_buffer_bpm_pattern_compile((gpu_bpm_peq_entry_t*)gpu_queries_entry,bpm_pattern);
  gpu_buffer_bpm_align->num_query_entries += num_entries;
}
void gpu_buffer_bpm_align_add_query_text(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    pattern_t* const pattern) {
  // Parameters
  const uint64_t key_length = pattern->key_length;
  const uint8_t* const key = pattern->key;
  const uint32_t query_buffer_offset = gpu_buffer_bpm_align->query_buffer_offset;
  // Add query text
  gpu_bpm_align_qry_entry_t* const gpu_query_buffer =
      gpu_bpm_align_buffer_get_queries_(gpu_buffer_bpm_align->buffer) + query_buffer_offset;
  memcpy(gpu_query_buffer,key,key_length);
  gpu_buffer_bpm_align->query_buffer_offset += key_length;
}
void gpu_buffer_bpm_align_add_pattern(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    pattern_t* const pattern) {
  // Stats
  if (gpu_buffer_bpm_align->current_candidates_added != 0) {
    gpu_buffer_bpm_align_record_candidates_per_tile(
        gpu_buffer_bpm_align,gpu_buffer_bpm_align->current_candidates_added);
    gpu_buffer_bpm_align->current_candidates_added = 0;
  }
  // Add Query (Metadata)
  gpu_buffer_bpm_align_add_query_info(gpu_buffer_bpm_align,pattern);
  // Add Query (Entries / PEQ pattern)
  gpu_buffer_bpm_align_add_query_entries(gpu_buffer_bpm_align,pattern);
  // Add Query Text
  gpu_buffer_bpm_align_add_query_text(gpu_buffer_bpm_align,pattern);
}
/*
 * Candidate accessor
 */
void gpu_buffer_bpm_align_add_candidate(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    const uint64_t tile_offset,
    uint8_t* const candidate_buffer,
    const uint64_t candidate_length,
    const bool left_gap_align) {
  // Insert candidate
#ifdef GEM_PROFILE
  gpu_bpm_align_qry_info_t* const gpu_queries_info =
      gpu_bpm_align_buffer_get_queries_info_(gpu_buffer_bpm_align->buffer) +
      (gpu_buffer_bpm_align->current_query_offset + tile_offset);
  PROF_INC_COUNTER(GP_GPU_BUFFER_BPM_ALIGN_NUM_QUERIES);
  PROF_ADD_COUNTER(GP_GPU_BUFFER_BPM_ALIGN_CANDIDATE_LENGTH,candidate_length);
  PROF_ADD_COUNTER(GP_GPU_BUFFER_BPM_ALIGN_CELLS,candidate_length*gpu_queries_info->size);
#endif
  // Add candidate info
  const uint64_t num_candidates = gpu_buffer_bpm_align->num_candidates;
  gpu_bpm_align_cand_info_t* const gpu_candidate_info =
      gpu_bpm_align_buffer_get_candidates_info_(gpu_buffer_bpm_align->buffer) + num_candidates;
  gpu_candidate_info->idQuery = gpu_buffer_bpm_align->current_query_offset + tile_offset;
  gpu_candidate_info->posEntryBase = gpu_buffer_bpm_align->current_query_offset;
  gpu_candidate_info->size = candidate_length;
  gpu_candidate_info->leftGapAlign = left_gap_align;
  ++(gpu_buffer_bpm_align->num_candidates);
  ++(gpu_buffer_bpm_align->current_candidates_added);
  // Add candidate text
  const uint64_t candidate_buffer_offset = gpu_buffer_bpm_align->candidate_buffer_offset;
  gpu_bpm_align_cand_entry_t* const gpu_candidate_buffer =
      gpu_bpm_align_buffer_get_candidates_(gpu_buffer_bpm_align->buffer) + candidate_buffer_offset;
  memcpy(gpu_candidate_buffer,candidate_buffer,candidate_length);
  gpu_buffer_bpm_align->candidate_buffer_offset += candidate_length;
  // Stats
  gpu_buffer_bpm_align_record_candidate_length(gpu_buffer_bpm_align,candidate_length);
}
/*
 * Alignment accessor (CIGAR)
 */
void gpu_buffer_bpm_align_retrieve_alignment(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    const uint64_t candidate_offset,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector) {
  // Retrieve alignment
  gpu_bpm_align_cigar_info_t* const gpu_alignment =
      gpu_bpm_align_buffer_get_cigars_info_(gpu_buffer_bpm_align->buffer) + candidate_offset;
  match_alignment->effective_length = gpu_alignment->matchEffLenght;
  match_alignment->match_position = 0;
  /* TODO  TODO  TODO
   * match_alignment->match_position =
   *   filtering_region->text_begin_position +
   *   tile_offset +
   *   gpu_alignment->initCood;
   * TODO  TODO  TODO
   */
  match_alignment->match_text_offset = gpu_alignment->initCood.x;
  // Allocate CIGAR string memory (worst case)
  match_alignment->cigar_offset = vector_get_used(cigar_vector); // Set CIGAR offset
  match_alignment->cigar_length = gpu_alignment->cigarLenght;
  vector_reserve_additional(cigar_vector,gpu_alignment->cigarLenght); // Reserve
  cigar_element_t* const cigar_buffer = vector_get_free_elm(cigar_vector,cigar_element_t); // Sentinel
  // Convert CIGAR
  gpu_bpm_align_cigar_entry_t* const gpu_alignment_cigar =
      gpu_bpm_align_buffer_get_cigars_(gpu_buffer_bpm_align->buffer) + gpu_alignment->cigarStartPos;
  uint64_t i;
  for (i=0;i<gpu_alignment->cigarLenght;++i) {
    switch (gpu_alignment_cigar[i].event) {
      case GPU_CIGAR_MATCH:
        cigar_buffer[i].type = cigar_match;
        cigar_buffer[i].length = gpu_alignment_cigar[i].occurrences;
        break;
      case GPU_CIGAR_MISSMATCH:
        cigar_buffer[i].type = cigar_mismatch;
        break;
      case GPU_CIGAR_INSERTION:
        cigar_buffer[i].type = cigar_ins;
        cigar_buffer[i].length = gpu_alignment_cigar[i].occurrences;
        break;
      case GPU_CIGAR_DELETION:
        cigar_buffer[i].type = cigar_del;
        cigar_buffer[i].length = gpu_alignment_cigar[i].occurrences;
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
  vector_add_used(cigar_vector,gpu_alignment->cigarLenght);
}
/*
 * Send/Receive
 */
void gpu_buffer_bpm_align_send(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align) {
  PROF_START(GP_GPU_BUFFER_BPM_ALIGN_SEND);
#ifdef GEM_PROFILE
  const uint64_t max_queries = gpu_buffer_bpm_align_get_max_queries(gpu_buffer_bpm_align);
  const uint64_t max_query_entries = gpu_buffer_bpm_align_get_max_query_entries(gpu_buffer_bpm_align);
  const uint64_t max_query_buffer_size = gpu_buffer_bpm_align_get_max_query_buffer_size(gpu_buffer_bpm_align);
  const uint64_t max_candidates = gpu_buffer_bpm_align_get_max_candidates(gpu_buffer_bpm_align);
  const uint64_t max_candidate_buffer_size = gpu_buffer_bpm_align_get_max_candidate_buffer_size(gpu_buffer_bpm_align);
  const uint64_t used_queries = gpu_buffer_bpm_align->num_queries;
  const uint64_t used_query_entries = gpu_buffer_bpm_align->num_query_entries;
  const uint64_t used_query_buffer = gpu_buffer_bpm_align->query_buffer_offset;
  const uint64_t used_candidates = gpu_buffer_bpm_align->num_candidates;
  const uint64_t used_candidates_buffer = gpu_buffer_bpm_align->candidate_buffer_offset;
  PROF_ADD_COUNTER(GP_GPU_BUFFER_BPM_ALIGN_USAGE_QUERIES,(100*used_queries)/max_queries);
  PROF_ADD_COUNTER(GP_GPU_BUFFER_BPM_ALIGN_USAGE_QUERY_ENTRIES,(100*used_query_entries)/max_query_entries);
  PROF_ADD_COUNTER(GP_GPU_BUFFER_BPM_ALIGN_USAGE_QUERY_BUFFER,(100*used_query_buffer)/max_query_buffer_size);
  PROF_ADD_COUNTER(GP_GPU_BUFFER_BPM_ALIGN_USAGE_CANDIDATES,(100*used_candidates)/max_candidates);
  PROF_ADD_COUNTER(GP_GPU_BUFFER_BPM_ALIGN_USAGE_CANDIDATE_BUFFER,(100*used_candidates_buffer)/max_candidate_buffer_size);
  TIMER_START(&gpu_buffer_bpm_align->timer);
#endif
  // Select computing device
  if (gpu_buffer_bpm_align->bpm_align_enabled) {
    if (gpu_buffer_bpm_align->num_candidates > 0) {
      gpu_bpm_align_send_buffer_(
          gpu_buffer_bpm_align->buffer,
          gpu_buffer_bpm_align->num_query_entries,
          gpu_buffer_bpm_align->query_buffer_offset,
          gpu_buffer_bpm_align->candidate_buffer_offset,
          gpu_buffer_bpm_align->num_queries,
          gpu_buffer_bpm_align->num_candidates,
          (gpu_buffer_bpm_align->query_same_length==UINT32_MAX) ?
              0 : gpu_buffer_bpm_align->query_same_length);
    }
  }
	PROF_STOP(GP_GPU_BUFFER_BPM_ALIGN_SEND);
}
void gpu_buffer_bpm_align_receive(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align) {
  PROF_START(GP_GPU_BUFFER_BPM_ALIGN_RECEIVE);
  if (gpu_buffer_bpm_align->bpm_align_enabled) {
    if (gpu_buffer_bpm_align->num_candidates > 0) {
      gpu_bpm_align_receive_buffer_(gpu_buffer_bpm_align->buffer);
    }
  }
  PROF_STOP(GP_GPU_BUFFER_BPM_ALIGN_RECEIVE);
#ifdef GEM_PROFILE
  TIMER_STOP(&gpu_buffer_bpm_align->timer);
  COUNTER_ADD(&PROF_GET_TIMER(GP_GPU_BUFFER_BPM_ALIGN_DUTY_CYCLE)->time_ns,
      gpu_buffer_bpm_align->timer.accumulated);
#endif
}
/*
 * CUDA NOT-Supported
 */
#else
/*
 * Setup
 */
gpu_buffer_bpm_align_t* gpu_buffer_bpm_align_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_no,
    const bool bpm_align_enabled) { return NULL; } //  GEM_CUDA_NOT_SUPPORTED();
void gpu_buffer_bpm_align_clear(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align) {  }
void gpu_buffer_bpm_align_delete(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align) {  }
/*
 * Limits
 */
uint64_t gpu_buffer_bpm_align_get_entry_length(void) { return 0; }
uint64_t gpu_buffer_bpm_align_get_max_queries(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align) { return 0; }
uint64_t gpu_buffer_bpm_align_get_max_query_entries(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align) { return 0; }
uint64_t gpu_buffer_bpm_align_get_max_query_buffer_size(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align) { return 0; }
uint64_t gpu_buffer_bpm_align_get_max_candidates(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align) { return 0; }
uint64_t gpu_buffer_bpm_align_get_max_candidate_buffer_size(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align) { return 0; }
uint64_t gpu_buffer_bpm_align_get_num_candidates(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align) { return 0; }
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
    uint64_t* const total_candidates,
    uint64_t* const total_candidates_length);
bool gpu_buffer_bpm_align_fits_in_buffer(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    const uint64_t total_queries,
    const uint64_t total_query_entries,
    const uint64_t total_query_length,
    const uint64_t total_candidates,
    const uint64_t total_candidates_length);
/*
 * Pattern accessor
 */
void gpu_buffer_bpm_align_add_pattern(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    pattern_t* const pattern) {  }
/*
 * Candidate accessor
 */
void gpu_buffer_bpm_align_add_candidate(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    const uint64_t tile_offset,
    uint8_t* const candidate_buffer,
    const uint64_t candidate_length,
    const bool left_gap_align) {  }
/*
 * Alignment accessor (CIGAR)
 */
void gpu_buffer_bpm_align_retrieve_alignment(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    const uint64_t candidate_offset,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector) {  }
/*
 * Send/Receive
 */
void gpu_buffer_bpm_align_send(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align) {  }
void gpu_buffer_bpm_align_receive(
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align) {  }
#endif /* HAVE_CUDA */

