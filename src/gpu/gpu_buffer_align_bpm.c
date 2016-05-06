/*
 * PROJECT: GEMMapper
 * FILE: gpu_buffer_align_bpm.c
 * DATE: 04/09/2014
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 *            Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "gpu/gpu_buffer_align_bpm.h"
#include "resources/gpu_modules/gpu_interface.h"
#include "align/align_bpm_distance.h"

/*
 * Errors
 */
#define GEM_ERROR_GPU_ALIGN_BPM_MAX_PATTERN_LENGTH "GPU.BPM.Align. Query pattern (%"PRIu64" entries) exceeds maximum buffer capacity (%"PRIu64" entries)"
#define GEM_ERROR_GPU_ALIGN_BPM_MAX_CANDIDATES "GPU.BPM.Align. Number of candidates (%"PRIu64") exceeds maximum buffer capacity (%"PRIu64" candidates)"
#define GEM_ERROR_GPU_ALIGN_BPM_MAX_QUERIES "GPU.BPM.Align. Number of queries (%"PRIu64") exceeds maximum buffer capacity (%"PRIu64" queries)"

/*
 * Constants :: GPU Align-BMP
 */
#define GPU_ALIGN_BPM_NUM_SUB_ENTRIES  GPU_BPM_PEQ_SUBENTRIES
#define GPU_ALIGN_BPM_ENTRY_LENGTH     GPU_BPM_PEQ_ENTRY_LENGTH
#define GPU_ALIGN_BPM_SUBENTRY_LENGTH  GPU_BPM_PEQ_SUBENTRY_LENGTH
#define GPU_ALIGN_BPM_ENTRY_SIZE       (GPU_ALIGN_BPM_ENTRY_LENGTH/UINT8_LENGTH)
#define GPU_ALIGN_BPM_ALPHABET_LENGTH  GPU_BPM_PEQ_ALPHABET_SIZE

/*
 * Constants :: Buffer Hints
 */
#define GPU_ALIGN_BPM_MIN_NUM_SAMPLES           1
#define GPU_ALIGN_BPM_AVERAGE_QUERY_LENGTH      150
#define GPU_ALIGN_BPM_CANDIDATES_PER_TILE       20

/*
 * CUDA Support
 */
#ifdef HAVE_CUDA
/*
 * Setup
 */
gpu_buffer_align_bpm_t* gpu_buffer_align_bpm_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_no,
    archive_text_t* const archive_text,
    text_collection_t* const text_collection,
    mm_stack_t* const mm_stack) {
  PROF_START(GP_GPU_BUFFER_ALIGN_BPM_ALLOC);
  // Alloc
  gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm = mm_alloc(gpu_buffer_align_bpm_t);
  // Dimensions Hints
  COUNTER_RESET(&gpu_buffer_align_bpm->query_length);
  COUNTER_RESET(&gpu_buffer_align_bpm->candidates_per_tile);
  gpu_buffer_align_bpm->query_same_length = 0;
  // Buffer state
  gpu_buffer_align_bpm->buffer = gpu_buffer_collection_get_buffer(gpu_buffer_collection,buffer_no);
  gpu_buffer_align_bpm->num_entries = 0;
  gpu_buffer_align_bpm->num_queries = 0;
  gpu_buffer_align_bpm->num_candidates = 0;
  gpu_buffer_align_bpm->current_query_offset = 0;
  // CPU Computation
  gpu_buffer_align_bpm->compute_cpu = false;
  gpu_buffer_align_bpm->archive_text = archive_text;
  gpu_buffer_align_bpm->text_collection = text_collection;
  gpu_buffer_align_bpm->mm_stack = mm_stack;
  TIMER_RESET(&gpu_buffer_align_bpm->timer);
  // Init buffer
  gpu_alloc_buffer_(gpu_buffer_align_bpm->buffer);
  gpu_bpm_init_buffer_(gpu_buffer_align_bpm->buffer,
      gpu_buffer_align_bpm_get_mean_query_length(gpu_buffer_align_bpm),
      gpu_buffer_align_bpm_get_mean_candidates_per_tile(gpu_buffer_align_bpm));
  PROF_STOP(GP_GPU_BUFFER_ALIGN_BPM_ALLOC);
  // Return
  return gpu_buffer_align_bpm;
}
void gpu_buffer_align_bpm_clear(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  // Init buffer
  gpu_bpm_init_buffer_(gpu_buffer_align_bpm->buffer,
      gpu_buffer_align_bpm_get_mean_query_length(gpu_buffer_align_bpm),
      gpu_buffer_align_bpm_get_mean_candidates_per_tile(gpu_buffer_align_bpm));
  // Dimensions Hints
  gpu_buffer_align_bpm->query_same_length = 0;
  gpu_buffer_align_bpm->num_entries = 0;
  gpu_buffer_align_bpm->num_queries = 0;
  gpu_buffer_align_bpm->num_candidates = 0;
  gpu_buffer_align_bpm->current_query_offset = 0;
}
void gpu_buffer_align_bpm_delete(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  mm_free(gpu_buffer_align_bpm);
}
/*
 * Computing Device
 */
void gpu_buffer_align_bpm_set_device_cpu(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  gpu_buffer_align_bpm->compute_cpu = true;
}
void gpu_buffer_align_bpm_set_device_gpu(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  gpu_buffer_align_bpm->compute_cpu = false;
}
/*
 * Occupancy & Limits
 */
uint64_t gpu_buffer_align_bpm_get_entry_length() {
  return GPU_ALIGN_BPM_ENTRY_LENGTH;
}
uint64_t gpu_buffer_align_bpm_get_max_candidates(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  return gpu_bpm_buffer_get_max_candidates_(gpu_buffer_align_bpm->buffer);
}
uint64_t gpu_buffer_align_bpm_get_max_queries(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  return gpu_bpm_buffer_get_max_queries_(gpu_buffer_align_bpm->buffer);
}
uint64_t gpu_buffer_align_bpm_get_max_entries(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  return gpu_bpm_buffer_get_max_peq_entries_(gpu_buffer_align_bpm->buffer);
}
uint64_t gpu_buffer_align_bpm_get_num_candidates(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  return gpu_buffer_align_bpm->num_candidates;
}
uint64_t gpu_buffer_align_bpm_get_num_queries(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  return gpu_buffer_align_bpm->num_queries;
}
void gpu_buffer_align_bpm_compute_dimensions(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    bpm_pattern_t* const bpm_pattern,
    bpm_pattern_t* const bpm_pattern_tiles,
    const uint64_t num_candidates,
    uint64_t* const total_entries,
    uint64_t* const total_queries,
    uint64_t* const total_candidates) {
  // Calculate dimensions
  const uint64_t num_entries = DIV_CEIL(bpm_pattern->pattern_length,GPU_ALIGN_BPM_ENTRY_LENGTH);
  const uint64_t num_tiles = bpm_pattern_tiles->num_pattern_tiles;
  *total_entries += num_entries;
  *total_queries += num_tiles;
  *total_candidates += num_tiles*num_candidates;
}
bool gpu_buffer_align_bpm_fits_in_buffer(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    const uint64_t total_entries,
    const uint64_t total_queries,
    const uint64_t total_candidates) {
  // Get Limits
  uint64_t max_PEQ_entries =  gpu_buffer_align_bpm_get_max_entries(gpu_buffer_align_bpm);
  uint64_t max_queries = gpu_buffer_align_bpm_get_max_queries(gpu_buffer_align_bpm);
  uint64_t max_candidates = gpu_buffer_align_bpm_get_max_candidates(gpu_buffer_align_bpm);
  // Check available space in buffer for the pattern
  if (gpu_buffer_align_bpm->num_queries+total_queries > max_queries ||
      gpu_buffer_align_bpm->num_entries+total_entries > max_PEQ_entries ||
      gpu_buffer_align_bpm->num_candidates+total_candidates > max_candidates) {
    // Check buffer occupancy
    if (gpu_buffer_align_bpm->num_queries > 0) {
      return false; // Leave it to the next fresh buffer
    }
    // Reallocate buffer
    gpu_bpm_init_and_realloc_buffer_(gpu_buffer_align_bpm->buffer,total_entries,total_candidates,total_queries);
    // Check reallocated buffer dimensions (error otherwise)
    max_queries = gpu_buffer_align_bpm_get_max_queries(gpu_buffer_align_bpm);
    gem_cond_fatal_error(total_queries > max_queries,GPU_ALIGN_BPM_MAX_QUERIES,total_queries,max_queries);
    max_PEQ_entries = gpu_buffer_align_bpm_get_max_entries(gpu_buffer_align_bpm);
    gem_cond_fatal_error(total_entries > max_PEQ_entries,GPU_ALIGN_BPM_MAX_PATTERN_LENGTH,total_entries,max_PEQ_entries);
    max_candidates = gpu_buffer_align_bpm_get_max_candidates(gpu_buffer_align_bpm);
    gem_cond_fatal_error(total_candidates > max_candidates,GPU_ALIGN_BPM_MAX_CANDIDATES,total_candidates,max_candidates);
    // Return OK (after reallocation)
    return true;
  }
  // Return OK (fits in buffer)
  return true;
}
/*
 * Pattern Compile/Decompile
 */
void gpu_buffer_align_bpm_pattern_decompile(
    gpu_bpm_qry_entry_t* const pattern_entry,
    const uint32_t pattern_length,
    bpm_pattern_t* const bpm_pattern,
    mm_stack_t* const mm_stack) {
  // Calculate dimensions
  const uint64_t pattern_num_words = DIV_CEIL(pattern_length,UINT64_LENGTH);
  const uint64_t pattern_mod = pattern_length%UINT64_LENGTH;
  // Init fields
  bpm_pattern->pattern_length = pattern_length;
  bpm_pattern->pattern_num_words64 = pattern_num_words;
  bpm_pattern->pattern_mod = pattern_mod;
  // Allocate memory
  const uint64_t gpu_pattern_num_words = DIV_CEIL(pattern_length,GPU_ALIGN_BPM_ENTRY_LENGTH);
  const uint64_t aux_vector_size = gpu_pattern_num_words*GPU_ALIGN_BPM_ENTRY_SIZE;
  const uint64_t PEQ_size = DNA__N_RANGE*aux_vector_size;
  const uint64_t score_size = pattern_num_words*UINT64_SIZE;
  uint64_t* PEQ = mm_stack_malloc(mm_stack,PEQ_size);
  bpm_pattern->PEQ = PEQ;
  bpm_pattern->P = mm_stack_malloc(mm_stack,aux_vector_size);
  bpm_pattern->M = mm_stack_malloc(mm_stack,aux_vector_size);
  bpm_pattern->level_mask = mm_stack_malloc(mm_stack,aux_vector_size);
  bpm_pattern->score = mm_stack_malloc(mm_stack,score_size);
  bpm_pattern->init_score = mm_stack_malloc(mm_stack,score_size);
  bpm_pattern->pattern_left = mm_stack_malloc(mm_stack,(pattern_num_words+1)*UINT64_SIZE);
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
void gpu_buffer_align_bpm_pattern_compile(
    gpu_bpm_qry_entry_t* const pattern_entry,
    const bpm_pattern_t* const bpm_pattern,
    const bpm_pattern_t* const bpm_pattern_tiles) {
  // Copy PEQ pattern
  const uint64_t bpm_pattern_num_words = bpm_pattern->pattern_num_words64;
  const uint64_t gpu_pattern_length = bpm_pattern_tiles->num_pattern_tiles*bpm_pattern_tiles->tile_length;
  const uint64_t gpu_pattern_num_words = gpu_pattern_length/UINT64_LENGTH;
  const uint32_t* PEQ = (uint32_t*) bpm_pattern->PEQ;
  uint64_t i, entry=0, subentry=0;
  for (i=0;i<bpm_pattern_num_words;++i) {
    // Update location
    if (subentry==GPU_ALIGN_BPM_NUM_SUB_ENTRIES) {
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
    if (subentry==GPU_ALIGN_BPM_NUM_SUB_ENTRIES) {
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
void gpu_buffer_align_bpm_add_pattern(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    bpm_pattern_t* const bpm_pattern,
    bpm_pattern_t* const bpm_pattern_tiles) {
  // Parameters (Buffer & Dimensions)
  const uint64_t pattern_length = bpm_pattern->pattern_length;
  const uint64_t buffer_num_entries = gpu_buffer_align_bpm->num_entries;
  void* const gpu_buffer = gpu_buffer_align_bpm->buffer;
  const uint64_t num_entries = DIV_CEIL(bpm_pattern->pattern_length,GPU_ALIGN_BPM_ENTRY_LENGTH);
  const uint64_t num_entries_per_tile = DIV_CEIL(bpm_pattern_tiles->tile_length,GPU_ALIGN_BPM_ENTRY_LENGTH);
  const uint64_t num_tiles = bpm_pattern_tiles->num_pattern_tiles;
  // Add query (Metadata)
  const uint32_t buffer_pattern_query_offset = gpu_buffer_align_bpm->num_queries;
  (gpu_buffer_align_bpm->num_queries) += num_tiles;
  gpu_buffer_align_bpm->current_query_offset = buffer_pattern_query_offset;
  // Add query (for all tiles)
  gpu_bpm_qry_info_t* buffer_pattern_query = gpu_bpm_buffer_get_peq_info_(gpu_buffer) + buffer_pattern_query_offset;
  const uint64_t tile_entry_length = num_entries_per_tile*GPU_ALIGN_BPM_ENTRY_LENGTH;
  uint64_t i, tile_length, remaining_pattern_length = pattern_length;
  for (i=0;i<num_tiles;++i,++buffer_pattern_query) {
    // Set entry & size
    buffer_pattern_query->posEntry = buffer_num_entries + i*num_entries_per_tile;
    tile_length = MIN(remaining_pattern_length,tile_entry_length);
    buffer_pattern_query->size = tile_length;
    remaining_pattern_length -= tile_length;
    // Check tile length
    if (gpu_buffer_align_bpm->query_same_length != UINT32_MAX) {
      const uint32_t tile_num_entries = DIV_CEIL(tile_length,GPU_ALIGN_BPM_ENTRY_LENGTH);
      if (gpu_buffer_align_bpm->query_same_length == 0) {
        gpu_buffer_align_bpm->query_same_length = tile_num_entries;
      } else {
        if (gpu_buffer_align_bpm->query_same_length != tile_num_entries) {
          gpu_buffer_align_bpm->query_same_length = UINT32_MAX; // Not the same length
        }
      }
    }
  }
  // [DTO] Compile PEQ pattern (Add pattern entries)
  gpu_bpm_qry_entry_t* const buffer_pattern_entry = gpu_bpm_buffer_get_peq_entries_(gpu_buffer) + buffer_num_entries;
  gpu_buffer_align_bpm->num_entries += num_entries;
  gpu_buffer_align_bpm_pattern_compile(buffer_pattern_entry,bpm_pattern,bpm_pattern_tiles);
}
void gpu_buffer_align_bpm_add_candidate(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    const uint64_t tile_offset,
    const uint64_t candidate_text_position,
    const uint64_t candidate_length) {
  // Insert candidate
#ifdef GEM_PROFILE
  gpu_bpm_qry_info_t* const pattern_query =
      gpu_bpm_buffer_get_peq_info_(gpu_buffer_align_bpm->buffer) +
      (gpu_buffer_align_bpm->current_query_offset + tile_offset);
  PROF_INC_COUNTER(GP_GPU_BUFFER_ALIGN_BPM_NUM_QUERIES);
  PROF_ADD_COUNTER(GP_GPU_BUFFER_ALIGN_BPM_CANDIDATE_LENGTH,candidate_length);
  PROF_ADD_COUNTER(GP_GPU_BUFFER_ALIGN_BPM_CELLS,candidate_length*pattern_query->size);
#endif
  const uint64_t candidate_offset = gpu_buffer_align_bpm->num_candidates;
  gpu_bpm_cand_info_t* const candidate =
      gpu_bpm_buffer_get_candidates_(gpu_buffer_align_bpm->buffer) + candidate_offset;
  candidate->query = gpu_buffer_align_bpm->current_query_offset + tile_offset;
  candidate->position = candidate_text_position;
  candidate->size = candidate_length;
  ++(gpu_buffer_align_bpm->num_candidates);
}
void gpu_buffer_align_bpm_get_candidate(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    const uint64_t candidate_offset,
    uint64_t* const candidate_text_position,
    uint32_t* const candidate_length) {
  // Get candidate
  gpu_bpm_cand_info_t* const candidate = gpu_bpm_buffer_get_candidates_(gpu_buffer_align_bpm->buffer) + candidate_offset;
  *candidate_text_position = candidate->position;
  *candidate_length = candidate->size;
}
void gpu_buffer_align_bpm_get_result(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    const uint64_t candidate_offset,
    uint32_t* const levenshtein_distance,
    uint32_t* const levenshtein_match_pos) {
  // Get candidate results
  gpu_bpm_alg_entry_t* const result = gpu_bpm_buffer_get_alignments_(gpu_buffer_align_bpm->buffer) + candidate_offset;
  *levenshtein_distance = result->score;
  *levenshtein_match_pos = result->column;
}
void gpu_buffer_align_bpm_retrieve_pattern(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    const uint64_t candidate_offset,
    bpm_pattern_t* const bpm_pattern,
    mm_stack_t* const mm_stack) {
  // Parameters
  void* const gpu_buffer = gpu_buffer_align_bpm->buffer;
  // Get candidate
  gpu_bpm_cand_info_t* const candidate = gpu_bpm_buffer_get_candidates_(gpu_buffer) + candidate_offset;
  // Get Query & Entry
  gpu_bpm_qry_info_t* const pattern_query = gpu_bpm_buffer_get_peq_info_(gpu_buffer) + candidate->query;
  gpu_bpm_qry_entry_t* const pattern_entry = gpu_bpm_buffer_get_peq_entries_(gpu_buffer) + pattern_query->posEntry;
  // Decompile
  gpu_buffer_align_bpm_pattern_decompile(pattern_entry,pattern_query->size,bpm_pattern,mm_stack);
}
/*
 * Hints
 */
void gpu_buffer_align_bpm_record_query_length(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    const uint64_t query_length) {
  COUNTER_ADD(&gpu_buffer_align_bpm->query_length,query_length);
}
void gpu_buffer_align_bpm_record_candidates_per_tile(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    const uint64_t candidates_per_tile) {
  PROF_ADD_COUNTER(GP_GPU_BUFFER_ALIGN_BPM_CANDIDATE_PER_TILE,candidates_per_tile);
  COUNTER_ADD(&gpu_buffer_align_bpm->candidates_per_tile,candidates_per_tile);
}
uint64_t gpu_buffer_align_bpm_get_mean_query_length(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  if (COUNTER_GET_NUM_SAMPLES(&gpu_buffer_align_bpm->query_length) >= GPU_ALIGN_BPM_MIN_NUM_SAMPLES) {
    return (uint64_t)ceil(COUNTER_GET_MEAN(&gpu_buffer_align_bpm->query_length));
  } else {
    return GPU_ALIGN_BPM_AVERAGE_QUERY_LENGTH;
  }
}
uint64_t gpu_buffer_align_bpm_get_mean_candidates_per_tile(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  if (COUNTER_GET_NUM_SAMPLES(&gpu_buffer_align_bpm->candidates_per_tile) >= GPU_ALIGN_BPM_MIN_NUM_SAMPLES) {
    return (uint64_t)ceil(COUNTER_GET_MEAN(&gpu_buffer_align_bpm->candidates_per_tile));
  } else {
    return GPU_ALIGN_BPM_CANDIDATES_PER_TILE;
  }
}
/*
 * CPU emulated
 */
void gpu_buffer_align_bpm_compute_cpu(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  // Parameters
  mm_stack_t* const mm_stack = gpu_buffer_align_bpm->mm_stack;
  text_collection_t* const text_collection = gpu_buffer_align_bpm->text_collection;
  void* const gpu_buffer = gpu_buffer_align_bpm->buffer;
  const uint64_t used_candidates = gpu_buffer_align_bpm->num_candidates; // Limits
  gpu_bpm_cand_info_t* buffer_candidates = gpu_bpm_buffer_get_candidates_(gpu_buffer);
  gpu_bpm_alg_entry_t* buffer_results = gpu_bpm_buffer_get_alignments_(gpu_buffer);
  // Traverse all candidates
  uint64_t candidate_pos;
  for (candidate_pos=0;candidate_pos<used_candidates;++candidate_pos) {
    // Get Pattern
    bpm_pattern_t bpm_pattern;
    mm_stack_push_state(mm_stack);
    gpu_buffer_align_bpm_retrieve_pattern(gpu_buffer_align_bpm,candidate_pos,&bpm_pattern,mm_stack);
    // Get Candidate Text
    const uint64_t text_length = buffer_candidates->size;
    const uint64_t text_trace_offset = archive_text_retrieve_collection(
        gpu_buffer_align_bpm->archive_text,text_collection,
        buffer_candidates->position,text_length,false,false,mm_stack); // Retrieve text(s)
    const text_trace_t* const text_trace = text_collection_get_trace(text_collection,text_trace_offset);
    const uint8_t* const text = text_trace->text; // Candidate
    // Align BPM & Set result
    uint64_t match_end_column, match_distance;
    bpm_compute_edit_distance(&bpm_pattern,text,text_length,
        &match_distance,&match_end_column,bpm_pattern.pattern_length,true);
    buffer_results->column = match_end_column;
    buffer_results->score = match_distance;
    // Next
    ++buffer_candidates;
    ++buffer_results;
    // Free
    mm_stack_pop_state(mm_stack);
    text_collection_clear(text_collection);
  }
}
/*
 * Send/Receive
 */
void gpu_buffer_align_bpm_send(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  PROF_START(GP_GPU_BUFFER_ALIGN_BPM_SEND);
#ifdef GEM_PROFILE
  const uint64_t max_candidates = gpu_buffer_align_bpm_get_max_candidates(gpu_buffer_align_bpm);
  const uint64_t max_queries = gpu_buffer_align_bpm_get_max_queries(gpu_buffer_align_bpm);
  const uint64_t max_entries = gpu_buffer_align_bpm_get_max_entries(gpu_buffer_align_bpm);
  const uint64_t used_candidates = gpu_buffer_align_bpm->num_candidates;
  const uint64_t used_queries = gpu_buffer_align_bpm->num_queries;
  const uint64_t used_entries = gpu_buffer_align_bpm->num_entries;
  PROF_ADD_COUNTER(GP_GPU_BUFFER_ALIGN_BPM_USAGE_CANDIDATES,(100*used_candidates)/max_candidates);
  PROF_ADD_COUNTER(GP_GPU_BUFFER_ALIGN_BPM_USAGE_QUERIES,(100*used_queries)/max_queries);
  PROF_ADD_COUNTER(GP_GPU_BUFFER_ALIGN_BPM_USAGE_PEQ_ENTRIES,(100*used_entries)/max_entries);
  TIMER_START(&gpu_buffer_align_bpm->timer);
#endif
  // Select computing device
  if (!gpu_buffer_align_bpm->compute_cpu) {
    if (gpu_buffer_align_bpm->num_candidates > 0) {
      gpu_bpm_send_buffer_(
          gpu_buffer_align_bpm->buffer,gpu_buffer_align_bpm->num_entries,
          gpu_buffer_align_bpm->num_queries,gpu_buffer_align_bpm->num_candidates,
          gpu_buffer_align_bpm->query_same_length);
    }
  }
	PROF_STOP(GP_GPU_BUFFER_ALIGN_BPM_SEND);
}
void gpu_buffer_align_bpm_receive(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  PROF_START(GP_GPU_BUFFER_ALIGN_BPM_RECEIVE);
  if (gpu_buffer_align_bpm->num_candidates > 0) {
    // Select computing device
    if (!gpu_buffer_align_bpm->compute_cpu) {
      gpu_bpm_receive_buffer_(gpu_buffer_align_bpm->buffer);
    } else {
      gpu_buffer_align_bpm_compute_cpu(gpu_buffer_align_bpm); // CPU emulated
    }
  }
  PROF_STOP(GP_GPU_BUFFER_ALIGN_BPM_RECEIVE);
#ifdef GEM_PROFILE
  TIMER_STOP(&gpu_buffer_align_bpm->timer);
  COUNTER_ADD(&PROF_GET_TIMER(GP_GPU_BUFFER_ALIGN_BPM_DUTY_CYCLE)->time_ns,gpu_buffer_align_bpm->timer.accumulated);
#endif
}
/*
 * CUDA NOT-Supported
 */
#else
/*
 * Setup
 */
gpu_buffer_align_bpm_t* gpu_buffer_align_bpm_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_no,
    archive_text_t* const archive_text,
    text_collection_t* const text_collection,
    mm_stack_t* const mm_stack) { GEM_CUDA_NOT_SUPPORTED(); return NULL; }
void gpu_buffer_align_bpm_clear(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) { GEM_CUDA_NOT_SUPPORTED(); }
void gpu_buffer_align_bpm_delete(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) { GEM_CUDA_NOT_SUPPORTED(); }
/*
 * Computing Device
 */
void gpu_buffer_align_bpm_set_device_cpu(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) { GEM_CUDA_NOT_SUPPORTED(); }
void gpu_buffer_align_bpm_set_device_gpu(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) { GEM_CUDA_NOT_SUPPORTED(); }
/*
 * Occupancy & Limits
 */
uint64_t gpu_buffer_align_bpm_get_entry_length() { GEM_CUDA_NOT_SUPPORTED(); return 0; }
uint64_t gpu_buffer_align_bpm_get_max_candidates(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
uint64_t gpu_buffer_align_bpm_get_max_queries(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
uint64_t gpu_buffer_align_bpm_get_num_candidates(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
uint64_t gpu_buffer_align_bpm_get_num_queries(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
void gpu_buffer_align_bpm_compute_dimensions(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    bpm_pattern_t* const bpm_pattern,
    bpm_pattern_t* const bpm_pattern_tiles,
    const uint64_t num_candidates,
    uint64_t* const total_entries,
    uint64_t* const total_queries,
    uint64_t* const total_candidates) { GEM_CUDA_NOT_SUPPORTED(); }
bool gpu_buffer_align_bpm_fits_in_buffer(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    const uint64_t total_entries,
    const uint64_t total_queries,
    const uint64_t total_candidates) { GEM_CUDA_NOT_SUPPORTED(); return false; }
/*
 * Accessors
 */
void gpu_buffer_align_bpm_add_pattern(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    bpm_pattern_t* const bpm_pattern,
    bpm_pattern_t* const bpm_pattern_tiles) { GEM_CUDA_NOT_SUPPORTED(); }
void gpu_buffer_align_bpm_add_candidate(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    const uint64_t tile_offset,
    const uint64_t candidate_text_position,
    const uint64_t candidate_length) { GEM_CUDA_NOT_SUPPORTED(); }
void gpu_buffer_align_bpm_get_candidate(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    const uint64_t candidate_offset,
    uint64_t* const candidate_text_position,
    uint32_t* const candidate_length) { GEM_CUDA_NOT_SUPPORTED(); }
void gpu_buffer_align_bpm_get_result(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    const uint64_t candidate_offset,
    uint32_t* const levenshtein_distance,
    uint32_t* const levenshtein_match_pos) { GEM_CUDA_NOT_SUPPORTED(); }
void gpu_buffer_align_bpm_retrieve_pattern(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    const uint64_t candidate_offset,
    bpm_pattern_t* const bpm_pattern,
    mm_stack_t* const mm_stack) { GEM_CUDA_NOT_SUPPORTED(); }
/*
 * Hints
 */
void gpu_buffer_align_bpm_record_query_length(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    const uint64_t query_length) { GEM_CUDA_NOT_SUPPORTED(); }
void gpu_buffer_align_bpm_record_candidates_per_tile(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    const uint64_t num_candidates) { GEM_CUDA_NOT_SUPPORTED(); }
uint64_t gpu_buffer_align_bpm_get_mean_query_length(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
uint64_t gpu_buffer_align_bpm_get_mean_candidates_per_tile(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
/*
 * Send/Receive
 */
void gpu_buffer_align_bpm_send(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) { GEM_CUDA_NOT_SUPPORTED(); }
void gpu_buffer_align_bpm_receive(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) { GEM_CUDA_NOT_SUPPORTED(); }
#endif /* HAVE_CUDA */

