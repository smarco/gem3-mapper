/*
 * PROJECT: GEMMapper
 * FILE: gpu_buffer_align_bpm.c
 * DATE: 04/09/2014
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 *            Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "gpu_buffer_align_bpm.h"
#include "../resources/gpu_modules/gpu_interface.h"

/*
 * Constants :: GPU Align-BMP
 */
#define GPU_ALIGN_BPM_NUM_SUB_ENTRIES  GPU_BPM_PEQ_SUBENTRIES
#define GPU_ALIGN_BPM_ENTRY_LENGTH     GPU_BPM_PEQ_ENTRY_LENGTH
#define GPU_ALIGN_BPM_SUBENTRY_LENGTH  GPU_BPM_PEQ_SUBENTRY_LENGTH
#define GPU_ALIGN_BPM_ENTRY_SIZE       (GPU_ALIGN_BPM_ENTRY_LENGTH/UINT8_SIZE)
#define GPU_ALIGN_BPM_ALPHABET_LENGTH  GPU_BPM_PEQ_ALPHABET_SIZE

/*
 * Constants :: Buffer Hints
 */
#define GPU_ALIGN_BPM_AVERAGE_QUERY_SIZE        150
#define GPU_ALIGN_BPM_CANDIDATES_PER_QUERY      20

/*
 * Pattern Setup
 */
GEM_INLINE void gpu_bpm_pattern_compile(bpm_pattern_t* const bpm_pattern,const uint64_t max_error) {
  // Init BPM-GPU Dimensions
  bpm_pattern->gpu_num_entries = DIV_CEIL(bpm_pattern->pattern_length,GPU_ALIGN_BPM_ENTRY_LENGTH);
  bpm_pattern->gpu_entries_per_chunk = DIV_CEIL(max_error,GPU_ALIGN_BPM_ENTRY_LENGTH);
  bpm_pattern->gpu_num_chunks = DIV_CEIL(bpm_pattern->gpu_num_entries,bpm_pattern->gpu_entries_per_chunk);
}
GEM_INLINE uint64_t gpu_bpm_pattern_get_entry_length() {
  return GPU_ALIGN_BPM_ENTRY_LENGTH;
}

/*
 * CUDA Support
 */
#ifdef HAVE_CUDA
/*
 * Setup
 */
GEM_INLINE gpu_buffer_align_bpm_t* gpu_buffer_align_bpm_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,const uint64_t buffer_no) {
  PROF_START(GP_GPU_BUFFER_ALIGN_BPM_ALLOC);
  // Alloc
  gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm = mm_alloc(gpu_buffer_align_bpm_t);
  // Init
  gpu_buffer_align_bpm->buffer = gpu_buffer_collection_get_buffer(gpu_buffer_collection,buffer_no);
  gpu_buffer_align_bpm->averageQuerySize = GPU_ALIGN_BPM_AVERAGE_QUERY_SIZE;
  gpu_buffer_align_bpm->candidatesPerQuery = GPU_ALIGN_BPM_CANDIDATES_PER_QUERY;
  gpu_buffer_align_bpm->candidates_same_length = 0;
  gpu_buffer_align_bpm->num_PEQ_entries = 0;
  gpu_buffer_align_bpm->num_queries = 0;
  gpu_buffer_align_bpm->num_candidates = 0;
  gpu_buffer_align_bpm->pattern_id = 0;
  gpu_buffer_align_bpm->compute_cpu = false;
  TIMER_RESET(&gpu_buffer_align_bpm->timer);
  // Init buffer
  gpu_alloc_buffer_(gpu_buffer_align_bpm->buffer);
  gpu_bpm_init_buffer_(gpu_buffer_align_bpm->buffer,
      gpu_buffer_align_bpm->averageQuerySize,gpu_buffer_align_bpm->candidatesPerQuery);
  PROF_STOP(GP_GPU_BUFFER_ALIGN_BPM_ALLOC);
  // Return
  return gpu_buffer_align_bpm;
}
GEM_INLINE void gpu_buffer_align_bpm_clear(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  // Init buffer
  gpu_bpm_init_buffer_(gpu_buffer_align_bpm->buffer,
      gpu_buffer_align_bpm->averageQuerySize,gpu_buffer_align_bpm->candidatesPerQuery);
  // Clear
  gpu_buffer_align_bpm->num_PEQ_entries = 0;
  gpu_buffer_align_bpm->num_queries = 0;
  gpu_buffer_align_bpm->num_candidates = 0;
  gpu_buffer_align_bpm->pattern_id = 0;
}
GEM_INLINE void gpu_buffer_align_bpm_delete(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  mm_free(gpu_buffer_align_bpm);
}
/*
 * Computing Device
 */
GEM_INLINE void gpu_buffer_align_bpm_set_device_cpu(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  gpu_buffer_align_bpm->compute_cpu = true;
}
GEM_INLINE void gpu_buffer_align_bpm_set_device_gpu(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  gpu_buffer_align_bpm->compute_cpu = false;
}
/*
 * Occupancy & Limits
 */
GEM_INLINE uint64_t gpu_buffer_align_bpm_get_max_candidates(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  return gpu_bpm_buffer_get_max_candidates_(gpu_buffer_align_bpm->buffer);
}
GEM_INLINE uint64_t gpu_buffer_align_bpm_get_max_queries(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  return gpu_bpm_buffer_get_max_queries_(gpu_buffer_align_bpm->buffer);
}
GEM_INLINE uint64_t gpu_buffer_align_bpm_get_max_entries(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  return gpu_bpm_buffer_get_max_peq_entries_(gpu_buffer_align_bpm->buffer);
}
GEM_INLINE uint64_t gpu_buffer_align_bpm_get_num_candidates(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  return gpu_buffer_align_bpm->num_candidates;
}
GEM_INLINE uint64_t gpu_buffer_align_bpm_get_num_queries(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  return gpu_buffer_align_bpm->num_queries;
}
GEM_INLINE void gpu_buffer_align_bpm_compute_dimensions(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,const pattern_t* const pattern,
    const uint64_t total_candidates,uint64_t* const total_entries,
    uint64_t* const total_query_chunks,uint64_t* const total_candidate_chunks) {
  const bpm_pattern_t* const bpm_pattern = &pattern->bpm_pattern;
  // Calculate dimensions
  const uint64_t pattern_num_entries = bpm_pattern->gpu_num_entries;
  const uint64_t pattern_num_chunks = bpm_pattern->gpu_num_chunks;
  *total_entries += pattern_num_entries;
  *total_query_chunks += pattern_num_chunks;
  *total_candidate_chunks += pattern_num_chunks*total_candidates;
}
GEM_INLINE bool gpu_buffer_align_bpm_fits_in_buffer(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,const uint64_t total_entries,
    const uint64_t total_query_chunks,const uint64_t total_candidate_chunks) {
  // Get Limits
  const uint64_t max_PEQ_entries =  gpu_buffer_align_bpm_get_max_entries(gpu_buffer_align_bpm);
  const uint64_t max_queries = gpu_buffer_align_bpm_get_max_queries(gpu_buffer_align_bpm);
  // Check available space in buffer for the pattern
  if (gpu_buffer_align_bpm->num_queries+total_query_chunks > max_queries ||
      gpu_buffer_align_bpm->num_PEQ_entries+total_entries > max_PEQ_entries) {
    // Check if the pattern can fit into an empty buffer
    gem_cond_fatal_error(total_entries > max_PEQ_entries,GPU_BPM_MAX_PATTERN_LENGTH,total_entries,max_PEQ_entries);
    return false;
  }
  // Check available space in buffer for the candidates
  const uint64_t max_candidates = gpu_buffer_align_bpm_get_max_candidates(gpu_buffer_align_bpm);
  if (gpu_buffer_align_bpm->num_candidates+total_candidate_chunks > max_candidates) {
    // Check if the pattern can fit into an empty buffer
    gem_cond_fatal_error(total_candidate_chunks > max_candidates,GPU_BPM_MAX_CANDIDATES,total_candidate_chunks,max_candidates);
    return false;
  }
  // Ok, go on
  return true;
}
/*
 * Accessors
 */
GEM_INLINE void gpu_buffer_align_bpm_pattern_decompile(
    gpu_bpm_qry_entry_t* const pattern_entry,const uint32_t size,
    const bpm_pattern_t* const bpm_pattern) {
//  /* BMP Pattern */
//  uint64_t* PEQ;              // Pattern equalities (Bit vector for Myers-DP)
//  uint64_t pattern_word_size; // Word size (in bytes)
//  uint64_t pattern_length;    // Length
//  uint64_t pattern_num_words; // ceil(Length / |w|)
//  uint64_t pattern_mod;       // Length % |w|
//  uint64_t PEQ_length;        // ceil(Length / |w|) * |w|
//  /* BPM Auxiliary data */
//  uint64_t* P;
//  uint64_t* M;
//  uint64_t* level_mask;
//  int64_t* score;
//  int64_t* init_score;
//  uint64_t* pattern_left;
//  /* BPM chunks (Pattern split in chunks) */
//  uint64_t words_per_chunk;
//  uint64_t num_pattern_chunks;
//  /* BPM-GPU Dimensions */
//  uint64_t gpu_num_entries;
//  uint64_t gpu_entries_per_chunk;
//  uint64_t gpu_num_chunks;
//  bpm_pattern_t* bpm_pattern_chunks;
}
GEM_INLINE void gpu_buffer_align_bpm_pattern_compile(
    gpu_bpm_qry_entry_t* const pattern_entry,const bpm_pattern_t* const bpm_pattern) {
  // Copy PEQ pattern
  const uint64_t bpm_pattern_num_words = bpm_pattern->pattern_num_words;
  const uint64_t gpu_buffer_align_bpm_pattern_length = bpm_pattern->gpu_num_entries*GPU_ALIGN_BPM_ENTRY_LENGTH;
  const uint64_t gpu_buffer_align_bpm_pattern_num_words = gpu_buffer_align_bpm_pattern_length/BPM_ALIGN_WORD_LENGTH;
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
  for (;i<gpu_buffer_align_bpm_pattern_num_words;++i) {
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
GEM_INLINE void gpu_buffer_align_bpm_add_pattern(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,pattern_t* const pattern) {
  // Parameters (Buffer & Dimensions)
  void* const gpu_buffer = gpu_buffer_align_bpm->buffer;
  const uint64_t num_PEQ_entries = gpu_buffer_align_bpm->num_PEQ_entries;
  const uint64_t key_length = pattern->key_length;
  bpm_pattern_t* bpm_pattern = &pattern->bpm_pattern;
  const uint64_t num_entries = bpm_pattern->gpu_num_entries;
  const uint64_t entries_per_chunk = bpm_pattern->gpu_entries_per_chunk;
  const uint64_t num_chunks = bpm_pattern->gpu_num_chunks;
  // Add query metadata (Add pattern ID(s))
  const uint32_t buffer_pattern_info_offset = gpu_buffer_align_bpm->num_queries;
  (gpu_buffer_align_bpm->num_queries) += num_chunks;
  gpu_buffer_align_bpm->pattern_id = buffer_pattern_info_offset;
  // Add query info (for all chunks in case of tiling)
  gpu_bpm_qry_info_t* buffer_pattern_info = gpu_bpm_buffer_get_peq_info_(gpu_buffer) + buffer_pattern_info_offset;
  uint64_t i, chunk_length, remaining_key_length = key_length;
  for (i=0;i<num_chunks;++i,++buffer_pattern_info) {
    buffer_pattern_info->posEntry = num_PEQ_entries + i*entries_per_chunk;
    chunk_length = MIN(remaining_key_length,GPU_ALIGN_BPM_ENTRY_LENGTH);
    buffer_pattern_info->size = chunk_length;
    remaining_key_length -= chunk_length;
  }
  // [DTO] Compile PEQ pattern (Add pattern entries)
  gpu_bpm_qry_entry_t* const buffer_pattern_entry = gpu_bpm_buffer_get_peq_entries_(gpu_buffer) + num_PEQ_entries;
  gpu_buffer_align_bpm->num_PEQ_entries += num_entries;
  gpu_buffer_align_bpm_pattern_compile(buffer_pattern_entry,bpm_pattern);
}
GEM_INLINE void gpu_buffer_align_bpm_add_candidate(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,const uint64_t candidate_text_position,
    const uint64_t candidate_length,const uint64_t pattern_chunk) {
  // Insert candidate
  PROF_INC_COUNTER(GP_GPU_BUFFER_ALIGN_BPM_NUM_QUERIES);
  PROF_ADD_COUNTER(GP_GPU_BUFFER_ALIGN_BPM_CANDIDATE_LENGTH,candidate_length);
  const uint64_t candidate_offset = gpu_buffer_align_bpm->num_candidates;
  gpu_bpm_cand_info_t* const candidate = gpu_bpm_buffer_get_candidates_(gpu_buffer_align_bpm->buffer) + candidate_offset;
  candidate->query = gpu_buffer_align_bpm->pattern_id + pattern_chunk;
  candidate->position = candidate_text_position;
  candidate->size = candidate_length;
  ++(gpu_buffer_align_bpm->num_candidates);
}
GEM_INLINE void gpu_buffer_align_bpm_get_candidate(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,const uint64_t buffer_pos,
    uint64_t* const candidate_text_position,uint32_t* const candidate_length) {
  // Get candidate
  gpu_bpm_cand_info_t* const candidate = gpu_bpm_buffer_get_candidates_(gpu_buffer_align_bpm->buffer) + buffer_pos;
  *candidate_text_position = candidate->position;
  *candidate_length = candidate->size;
}
GEM_INLINE void gpu_buffer_align_bpm_get_result(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,const uint64_t buffer_pos,
    uint32_t* const levenshtein_distance,uint32_t* const levenshtein_match_pos) {
  // Get candidate results
  gpu_bpm_alg_entry_t* const result = gpu_bpm_buffer_get_alignments_(gpu_buffer_align_bpm->buffer) + buffer_pos;
  *levenshtein_distance = result->score;
  *levenshtein_match_pos = result->column;
}
/*
 * CPU emulated
 */
GEM_INLINE void gpu_buffer_align_bpm_compute_cpu(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  GEM_NOT_IMPLEMENTED();
//  // Parameters
//  void* const gpu_buffer = gpu_buffer_align_bpm->buffer;
//  gpu_bpm_qry_info_t* const buffer_pattern_info = gpu_bpm_buffer_get_peq_info_(gpu_buffer);
//  gpu_bpm_qry_entry_t* const buffer_pattern_entry = gpu_bpm_buffer_get_peq_entries_(gpu_buffer);
//  gpu_bpm_cand_info_t* buffer_candidates = gpu_bpm_buffer_get_candidates_(gpu_buffer);
//  gpu_bpm_alg_entry_t* buffer_results = gpu_bpm_buffer_get_alignments_(gpu_buffer);
//  // Limits
//  const uint64_t used_candidates = gpu_buffer_align_bpm->num_candidates;
//  const uint64_t used_queries = gpu_buffer_align_bpm->num_queries;
//  // Traverse all candidates
//  uint64_t i;
//  for (i=0;i<used_candidates;++i,++buffer_candidates,++buffer_results) {
//    // Get Pattern
//    gpu_bpm_qry_info_t* const pattern_info = buffer_pattern_info + buffer_candidates->query;
//    gpu_bpm_qry_entry_t* const pattern_entry = buffer_pattern_entry+pattern_info->posEntry;
//    bpm_pattern_t bpm_pattern;
//    gpu_buffer_align_bpm_pattern_decompile(pattern_entry,pattern_info->size,&bpm_pattern);
//    // Get candidate text
//    const uint64_t text_length = buffer_candidates->size;
//    const uint64_t text_trace_offset = archive_text_retrieve(archive_text,text_collection,
//        buffer_candidates->position,text_length,false,mm_stack); // Retrieve text(s)
//    const text_trace_t* const text_trace = text_collection_get_trace(text_collection,text_trace_offset);
//    const uint8_t* const text = text_trace->text; // Candidate
//    // Align BPM
//    uint64_t match_end_column, match_distance;
//    bpm_compute_edit_distance_cutoff(&bpm_pattern,text,text_length,
//        &match_end_column,&match_distance,max_error,true);
//    // Set result
//    buffer_results->column = match_end_column;
//    buffer_results->score = match_distance;
//  }
}
/*
 * Send/Receive
 */
GEM_INLINE void gpu_buffer_align_bpm_send(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  PROF_START(GP_GPU_BUFFER_ALIGN_BPM_SEND);
#ifdef GEM_PROFILE
  const uint64_t max_candidates = gpu_buffer_align_bpm_get_max_candidates(gpu_buffer_align_bpm);
  const uint64_t max_queries = gpu_buffer_align_bpm_get_max_queries(gpu_buffer_align_bpm);
  const uint64_t max_peq_entries = gpu_buffer_align_bpm_get_max_entries(gpu_buffer_align_bpm);
  const uint64_t used_candidates = gpu_buffer_align_bpm->num_candidates;
  const uint64_t used_queries = gpu_buffer_align_bpm->num_queries;
  const uint64_t used_peq_entries = gpu_buffer_align_bpm->num_PEQ_entries;
  PROF_ADD_COUNTER(GP_GPU_BUFFER_ALIGN_BPM_USAGE_CANDIDATES,(100*used_candidates)/max_candidates);
  PROF_ADD_COUNTER(GP_GPU_BUFFER_ALIGN_BPM_USAGE_QUERIES,(100*used_queries)/max_queries);
  PROF_ADD_COUNTER(GP_GPU_BUFFER_ALIGN_BPM_USAGE_PEQ_ENTRIES,(100*used_peq_entries)/max_peq_entries);
  TIMER_START(&gpu_buffer_align_bpm->timer);
#endif
  // Select computing device
  if (!gpu_buffer_align_bpm->compute_cpu) {
    gpu_bpm_send_buffer_(
        gpu_buffer_align_bpm->buffer,gpu_buffer_align_bpm->num_PEQ_entries,
        gpu_buffer_align_bpm->num_queries,gpu_buffer_align_bpm->num_candidates,
        gpu_buffer_align_bpm->candidates_same_length);
  }
	PROF_STOP(GP_GPU_BUFFER_ALIGN_BPM_SEND);
}
GEM_INLINE void gpu_buffer_align_bpm_receive(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  PROF_START(GP_GPU_BUFFER_ALIGN_BPM_RECEIVE);
  // Select computing device
  if (!gpu_buffer_align_bpm->compute_cpu) {
    gpu_bpm_receive_buffer_(gpu_buffer_align_bpm->buffer);
  } else {
    // CPU emulated
    gpu_buffer_align_bpm_compute_cpu(gpu_buffer_align_bpm);
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
GEM_INLINE gpu_buffer_align_bpm_t* gpu_buffer_align_bpm_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,const uint64_t buffer_no) { GEM_CUDA_NOT_SUPPORTED(); return NULL; }
GEM_INLINE void gpu_buffer_align_bpm_clear(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) { GEM_CUDA_NOT_SUPPORTED(); }
GEM_INLINE void gpu_buffer_align_bpm_delete(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) { GEM_CUDA_NOT_SUPPORTED(); }
/*
 * Computing Device
 */
GEM_INLINE void gpu_buffer_align_bpm_set_device_cpu(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) { GEM_CUDA_NOT_SUPPORTED(); }
GEM_INLINE void gpu_buffer_align_bpm_set_device_gpu(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) { GEM_CUDA_NOT_SUPPORTED(); }
/*
 * Occupancy & Limits
 */
GEM_INLINE uint64_t gpu_buffer_align_bpm_get_max_candidates(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
GEM_INLINE uint64_t gpu_buffer_align_bpm_get_max_queries(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
GEM_INLINE uint64_t gpu_buffer_align_bpm_get_num_candidates(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
GEM_INLINE uint64_t gpu_buffer_align_bpm_get_num_queries(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
GEM_INLINE void gpu_buffer_align_bpm_compute_dimensions(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,const pattern_t* const pattern,
    const uint64_t total_candidates,uint64_t* const total_entries,
    uint64_t* const total_query_chunks,uint64_t* const total_candidate_chunks) { GEM_CUDA_NOT_SUPPORTED(); }
GEM_INLINE bool gpu_buffer_align_bpm_fits_in_buffer(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,const uint64_t total_entries,
    const uint64_t total_query_chunks,const uint64_t total_candidate_chunks) { GEM_CUDA_NOT_SUPPORTED(); return false; }
/*
 * Accessors
 */
GEM_INLINE void gpu_buffer_align_bpm_add_pattern(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,pattern_t* const pattern) { GEM_CUDA_NOT_SUPPORTED(); }
GEM_INLINE void gpu_buffer_align_bpm_add_candidate(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,const uint64_t candidate_text_position,
        const uint64_t candidate_length,const uint64_t pattern_chunk) { GEM_CUDA_NOT_SUPPORTED(); }
GEM_INLINE void gpu_buffer_align_bpm_get_candidate(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,const uint64_t buffer_pos,
    uint64_t* const candidate_text_position,uint32_t* const candidate_length) { GEM_CUDA_NOT_SUPPORTED(); }
GEM_INLINE void gpu_buffer_align_bpm_get_result(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,const uint64_t buffer_pos,
    uint32_t* const levenshtein_distance,uint32_t* const levenshtein_match_pos) { GEM_CUDA_NOT_SUPPORTED(); }
/*
 * Send/Receive
 */
GEM_INLINE void gpu_buffer_align_bpm_send(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) { GEM_CUDA_NOT_SUPPORTED(); }
GEM_INLINE void gpu_buffer_align_bpm_receive(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) { GEM_CUDA_NOT_SUPPORTED(); }
#endif /* HAVE_CUDA */

