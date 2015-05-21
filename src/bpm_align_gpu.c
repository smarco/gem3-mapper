/*
 * PROJECT: GEMMapper
 * FILE: bpm_align_gpu.c
 * DATE: 04/09/2014
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 *            Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "bpm_align_gpu.h"

/*
 * No-CUDA Support
 */
#ifndef HAVE_CUDA
  // BPM_GPU Setup
  GEM_INLINE bpm_gpu_buffer_collection_t* bpm_gpu_init(
      archive_text_t* const archive_text,const uint32_t num_buffers,const uint32_t buffer_size,
              const int32_t average_query_size,const int32_t candidates_per_query,const bool verbose) {
    GEM_CUDA_NOT_SUPPORTED();
    return NULL;
  }
  GEM_INLINE void bpm_gpu_destroy(bpm_gpu_buffer_collection_t* const buffer_collection) { GEM_CUDA_NOT_SUPPORTED(); }
  GEM_INLINE bool bpm_gpu_support() { return false; }
  // Buffer Accessors
  GEM_INLINE void bpm_gpu_buffer_clear(bpm_gpu_buffer_t* const bpm_gpu_buffer) { GEM_CUDA_NOT_SUPPORTED(); }
  GEM_INLINE uint64_t bpm_gpu_buffer_get_max_candidates(bpm_gpu_buffer_t* const bpm_gpu_buffer) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
  GEM_INLINE uint64_t bpm_gpu_buffer_get_max_queries(bpm_gpu_buffer_t* const bpm_gpu_buffer) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
  GEM_INLINE uint64_t bpm_gpu_buffer_get_num_candidates(bpm_gpu_buffer_t* const bpm_gpu_buffer) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
  GEM_INLINE uint64_t bpm_gpu_buffer_get_num_queries(bpm_gpu_buffer_t* const bpm_gpu_buffer) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
  GEM_INLINE void bpm_gpu_buffer_compute_dimensions(
      bpm_gpu_buffer_t* const bpm_gpu_buffer,const pattern_t* const pattern,
      const uint64_t total_candidates,uint64_t* const total_entries,
      uint64_t* const total_query_chunks,uint64_t* const total_candidate_chunks) { GEM_CUDA_NOT_SUPPORTED(); }
  GEM_INLINE bool bpm_gpu_buffer_fits_in_buffer(
      bpm_gpu_buffer_t* const bpm_gpu_buffer,const uint64_t total_entries,
      const uint64_t total_query_chunks,const uint64_t total_candidate_chunks) { GEM_CUDA_NOT_SUPPORTED(); return false; }
  GEM_INLINE void bpm_gpu_buffer_put_pattern(
      bpm_gpu_buffer_t* const bpm_gpu_buffer,pattern_t* const pattern) { GEM_CUDA_NOT_SUPPORTED(); }
  GEM_INLINE void bpm_gpu_buffer_get_candidate(
      bpm_gpu_buffer_t* const bpm_gpu_buffer,const uint64_t position,
      uint64_t* const candidate_text_position,uint32_t* const candidate_length) { GEM_CUDA_NOT_SUPPORTED(); }
  GEM_INLINE void bpm_gpu_buffer_get_candidate_result(
      bpm_gpu_buffer_t* const bpm_gpu_buffer,const uint64_t position,
      uint32_t* const levenshtein_distance,uint32_t* const levenshtein_match_pos) { GEM_CUDA_NOT_SUPPORTED(); }
  GEM_INLINE void bpm_gpu_buffer_put_candidate(
      bpm_gpu_buffer_t* const bpm_gpu_buffer,const uint64_t candidate_text_position,
          const uint64_t candidate_length,const uint64_t pattern_chunk) { GEM_CUDA_NOT_SUPPORTED(); }
  // Init Buffer
  GEM_INLINE void bpm_gpu_init_buffer(bpm_gpu_buffer_t* const bpm_gpu_buffer) { GEM_CUDA_NOT_SUPPORTED(); }
  // Send/Receive Buffer
  GEM_INLINE void bpm_gpu_buffer_send(bpm_gpu_buffer_t* const bpm_gpu_buffer) { GEM_CUDA_NOT_SUPPORTED(); }
  GEM_INLINE void bpm_gpu_buffer_receive(bpm_gpu_buffer_t* const bpm_gpu_buffer) { GEM_CUDA_NOT_SUPPORTED(); }
#else
/*
 * BPM-GPU Setup
 */
GEM_INLINE bpm_gpu_buffer_collection_t* bpm_gpu_init(
    archive_text_t* const archive_text,const uint32_t num_buffers,const uint32_t buffer_size,
    const int32_t average_query_size,const int32_t candidates_per_query,const bool verbose) {
  GEM_CHECK_POSITIVE(average_query_size);
  GEM_CHECK_POSITIVE(candidates_per_query);
  PROF_START(GP_BPM_GPU_INIT);
  // Allocate Buffer Collection
  bpm_gpu_buffer_collection_t* const buffer_collection = mm_alloc(bpm_gpu_buffer_collection_t);
  buffer_collection->bpm_gpu_buffers = mm_calloc(num_buffers,bpm_gpu_buffer_t,true);
  buffer_collection->num_buffers = num_buffers;
  // Initialize Myers
  const char* const text = (const char* const) dna_text_get_text(archive_text->enc_text);
  const uint64_t text_length = dna_text_get_length(archive_text->enc_text);
  const bpm_gpu_ref_coding_t reference_encoding = (archive_text->explicit_complement) ? GEM_FULL : GEM_ONLY_FORWARD;
  bpm_gpu_init_(&buffer_collection->internal_buffers,num_buffers,CONVERT_B_TO_MB(buffer_size),text,reference_encoding,
      text_length,average_query_size,candidates_per_query,ARCH_SUPPORTED,LOCAL_OR_REMOTE_REFERENCE,verbose);
  // Initialize Buffers
  uint64_t i;
  for (i=0;i<num_buffers;++i) {
    bpm_gpu_buffer_t* const bpm_gpu_buffer = buffer_collection->bpm_gpu_buffers + i;
    bpm_gpu_buffer->buffer = buffer_collection->internal_buffers[i];
    bpm_gpu_buffer->num_PEQ_entries = 0;
    bpm_gpu_buffer->num_queries = 0;
    bpm_gpu_buffer->num_candidates = 0;
    bpm_gpu_buffer->pattern_id = 0;
    TIMER_RESET(&bpm_gpu_buffer->timer);
#ifdef BPM_GPU_PATTERN_DEBUG
    bpm_gpu_buffer->enc_text = enc_text;
#endif
  }
  // Display
  if (verbose) {
    const uint64_t max_queries = bpm_gpu_buffer_get_max_queries_(buffer_collection->bpm_gpu_buffers->buffer);
    const uint64_t max_PEQ_entries = bpm_gpu_buffer_get_max_peq_entries_(buffer_collection->bpm_gpu_buffers->buffer);
    const uint64_t max_candidates = bpm_gpu_buffer_get_max_candidates_(buffer_collection->bpm_gpu_buffers->buffer);
    gem_log("[BPM-GPU Init] Total %lu buffers allocated (Each %lu MB {%lu queries,%lu PeqEntries,%lu candidates})",
        num_buffers,CONVERT_B_TO_MB(buffer_size),max_queries,max_PEQ_entries,max_candidates);
  }
  // Return
  PROF_STOP(GP_BPM_GPU_INIT);
  return buffer_collection;
}
GEM_INLINE void bpm_gpu_destroy(bpm_gpu_buffer_collection_t* const buffer_collection) {
  // Destroy (Myers)
  bpm_gpu_destroy_(&buffer_collection->internal_buffers);
  // Free HUB
  mm_free(buffer_collection->bpm_gpu_buffers);
  mm_free(buffer_collection);
}
GEM_INLINE bool bpm_gpu_support() {
  return bpm_gpu_get_num_supported_devices_()>0;
}
/*
 * Buffer Accessors
 */
GEM_INLINE void bpm_gpu_buffer_clear(bpm_gpu_buffer_t* const bpm_gpu_buffer) {
  bpm_gpu_buffer->num_PEQ_entries = 0;
  bpm_gpu_buffer->num_queries = 0;
  bpm_gpu_buffer->num_candidates = 0;
  bpm_gpu_buffer->pattern_id = 0;
}
GEM_INLINE uint64_t bpm_gpu_buffer_get_max_candidates(bpm_gpu_buffer_t* const bpm_gpu_buffer) {
  return bpm_gpu_buffer_get_max_candidates_(bpm_gpu_buffer->buffer);
}
GEM_INLINE uint64_t bpm_gpu_buffer_get_max_queries(bpm_gpu_buffer_t* const bpm_gpu_buffer) {
  return bpm_gpu_buffer_get_max_queries_(bpm_gpu_buffer->buffer);
}
GEM_INLINE uint64_t bpm_gpu_buffer_get_num_candidates(bpm_gpu_buffer_t* const bpm_gpu_buffer) {
  return bpm_gpu_buffer->num_candidates;
}
GEM_INLINE uint64_t bpm_gpu_buffer_get_num_queries(bpm_gpu_buffer_t* const bpm_gpu_buffer) {
  return bpm_gpu_buffer->num_queries;
}
GEM_INLINE void bpm_gpu_buffer_compute_dimensions(
    bpm_gpu_buffer_t* const bpm_gpu_buffer,const pattern_t* const pattern,
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
GEM_INLINE bool bpm_gpu_buffer_fits_in_buffer(
    bpm_gpu_buffer_t* const bpm_gpu_buffer,const uint64_t total_entries,
    const uint64_t total_query_chunks,const uint64_t total_candidate_chunks) {
  // Get Limits
  const uint64_t max_PEQ_entries =  bpm_gpu_buffer_get_max_peq_entries_(bpm_gpu_buffer->buffer);
  const uint64_t max_queries = bpm_gpu_buffer_get_max_queries_(bpm_gpu_buffer->buffer);
  // Check available space in buffer for the pattern
  if (bpm_gpu_buffer->num_queries+total_query_chunks > max_queries ||
      bpm_gpu_buffer->num_PEQ_entries+total_entries > max_PEQ_entries) {
    // Check if the pattern can fit into an empty buffer
    gem_cond_fatal_error(total_entries > max_PEQ_entries,BPM_GPU_MAX_PATTERN_LENGTH,total_entries,max_PEQ_entries);
    return false;
  }
  // Check available space in buffer for the candidates
  const uint64_t max_candidates = bpm_gpu_buffer_get_max_candidates_(bpm_gpu_buffer->buffer);
  if (bpm_gpu_buffer->num_candidates+total_candidate_chunks > max_candidates) {
    // Check if the pattern can fit into an empty buffer
    gem_cond_fatal_error(total_candidate_chunks > max_candidates,BPM_GPU_MAX_CANDIDATES,total_candidate_chunks,max_candidates);
    return false;
  }
  // Ok, go on
  return true;
}
GEM_INLINE void bpm_gpu_buffer_put_candidate(
    bpm_gpu_buffer_t* const bpm_gpu_buffer,const uint64_t candidate_text_position,
    const uint64_t candidate_length,const uint64_t pattern_chunk) {
  // Insert candidate
  PROF_INC_COUNTER(GP_BPM_GPU_BUFFER_NUM_CANDIDATES);
  PROF_ADD_COUNTER(GP_BPM_GPU_BUFFER_CANDIDATES_LENGTH,candidate_length);
  const uint64_t candidate_offset = bpm_gpu_buffer->num_candidates;
  bpm_gpu_cand_info_t* const query_candidate = bpm_gpu_buffer_get_candidates_(bpm_gpu_buffer->buffer) + candidate_offset;
  query_candidate->query = bpm_gpu_buffer->pattern_id + pattern_chunk;
  query_candidate->position = candidate_text_position;
  query_candidate->size = candidate_length;
  ++(bpm_gpu_buffer->num_candidates);
}
GEM_INLINE void bpm_gpu_buffer_get_candidate(
    bpm_gpu_buffer_t* const bpm_gpu_buffer,const uint64_t position,
    uint64_t* const candidate_text_position,uint32_t* const candidate_length) {
  // Get candidate
  bpm_gpu_cand_info_t* const query_candidate =bpm_gpu_buffer_get_candidates_(bpm_gpu_buffer->buffer) + position;
  *candidate_text_position = query_candidate->position;
  *candidate_length = query_candidate->size;
}
GEM_INLINE void bpm_gpu_buffer_get_candidate_result(
    bpm_gpu_buffer_t* const bpm_gpu_buffer,const uint64_t position,
    uint32_t* const levenshtein_distance,uint32_t* const levenshtein_match_pos) {
  // Get candidate results
  bpm_gpu_res_entry_t* const results_buffer = bpm_gpu_buffer_get_results_(bpm_gpu_buffer->buffer) + position;
  *levenshtein_distance = results_buffer->score;
  *levenshtein_match_pos = results_buffer->column;
}
GEM_INLINE void bpm_gpu_buffer_put_pattern(bpm_gpu_buffer_t* const bpm_gpu_buffer,pattern_t* const pattern) {
  // Fetch dimensions
  const uint64_t key_length = pattern->key_length;
  const bpm_pattern_t* const bpm_pattern = &pattern->bpm_pattern;
  const uint64_t num_entries = bpm_pattern->gpu_num_entries;
  const uint64_t entries_per_chunk = bpm_pattern->gpu_entries_per_chunk;
  const uint64_t num_chunks = bpm_pattern->gpu_num_chunks;
  /*
   * Insert query metadata
   */
  // Insert pattern ID(s)
  bpm_gpu_buffer->pattern_id = bpm_gpu_buffer->num_queries;
  (bpm_gpu_buffer->num_queries) += num_chunks;
  bpm_gpu_qry_info_t* query_info = bpm_gpu_buffer_get_peq_info_(bpm_gpu_buffer->buffer) + bpm_gpu_buffer->pattern_id;
  uint64_t i, remaining_key_length = key_length;
  for (i=0;i<num_chunks;++i,++query_info) {
    query_info->posEntry = bpm_gpu_buffer->num_PEQ_entries + i*entries_per_chunk;
    const uint64_t chunk_length = (remaining_key_length > BPM_GPU_PATTERN_ENTRY_LENGTH) ?
        BPM_GPU_PATTERN_ENTRY_LENGTH : remaining_key_length;
    query_info->size = chunk_length;
    remaining_key_length -= chunk_length;
  }
  /*
   * [DTO] Compile PEQ pattern
   */
  // Insert pattern query/queries
  const uint64_t PEQ_entry_offset = bpm_gpu_buffer->num_PEQ_entries;
  bpm_gpu_qry_entry_t* const query_pattern = bpm_gpu_buffer_get_peq_entries_(bpm_gpu_buffer->buffer) + PEQ_entry_offset;
  bpm_gpu_buffer->num_PEQ_entries += num_entries;
  // Copy PEQ pattern
  const uint64_t bpm_pattern_num_words = pattern->bpm_pattern.pattern_num_words;
  const uint64_t bpm_gpu_pattern_length = num_entries*BPM_GPU_PATTERN_ENTRY_LENGTH;
  const uint64_t bpm_gpu_pattern_num_words = bpm_gpu_pattern_length/BPM_ALIGN_WORD_LENGTH;
  const uint32_t* PEQ = (uint32_t*) pattern->bpm_pattern.PEQ;
  uint64_t entry=0, subentry=0;
  for (i=0;i<bpm_pattern_num_words;++i) {
    // Update location
    if (subentry==BPM_GPU_PATTERN_NUM_SUB_ENTRIES) {
      subentry = 0; ++entry;
    }
    // Copy pattern
    uint8_t enc_char;
    for (enc_char=0;enc_char<DNA__N_RANGE;++enc_char) {
      query_pattern[entry].bitmap[enc_char][subentry] = *PEQ; ++PEQ;
      query_pattern[entry].bitmap[enc_char][subentry+1] = *PEQ; ++PEQ;
    }
    subentry += 2;
  }
  for (;i<bpm_gpu_pattern_num_words;++i) {
    // Update location
    if (subentry==BPM_GPU_PATTERN_NUM_SUB_ENTRIES) {
      subentry = 0; ++entry;
    }
    // Copy pattern
    uint8_t enc_char;
    for (enc_char=0;enc_char<DNA__N_RANGE;++enc_char) {
      query_pattern[entry].bitmap[enc_char][subentry] = UINT32_ONES;
      query_pattern[entry].bitmap[enc_char][subentry+1] = UINT32_ONES;
    }
    subentry += 2;
  }
}
/*
 * Init the local thread BMP Buffers
 */
GEM_INLINE void bpm_gpu_init_buffer(bpm_gpu_buffer_t* const bpm_gpu_buffer) {
  PROF_START(GP_BPM_GPU_BUFFER_INIT);
  bpm_gpu_init_buffer_(bpm_gpu_buffer->buffer);
  PROF_STOP(GP_BPM_GPU_BUFFER_INIT);
}
/*
 * Send/Receive Buffer
 */
GEM_INLINE void bpm_gpu_buffer_send(bpm_gpu_buffer_t* const bpm_gpu_buffer) {
  PROF_START_TIMER(GP_BPM_GPU_BUFFER_SEND);
  PROF_BLOCK() {
    #ifndef GEM_NOPROFILE
    const uint64_t max_candidates = bpm_gpu_buffer_get_max_candidates_(bpm_gpu_buffer->buffer);
    const uint64_t max_queries = bpm_gpu_buffer_get_max_queries_(bpm_gpu_buffer->buffer);
    const uint64_t max_peq_entries = bpm_gpu_buffer_get_max_peq_entries_(bpm_gpu_buffer->buffer);
    const uint64_t used_candidates = bpm_gpu_buffer->num_candidates;
    const uint64_t used_queries = bpm_gpu_buffer->num_queries;
    const uint64_t used_peq_entries = bpm_gpu_buffer->num_PEQ_entries;
    PROF_ADD_COUNTER(GP_BPM_GPU_BUFFER_USAGE_CANDIDATES,(100*used_candidates)/max_candidates);
    PROF_ADD_COUNTER(GP_BPM_GPU_BUFFER_USAGE_QUERIES,(100*used_queries)/max_queries);
    PROF_ADD_COUNTER(GP_BPM_GPU_BUFFER_USAGE_PEQ_ENTRIES,(100*used_peq_entries)/max_peq_entries);
    TIMER_START(&bpm_gpu_buffer->timer);
    #endif
  }
	bpm_gpu_send_buffer_(bpm_gpu_buffer->buffer,
      bpm_gpu_buffer->num_PEQ_entries,bpm_gpu_buffer->num_queries,bpm_gpu_buffer->num_candidates);
	PROF_STOP_TIMER(GP_BPM_GPU_BUFFER_SEND);
}
GEM_INLINE void bpm_gpu_buffer_receive(bpm_gpu_buffer_t* const bpm_gpu_buffer) {
	bpm_gpu_receive_buffer_(bpm_gpu_buffer->buffer);
  PROF_BLOCK() {
    #ifndef GEM_NOPROFILE
    TIMER_STOP(&bpm_gpu_buffer->timer);
    COUNTER_ADD(&PROF_GET_TIMER(GP_BPM_GPU_BUFFER_CHECK_TIME)->time_ns,bpm_gpu_buffer->timer.accumulated);
    #endif
  }
}
#endif /* HAVE_CUDA */

