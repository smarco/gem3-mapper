/*
 * PROJECT: GEMMapper
 * FILE: bpm_align_gpu.c
 * DATE: 04/09/2014
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 *            Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "bpm_align_gpu.h"
#include "../resources/myers_gpu/myers-interface.h"

/*
 * Constants
 */
#define BPM_GPU_PATTERN_NUM_SUB_ENTRIES  BMP_GPU_PEQ_SUBENTRIES
#define BPM_GPU_PATTERN_ENTRY_LENGTH     BMP_GPU_PEQ_ENTRY_LENGTH
#define BPM_GPU_PATTERN_SUBENTRY_LENGTH  BMP_GPU_PEQ_SUBENTRY_LENGTH
#define BPM_GPU_PATTERN_ENTRY_SIZE       (BPM_GPU_PATTERN_ENTRY_LENGTH/UINT8_SIZE)
#define BPM_GPU_PATTERN_ALPHABET_LENGTH  BMP_GPU_PEQ_ALPHABET_SIZE
#define BPM_GPU_BUFFER_SIZE              BUFFER_SIZE_32M

/*
 * No-CUDA Support
 */
#ifndef HAVE_CUDA
  // BPM_GPU Setup
  GEM_INLINE bpm_gpu_buffer_collection_t* bpm_gpu_init(
      dna_text_t* const enc_text,const uint32_t num_buffers,
      const int32_t average_query_size,const int32_t candidates_per_query) {
    GEM_CUDA_NOT_SUPPORTED();
    return NULL;
  }
  GEM_INLINE void bpm_gpu_destroy(bpm_gpu_buffer_collection_t* const buffer_collection) { GEM_CUDA_NOT_SUPPORTED(); }
  GEM_INLINE bool bpm_gpu_support() { return false; }
  // Buffer Accessors
  GEM_INLINE void bpm_buffer_clear(bpm_gpu_buffer_t* const bpm_gpu_buffer) { GEM_CUDA_NOT_SUPPORTED(); }
  GEM_INLINE uint64_t bpm_gpu_buffer_get_max_candidates(bpm_gpu_buffer_t* const bpm_gpu_buffer) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
  GEM_INLINE uint64_t bpm_gpu_buffer_get_max_queries(bpm_gpu_buffer_t* const bpm_gpu_buffer) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
  GEM_INLINE uint64_t bpm_gpu_buffer_get_num_candidates(bpm_gpu_buffer_t* const bpm_gpu_buffer) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
  GEM_INLINE uint64_t bpm_gpu_buffer_get_num_queries(bpm_gpu_buffer_t* const bpm_gpu_buffer) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
  GEM_INLINE bool bpm_gpu_buffer_fits_in_buffer(
      bpm_gpu_buffer_t* const bpm_gpu_buffer,
      const uint64_t num_patterns,const uint64_t total_pattern_length,const uint64_t total_candidates) { GEM_CUDA_NOT_SUPPORTED(); return false; }
  GEM_INLINE void bpm_gpu_buffer_put_pattern(
      bpm_gpu_buffer_t* const bpm_gpu_buffer,pattern_t* const pattern) { GEM_CUDA_NOT_SUPPORTED(); }
  GEM_INLINE void bpm_gpu_buffer_get_candidate(
      bpm_gpu_buffer_t* const bpm_gpu_buffer,const uint64_t position,
      uint32_t* const candidate_text_position,uint32_t* const candidate_length) { GEM_CUDA_NOT_SUPPORTED(); }
  GEM_INLINE void bpm_gpu_buffer_get_candidate_result(
      bpm_gpu_buffer_t* const bpm_gpu_buffer,const uint64_t position,
      uint32_t* const levenshtein_distance,uint32_t* const levenshtein_match_pos) { GEM_CUDA_NOT_SUPPORTED(); }
  GEM_INLINE void bpm_gpu_buffer_put_candidate(
      bpm_gpu_buffer_t* const bpm_gpu_buffer,
      const uint64_t candidate_text_position,const uint64_t candidate_length) { GEM_CUDA_NOT_SUPPORTED(); }
  // Send/Receive Buffer
  GEM_INLINE void bpm_gpu_buffer_send(bpm_gpu_buffer_t* const bpm_gpu_buffer) { GEM_CUDA_NOT_SUPPORTED(); }
  GEM_INLINE void bpm_gpu_buffer_receive(bpm_gpu_buffer_t* const bpm_gpu_buffer) { GEM_CUDA_NOT_SUPPORTED(); }
#else
/*
 * BPM-GPU Setup
 */
GEM_INLINE bpm_gpu_buffer_collection_t* bpm_gpu_init(
    dna_text_t* const enc_text,const uint32_t num_buffers,
    const int32_t average_query_size,const int32_t candidates_per_query) {
  GEM_CHECK_POSITIVE(average_query_size); // TODO More checkers here
  GEM_CHECK_POSITIVE(candidates_per_query);
  // Allocate Buffer Collection
  bpm_gpu_buffer_collection_t* const buffer_collection = mm_alloc(bpm_gpu_buffer_collection_t);
  buffer_collection->bpm_gpu_buffers = mm_calloc(num_buffers,bpm_gpu_buffer_t,true);
  buffer_collection->num_buffers = num_buffers;
  // Initialize Myers
  const char* const text = (const char* const) dna_text_get_buffer(enc_text);
  const uint64_t text_length = dna_text_get_length(enc_text);
  bpm_gpu_init_(&buffer_collection->internal_buffers,num_buffers,CONVERT_B_TO_MB(BPM_GPU_BUFFER_SIZE),
      text,GEM,text_length,average_query_size,candidates_per_query,ARCH_SUPPORTED);
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
  // Return
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
GEM_INLINE void bpm_buffer_clear(bpm_gpu_buffer_t* const bpm_gpu_buffer) {
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
GEM_INLINE bool bpm_gpu_buffer_fits_in_buffer(
    bpm_gpu_buffer_t* const bpm_gpu_buffer,
    const uint64_t num_patterns,const uint64_t total_pattern_length,const uint64_t total_candidates) {
  // Calculate dimensions
  const uint64_t pattern_PEQ_num_entries = DIV_CEIL(total_pattern_length,BPM_GPU_PATTERN_ENTRY_LENGTH);
  // Get Limits
  const uint64_t max_PEQ_entries =  bpm_gpu_buffer_get_max_peq_entries_(bpm_gpu_buffer->buffer);
  const uint64_t max_queries = bpm_gpu_buffer_get_max_queries_(bpm_gpu_buffer->buffer);
  // Check available space in buffer for the pattern
  if (bpm_gpu_buffer->num_queries+num_patterns > max_queries ||
      bpm_gpu_buffer->num_PEQ_entries+pattern_PEQ_num_entries > max_PEQ_entries) {
    // Check if the pattern can fit into an empty buffer
    gem_cond_fatal_error(pattern_PEQ_num_entries > max_PEQ_entries,
        BPM_GPU_MAX_PATTERN_LENGTH,pattern_PEQ_num_entries,max_PEQ_entries);
    return false;
  }
  // Check available space in buffer for the candidates
  const uint64_t max_candidates = bpm_gpu_buffer_get_max_candidates_(bpm_gpu_buffer->buffer);
  if (bpm_gpu_buffer->num_candidates+total_candidates > max_candidates) {
    // Check if the pattern can fit into an empty buffer
    gem_cond_fatal_error(total_candidates > max_candidates,
        BPM_GPU_MAX_CANDIDATES,total_candidates,max_candidates);
    return false;
  }
  // Ok, go on
  return true;
}
GEM_INLINE void bpm_gpu_buffer_get_candidate(
    bpm_gpu_buffer_t* const bpm_gpu_buffer,const uint64_t position,
    uint32_t* const candidate_text_position,uint32_t* const candidate_length) {
  // Get candidate
  bpm_gpu_cand_info_t* const query_candidate =bpm_gpu_buffer_get_candidates_(bpm_gpu_buffer->buffer) + position;
  *candidate_text_position = query_candidate->position;
  *candidate_length = query_candidate->size;
}
GEM_INLINE void bpm_gpu_buffer_get_candidate_result(
    bpm_gpu_buffer_t* const bpm_gpu_buffer,const uint64_t position,
    uint32_t* const levenshtein_distance,uint32_t* const levenshtein_match_pos) {
  // Get candidate results
  bpm_gpu_res_entry_t* const results_buffer =
      bpm_gpu_buffer_get_results_(bpm_gpu_buffer->buffer) + position;
  *levenshtein_distance = results_buffer->score;
  *levenshtein_match_pos = results_buffer->column;
}
GEM_INLINE void bpm_gpu_buffer_put_candidate(
    bpm_gpu_buffer_t* const bpm_gpu_buffer,
    const uint64_t candidate_text_position,const uint64_t candidate_length) {
  // Insert candidate
  bpm_gpu_cand_info_t* const query_candidate =
      bpm_gpu_buffer_get_candidates_(bpm_gpu_buffer->buffer) + bpm_gpu_buffer->num_candidates;
  query_candidate->query = bpm_gpu_buffer->pattern_id;
  query_candidate->position = candidate_text_position;
  query_candidate->size = candidate_length;
  ++(bpm_gpu_buffer->num_candidates);
}
// Debug Guard
# ifndef BPM_GPU_PATTERN_DEBUG
GEM_INLINE void bpm_gpu_buffer_put_pattern(
    bpm_gpu_buffer_t* const bpm_gpu_buffer,pattern_t* const pattern) {
  const uint8_t* const key = pattern->key;
  const uint64_t key_length = pattern->key_length;
  // Calculate PEQ dimensions
  const uint64_t pattern_PEQ_num_entries = DIV_CEIL(key_length,BPM_GPU_PATTERN_ENTRY_LENGTH);
  const uint64_t pattern_PEQ_length = pattern_PEQ_num_entries*BPM_GPU_PATTERN_ENTRY_LENGTH;
  // Insert query metadata
  bpm_gpu_buffer->pattern_id = bpm_gpu_buffer->num_queries;
  (bpm_gpu_buffer->num_queries)++;
  bpm_gpu_qry_info_t* const query_info =
      bpm_gpu_buffer_get_peq_info_(bpm_gpu_buffer->buffer) + bpm_gpu_buffer->pattern_id;
  query_info->posEntry = bpm_gpu_buffer->num_PEQ_entries;
  query_info->size = key_length;
  // Insert query pattern
  bpm_gpu_qry_entry_t* const query_pattern =
      bpm_gpu_buffer_get_peq_entries_(bpm_gpu_buffer->buffer) + bpm_gpu_buffer->num_PEQ_entries;
  bpm_gpu_buffer->num_PEQ_entries += pattern_PEQ_num_entries;
  // Clear PEQ pattern
  uint64_t entry, subentry, alpha;
  for (entry=0;entry<pattern_PEQ_num_entries;++entry) {
    for (alpha=0;alpha<BPM_GPU_PATTERN_ALPHABET_LENGTH;++alpha) {
      for (subentry=0;subentry<BPM_GPU_PATTERN_NUM_SUB_ENTRIES;++subentry) {
        query_pattern[entry].bitmap[alpha][subentry] = 0;
      }
    }
  }
  // Compile PEQ pattern
  uint64_t i, offset;
  for (i=0,offset=0,entry=0,subentry=0;i<key_length;++i,++offset) {
    // Update location
    if (offset==BPM_GPU_PATTERN_SUBENTRY_LENGTH) {
      offset = 0; ++subentry;
      if (subentry==BPM_GPU_PATTERN_NUM_SUB_ENTRIES) {
        subentry = 0; ++entry;
      }
    }
    // Encode pattern char
    const uint8_t enc_char = key[i];
    const uint32_t mask = (UINT32_ONE_MASK<<offset);
    query_pattern[entry].bitmap[enc_char][subentry] |= mask;
  }
  for (;i<pattern_PEQ_length;++i) {
    // Update location
    if (offset==BPM_GPU_PATTERN_SUBENTRY_LENGTH) {
      offset = 0; ++subentry;
      if (subentry==BPM_GPU_PATTERN_NUM_SUB_ENTRIES) {
        subentry = 0; ++entry;
      }
    }
    // Encode remaining chars
    uint8_t enc_char;
    const uint32_t mask = (UINT32_ONE_MASK<<offset);
    for (enc_char=0;enc_char<DNA__N_RANGE;++enc_char) {
      query_pattern[entry].bitmap[enc_char][subentry] |= mask;
    }
  }
}
/*
 * Send/Receive Buffer
 */
GEM_INLINE void bpm_gpu_buffer_send(bpm_gpu_buffer_t* const bpm_gpu_buffer) {
  PROF_START_TIMER(GP_BPM_GPU_BUFFER_SEND);
  PROF_BLOCK() {
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
  }
	bpm_gpu_send_buffer_(bpm_gpu_buffer->buffer,
      bpm_gpu_buffer->num_PEQ_entries,bpm_gpu_buffer->num_queries,bpm_gpu_buffer->num_candidates);
	PROF_STOP_TIMER(GP_BPM_GPU_BUFFER_SEND);
}
GEM_INLINE void bpm_gpu_buffer_receive(bpm_gpu_buffer_t* const bpm_gpu_buffer) {
	bpm_gpu_receive_buffer_(bpm_gpu_buffer->buffer);
  PROF_BLOCK() {
    TIMER_STOP(&bpm_gpu_buffer->timer);
    PROF_ADD_COUNTER(GP_BPM_GPU_BUFFER_CHECK_TIME,bpm_gpu_buffer->timer.accumulated);
  }
}
# else /* BPM_GPU_PATTERN_DEBUG */
GEM_INLINE void bpm_gpu_buffer_put_pattern(
    bpm_gpu_buffer_t* const bpm_gpu_buffer,pattern_t* const pattern) {
  bpm_pattern_t* const bpm_pattern = &pattern->bpm_pattern;
  // Calculate PEQ dimensions
  const uint64_t pattern_PEQ_length = bpm_pattern->PEQ_length;
  const uint64_t pattern_PEQ_num_entries = DIV_CEIL(pattern_PEQ_length,BPM_GPU_PATTERN_ENTRY_LENGTH);
  // Insert query metadata
  bpm_gpu_buffer->pattern_id = bpm_gpu_buffer->num_queries;
  (bpm_gpu_buffer->num_queries)++;
  // Insert query pattern
  bpm_gpu_buffer->num_PEQ_entries += pattern_PEQ_num_entries;
  // [DEBUG] Insert pointer to pattern
  pattern_t** const query_pattern =
      (pattern_t**)bpm_gpu_buffer_get_peq_entries_(bpm_gpu_buffer->buffer) + bpm_gpu_buffer->pattern_id;
  *query_pattern = pattern;
}
/*
 * Send/Receive Buffer
 */
GEM_INLINE void bpm_gpu_buffer_send(bpm_gpu_buffer_t* const bpm_gpu_buffer) {
  // [DEBUG] Solve all queries using CPU BPM-align
  PROF_START_TIMER(GP_BPM_GPU_BUFFER_SEND);
  PROF_START_TIMER(GP_BPM_GPU_BUFFER_CHECK_TIME);
  PROF_BLOCK() {
    const uint64_t max_candidates = bpm_gpu_buffer_get_max_candidates_(bpm_gpu_buffer->buffer);
    const uint64_t max_queries = bpm_gpu_buffer_get_max_queries_(bpm_gpu_buffer->buffer);
    const uint64_t max_peq_entries = bpm_gpu_buffer_get_max_peq_entries_(bpm_gpu_buffer->buffer);
    const uint64_t used_candidates = bpm_gpu_buffer->num_candidates;
    const uint64_t used_queries = bpm_gpu_buffer->num_queries;
    const uint64_t used_peq_entries = bpm_gpu_buffer->num_PEQ_entries;
    PROF_ADD_COUNTER(GP_BPM_GPU_BUFFER_USAGE_CANDIDATES,(100*used_candidates)/max_candidates);
    PROF_ADD_COUNTER(GP_BPM_GPU_BUFFER_USAGE_QUERIES,(100*used_queries)/max_queries);
    PROF_ADD_COUNTER(GP_BPM_GPU_BUFFER_USAGE_PEQ_ENTRIES,(100*used_peq_entries)/max_peq_entries);
  }
  const uint64_t num_candidates = bpm_gpu_buffer->num_candidates;
  bpm_gpu_cand_info_t* query_candidate = bpm_gpu_buffer_get_candidates_(bpm_gpu_buffer->buffer);
  bpm_gpu_res_entry_t* query_result = bpm_gpu_buffer_get_results_(bpm_gpu_buffer->buffer);
  const uint64_t num_queries = bpm_gpu_buffer->num_queries;
  uint64_t query_pos, candidate_pos;
  for (query_pos=0,candidate_pos=0;query_pos<num_queries;++query_pos) {
    // Get pattern
    pattern_t* const query_pattern = *((pattern_t**)bpm_gpu_buffer_get_peq_entries_(bpm_gpu_buffer->buffer) + query_pos);
#ifdef BPM_GPU_GENERATE_CANDIDATES_PROFILE
    uint64_t i;
    for (i=0;i<query_pattern->key_length;++i) {
      fprintf(stdout,"%c",dna_decode(query_pattern->key[i]));
    }
    fprintf(stdout,"\t");
#endif
    // Traverse all candidates
    while (candidate_pos < num_candidates && query_candidate->query==query_pos) {
      // Run BPM
      const uint8_t* const sequence = dna_text_get_buffer(bpm_gpu_buffer->enc_text) + query_candidate->position;
      uint64_t position, distance;

#ifdef BPM_GPU_GENERATE_CANDIDATES_PROFILE
    bpm_get_distance(&query_pattern->bpm_pattern,
      sequence,query_candidate->size,&position,&distance);
    fprintf(stderr,"B p=%lu e=%lu c=%lu\n",query_candidate->position,distance,position);
    fprintf(stdout,"\tchr1:+:%lu:%ld",query_candidate->position,distance);
#else
    bpm_get_distance__cutoff(&query_pattern->bpm_pattern,
        sequence,query_candidate->size,&position,&distance,query_pattern->max_effective_filtering_error);
#endif

      // Copy results
      query_result->column = position;
      query_result->score = distance;
      // Next
      ++candidate_pos;
      ++query_candidate;
      ++query_result;
    }
#ifdef BPM_GPU_GENERATE_CANDIDATES_PROFILE
    fprintf(stdout,"\n");
#endif
  }
  PROF_STOP_TIMER(GP_BPM_GPU_BUFFER_CHECK_TIME);
  PROF_STOP_TIMER(GP_BPM_GPU_BUFFER_SEND);
}
GEM_INLINE void bpm_gpu_buffer_receive(bpm_gpu_buffer_t* const bpm_gpu_buffer) { /* [DEBUG] Do nothing */ }
# endif /* BPM_GPU_PATTERN_DEBUG */
#endif /* HAVE_CUDA */

