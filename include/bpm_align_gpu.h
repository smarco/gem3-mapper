/*
 * PROJECT: GEMMapper
 * FILE: bpm_align_gpu.h
 * DATE: 06/06/2012
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 *            Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#ifndef BPM_ALIGN_GPU_H_
#define BPM_ALIGN_GPU_H_

#include "essentials.h"
#include "archive.h"
#include "pattern.h"
#include "bpm_align.h"

/*
 * BMP-GPU Buffer & Collection
 */
typedef struct {
  /* BMP-GPU Buffer*/
  void* buffer;
  /* Buffer state */
  uint32_t num_PEQ_entries;
  uint32_t num_queries;
  uint32_t num_candidates;
  /* Pattern state */
  uint32_t pattern_id;
  /* Misc */
  dna_text_t* enc_text; /* BPM_GPU_PATTERN_DEBUG */
  gem_timer_t timer;    /* !GEM_NOPROFILE */
} bpm_gpu_buffer_t;
typedef struct {
  void** internal_buffers;            // Internal Buffers
  bpm_gpu_buffer_t* bpm_gpu_buffers;  // Wrapped Buffers
  uint64_t num_buffers;               // Total number of buffers allocated
} bpm_gpu_buffer_collection_t;

/*
 * BPM_GPU Setup
 */
GEM_INLINE bpm_gpu_buffer_collection_t* bpm_gpu_init(
    dna_text_t* const enc_text,const uint32_t num_buffers,const uint32_t buffer_size,
    const int32_t average_query_size,const int32_t candidates_per_query,const bool verbose);
GEM_INLINE void bpm_gpu_destroy(bpm_gpu_buffer_collection_t* const buffer_collection);
GEM_INLINE bool bpm_gpu_support();

/*
 * Buffer Accessors
 */
GEM_INLINE void bpm_gpu_buffer_clear(bpm_gpu_buffer_t* const bpm_gpu_buffer);

GEM_INLINE uint64_t bpm_gpu_buffer_get_max_candidates(bpm_gpu_buffer_t* const bpm_gpu_buffer);
GEM_INLINE uint64_t bpm_gpu_buffer_get_max_queries(bpm_gpu_buffer_t* const bpm_gpu_buffer);

GEM_INLINE uint64_t bpm_gpu_buffer_get_num_candidates(bpm_gpu_buffer_t* const bpm_gpu_buffer);
GEM_INLINE uint64_t bpm_gpu_buffer_get_num_queries(bpm_gpu_buffer_t* const bpm_gpu_buffer);

GEM_INLINE bool bpm_gpu_buffer_fits_in_buffer(
    bpm_gpu_buffer_t* const bpm_gpu_buffer,const pattern_t* const pattern,
    const uint64_t num_pattern_dups,const uint64_t total_candidates);

GEM_INLINE void bpm_gpu_buffer_put_pattern(
    bpm_gpu_buffer_t* const bpm_gpu_buffer,pattern_t* const pattern);
GEM_INLINE void bpm_gpu_buffer_put_candidate(
    bpm_gpu_buffer_t* const bpm_gpu_buffer,const uint64_t candidate_text_position,
    const uint64_t candidate_length,const uint64_t pattern_chunk);
GEM_INLINE void bpm_gpu_buffer_get_candidate(
    bpm_gpu_buffer_t* const bpm_gpu_buffer,const uint64_t position,
    uint64_t* const candidate_text_position,uint32_t* const candidate_length);
GEM_INLINE void bpm_gpu_buffer_get_candidate_result(
    bpm_gpu_buffer_t* const bpm_gpu_buffer,const uint64_t position,
    uint32_t* const levenshtein_distance,uint32_t* const levenshtein_match_pos);

/*
 * Init the local thread BMP Buffers
 */
GEM_INLINE void bpm_gpu_init_buffer(bpm_gpu_buffer_t* const bpm_gpu_buffer);

/*
 * Send/Receive Buffer
 */
GEM_INLINE void bpm_gpu_buffer_send(bpm_gpu_buffer_t* const bpm_gpu_buffer);
GEM_INLINE void bpm_gpu_buffer_receive(bpm_gpu_buffer_t* const bpm_gpu_buffer);

/*
 * Errors
 */
#define GEM_ERROR_BPM_GPU_MAX_PATTERN_LENGTH "BPM-GPU. Query pattern (%lu entries) exceeds maximum buffer capacity (%lu entries)"
#define GEM_ERROR_BPM_GPU_MAX_CANDIDATES "BPM-GPU. Number of candidates (%lu) exceeds maximum buffer capacity (%lu candidates)"

#endif /* BPM_ALIGN_GPU_H_ */
