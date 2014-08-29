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
#include "bpm_align.h"

typedef struct {
  /* BMP-GPU Buffer*/
  void* buffer;
  /* Buffer state */
  uint32_t num_PEQ_entries;
  uint32_t num_queries;
  uint32_t num_candidates;
  /* Pattern ID generator */
  uint32_t pattern_id;
  uint64_t candidates_left;           // Number of candidates left to put in buffer // FIXME maybe remove
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
    const dna_text_t* const enc_text,const uint32_t num_buffers,
    const int32_t average_query_size,const int32_t candidates_per_query);
GEM_INLINE void bpm_gpu_destroy(bpm_gpu_buffer_collection_t* const buffer_collection);

/*
 * Buffer Accessors
 */
GEM_INLINE uint64_t bpm_gpu_buffer_get_num_candidates(bpm_gpu_buffer_t* const bpm_gpu_buffer);
GEM_INLINE uint64_t bpm_gpu_buffer_get_num_queries(bpm_gpu_buffer_t* const bpm_gpu_buffer);

GEM_INLINE bool bpm_gpu_buffer_fits_in_buffer(
    bpm_gpu_buffer_t* const bpm_gpu_buffer,
    const uint64_t num_patterns,const uint64_t total_pattern_length,const uint64_t total_candidates);

GEM_INLINE void bpm_gpu_buffer_put_pattern(
    bpm_gpu_buffer_t* const bpm_gpu_buffer,const bpm_pattern_t* const bpm_pattern);
GEM_INLINE void bpm_gpu_buffer_put_candidate(
    bpm_gpu_buffer_t* const bpm_gpu_buffer,
    const uint64_t candidate_text_position,const uint64_t candidate_length);

/*
 * Errors
 */
#define GEM_ERROR_BPM_GPU_MAX_PATTERN_LENGTH "BPM-GPU. Query pattern (%lu entries) exceeds maximum buffer capacity (%lu entries)"
#define GEM_ERROR_BPM_GPU_MAX_CANDIDATES "BPM-GPU. Number of candidates (%lu) exceeds maximum buffer capacity (%lu candidates)"

#endif /* BPM_ALIGN_GPU_H_ */
