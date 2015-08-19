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
bpm_gpu_buffer_collection_t* bpm_gpu_init(
    archive_text_t* const archive_text,const uint32_t num_buffers,const uint32_t buffer_size,
    const int32_t average_query_size,const int32_t candidates_per_query,const bool verbose);
void bpm_gpu_destroy(bpm_gpu_buffer_collection_t* const buffer_collection);
bool bpm_gpu_support();

/*
 * Buffer Accessors
 */
void bpm_gpu_buffer_clear(bpm_gpu_buffer_t* const bpm_gpu_buffer);

uint64_t bpm_gpu_buffer_get_max_candidates(bpm_gpu_buffer_t* const bpm_gpu_buffer);
uint64_t bpm_gpu_buffer_get_max_queries(bpm_gpu_buffer_t* const bpm_gpu_buffer);

uint64_t bpm_gpu_buffer_get_num_candidates(bpm_gpu_buffer_t* const bpm_gpu_buffer);
uint64_t bpm_gpu_buffer_get_num_queries(bpm_gpu_buffer_t* const bpm_gpu_buffer);

void bpm_gpu_buffer_compute_dimensions(
    bpm_gpu_buffer_t* const bpm_gpu_buffer,const pattern_t* const pattern,
    const uint64_t total_candidates,uint64_t* const total_entries,
    uint64_t* const total_query_chunks,uint64_t* const total_candidate_chunks);
bool bpm_gpu_buffer_fits_in_buffer(
    bpm_gpu_buffer_t* const bpm_gpu_buffer,const uint64_t total_entries,
    const uint64_t total_query_chunks,const uint64_t total_candidate_chunks);

void bpm_gpu_buffer_put_pattern(
    bpm_gpu_buffer_t* const bpm_gpu_buffer,pattern_t* const pattern);
void bpm_gpu_buffer_put_candidate(
    bpm_gpu_buffer_t* const bpm_gpu_buffer,const uint64_t candidate_text_position,
    const uint64_t candidate_length,const uint64_t pattern_chunk);
void bpm_gpu_buffer_get_candidate(
    bpm_gpu_buffer_t* const bpm_gpu_buffer,const uint64_t position,
    uint64_t* const candidate_text_position,uint32_t* const candidate_length);
void bpm_gpu_buffer_get_candidate_result(
    bpm_gpu_buffer_t* const bpm_gpu_buffer,const uint64_t position,
    uint32_t* const levenshtein_distance,uint32_t* const levenshtein_match_pos);

/*
 * Init the local thread BMP Buffers
 */
void bpm_gpu_init_buffer(bpm_gpu_buffer_t* const bpm_gpu_buffer);

/*
 * Send/Receive Buffer
 */
void bpm_gpu_buffer_send(bpm_gpu_buffer_t* const bpm_gpu_buffer);
void bpm_gpu_buffer_receive(bpm_gpu_buffer_t* const bpm_gpu_buffer);

/*
 * Errors
 */
#define GEM_ERROR_BPM_GPU_MAX_PATTERN_LENGTH "BPM-GPU. Query pattern (%"PRIu64" entries) exceeds maximum buffer capacity (%"PRIu64" entries)"
#define GEM_ERROR_BPM_GPU_MAX_CANDIDATES "BPM-GPU. Number of candidates (%"PRIu64") exceeds maximum buffer capacity (%"PRIu64" candidates)"

#endif /* BPM_ALIGN_GPU_H_ */
