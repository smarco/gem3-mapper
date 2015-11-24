/*
 * PROJECT: GEMMapper
 * FILE: gpu_buffer_align_bpm.h
 * DATE: 06/06/2012
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 *            Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#ifndef GPU_BUFFER_ALIGN_BPM_H_
#define GPU_BUFFER_ALIGN_BPM_H_

#include "essentials.h"
#include "profiler_timer.h"
#include "pattern.h"
#include "gpu_buffer_collection.h"

/*
 * GPU BMP Buffer
 */
typedef struct {
  /* GPU Generic Buffer*/
  void* buffer;             // GPU Generic Buffer
  /* Dimensions Hints */
  uint32_t averageNumPEQEntries;
  uint32_t candidatesPerQuery;
  uint32_t candidates_same_length; // TODO
  /* Buffer state */
  uint32_t pattern_id;      // Pattern ID (Generator)
  uint32_t num_PEQ_entries;
  uint32_t num_queries;
  uint32_t num_candidates;
  /* CPU Computation */
  bool compute_cpu;         // Computing Using CPU (disable GPU)
  /* Profile */
  gem_timer_t timer;
} gpu_buffer_align_bpm_t;

/*
 * Pattern Setup
 */
void gpu_bpm_pattern_compile(bpm_pattern_t* const bpm_pattern,const uint64_t max_error);
uint64_t gpu_bpm_pattern_get_entry_length();

/*
 * Setup
 */
gpu_buffer_align_bpm_t* gpu_buffer_align_bpm_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,const uint64_t buffer_no);
void gpu_buffer_align_bpm_clear(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm);
void gpu_buffer_align_bpm_delete(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm);

/*
 * Computing Device
 */
void gpu_buffer_align_bpm_set_device_cpu(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm);
void gpu_buffer_align_bpm_set_device_gpu(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm);

/*
 * Occupancy & Limits
 */
uint64_t gpu_buffer_align_bpm_get_max_candidates(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm);
uint64_t gpu_buffer_align_bpm_get_max_queries(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm);

uint64_t gpu_buffer_align_bpm_get_num_candidates(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm);
uint64_t gpu_buffer_align_bpm_get_num_queries(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm);

void gpu_buffer_align_bpm_compute_dimensions(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,const pattern_t* const pattern,
    const uint64_t total_candidates,uint64_t* const total_entries,
    uint64_t* const total_query_chunks,uint64_t* const total_candidate_chunks);
bool gpu_buffer_align_bpm_fits_in_buffer(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,const uint64_t total_entries,
    const uint64_t total_query_chunks,const uint64_t total_candidate_chunks);

/*
 * Accessors
 */
void gpu_buffer_align_bpm_add_pattern(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,pattern_t* const pattern);
void gpu_buffer_align_bpm_add_candidate(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,const uint64_t candidate_text_position,
    const uint64_t candidate_length,const uint64_t pattern_chunk);
void gpu_buffer_align_bpm_get_candidate(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,const uint64_t buffer_pos,
    uint64_t* const candidate_text_position,uint32_t* const candidate_length);
void gpu_buffer_align_bpm_get_result(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,const uint64_t buffer_pos,
    uint32_t* const levenshtein_distance,uint32_t* const levenshtein_match_pos);

/*
 * Send/Receive
 */
void gpu_buffer_align_bpm_send(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm);
void gpu_buffer_align_bpm_receive(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm);

/*
 * Errors
 */
#define GEM_ERROR_GPU_BPM_MAX_PATTERN_LENGTH "BPM-GPU. Query pattern (%"PRIu64" entries) exceeds maximum buffer capacity (%"PRIu64" entries)"
#define GEM_ERROR_GPU_BPM_MAX_CANDIDATES "BPM-GPU. Number of candidates (%"PRIu64") exceeds maximum buffer capacity (%"PRIu64" candidates)"

#endif /* GPU_BUFFER_ALIGN_BPM_H_ */
