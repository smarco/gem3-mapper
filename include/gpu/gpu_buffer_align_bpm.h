/*
 * PROJECT: GEMMapper
 * FILE: gpu_buffer_align_bpm.h
 * DATE: 06/06/2012
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 *            Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#ifndef GPU_BUFFER_ALIGN_BPM_H_
#define GPU_BUFFER_ALIGN_BPM_H_

#include "utils/essentials.h"
#include "system/profiler_timer.h"
#include "data_structures/pattern.h"
#include "gpu/gpu_buffer_collection.h"

/*
 * GPU BMP Buffer
 */
typedef struct {
  /* GPU Generic Buffer*/
  void* buffer;                       // GPU Generic Buffer
  /* Dimensions Hints */
  gem_counter_t query_length;         // Tracks queries' length
  gem_counter_t candidates_per_tile;  // Tracks number of candidates per tile
  uint32_t query_same_length;         // Tracks same-read-length buffers
  /* Buffer state */
  uint32_t current_query_offset;
  uint32_t num_entries;
  uint32_t num_queries;
  uint32_t num_candidates;
  /* CPU Computation */
  bool compute_cpu;                   // Computing Using CPU (disable GPU)
  archive_text_t* archive_text;
  text_collection_t* text_collection;
  mm_stack_t* mm_stack;
  /* Profile */
  gem_timer_t timer;
} gpu_buffer_align_bpm_t;

/*
 * Pattern Setup
 */
void gpu_bpm_pattern_compile(
    bpm_pattern_t* const bpm_pattern,
    const uint64_t words128_per_tile,
    const uint64_t max_error);
uint64_t gpu_bpm_pattern_get_entry_length();

/*
 * Setup
 */
gpu_buffer_align_bpm_t* gpu_buffer_align_bpm_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_no,
    archive_text_t* const archive_text,
    text_collection_t* const text_collection,
    mm_stack_t* const mm_stack);
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
uint64_t gpu_buffer_align_bpm_get_entry_length();
uint64_t gpu_buffer_align_bpm_get_max_candidates(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm);
uint64_t gpu_buffer_align_bpm_get_max_queries(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm);

uint64_t gpu_buffer_align_bpm_get_num_candidates(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm);
uint64_t gpu_buffer_align_bpm_get_num_queries(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm);

void gpu_buffer_align_bpm_compute_dimensions(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    bpm_pattern_t* const bpm_pattern,
    bpm_pattern_t* const bpm_pattern_tiles,
    const uint64_t num_candidates,
    uint64_t* const total_entries,
    uint64_t* const total_queries,
    uint64_t* const total_candidates);
bool gpu_buffer_align_bpm_fits_in_buffer(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    const uint64_t total_entries,
    const uint64_t total_queries,
    const uint64_t total_candidates);

/*
 * Accessors
 */
void gpu_buffer_align_bpm_add_pattern(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    bpm_pattern_t* const bpm_pattern,
    bpm_pattern_t* const bpm_pattern_tiles);
void gpu_buffer_align_bpm_add_candidate(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    const uint64_t tile_offset,
    const uint64_t candidate_text_position,
    const uint64_t candidate_length);

void gpu_buffer_align_bpm_get_candidate(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    const uint64_t candidate_offset,
    uint64_t* const candidate_text_position,
    uint32_t* const candidate_length);
void gpu_buffer_align_bpm_get_result(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    const uint64_t candidate_offset,
    uint32_t* const levenshtein_distance,
    uint32_t* const levenshtein_match_pos);
void gpu_buffer_align_bpm_retrieve_pattern(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    const uint64_t candidate_offset,
    bpm_pattern_t* const bpm_pattern,
    mm_stack_t* const mm_stack);

/*
 * Hints
 */
void gpu_buffer_align_bpm_record_query_length(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    const uint64_t query_length);
void gpu_buffer_align_bpm_record_candidates_per_tile(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    const uint64_t num_candidates);
uint64_t gpu_buffer_align_bpm_get_mean_query_length(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm);
uint64_t gpu_buffer_align_bpm_get_mean_candidates_per_tile(
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm);

/*
 * Send/Receive
 */
void gpu_buffer_align_bpm_send(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm);
void gpu_buffer_align_bpm_receive(gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm);

#endif /* GPU_BUFFER_ALIGN_BPM_H_ */
