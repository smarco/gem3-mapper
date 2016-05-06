/*
 * PROJECT: GEMMapper
 * FILE: mapper_cuda.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MAPPER_CUDA_H_
#define MAPPER_CUDA_H_

#include "mapper/mapper.h"
#include "search_pipeline/search_pipeline.h"
#include "stats/report_stats_mstats.h"

/*
 * Define to use CPU adaptive-region-profile & disable GPU static-region-profile
 */
#define MAPPER_CUDA_ADAPTIVE_REGION_PROFILE

/*
 * Mapper-CUDA Search
 */
typedef struct {
  /* Thread Info */
  uint64_t thread_id;
  pthread_t* thread_data;
  /* Parameters */
  mapper_parameters_t* mapper_parameters;
  /* I/O */
  buffered_input_file_t* buffered_fasta_input_end1;
  buffered_input_file_t* buffered_fasta_input_end2;
  buffered_output_file_t* buffered_output_file;
  /* GPU Buffers */
  gpu_buffer_collection_t* gpu_buffer_collection;
  uint64_t gpu_buffers_offset;
  /* Search Pipeline (and auxiliary variables) */
  search_pipeline_t* search_pipeline;
  archive_search_t* pending_search_region_profile_end1;
  archive_search_t* pending_search_region_profile_end2;
  archive_search_t* pending_search_decode_candidates_end1;
  archive_search_t* pending_search_decode_candidates_end2;
  archive_search_t* pending_search_verify_candidates_end1;
  archive_search_t* pending_search_verify_candidates_end2;
  /* Stats */
  mapping_stats_t* mapping_stats; // Per thread stats report structures
  /* Progress */
  ticker_t* ticker;
  uint64_t reads_processed;
} mapper_cuda_search_t;

/*
 * Setup
 */
void mapper_cuda_search_init(mapper_cuda_search_t* const mapper_cuda_search);

/*
 * SE-CUDA Mapper
 */
void mapper_cuda_se_run(mapper_parameters_t* const mapper_parameters);

/*
 * PE-CUDA Mapper
 */
void mapper_cuda_pe_run(mapper_parameters_t* const mapper_parameters);

/*
 * Error Messages
 */
#define GEM_ERROR_MAPPER_CUDA_ERROR_PARSING "Mapper-CUDA. Error parsing sequence"

#endif /* MAPPER_CUDA_H_ */
