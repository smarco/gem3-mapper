/*
 * PROJECT: GEMMapper
 * FILE: mapper_cuda.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MAPPER_CUDA_H_
#define MAPPER_CUDA_H_

#include "mapper.h"
#include "bpm_align_gpu.h"
#include "archive_search_group.h"

/*
 * Mapper-CUDA Parameters
 */
typedef struct {
  /* BPM Buffers */
  uint64_t num_search_groups;      // Total number of search-groups deployed
  uint64_t average_query_size;     // HINT on the average query size
  uint64_t candidates_per_query;   // HINT on the number of candidates per query
  /* System */
  uint64_t num_generating_threads; // Total number of threads generating candidates
  uint64_t num_selecting_threads;  // Total number of threads selecting candidates
} mapper_cuda_parameters_t;

/*
 * CUDA Mapper parameters
 */
GEM_INLINE void mapper_cuda_parameters_set_defaults(mapper_cuda_parameters_t* const mapper_cuda_parameters);

/*
 * SE-CUDA Mapper
 */
GEM_INLINE void mapper_SE_CUDA_run(
    mapper_parameters_t* const mapper_parameters,mapper_cuda_parameters_t* const cuda_parameters);

/*
 * Error Messages
 */
#define GEM_ERROR_MAPPER_CUDA_ERROR_PARSING "Error parsing sequence"

#endif /* MAPPER_CUDA_H_ */
