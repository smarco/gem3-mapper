/*
 * PROJECT: GEMMapper
 * FILE: mapper_cuda.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

/*
 * TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
 *
 * 1.- Add Macro call to NOT_SUPPORTED_CUDA
 *
 *
 * TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
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
  /* I/O */
  output_file_type output_file_type;
  /* Single-end Alignment */
//  mapping_mode_t mapping_mode;
  /* Paired-end Alignment */
  /* Reporting */
  /* System */
  uint64_t num_generating_threads; // Total number of threads generating candidates
  uint64_t num_selecting_threads;  // TOtal number of threads selecting candidates
  /* Miscellaneous */
  /* Extras */
} mapper_cuda_parameters_t;

/*
 * CUDA Mapper parameters
 */
GEM_INLINE void mapper_cuda_parameters_set_defaults(mapper_cuda_parameters_t* const mapper_cuda_parameters);

/*
 * SE-CUDA Mapper
 */
GEM_INLINE void mapper_SE_CUDA_run(
    mapper_parameters_t* const mapper_parameters,
    mapper_cuda_parameters_t* const mapper_cuda_parameters);

#endif /* MAPPER_CUDA_H_ */
