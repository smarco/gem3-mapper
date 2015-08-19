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
 * SE-CUDA Mapper
 */
void mapper_SE_CUDA_run(mapper_parameters_t* const mapper_parameters);

/*
 * PE-CUDA Mapper
 */
void mapper_PE_CUDA_run(mapper_parameters_t* const mapper_parameters);

/*
 * Error Messages
 */
#define GEM_ERROR_MAPPER_CUDA_ERROR_PARSING "Mapper-CUDA. Error parsing sequence"

#endif /* MAPPER_CUDA_H_ */
