/*
 * PROJECT: GEMMapper
 * FILE: mapper_cuda_se.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MAPPER_CUDA_SE_H_
#define MAPPER_CUDA_SE_H_

#include "mapper/mapper_cuda.h"

/*
 * Mapper SE-CUDA
 */
void* mapper_cuda_se_thread(mapper_cuda_search_t* const mapper_search);

#endif /* MAPPER_CUDA_SE_H_ */
