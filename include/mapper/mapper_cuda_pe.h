/*
 * PROJECT: GEMMapper
 * FILE: mapper_cuda_pe.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MAPPER_CUDA_PE_H_
#define MAPPER_CUDA_PE_H_

#include "mapper/mapper_cuda.h"

/*
 * PE CUDA Mapper
 */
void* mapper_cuda_pe_thread(mapper_cuda_search_t* const mapper_search);

#endif /* MAPPER_CUDA_PE_H_ */
