/*
 * PROJECT: GEMMapper
 * FILE: profiler_cuda.h
 * DATE: 06/06/2012
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 *            Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: CUDA profile module
 */
#ifndef PROFILE_CUDA_H_
#define PROFILE_CUDA_H_

/*
 * Include
 */
#include "system/commons.h"

/*
 * Tags Colors
 */
#define PROFILER_CUDA_TAGS_NUM_COLORS 7
extern const uint32_t profiler_cuda_tags_colors[];

/*
 * Profile Start/Stop
 */
void PROFILE_CUDA_START(char* const name,const uint64_t cid);
void PROFILE_CUDA_STOP();

#endif /* PROFILE_CUDA_H_ */
