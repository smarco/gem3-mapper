/*
 * PROJECT: GEMMapper
 * FILE: gpu_config.h
 * DATE: 06/06/2012
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 *            Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#ifndef GPU_CONFIG_H_
#define GPU_CONFIG_H_

#include "utils/essentials.h"

/*
 * Tiling
 */
#define GPU_WORDS128_PER_TILE 2

/*
 * Benchmark Generation
 */
//#define CUDA_BENCHMARK_GENERATE_REGION_PROFILE
//#define CUDA_BENCHMARK_GENERATE_DECODE_CANDIDATES
//#define CUDA_BENCHMARK_GENERATE_VERIFY_CANDIDATES

/*
 * Debug GPU results (from GPU-Buffers at each stage)
 */
//#define CUDA_CHECK_BUFFERED_REGION_PROFILE
//#define CUDA_CHECK_BUFFERED_DECODE_POSITIONS
//#define CUDA_CHECK_BUFFERED_VERIFY_CANDIDATES

/*
 * GPU Support
 */
bool gpu_supported();

#endif /* GPU_CONFIG_H_ */
