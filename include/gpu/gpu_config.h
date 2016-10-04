/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *  Copyright (c) 2011-2017 by Alejandro Chacon <alejandro.chacon@uab.es>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   GPU-adaptor module provides support functions to configure the GPU workflow
 */

#ifndef GPU_CONFIG_H_
#define GPU_CONFIG_H_

#include "utils/essentials.h"

/*
 * Region-Profile Adaptive (if not define, static seeds are generated)
 */
#define GPU_REGION_PROFILE_ADAPTIVE

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
bool gpu_supported(void);

#endif /* GPU_CONFIG_H_ */
