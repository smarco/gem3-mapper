/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *  Copyright (c) 2013-2017 by Alejandro Chacon <alejandro.chacond@gmail.com>
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
 *            Alejandro Chacon <alejandro.chacond@gmail.com>
 * DESCRIPTION:
 *   GPU-adaptor module provides support functions to configure the GPU workflow
 */

#ifndef GPU_CONFIG_H_
#define GPU_CONFIG_H_

#include "utils/essentials.h"

#ifdef HAVE_CUDA
#include "resources/gem-cutter/gpu_interface.h"
#else
#define GPU_BPM_ALIGN_PEQ_ALPHABET_SIZE     5
#define GPU_BPM_ALIGN_PEQ_ENTRY_LENGTH      128
#define GPU_BPM_ALIGN_PEQ_SUBENTRY_LENGTH   32
#define GPU_BPM_ALIGN_PEQ_SUBENTRIES        (GPU_BPM_ALIGN_PEQ_ENTRY_LENGTH / UINT32_LENGTH)

#define GPU_BPM_FILTER_PEQ_ALPHABET_SIZE     5
#define GPU_BPM_FILTER_PEQ_ENTRY_LENGTH      128
#define GPU_BPM_FILTER_PEQ_SUBENTRY_LENGTH   32
#define GPU_BPM_FILTER_PEQ_SUBENTRIES        (GPU_BPM_FILTER_PEQ_ENTRY_LENGTH / UINT32_LENGTH)
#endif

/*
 * Region-Profile Adaptive (if not define, static seeds are generated)
 */
#define GPU_REGION_PROFILE_ADAPTIVE

/*
 * Benchmark Generation
 */
//#define CUDA_BENCHMARK_GENERATE_REGION_PROFILE
//#define CUDA_BENCHMARK_GENERATE_DECODE_CANDIDATES

/*
 * Debug GPU results (from GPU-Buffers at each stage)
 */
//#define GPU_CHECK_REGION_PROFILE
//#define GPU_CHECK_DECODE_POSITIONS
//#define GPU_CHECK_KMER_FILTER
//#define GPU_CHECK_BPM_DISTANCE
//#define GPU_CHECK_BPM_ALIGN

/*
 * GPU Support
 */
bool gpu_supported(void);

#endif /* GPU_CONFIG_H_ */
