/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_KMER_CORE_H_
#define GPU_KMER_CORE_H_

extern "C" {
#include "gpu_commons.h"
#include "gpu_reference.h"
#include "gpu_buffer.h"
}
#include "gpu_resources.h"

// Constants
#define GPU_KMER_FILTER_COUNTING_MASK_3                      0x0000003Fu
#define GPU_KMER_FILTER_COUNTING_MASK_4                      0x000000FFu
#define GPU_KMER_FILTER_COUNTING_MASK_5                      0x000003FFu
#define GPU_KMER_FILTER_COUNTING_MASK_6                      0x00000FFFu
#define GPU_KMER_FILTER_COUNTING_MASK_7                      0x00003FFFu

#define GPU_KMER_FILTER_COUNTING_LENGTH                      6
#define GPU_KMER_FILTER_COUNTING_MASK                        GPU_KMER_FILTER_COUNTING_MASK_6

#define GPU_KMER_FILTER_BASE_QUERY_LENGTH                    8
#define GPU_KMER_FILTER_BASE_QUERY_MASK                      (~(GPU_UINT64_ONES << GPU_KMER_FILTER_BASE_QUERY_LENGTH))
#define GPU_KMER_FILTER_BASES_PER_QUERY_ENTRY                (GPU_UINT64_LENGTH / GPU_KMER_FILTER_BASE_QUERY_LENGTH)

#define GPU_KMER_FILTER_BASE_CANDIDATE_LENGTH                2
#define GPU_KMER_FILTER_BASE_CANDIDATE_MASK                  (~(GPU_UINT64_ONES << GPU_KMER_FILTER_BASE_CANDIDATE_LENGTH))
#define GPU_KMER_FILTER_BASES_PER_CANDIDATE_ENTRY            (GPU_UINT64_LENGTH / GPU_KMER_FILTER_BASE_CANDIDATE_LENGTH)

#define GPU_ALIGN_DISTANCE_ZERO                              GPU_UINT32_MASK_ONE_HIGH
#define GPU_ALIGN_DISTANCE_INF                               GPU_UINT32_ONES

#define GPU_KMER_FILTER_COUNTING_NUM_KMERS                   GPU_POW4(GPU_KMER_FILTER_COUNTING_LENGTH)
#define GPU_KMER_FILTER_COUNTING_ADD_INDEX(kmerIdx, encBase) kmerIdx = (((kmerIdx << GPU_REFERENCE_CHAR_LENGTH) | (encBase & GPU_REFERENCE_UINT32_MASK_BASE)) & GPU_KMER_FILTER_COUNTING_MASK)

#endif /* GPU_KMER_CORE_H_ */

