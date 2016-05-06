/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *                2013-2016 by Santiago Marco-Sola <santiagomsola@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_BPM_INTERFACE_H_
#define GPU_BPM_INTERFACE_H_

/*
 * Constants
 */
#define GPU_BPM_PEQ_ALPHABET_SIZE     5
#define GPU_BPM_PEQ_ENTRY_LENGTH      128
#define GPU_BPM_PEQ_SUBENTRY_LENGTH   32
#define GPU_BPM_PEQ_SUBENTRIES        (GPU_BPM_PEQ_ENTRY_LENGTH / GPU_UINT32_LENGTH)

/*
 * Enum types for Device & Host
 */
typedef enum
{
  GPU_REF_NONE,
  GPU_REF_MFASTA_FILE,
  GPU_REF_PROFILE_FILE,
  GPU_REF_ASCII,
  GPU_REF_GEM_FILE,
  GPU_REF_GEM_FULL,
  GPU_REF_GEM_ONLY_FORWARD
} gpu_ref_coding_t;


/*
 * Common types for Device & Host
 */
typedef struct { /* each row 1 PEQ Entry (128bits) */
  uint32_t bitmap[GPU_BPM_PEQ_ALPHABET_SIZE][GPU_BPM_PEQ_SUBENTRIES];
} gpu_bpm_qry_entry_t;

typedef struct {
  uint32_t column;
  uint32_t score;
} gpu_bpm_alg_entry_t;

typedef struct {
  uint64_t position;
  uint32_t query;
  uint32_t size;
} gpu_bpm_cand_info_t;

typedef struct {
  uint32_t posEntry;
  uint32_t size;
} gpu_bpm_qry_info_t;

typedef struct {
  char*            reference;   // GEM reference using 8 bases
  gpu_ref_coding_t ref_coding;  // GEM reference coding (F or FR)
  uint64_t         ref_length;  // Length of the reference
} gpu_gem_ref_dto_t;

typedef struct {
  char              *reference;
  gpu_ref_coding_t  refCoding;
  uint64_t          refSize;
} gpu_reference_dto_t;

/*
 * Obtain Buffers
 */
gpu_bpm_qry_entry_t* gpu_bpm_buffer_get_peq_entries_(const void* const bpmBuffer);
gpu_bpm_cand_info_t* gpu_bpm_buffer_get_candidates_(const void* const bpmBuffer);
gpu_bpm_qry_info_t*  gpu_bpm_buffer_get_peq_info_(const void* const bpmBuffer);
gpu_bpm_alg_entry_t* gpu_bpm_buffer_get_alignments_(const void* const bpmBuffer);


/*
 * Get elements
 */
uint32_t gpu_bpm_buffer_get_max_peq_entries_(const void* const bpmBuffer);
uint32_t gpu_bpm_buffer_get_max_candidates_(const void* const bpmBuffer);
uint32_t gpu_bpm_buffer_get_max_queries_(const void* const bpmBuffer);

/*
 * Main functions
 */
void gpu_bpm_init_buffer_(void* const bpmBuffer, const uint32_t averageQuerySize, const uint32_t candidatesPerQuery);
void gpu_bpm_send_buffer_(void* const bpmBuffer, const uint32_t numPEQEntries, const uint32_t numQueries, const uint32_t numCandidates, const uint32_t sizeCandidates);
void gpu_bpm_receive_buffer_(void* const bpmBuffer);
void gpu_bpm_init_and_realloc_buffer_(void *bpmBuffer, const uint32_t totalPEQEntries, const uint32_t totalCandidates, const uint32_t totalQueries);

#endif /* GPU_BPM_INTERFACE_H_ */
