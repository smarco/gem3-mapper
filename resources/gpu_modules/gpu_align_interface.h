/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *                2013-2016 by Santiago Marco-Sola <santiagomsola@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_ALIGN_INTERFACE_H_
#define GPU_ALIGN_INTERFACE_H_

/*
 * Constants
 */
#define GPU_BPM_ALIGN_PEQ_ALPHABET_SIZE     5
#define GPU_BPM_ALIGN_PEQ_ENTRY_LENGTH      128
#define GPU_BPM_ALIGN_PEQ_SUBENTRY_LENGTH   32
#define GPU_BPM_ALIGN_PEQ_SUBENTRIES        (GPU_BPM_ALIGN_PEQ_ENTRY_LENGTH / GPU_UINT32_LENGTH)

/*
 * Common types for Device & Host
 */

/* BPM align data structures */

#define GPU_CIGAR_NULL        0
#define GPU_CIGAR_MATCH       1
#define GPU_CIGAR_MISSMATCH   2
#define GPU_CIGAR_INSERTION   3
#define GPU_CIGAR_DELETION    4

typedef char  gpu_bpm_align_cigar_event_t;
typedef char  gpu_bpm_align_cand_entry_t;
typedef char  gpu_bpm_align_qry_entry_t;

typedef struct{
  char      event;
  uint32_t	occurrences;
}gpu_bpm_align_cigar_entry_t;

typedef struct{
  uint32_t x;
  uint32_t y;
}gpu_bpm_align_coord_t;

typedef struct { /* each row 1 PEQ Entry (128bits) */
  uint32_t bitmap[GPU_BPM_ALIGN_PEQ_ALPHABET_SIZE][GPU_BPM_ALIGN_PEQ_SUBENTRIES];
} gpu_bpm_align_peq_entry_t;

typedef struct {
  uint32_t posEntryBase;
  uint32_t idQuery;
  uint32_t size;
  bool     leftGapAlign;
} gpu_bpm_align_cand_info_t;

typedef struct {
  uint32_t posEntryPEQ;
  uint32_t posEntryBase;
  uint32_t size;
} gpu_bpm_align_qry_info_t;

typedef struct {
  // Max allocation parameters for internal cigar results
  uint32_t                    offsetCigarStart;
  // Return Cigar results
  gpu_bpm_align_coord_t       initCood;
  gpu_bpm_align_coord_t       endCood;
  uint32_t 				      matchEffLenght;
  uint32_t 				      cigarStartPos;
  uint32_t				      cigarLenght;
} gpu_bpm_align_cigar_info_t;

/*
 * Obtain Buffers
 */
/* BPM align get primitives: internal buffers*/
gpu_bpm_align_qry_entry_t*   gpu_bpm_align_buffer_get_queries_(const void* const bpmBuffer);
gpu_bpm_align_peq_entry_t*   gpu_bpm_align_buffer_get_peq_entries_(const void* const bpmBuffer);
gpu_bpm_align_qry_info_t*    gpu_bpm_align_buffer_get_queries_info_(const void* const bpmBuffer);
gpu_bpm_align_cand_entry_t*  gpu_bpm_align_buffer_get_candidates_(const void* const bpmBuffer);
gpu_bpm_align_cand_info_t*   gpu_bpm_align_buffer_get_candidates_info_(const void* const bpmBuffer);
gpu_bpm_align_cigar_entry_t* gpu_bpm_align_buffer_get_cigars_(const void* const bpmBuffer);
gpu_bpm_align_cigar_info_t*  gpu_bpm_align_buffer_get_cigars_info_(const void* const bpmBuffer);

/*
 * Get elements
 */
/* BPM align get primitives: maximum allocatable entries*/
uint32_t gpu_bpm_align_buffer_get_max_peq_entries_(const void* const bpmBuffer);
uint32_t gpu_bpm_align_buffer_get_max_candidates_(const void* const bpmBuffer);
uint32_t gpu_bpm_align_buffer_get_max_candidate_size_(const void* const bpmBuffer);
uint32_t gpu_bpm_align_buffer_get_max_queries_(const void* const bpmBuffer);
uint32_t gpu_bpm_align_buffer_get_max_query_bases_(const void* const bpmBuffer);
uint32_t gpu_bpm_align_buffer_get_max_candidate_bases_(const void* const bpmBuffer);

/*
 * Main functions
 */
/* BPM align buffer primitives */
void gpu_bpm_align_init_buffer_(void* const bpmBuffer, const uint32_t averageQuerySize, const uint32_t averageCandidateSize, const uint32_t candidatesPerQuery);
void gpu_bpm_align_send_buffer_(void* const bpmBuffer, const uint32_t numPEQEntries, const uint32_t numQueryBases, const uint32_t numCandidateBases,
                                const uint32_t numQueries, const uint32_t numCandidates, const uint32_t queryBinSize);
void gpu_bpm_align_receive_buffer_(void* const bpmBuffer);
void gpu_bpm_align_init_and_realloc_buffer_(void *bpmBuffer, const uint32_t numPEQEntries, const uint32_t numQueryBases, const uint32_t numCandidateBases,
                                            const uint32_t numQueries, const uint32_t numCandidates);
#endif /* GPU_BPM_ALIGN_INTERFACE_H_ */
