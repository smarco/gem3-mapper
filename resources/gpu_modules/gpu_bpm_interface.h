/*
 * PROJECT: Bit-Parallel Myers on GPU
 * FILE: myers-interface.h
 * DATE: 4/7/2014
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 * DESCRIPTION: Interface for BPM on GPU
 */

#ifndef GPU_BPM_INTERFACE_H_
#define GPU_BPM_INTERFACE_H_

/*
 * Constants
 */
#define GPU_BPM_PEQ_ALPHABET_SIZE    5
#define GPU_BPM_PEQ_ENTRY_LENGTH     128
#define GPU_BPM_PEQ_SUBENTRY_LENGTH  32
#define GPU_BPM_PEQ_SUBENTRIES       (GPU_BPM_PEQ_ENTRY_LENGTH / GPU_UINT32_LENGTH)

/*
 * Enum types for Device & Host
 */
typedef enum
{
	GPU_REF_MFASTA_FILE,
	GPU_REF_PROFILE_FILE,
	GPU_REF_ASCII,
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

/*
 * Obtain Buffers
 */
inline gpu_bpm_qry_entry_t* gpu_bpm_buffer_get_peq_entries_(void* bpmBuffer);
inline gpu_bpm_cand_info_t* gpu_bpm_buffer_get_candidates_(void* bpmBuffer);
inline gpu_bpm_qry_info_t*  gpu_bpm_buffer_get_peq_info_(void* bpmBuffer);
inline gpu_bpm_alg_entry_t* gpu_bpm_buffer_get_alignments_(void* bpmBuffer);


/*
 * Get elements
 */
inline uint32_t gpu_bpm_buffer_get_max_peq_entries_(void* bpmBuffer);
inline uint32_t gpu_bpm_buffer_get_max_candidates_(void* bpmBuffer);
inline uint32_t gpu_bpm_buffer_get_max_queries_(void* bpmBuffer);

/*
 * Main functions
 */
inline void gpu_bpm_init_buffer_(void* bpmBuffer, const uint32_t averageNumPEQEntries, const uint32_t candidatesPerQuery);
inline void gpu_bpm_send_buffer_(void* bpmBuffer, const uint32_t numPEQEntries, const uint32_t numQueries, const uint32_t numCandidates, const sizeCandidates);
inline void gpu_bpm_receive_buffer_(void* bpmBuffer);

#endif /* GPU_BPM_INTERFACE_H_ */
