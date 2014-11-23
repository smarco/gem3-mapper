/*
 * PROJECT: Bit-Parallel Myers on GPU
 * FILE: myers-interface.h
 * DATE: 4/7/2014
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 * DESCRIPTION: Interface for BPM on GPU
 */

#include <stdbool.h>
#include <stdint.h>

/*
 * Constants
 */
#define BMP_GPU_PEQ_ALPHABET_SIZE    5
#define BMP_GPU_PEQ_ENTRY_LENGTH     128
#define BMP_GPU_PEQ_SUBENTRY_LENGTH  32
#define BMP_GPU_UINT32_LENGTH        32
#define BMP_GPU_UINT64_LENGTH        64
#define BMP_GPU_PEQ_SUBENTRIES       (BMP_GPU_PEQ_ENTRY_LENGTH / BMP_GPU_UINT32_LENGTH)
#define BMP_GPU_UINT32_ONE_MASK      0x00000001u

/*
 * Enum types for Device & Host
 */
typedef enum
{
 MFASTA_FILE,
 PROFILE_REFERENCE_FILE,
 ASCII,
 GEM
} bpm_gpu_ref_coding_t;

typedef enum
{
	LOCAL_REFERENCE,    /* Default */
	REMOTE_REFERENCE,
	LOCAL_OR_REMOTE_REFERENCE
} bpm_gpu_ref_location_t;

typedef enum
{
  ARCH_TESLA      = BMP_GPU_UINT32_ONE_MASK << 0,
  ARCH_FERMI_1G   = BMP_GPU_UINT32_ONE_MASK << 1,
  ARCH_FERMI_2G   = BMP_GPU_UINT32_ONE_MASK << 2,
  ARCH_KEPLER_1G  = BMP_GPU_UINT32_ONE_MASK << 3,
  ARCH_KEPLER_2G  = BMP_GPU_UINT32_ONE_MASK << 4,
  ARCH_MAXWELL    = BMP_GPU_UINT32_ONE_MASK << 5,
  ARCH_NEWGEN     = BMP_GPU_UINT32_ONE_MASK << 31,
  ARCH_SUPPORTED  = ARCH_FERMI_1G | ARCH_FERMI_2G | ARCH_KEPLER_1G | ARCH_KEPLER_2G | ARCH_MAXWELL | ARCH_NEWGEN
} bpm_gpu_dev_arch_t;

/*
 * Common types for Device & Host
 */
typedef struct { /* each row 1 PEQ Entry (128bits) */
 uint32_t bitmap[BMP_GPU_PEQ_ALPHABET_SIZE][BMP_GPU_PEQ_SUBENTRIES];
} bpm_gpu_qry_entry_t;

typedef struct {
 uint32_t column;
 uint32_t score;
} bpm_gpu_res_entry_t;

typedef struct {
 uint64_t position;
 uint32_t query;
 uint32_t size;
} bpm_gpu_cand_info_t;

typedef struct {
 uint32_t posEntry;
 uint32_t size;
} bpm_gpu_qry_info_t;

/*
 * Obtain Buffers
 */
inline bpm_gpu_qry_entry_t* bpm_gpu_buffer_get_peq_entries_(void* myersBuffer);
inline bpm_gpu_cand_info_t* bpm_gpu_buffer_get_candidates_(void* myersBuffer);
inline bpm_gpu_qry_info_t*  bpm_gpu_buffer_get_peq_info_(void* myersBuffer);
inline bpm_gpu_res_entry_t* bpm_gpu_buffer_get_results_(void* myersBuffer);

/*
 * Get elements
 */
inline uint32_t bpm_gpu_get_num_supported_devices_();
inline uint32_t bpm_gpu_buffer_get_max_peq_entries_(void* myersBuffer);
inline uint32_t bpm_gpu_buffer_get_max_candidates_(void* myersBuffer);
inline uint32_t bpm_gpu_buffer_get_max_queries_(void* myersBuffer);
inline uint32_t bpm_gpu_buffer_get_id_device_(void* myersBuffer);

/*
 * Main functions
 */
void bpm_gpu_init_(void*** myersBuffer,uint32_t numBuffers,uint32_t maxMbPerBuffer,
    const char* referenceRaw,bpm_gpu_ref_coding_t refCoding,const uint64_t refSize,
    uint32_t averageQuerySize,uint32_t candidatesPerQuery,
	bpm_gpu_dev_arch_t selectedArchitectures, bpm_gpu_ref_location_t userReferenceAllocOption,
	const bool verbose);
void bpm_gpu_init_buffer_(void *myersBuffer);
void bpm_gpu_send_buffer_(void* myersBuffer,uint32_t numPEQEntries,uint32_t numQueries,uint32_t numCandidates);
void bpm_gpu_receive_buffer_(void* myersBuffer);
void bpm_gpu_destroy_(void*** myersBuffer);

