/*
 * PROJECT: Bit-Parallel Myers on GPU
 * FILE: myers-interface.h
 * DATE: 4/7/2014
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 * DESCRIPTION: Interface for BPM on GPU
 */

#include "gpu_bpm_interface.h"
#include "gpu_fmi_interface.h"

/*
 * Constants
 */
#define GPU_UINT32_ONE_MASK     0x00000001u
#define	GPU_UINT32_LENGTH		32

/*
 * Enum types for Device & Host
 */
typedef enum
{
	GPU_LOCAL_DATA,    /* Default */
	GPU_REMOTE_DATA,
	GPU_LOCAL_OR_REMOTE_DATA,
	GPU_NONE_DATA
} gpu_data_location_t;

typedef enum
{
	GPU_FMI_EXACT_SEARCH	= GPU_UINT32_ONE_MASK << 0,
	GPU_FMI_DECODE_POS		= GPU_UINT32_ONE_MASK << 1,
	GPU_BPM					= GPU_UINT32_ONE_MASK << 2,

	GPU_NONE_MODULES		= 0,
	GPU_ALL_MODULES 		= GPU_FMI_EXACT_SEARCH | GPU_FMI_DECODE_POS | GPU_BPM
} gpu_module_t;

typedef enum
{
	GPU_ARCH_TESLA      = GPU_UINT32_ONE_MASK << 0,
	GPU_ARCH_FERMI_1G   = GPU_UINT32_ONE_MASK << 1,
	GPU_ARCH_FERMI_2G   = GPU_UINT32_ONE_MASK << 2,
	GPU_ARCH_KEPLER_1G  = GPU_UINT32_ONE_MASK << 3,
	GPU_ARCH_KEPLER_2G  = GPU_UINT32_ONE_MASK << 4,
	GPU_ARCH_MAXWELL_1G = GPU_UINT32_ONE_MASK << 5,
	GPU_ARCH_MAXWELL_2G = GPU_UINT32_ONE_MASK << 6,
	GPU_ARCH_PASCAL_1G  = GPU_UINT32_ONE_MASK << 7,
	GPU_ARCH_PASCAL_2G  = GPU_UINT32_ONE_MASK << 8,

	GPU_ARCH_FERMI      = GPU_ARCH_FERMI_1G   | GPU_ARCH_FERMI_2G,
	GPU_ARCH_KEPLER     = GPU_ARCH_KEPLER_1G  | GPU_ARCH_KEPLER_2G,
	GPU_ARCH_MAXWELL    = GPU_ARCH_MAXWELL_1G | GPU_ARCH_MAXWELL_2G,
	GPU_ARCH_PASCAL     = GPU_ARCH_PASCAL_1G  | GPU_ARCH_PASCAL_2G,

	GPU_ARCH_NEWGEN     = GPU_UINT32_ONE_MASK << 31,
	GPU_ARCH_SUPPORTED  = GPU_ARCH_FERMI | GPU_ARCH_KEPLER | GPU_ARCH_MAXWELL | GPU_ARCH_NEWGEN
} gpu_dev_arch_t;


/*
 * Get elements
 */
uint32_t gpu_get_num_supported_devices_(gpu_dev_arch_t selectedArchitectures);
uint32_t gpu_buffer_get_id_device_(void* gpu_buffer);


/*
 * Main functions
 */

inline void gpu_init_buffers_(void ***gpuBuffer, uint32_t numBuffers, float maxMbPerBuffer,
						  	  const char *referenceRaw, gpu_ref_coding_t refCoding, const uint64_t refSize,
						  	  void *indexRaw, gpu_index_coding_t indexCoding, const uint64_t bwtSize,
						  	  const gpu_module_t activeModules, gpu_dev_arch_t selectedArchitectures,
						  	  gpu_data_location_t userAllocOption, const bool verbose);

inline void gpu_alloc_buffer_(void *gpuBuffer);
inline void gpu_destroy_buffers_(void*** gpuBuffer);

