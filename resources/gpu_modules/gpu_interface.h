/*
 * PROJECT: Bit-Parallel Myers on GPU
 * FILE: myers-interface.h
 * DATE: 4/7/2014
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 * DESCRIPTION: Interface for BPM on GPU
 */

#include <stdint.h>
#include <stdbool.h>

/*
 * Constants
 */
#define GPU_UINT32_ONE_MASK   0x00000001u
#define	GPU_UINT32_LENGTH     32

#include "gpu_filter_interface.h"
#include "gpu_index_interface.h"

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
  /* GPU modules */
  GPU_FMI_EXACT_SEARCH	= GPU_UINT32_ONE_MASK << 0,
  GPU_FMI_DECODE_POS    = GPU_UINT32_ONE_MASK << 1,
  GPU_SA_DECODE_POS     = GPU_UINT32_ONE_MASK << 2,
  GPU_BPM               = GPU_UINT32_ONE_MASK << 3,
  /* GPU data structures */
  GPU_FMI               = GPU_FMI_EXACT_SEARCH | GPU_FMI_DECODE_POS | GPU_SA_DECODE_POS,
  GPU_SA                = GPU_SA_DECODE_POS,
  GPU_INDEX             = GPU_FMI | GPU_SA,
  GPU_REFERENCE         = GPU_BPM,
  /* GPU stages          */
  GPU_SEEDING           = GPU_INDEX,
  GPU_FILTERING         = GPU_BPM,

  GPU_NONE_MODULES      = 0,
  GPU_ALL_MODULES       = GPU_SEEDING | GPU_FILTERING
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

typedef struct {
  gpu_dev_arch_t      selectedArchitectures;
  gpu_data_location_t userAllocOption;
  gpu_module_t        activeModules;
  gpu_module_t        remoteStructures;
} gpu_info_dto_t;

typedef struct {
  void                **buffer;
  uint32_t            numBuffers;
  float               maxMbPerBuffer;
  gpu_module_t        activeModules;
} gpu_buffers_dto_t;


/*
 * Get elements
 */
uint32_t gpu_get_num_supported_devices_(const gpu_dev_arch_t selectedArchitectures);
uint32_t gpu_buffer_get_id_device_(const void* const gpu_buffer);
uint32_t gpu_buffer_get_id_supported_device_(const void* const gpuBuffer);


/*
 * Main functions
 */
void gpu_save_indexed_structures_GEM_(const char* const fileName, const gpu_gem_fmi_dto_t* const gemFMindex, const gpu_gem_ref_dto_t* const gemRef, const gpu_gem_sa_dto_t* const gemSAindex, const gpu_module_t activeModules);
void gpu_init_buffers_(gpu_buffers_dto_t* const buff, gpu_index_dto_t* const rawIndex, gpu_reference_dto_t* const rawRef, gpu_info_dto_t* const sys, const bool verbose);
void gpu_alloc_buffer_(void* const gpuBuffer);
void gpu_realloc_buffer_(void* const gpuBuffer, const float maxMbPerBuffer);
void gpu_destroy_buffers_(gpu_buffers_dto_t* buff);


