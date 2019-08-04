/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2011-2018 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *                2013-2016 by Santiago Marco-Sola <santiagomsola@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#include <stdint.h>
#include <stdbool.h>

/*
 * Constants
 */
#define GPU_UINT32_ONE_MASK   0x00000001u
#define	GPU_UINT32_LENGTH     32

#include "gpu_filter_interface.h"
#include "gpu_align_interface.h"
#include "gpu_index_interface.h"

/*
 * Enum types for Device & Host
 */
typedef enum
{
  GPU_LOCAL_DATA,    /* Default */
  GPU_REMOTE_DATA,
  GPU_LOCAL_OR_REMOTE_DATA,
  GPU_GEM_POLICY,
  GPU_NONE_DATA             = 0
} gpu_data_location_t;

typedef enum
{
  /* GPU modules */
  GPU_FMI_ADAPT_SEARCH  = GPU_UINT32_ONE_MASK << 0,
  GPU_FMI_EXACT_SEARCH	= GPU_UINT32_ONE_MASK << 1,
  GPU_FMI_DECODE_POS    = GPU_UINT32_ONE_MASK << 2,
  GPU_SA_DECODE_POS     = GPU_UINT32_ONE_MASK << 3,
  GPU_BPM_FILTER        = GPU_UINT32_ONE_MASK << 4,
  GPU_KMER_FILTER       = GPU_UINT32_ONE_MASK << 5,
  GPU_BPM_ALIGN         = GPU_UINT32_ONE_MASK << 6,
  GPU_SWG_ALIGN         = GPU_UINT32_ONE_MASK << 7,
  /* GPU data structures */
  GPU_FMI               = GPU_FMI_ADAPT_SEARCH | GPU_FMI_EXACT_SEARCH | GPU_FMI_DECODE_POS,
  GPU_SA                = GPU_SA_DECODE_POS,
  GPU_INDEX             = GPU_FMI | GPU_SA,
  GPU_REFERENCE_MASKED	= GPU_BPM_ALIGN,
  GPU_REFERENCE_PLAIN   = GPU_BPM_FILTER | GPU_KMER_FILTER,
  GPU_REFERENCE         = GPU_REFERENCE_PLAIN | GPU_REFERENCE_MASKED,
  /* GPU stages          */
  GPU_SEEDING           = GPU_INDEX,
  GPU_FILTERING         = GPU_REFERENCE_PLAIN,
  GPU_ALIGNMENT         = GPU_REFERENCE,
  /* General setups      */
  GPU_NONE_MODULES      = 0,
  GPU_ALL_MODULES       = GPU_SEEDING | GPU_FILTERING | GPU_ALIGNMENT
} gpu_module_t;

typedef enum
{
  /* Non-supported GPU architectures */
  GPU_ARCH_TESLA      = GPU_UINT32_ONE_MASK << 0,
  /* Supported GPU architectures     */
  GPU_ARCH_FERMI_1G   = GPU_UINT32_ONE_MASK << 1,
  GPU_ARCH_FERMI_2G   = GPU_UINT32_ONE_MASK << 2,
  GPU_ARCH_KEPLER_1G  = GPU_UINT32_ONE_MASK << 3,
  GPU_ARCH_KEPLER_2G  = GPU_UINT32_ONE_MASK << 4,
  GPU_ARCH_MAXWELL_1G = GPU_UINT32_ONE_MASK << 5,
  GPU_ARCH_MAXWELL_2G = GPU_UINT32_ONE_MASK << 6,
  GPU_ARCH_PASCAL_1G  = GPU_UINT32_ONE_MASK << 7,
  GPU_ARCH_PASCAL_2G  = GPU_UINT32_ONE_MASK << 8,
  GPU_ARCH_VOLTA_1G   = GPU_UINT32_ONE_MASK << 9,
  GPU_ARCH_VOLTA_2G   = GPU_UINT32_ONE_MASK << 10,
  /* Main GPU Architectures          */
  GPU_ARCH_FERMI      = GPU_ARCH_FERMI_1G   | GPU_ARCH_FERMI_2G,
  GPU_ARCH_KEPLER     = GPU_ARCH_KEPLER_1G  | GPU_ARCH_KEPLER_2G,
  GPU_ARCH_MAXWELL    = GPU_ARCH_MAXWELL_1G | GPU_ARCH_MAXWELL_2G,
  GPU_ARCH_PASCAL     = GPU_ARCH_PASCAL_1G  | GPU_ARCH_PASCAL_2G,
  GPU_ARCH_VOLTA      = GPU_ARCH_VOLTA_1G   | GPU_ARCH_VOLTA_2G,
  /* General setups                  */
  GPU_ARCH_NEWGEN     = GPU_UINT32_ONE_MASK << 31,
  GPU_ARCH_SUPPORTED  = GPU_ARCH_FERMI | GPU_ARCH_KEPLER | GPU_ARCH_MAXWELL | GPU_ARCH_PASCAL | GPU_ARCH_VOLTA | GPU_ARCH_NEWGEN
} gpu_dev_arch_t;

typedef struct {
  gpu_dev_arch_t      selectedArchitectures;
  gpu_data_location_t userAllocOption;
  gpu_module_t        activatedModules;
  gpu_module_t        allocatedStructures;
  bool                verbose;
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
void gpu_io_save_indexed_structures_GEM_(const char* const fileName, const gpu_gem_fmi_dto_t* const gemFMindex, const gpu_gem_ref_dto_t* const gemRef, const gpu_gem_sa_dto_t* const gemSAindex, const gpu_module_t activeModules);
void gpu_init_buffers_(gpu_buffers_dto_t* const buff, gpu_index_dto_t* const rawIndex, gpu_reference_dto_t* const rawRef, gpu_info_dto_t* const sys);
void gpu_alloc_buffer_(void* const gpuBuffer, const uint64_t idThread);
void gpu_realloc_buffer_(void* const gpuBuffer, const float maxMbPerBuffer);
void gpu_destroy_buffers_(gpu_buffers_dto_t* buff);


