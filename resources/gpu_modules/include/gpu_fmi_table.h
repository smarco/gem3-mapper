/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_FMI_TABLE_H_
#define GPU_FMI_TABLE_H_

#include "gpu_commons.h"
#include "gpu_devices.h"

#define GPU_FMI_TABLE_ALPHABET_SIZE       4
#define GPU_FMI_TABLE_MIN_ELEMENTS        2
#define GPU_FMI_TABLE_MIN_LEVELS          1
#define GPU_FMI_TABLE_DEFAULT_LEVELS      11
#define GPU_FMI_TABLE_DEFAULT_SKIP_LEVELS 0

#define GPU_FMI_TABLE_ENTRY_LENGTH        GPU_UINT64_LENGTH
#define GPU_FMI_TABLE_LINK_LENGTH         4
#define GPU_FMI_TABLE_FIELD_LENGTH        (GPU_FMI_TABLE_ENTRY_LENGTH - GPU_FMI_TABLE_LINK_LENGTH)
#define GPU_FMI_TABLE_LINK_MASK           (GPU_UINT64_ONES << GPU_FMI_TABLE_FIELD_LENGTH)
#define GPU_FMI_TABLE_FIELD_MASK          (GPU_UINT64_ONES >> GPU_FMI_TABLE_LINK_LENGTH)

#define GPU_FMI_TABLE_KEY_LENGTH          2
#define GPU_FMI_TABLE_KEY_MASK            (GPU_UINT32_ONES >> (GPU_UINT32_LENGTH - GPU_FMI_TABLE_KEY_LENGTH))

typedef struct{
  uint32_t init;
  uint32_t top;
} offset_table_t;

typedef enum
{
  GPU_FMI_TABLE_MULTILEVEL,
  GPU_FMI_TABLE_MULTILEVEL_LINKED,  /* Default */
  GPU_FMI_TABLE_DISABLED
} gpu_fmi_table_format_t;

typedef struct{
  gpu_fmi_table_format_t formatTableLUT;
  uint32_t               maxLevelsTableLUT;
  uint32_t               skipLevelsTableLUT;
  uint32_t               totalElemTableLUT;
  uint32_t               occThresholdTableLUT;
  offset_table_t*        h_offsetsTableLUT;
  offset_table_t**       d_offsetsTableLUT;
  gpu_sa_entry_t*        h_fmiTableLUT;
  gpu_sa_entry_t**       d_fmiTableLUT;
  memory_stats_t         hostAllocStats;
  memory_alloc_t*        memorySpace;
} gpu_fmi_table_t;

/* High level functions */
gpu_error_t gpu_fmi_table_get_num_elements(const uint32_t numLevels, uint32_t* const totalElemTableLUT);
gpu_error_t gpu_fmi_table_init_dto(gpu_fmi_table_t* const fmiTable);
gpu_error_t gpu_fmi_table_init(gpu_fmi_table_t* const fmiTable, const uint32_t numLevels, const uint32_t numSupportedDevices);
gpu_error_t gpu_fmi_table_allocate(gpu_fmi_table_t* const fmiTable);
gpu_error_t gpu_fmi_table_read_specs(int fp, gpu_fmi_table_t* const fmiTable);
gpu_error_t gpu_fmi_table_load_default_specs(gpu_fmi_table_t* const fmiTable);
gpu_error_t gpu_fmi_table_read(int fp, gpu_fmi_table_t* const fmiTable);
gpu_error_t gpu_fmi_table_write_specs(int fp, const gpu_fmi_table_t* const fmiTable);
gpu_error_t gpu_fmi_table_write(int fp, const gpu_fmi_table_t* const fmiTable);
gpu_error_t gpu_fmi_table_transfer_CPU_to_GPUs(gpu_fmi_table_t* const fmiTable, gpu_device_info_t** const devices);
gpu_error_t gpu_fmi_table_construction(gpu_fmi_table_t* const fmiTable, const gpu_fmi_entry_t* const h_fmi, const uint64_t bwtSize);
gpu_error_t gpu_fmi_table_get_size(const gpu_fmi_table_t* const fmiTable, size_t* const bytesPerFmiTable);
gpu_error_t gpu_fmi_table_free_host(gpu_fmi_table_t* const fmiTable);
gpu_error_t gpu_fmi_table_free_unused_host(gpu_fmi_table_t* const fmiTable, gpu_device_info_t** const devices);
gpu_error_t gpu_fmi_table_free_device(gpu_fmi_table_t* const fmiTable, gpu_device_info_t** const devices);
gpu_error_t gpu_fmi_table_free_metainfo(gpu_fmi_table_t* const fmiTable);

/* Low level functions */
gpu_error_t gpu_fmi_table_print_links(const offset_table_t* const offsetsTableLUT, gpu_sa_entry_t* const fmiTableLUT);
gpu_error_t gpu_fmi_table_process_links(const uint32_t numLevels, const uint64_t occThreshold, const offset_table_t* const offsetsTableLUT, gpu_sa_entry_t* const fmiTableLUT);
gpu_error_t gpu_fmi_table_build(gpu_fmi_table_t* const fmiTable, const gpu_fmi_entry_t* const h_fmi, const uint64_t bwtSize);
gpu_error_t gpu_fmi_table_init_offsets(offset_table_t* const offsetsTableLUT, const uint32_t numLevels);
void        gpu_fmi_table_process_entry(const gpu_fmi_entry_t* const h_fmi, gpu_sa_entry_t* const currentIntervals, const gpu_sa_entry_t L, const gpu_sa_entry_t R);
gpu_error_t gpu_fmi_table_process_forward_links(const uint32_t idLevel, const uint64_t occThreshold, const offset_table_t* const offsetsTableLUT, gpu_sa_entry_t* const fmiTableLUT);
gpu_error_t gpu_fmi_table_process_backward_links(const uint32_t idLevel, const uint64_t occThreshold, const offset_table_t* const offsetsTableLUT, gpu_sa_entry_t* const fmiTableLUT);
gpu_error_t gpu_fmi_table_process_forward_level(const gpu_fmi_entry_t* const h_fmi, const uint32_t idLevel, offset_table_t* const offsetsTableLUT, gpu_sa_entry_t* const fmiTableLUT);
gpu_error_t gpu_fmi_table_process_backward_level(const gpu_fmi_entry_t* const h_fmi, const uint32_t idLevel, offset_table_t* const offsetsTableLUT, gpu_sa_entry_t* const fmiTableLUT);
gpu_error_t gpu_fmi_table_get_positions(const uint32_t idLevel, const uint32_t idEntry, const offset_table_t* const offsetsTableLUT, uint32_t* const idGlobalL, uint32_t* const idGlobalR);

#endif /* GPU_FMI_TABLE_H_ */
