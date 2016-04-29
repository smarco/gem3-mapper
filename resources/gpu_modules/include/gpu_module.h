/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */


#ifndef GPU_MODULE_H_
#define GPU_MODULE_H_

/* Include the required objects */
#include "gpu_devices.h"
#include "gpu_reference.h"
#include "gpu_index.h"
#include "gpu_buffer.h"

/* Primitives to get information from modules */
uint32_t      gpu_module_get_num_allocated(const gpu_module_t activeModules);
uint32_t      gpu_module_get_num_structures();
gpu_module_t* gpu_module_get_list_structures();
gpu_error_t   gpu_module_get_min_memory(gpu_reference_buffer_t* const __restrict__ reference, const gpu_index_buffer_t* const __restrict__ index,
                                        const uint32_t numBuffers, const gpu_module_t activeModules, size_t *minimumMemorySize);
gpu_error_t   gpu_module_set_device_allocation(gpu_reference_buffer_t* const __restrict__ reference, gpu_index_buffer_t* const __restrict__ index,
                                               uint32_t idSupDevice, gpu_module_t allocatedModules);
gpu_error_t   gpu_module_get_device_allocation(gpu_reference_buffer_t* const __restrict__ reference, gpu_index_buffer_t* const __restrict__ index,
                                               uint32_t idSupDevice, gpu_module_t* allocatedModules);

/* Primitives initialize modules */
gpu_error_t   gpu_module_allocator_per_device(gpu_reference_buffer_t* const __restrict__ reference, gpu_index_buffer_t* const __restrict__ index, const uint32_t idDevice,
                                              const uint32_t numBuffers, const gpu_module_t activeModules, gpu_module_t* const __restrict__ allocatedModules);
gpu_error_t   gpu_module_manager_per_device(gpu_reference_buffer_t* const __restrict__ reference, gpu_index_buffer_t* const __restrict__ index,
                                            const uint32_t idDevice, const uint32_t numBuffers, const gpu_data_location_t userAllocOption,
                                            gpu_module_t* const __restrict__ activatedModules, gpu_module_t* const __restrict__ allocatedModules);
gpu_error_t   gpu_module_manager_all_system(gpu_reference_buffer_t* const __restrict__ reference, gpu_index_buffer_t* const __restrict__ index,
                                            const uint32_t numBuffers, const gpu_dev_arch_t selectedArchitectures,
                                            const gpu_data_location_t userAllocOption, gpu_module_t* const __restrict__ globalModules,
                                            gpu_module_t* const __restrict__ globalStructures);

/* Primitives configure modules */
gpu_error_t   gpu_module_memory_requirements_per_device(gpu_reference_buffer_t* const __restrict__ reference, gpu_index_buffer_t* const __restrict__ index,
                                                        const uint32_t idDevice, const uint32_t numBuffers, const gpu_data_location_t userAllocOption,
                                                        gpu_module_t* const __restrict__ maxAllocatedModules, bool* const __restrict__ maskedDevice);
gpu_error_t   gpu_module_configure_system(gpu_reference_buffer_t* const __restrict__ reference, gpu_index_buffer_t* const __restrict__ index,
                                          gpu_device_info_t ***devices, const uint32_t numBuffers,
                                          const gpu_dev_arch_t selectedArchitectures, const gpu_data_location_t userAllocOption,
                                          gpu_module_t* const __restrict__ activatedModules, gpu_module_t* const __restrict__ allocatedStructures);
gpu_error_t   gpu_module_manager_memory(gpu_reference_buffer_t* const __restrict__ reference, gpu_index_buffer_t* const __restrict__ index,
                                        const uint32_t numBuffers, const gpu_dev_arch_t selectedArchitectures,
                                        const gpu_data_location_t userAllocOption,
                                        size_t* const __restrict__ recomendedMemorySize, size_t* const __restrict__ requiredMemorySize,
                                        gpu_module_t* const __restrict__ modules, gpu_module_t* const __restrict__ structures);
/* Primitives search module configurations */
gpu_error_t   gpu_module_search_active(gpu_module_t* const __restrict__ allocatedModulesPerDevice, const uint32_t numSupportedDevices,
                                       gpu_module_t* const __restrict__ activatedModules);
gpu_error_t   gpu_module_search_structures(gpu_module_t* const __restrict__ allocatedModulesPerDevice, gpu_module_t* const __restrict__ allocatedStructuresPerDevice,
                                           const uint32_t numSupportedDevices, const gpu_module_t activatedModules,
                                           gpu_module_t* const __restrict__ allocatedStructures);


#endif /* GPU_MODULE_H_ */
