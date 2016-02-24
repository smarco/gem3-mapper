/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_MODULE_C_
#define GPU_MODULE_C_

#include "../include/gpu_module.h"

// Sorted by allocation priority
gpu_module_t structureList[] = {GPU_REFERENCE, GPU_FMI, GPU_SA};

uint32_t gpu_module_get_num_allocated(const gpu_module_t activeModules)
{
  return (gpu_count_active_bits(activeModules));
}

uint32_t gpu_module_get_num_structures()
{
  return ((uint32_t)(sizeof(structureList) / sizeof(structureList[0])));
}

gpu_module_t* gpu_module_get_list_structures()
{
  return(structureList);
}

gpu_error_t gpu_module_get_min_memory(gpu_reference_buffer_t* const reference, const gpu_index_buffer_t* const index,
                                      const uint32_t numBuffers, const gpu_module_t activeModules, size_t *minimumMemorySize)
{
  size_t memorySize;
  size_t bytesPerReference = 0, bytesPerFMIndex = 0, bytesPerSAIndex = 0, bytesPerBuffer = 0;

  GPU_ERROR(gpu_buffer_get_min_memory_size(&bytesPerBuffer));
  GPU_ERROR(gpu_reference_get_size(reference, &bytesPerReference));
  GPU_ERROR(gpu_index_get_size(index, &bytesPerFMIndex, GPU_FMI));
  GPU_ERROR(gpu_index_get_size(index, &bytesPerSAIndex, GPU_SA));

  memorySize = numBuffers * bytesPerBuffer;
  if(activeModules & GPU_REFERENCE) memorySize += bytesPerReference;
  if(activeModules & GPU_FMI) memorySize += bytesPerFMIndex;
  if(activeModules & GPU_SA) memorySize += bytesPerSAIndex;

  (* minimumMemorySize) = memorySize;
  return (SUCCESS);
}

gpu_error_t gpu_module_set_device_allocation(gpu_reference_buffer_t* const reference, gpu_index_buffer_t* const index,
                                            uint32_t idSupDevice, gpu_module_t allocatedModules)
{
  memory_alloc_t moduleMemorySpace;

  moduleMemorySpace = GPU_HOST_MAPPED;
  if(allocatedModules & GPU_REFERENCE) moduleMemorySpace = GPU_DEVICE_MAPPED;
  reference->memorySpace[idSupDevice] = moduleMemorySpace;

  moduleMemorySpace = GPU_HOST_MAPPED;
  if(allocatedModules & GPU_FMI) moduleMemorySpace = GPU_DEVICE_MAPPED;
  index->fmi.memorySpace[idSupDevice] = moduleMemorySpace;

  moduleMemorySpace = GPU_HOST_MAPPED;
  if(allocatedModules & GPU_SA) moduleMemorySpace = GPU_DEVICE_MAPPED;
  index->sa.memorySpace[idSupDevice] = moduleMemorySpace;

  return (SUCCESS);
}

gpu_error_t gpu_module_get_device_allocation(gpu_reference_buffer_t* const reference, gpu_index_buffer_t* const index,
                                            uint32_t idSupDevice, gpu_module_t* allocatedModules)
{
  (* allocatedModules) = GPU_NONE_MODULES;
  if (reference->memorySpace[idSupDevice] == GPU_DEVICE_MAPPED) (* allocatedModules) |= GPU_REFERENCE;
  if (index->fmi.memorySpace[idSupDevice] == GPU_DEVICE_MAPPED) (* allocatedModules) |= GPU_FMI;
  if (index->sa.memorySpace[idSupDevice]  == GPU_DEVICE_MAPPED) (* allocatedModules) |= GPU_SA;
  return (SUCCESS);
}

gpu_error_t gpu_module_allocator_per_device(gpu_reference_buffer_t* const reference, gpu_index_buffer_t* const index, const uint32_t idDevice,
                                            const uint32_t numBuffers, const gpu_module_t activeModules, gpu_module_t* const allocatedModules)
{
  const gpu_module_t* const moduleList = gpu_module_get_list_structures();
  size_t memoryFree  = gpu_device_get_free_memory(idDevice), currentMemorySize = 0;
  gpu_module_t currentModule, totalActiveModules = GPU_NONE_MODULES;
  uint32_t idModule, numModules = gpu_module_get_num_structures();

  for(idModule = 0; idModule < numModules; ++idModule){
    currentModule = moduleList[idModule];
    if(activeModules & currentModule){
      GPU_ERROR(gpu_module_get_min_memory(reference, index, numBuffers, totalActiveModules|currentModule, &currentMemorySize));
      // allocate the module in the device
      if (currentMemorySize < memoryFree) totalActiveModules |= currentModule;
    }
  }

  (* allocatedModules) = totalActiveModules;
  return(SUCCESS);
}

gpu_error_t gpu_module_manager_per_device(gpu_reference_buffer_t* const reference, gpu_index_buffer_t* const index,
                                          const uint32_t idDevice, const uint32_t numBuffers, const gpu_data_location_t userAllocOption,
                                          gpu_module_t* const activatedModules, gpu_module_t* const allocatedModules)
{
  const gpu_module_t userSelectedModules = reference->activeModules | index->activeModules;
  gpu_module_t tmpAllocatedModules = GPU_NONE_MODULES;

  switch (userAllocOption){
    case GPU_REMOTE_DATA:          // (Force to allocate all structures in HOST)
      (* activatedModules) = userSelectedModules;
      (* allocatedModules) = GPU_NONE_MODULES;
      break;
    case GPU_LOCAL_DATA:           // (Force to allocate all structures in DEVICE)
      (* activatedModules) = userSelectedModules;
      (* allocatedModules) = userSelectedModules;
      break;
    case GPU_LOCAL_OR_REMOTE_DATA: // (Best effort allocating structures in DEVICE and the rest in HOST)
      GPU_ERROR(gpu_module_allocator_per_device(reference, index, idDevice, numBuffers, userSelectedModules, &tmpAllocatedModules));
      (* activatedModules) = userSelectedModules;
      (* allocatedModules) = tmpAllocatedModules;
      break;
    case GPU_GEM_POLICY:           // (GEM Default Allocation)
      GPU_ERROR(gpu_module_allocator_per_device(reference, index, idDevice, numBuffers, userSelectedModules, &tmpAllocatedModules));
      (* activatedModules) = tmpAllocatedModules | GPU_REFERENCE; // Reference is the minimum required module to be execute
      (* allocatedModules) = tmpAllocatedModules;
      break;
    default:
      return(E_DATA_NOT_ALLOCATED);
  }
  return(SUCCESS);
}

gpu_error_t gpu_module_memory_requirements_per_device(gpu_reference_buffer_t* const reference, gpu_index_buffer_t* const index,
                                                      const uint32_t idDevice, const uint32_t numBuffers, const gpu_data_location_t userAllocOption,
                                                      gpu_module_t* const maxAllocatedModules, bool* const maskedDevice)
{
  const gpu_module_t userSelectedModules = reference->activeModules | index->activeModules;
  size_t memoryFree  = gpu_device_get_free_memory(idDevice);
  size_t minimumMemorySize = 0;

  switch (userAllocOption){
    case GPU_REMOTE_DATA:          // (Force to allocate all structures in HOST)
    case GPU_LOCAL_OR_REMOTE_DATA: // (Best effort allocating structures in DEVICE and the rest in HOST)
      GPU_ERROR(gpu_module_get_min_memory(reference, index, numBuffers, GPU_NONE_MODULES, &minimumMemorySize));
      break;
    case GPU_LOCAL_DATA:           // (Force to allocate all structures in DEVICE)
      GPU_ERROR(gpu_module_get_min_memory(reference, index, numBuffers, userSelectedModules, &minimumMemorySize));
      break;
    case GPU_GEM_POLICY:           // (GEM Default Allocation)
      GPU_ERROR(gpu_module_get_min_memory(reference, index, numBuffers, GPU_NONE_MODULES, &minimumMemorySize));
      break;
    default:
      return(E_DATA_NOT_ALLOCATED);
  }

  (* maskedDevice)    = memoryFree < minimumMemorySize;
  //(* memoryAllocated) = minimumMemorySize;
  return(SUCCESS);
}

gpu_error_t gpu_module_manager_all_system(gpu_reference_buffer_t* const reference, gpu_index_buffer_t* const index,
                                          const uint32_t numBuffers, const gpu_dev_arch_t selectedArchitectures,
                                          const gpu_data_location_t userAllocOption, gpu_module_t* const globalModules,
                                          gpu_module_t* const globalStructures)
{
  const uint32_t numDevices = gpu_device_get_num_all();
  const uint32_t numSupportedDevices = gpu_get_num_supported_devices_(selectedArchitectures);

  uint32_t idDevice, idSupportedDevice;
  uint32_t numActiveModules, maxNumActiveModules;
  gpu_module_t maxActiveModules, maxAllocStructures, activatedModules, allocatedStructure;
  gpu_module_t allocatedModulesPerDevice[numSupportedDevices], allocatedStructuresPerDevice[numSupportedDevices];

  // Analyze all the devices on the system and annotate the module requirements
  for(idDevice = 0, idSupportedDevice = 0; idDevice < numDevices; ++idDevice){
    const bool deviceArchSupported = gpu_device_get_architecture(idDevice) & selectedArchitectures;
    if(deviceArchSupported){
      GPU_ERROR(gpu_module_manager_per_device(reference, index, idDevice, numBuffers, userAllocOption,
                                              &activatedModules, &allocatedStructure));
      allocatedModulesPerDevice[idSupportedDevice]    = activatedModules;
      allocatedStructuresPerDevice[idSupportedDevice] = allocatedStructure;
      idSupportedDevice++;
    }
  }

  // Initialization for the best fitting modules exploration
  maxActiveModules    = GPU_NONE_MODULES;
  maxAllocStructures  = GPU_NONE_MODULES;
  numActiveModules    = 0;
  maxNumActiveModules = 0;

  // Explore all the system requirements to fit the best global module configuration
  for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
    numActiveModules = gpu_module_get_num_allocated(allocatedModulesPerDevice[idSupportedDevice]);
    if(numActiveModules > maxNumActiveModules){
      maxNumActiveModules = numActiveModules;
      maxActiveModules    = allocatedModulesPerDevice[idSupportedDevice];
      maxAllocStructures  = allocatedStructuresPerDevice[idSupportedDevice];
    }
  }

  (* globalModules)    = maxActiveModules;
  (* globalStructures) = maxAllocStructures;
  return(SUCCESS);
}

gpu_error_t gpu_module_configure_system(gpu_reference_buffer_t* const reference, gpu_index_buffer_t* const index,
                                        gpu_device_info_t ***devices, const uint32_t numBuffers,
                                        const gpu_dev_arch_t selectedArchitectures, const gpu_data_location_t userAllocOption,
                                        gpu_module_t* const activatedModules, gpu_module_t* const allocatedStructures)
{
  const uint32_t numDevices = gpu_device_get_num_all();
  const uint32_t numSupportedDevices = gpu_get_num_supported_devices_(selectedArchitectures);
  gpu_device_info_t **dev = (gpu_device_info_t **) malloc(numSupportedDevices * sizeof(gpu_device_info_t *));

  const gpu_module_t userRequestedModules = reference->activeModules | index->activeModules;
  size_t recomendedMemorySize, requiredMemorySize;
  uint32_t idDevice, idSupportedDevice = 0;
  gpu_module_t globalModules = GPU_NONE_MODULES, globalStructures = GPU_NONE_MODULES;

  // Explores the best combination of active modules and stored structures for all the devices
  GPU_ERROR(gpu_module_manager_all_system(reference, index, numBuffers, selectedArchitectures, userAllocOption, &globalModules, &globalStructures));
  // Calculates the device memory constrains for the above configuration
  GPU_ERROR(gpu_module_get_min_memory(reference, index, numBuffers, GPU_NONE_MODULES, &requiredMemorySize));
  GPU_ERROR(gpu_module_get_min_memory(reference, index, numBuffers, userRequestedModules, &recomendedMemorySize));

  // Activates the devices with enough requirements in the system
  for(idDevice = 0, idSupportedDevice = 0; idDevice < numDevices; ++idDevice){
    size_t memoryFree               = gpu_device_get_free_memory(idDevice);
    const bool deviceArchSupported  = gpu_device_get_architecture(idDevice) & selectedArchitectures;
    const bool dataFitsMemoryDevice = requiredMemorySize < memoryFree;
    gpu_device_screen_status(idDevice, deviceArchSupported, recomendedMemorySize, requiredMemorySize);
    if(deviceArchSupported && dataFitsMemoryDevice){ //Data fits on memory device
      GPU_ERROR(gpu_device_init(&dev[idDevice], idDevice, idSupportedDevice, selectedArchitectures));
      GPU_ERROR(gpu_module_set_device_allocation(reference, index, idSupportedDevice, globalModules)); // TODO: recalculate and use globalStructures
      idSupportedDevice++;
    }
  }

  GPU_ERROR(gpu_device_characterize_all(dev, idSupportedDevice));

  reference->activeModules = globalModules & GPU_REFERENCE;
  index->activeModules     = globalModules & GPU_INDEX;
  (*activatedModules)      = globalModules;
  (*allocatedStructures)   = globalStructures;

  (* devices)              = dev;
  return(SUCCESS);
}

#endif /* GPU_MODULE_C_ */
