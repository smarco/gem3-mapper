/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2011-2018 by Alejandro Chacon    <alejandro.chacond@gmail.com>
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
  GPU_ERROR(gpu_reference_get_size(reference, &bytesPerReference, GPU_REFERENCE));
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
  index->fmi.table.memorySpace[idSupDevice] = moduleMemorySpace;

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
      (* activatedModules) = tmpAllocatedModules | GPU_REFERENCE | GPU_BPM_ALIGN; // Reference is the minimum required module to be execute
      (* allocatedModules) = tmpAllocatedModules;
      break;
    default:
      return(E_NOT_SUPPORTED_ALLOC_POLICY);
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
      return(E_NOT_SUPPORTED_ALLOC_POLICY);
  }

  (* maskedDevice)    = memoryFree < minimumMemorySize;
  return(SUCCESS);
}

gpu_error_t gpu_module_search_structures(gpu_module_t* const allocatedModulesPerDevice, gpu_module_t* const allocatedStructuresPerDevice,
                                         const uint32_t numSupportedDevices, const gpu_module_t activatedModules,
                                         gpu_module_t* const allocatedStructures)
{
  // Initialization for the best fitting structures exploration
  const uint32_t maxNumActiveModules            = gpu_module_get_num_allocated(activatedModules);
  gpu_module_t   minAllocatedStructures         = GPU_ALL_MODULES;
  uint32_t       minNumAllocatedStructures      = gpu_module_get_num_allocated(minAllocatedStructures);
  uint32_t       currentNumAllocatedStructures  = 0;
  uint32_t       currentNumActiveModules        = 0;
  // Initialization for the device exploration
  uint32_t       idSupportedDevice              = 0;

  // Discover which devices match with the maximum active Modules and sets the local/remote structures
  for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
    currentNumActiveModules       = gpu_module_get_num_allocated(allocatedModulesPerDevice[idSupportedDevice]);
    currentNumAllocatedStructures = gpu_module_get_num_allocated(allocatedStructuresPerDevice[idSupportedDevice]);
    if((currentNumActiveModules == maxNumActiveModules) && (currentNumAllocatedStructures < minNumAllocatedStructures)){
      minAllocatedStructures    = allocatedStructuresPerDevice[idSupportedDevice];
      minNumAllocatedStructures = gpu_module_get_num_allocated(minAllocatedStructures);
    }
  }

  (* allocatedStructures) = minAllocatedStructures;
  return(SUCCESS);
}

gpu_error_t gpu_module_search_active(gpu_module_t* const allocatedModulesPerDevice, const uint32_t numSupportedDevices,
                                     gpu_module_t* const activatedModules)
{
  // Initialization for the best fitting modules exploration
  gpu_module_t maxActiveModules    = GPU_NONE_MODULES;
  uint32_t     maxNumActiveModules = gpu_module_get_num_allocated(maxActiveModules);
  uint32_t     numActiveModules    = 0;
  // Initialization for the device exploration
  uint32_t     idSupportedDevice   = 0;

  // Explore all the system requirements to fit the best global module configuration
  for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
    numActiveModules      = gpu_module_get_num_allocated(allocatedModulesPerDevice[idSupportedDevice]);
    if(numActiveModules > maxNumActiveModules){
      maxActiveModules    = allocatedModulesPerDevice[idSupportedDevice];
      maxNumActiveModules = gpu_module_get_num_allocated(maxActiveModules);
    }
  }

  (* activatedModules) = maxActiveModules;
  return(SUCCESS);
}

gpu_error_t gpu_module_manager_all_system(gpu_reference_buffer_t* const reference, gpu_index_buffer_t* const index,
                                          const uint32_t numBuffers, const gpu_dev_arch_t selectedArchitectures,
                                          const gpu_data_location_t userAllocOption, gpu_module_t* const globalModules,
                                          gpu_module_t* const globalStructures)
{
  // Device initialization variables
  const uint32_t numDevices = gpu_device_get_num_all();
  const uint32_t numSupportedDevices = gpu_get_num_supported_devices_(selectedArchitectures);
  uint32_t idDevice, idSupportedDevice;

  // Lists initialization for the best fitting modules exploration
  gpu_module_t allocatedModulesPerDevice[numSupportedDevices];
  gpu_module_t allocatedStructuresPerDevice[numSupportedDevices];
  gpu_module_t activatedModules = GPU_NONE_MODULES, allocatedStructures = GPU_NONE_MODULES;

  // Module lists initialization
  for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
    allocatedModulesPerDevice[idSupportedDevice]    = activatedModules;
    allocatedStructuresPerDevice[idSupportedDevice] = allocatedStructures;
  }

  // Analyze all the devices on the system and annotate the module requirements
  for(idDevice = 0, idSupportedDevice = 0; idDevice < numDevices; ++idDevice){
    const bool deviceArchSupported = gpu_device_get_architecture(idDevice) & selectedArchitectures;
    activatedModules = GPU_NONE_MODULES; allocatedStructures = GPU_NONE_MODULES;
    if(deviceArchSupported){
      GPU_ERROR(gpu_module_manager_per_device(reference, index, idDevice, numBuffers, userAllocOption, &activatedModules, &allocatedStructures));
      allocatedModulesPerDevice[idSupportedDevice]    = activatedModules;
      allocatedStructuresPerDevice[idSupportedDevice] = allocatedStructures;
      idSupportedDevice++;
    }
  }

  // Module exploration to define the module and structures configuration for all the system
  GPU_ERROR(gpu_module_search_active(allocatedModulesPerDevice, numSupportedDevices, &activatedModules));
  GPU_ERROR(gpu_module_search_structures(allocatedModulesPerDevice, allocatedStructuresPerDevice, numSupportedDevices,
                                         activatedModules, &allocatedStructures));

  (* globalModules)    = activatedModules;
  (* globalStructures) = allocatedStructures;
  return(SUCCESS);
}

gpu_error_t gpu_module_manager_memory_policies(const gpu_data_location_t userAllocOption, const gpu_module_t userRequestedModules,
                                               gpu_module_t* const requiredModules, gpu_module_t* const recomendedModules,
                                               gpu_module_t* const globalStructures)
{
  switch (userAllocOption){
    case GPU_REMOTE_DATA:          // (Force to allocate all structures in HOST)
      (* globalStructures)  = GPU_NONE_MODULES;
      (* requiredModules)   = GPU_NONE_MODULES;
      (* recomendedModules) = userRequestedModules;
      break;
    case GPU_LOCAL_DATA:           // (Force to allocate all structures in DEVICE)
      (* requiredModules)   = userRequestedModules;
      (* recomendedModules) = userRequestedModules;
      break;
    case GPU_LOCAL_OR_REMOTE_DATA: // (Best effort allocating structures in DEVICE and the rest in HOST)
      (* requiredModules)   = GPU_NONE_MODULES;
      (* recomendedModules) = userRequestedModules;
      break;
    case GPU_GEM_POLICY:           // (GEM Default Allocation)
      (* requiredModules)   = (* globalStructures); // Reference is the minimum required module to be execute
      (* recomendedModules) = userRequestedModules;
      break;
    default:
      return(E_NOT_SUPPORTED_ALLOC_POLICY);
  }
  return(SUCCESS);
}

gpu_error_t gpu_module_manager_memory(gpu_reference_buffer_t* const reference, gpu_index_buffer_t* const index,
                                      const uint32_t numBuffers, const gpu_dev_arch_t selectedArchitectures,
                                      const gpu_data_location_t userAllocOption,
                                      size_t* const recomendedMemorySize, size_t* const requiredMemorySize,
                                      gpu_module_t* const modules, gpu_module_t* const structures)
{
  const gpu_module_t userRequestedModules = index->activeModules | reference->activeModules;
  // Defines to obtain the module requirements
  gpu_module_t requiredModules = GPU_NONE_MODULES, recomendedModules = GPU_NONE_MODULES;
  gpu_module_t localModules    = GPU_NONE_MODULES, localStructures   = GPU_NONE_MODULES;
  size_t localRecomendedMemorySize, localRequiredMemorySize;

  // Explores the best combination of active modules and stored structures for all the devices
  GPU_ERROR(gpu_module_manager_all_system(reference, index, numBuffers, selectedArchitectures, userAllocOption, &localModules, &localStructures));
  // Calculates the device memory constrains for the above configuration
  GPU_ERROR(gpu_module_manager_memory_policies(userAllocOption, userRequestedModules, &requiredModules, &recomendedModules, &localStructures));
  GPU_ERROR(gpu_module_get_min_memory(reference, index, numBuffers, requiredModules, &localRequiredMemorySize));
  GPU_ERROR(gpu_module_get_min_memory(reference, index, numBuffers, recomendedModules, &localRecomendedMemorySize));

  (* modules)              = localModules;
  (* structures)           = localStructures;
  (* recomendedMemorySize) = localRecomendedMemorySize;
  (* requiredMemorySize)   = localRequiredMemorySize;

  return(SUCCESS);
}

gpu_error_t gpu_module_configure_system(gpu_reference_buffer_t* const reference, gpu_index_buffer_t* const index,
                                        gpu_device_info_t ***devices, const uint32_t numBuffers,
                                        const gpu_dev_arch_t selectedArchitectures, const gpu_data_location_t userAllocOption,
                                        gpu_module_t* const activatedModules, gpu_module_t* const allocatedStructures)
{
  // Defines for the system devices
  const uint32_t numDevices = gpu_device_get_num_all();
  const uint32_t numSupportedDevices = gpu_get_num_supported_devices_(selectedArchitectures);
  gpu_device_info_t **dev = (gpu_device_info_t **) malloc(numSupportedDevices * sizeof(gpu_device_info_t *));
  uint32_t idDevice, idSupportedDevice;
  // Defines for the module requirements
  gpu_module_t globalModules, globalStructures;
  size_t recomendedMemorySize, requiredMemorySize;

  // Choose the best modules-structures and the memory requirements for all the system
  GPU_ERROR(gpu_module_manager_memory(reference, index, numBuffers, selectedArchitectures, userAllocOption,
                                      &recomendedMemorySize, &requiredMemorySize, &globalModules, &globalStructures));

  // Activates the system devices with enough arch and memory requirements
  for(idDevice = 0, idSupportedDevice = 0; idDevice < numDevices; ++idDevice){
    size_t memoryFree               = gpu_device_get_free_memory(idDevice);
    const bool deviceArchSupported  = gpu_device_get_architecture(idDevice) & selectedArchitectures;
    const bool dataFitsMemoryDevice = requiredMemorySize < memoryFree;
    gpu_device_screen_status(idDevice, deviceArchSupported, recomendedMemorySize, requiredMemorySize);
    if(deviceArchSupported && dataFitsMemoryDevice){ //Data fits on memory device
      GPU_ERROR(gpu_device_init(&dev[idDevice], idDevice, idSupportedDevice, selectedArchitectures));
      GPU_ERROR(gpu_module_set_device_allocation(reference, index, idSupportedDevice, globalStructures));
      idSupportedDevice++;
    }
  }

  // Analyze and record the characteristic of each supported device
  GPU_ERROR(gpu_device_characterize_all(dev, idSupportedDevice));

  reference->activeModules = globalModules & GPU_REFERENCE;
  index->activeModules     = globalModules & GPU_INDEX;
  (*activatedModules)      = globalModules;
  (*allocatedStructures)   = globalStructures;

  (* devices)              = dev;
  return(SUCCESS);
}

#endif /* GPU_MODULE_C_ */
