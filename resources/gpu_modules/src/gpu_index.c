#ifndef GPU_INDEX_C_
#define GPU_INDEX_C_

#include "../include/gpu_index.h"
#include "../include/gpu_io.h"

/************************************************************
INPUT / OUPUT Functions
************************************************************/

gpu_error_t gpu_read_index(FILE* fp, gpu_index_buffer_t* const index, const gpu_module_t activeModules)
{
  if(activeModules & GPU_FMI) GPU_ERROR(gpu_read_fmi_index(fp, &index->fmi));
  if(activeModules & GPU_SA)  GPU_ERROR(gpu_read_sa_index(fp, &index->sa));
  return (SUCCESS);
}

gpu_error_t gpu_write_index(FILE* fp, const gpu_index_buffer_t* const index, const gpu_module_t activeModules)
{
  if(activeModules & GPU_FMI) GPU_ERROR(gpu_write_fmi_index(fp, &index->fmi));
  if(activeModules & GPU_SA)  GPU_ERROR(gpu_write_sa_index(fp, &index->sa));
  return(SUCCESS);
}


/************************************************************
Functions to transference the index (HOST <-> DEVICES)
************************************************************/

gpu_error_t gpu_transfer_index_CPU_to_GPUs(gpu_index_buffer_t* const index, gpu_device_info_t** const devices, const gpu_module_t activeModules)
{
  if(activeModules & GPU_FMI) GPU_ERROR(gpu_transfer_fmi_index_CPU_to_GPUs(&index->fmi, devices));
  if(activeModules & GPU_SA)  GPU_ERROR(gpu_transfer_sa_index_CPU_to_GPUs(&index->sa, devices));
  return (SUCCESS);
}


/************************************************************
Functions to transform the index
************************************************************/

gpu_error_t gpu_transform_index_ASCII(const char* const textRaw, gpu_index_buffer_t* const index, const gpu_module_t activeModules)
{
  if(activeModules & GPU_FMI) GPU_ERROR(gpu_transform_fmi_index_ASCII(textRaw, &index->fmi));
  if(activeModules & GPU_SA)  GPU_ERROR(gpu_transform_sa_index_ASCII(textRaw, &index->sa));
  return (SUCCESS);
}

gpu_error_t gpu_transform_index_GEM_FULL(const char* const indexRaw, gpu_index_buffer_t* const index, const gpu_module_t activeModules)
{
  if(activeModules & GPU_FMI) GPU_ERROR(gpu_transform_fmi_index_GEM_FULL((gpu_gem_fmi_dto_t*)indexRaw, &index->fmi));
  if(activeModules & GPU_SA)  GPU_ERROR(gpu_transform_sa_index_GEM_FULL((gpu_gem_sa_dto_t*)indexRaw, &index->sa));
  return (SUCCESS);
}


gpu_error_t gpu_transform_index_MFASTA_FULL(const char* const indexRaw, gpu_index_buffer_t* const index, const gpu_module_t activeModules)
{
  if(activeModules & GPU_FMI) GPU_ERROR(gpu_transform_fmi_index_MFASTA_FULL(indexRaw, &index->fmi));
  if(activeModules & GPU_SA)  GPU_ERROR(gpu_transform_sa_index_MFASTA_FULL(indexRaw, &index->sa));
  return (SUCCESS);
}

gpu_error_t gpu_transform_index(const char* const indexRaw, gpu_index_buffer_t* const index,
                                const gpu_index_coding_t indexCoding, const gpu_module_t activeModules)
{
  switch(indexCoding){
    case GPU_INDEX_ASCII:
      GPU_ERROR(gpu_transform_index_ASCII(indexRaw, index, activeModules));
      break;
    case GPU_INDEX_GEM_FULL:
      GPU_ERROR(gpu_transform_index_GEM_FULL(indexRaw, index, activeModules));
      break;
    case GPU_INDEX_GEM_FILE:
      GPU_ERROR(gpu_load_index_GEM_FULL(indexRaw, index, activeModules));
      //GPU_ERROR(gpu_save_index_PROFILE("internalIndexGEM", index, activeModules)); //DEBUG: backup the index
      break;
    case GPU_INDEX_MFASTA_FILE:
      GPU_ERROR(gpu_transform_index_MFASTA_FULL(indexRaw, index, activeModules));
      break;
    case GPU_INDEX_PROFILE_FILE:
      GPU_ERROR(gpu_load_index_PROFILE(indexRaw, index, activeModules));
      break;
    default:
      GPU_ERROR(E_INDEX_CODING);
    break;
  }

  return(SUCCESS);
}


/************************************************************
Index initialization functions
************************************************************/

gpu_error_t gpu_init_index_dto(gpu_index_buffer_t* const index, const gpu_module_t activeModules)
{
  //Initialize the the active index modules
  index->activeModules = GPU_NONE_MODULES;
  if(activeModules & GPU_FMI) GPU_ERROR(gpu_init_fmi_index_dto(&index->fmi));
  if(activeModules & GPU_SA)  GPU_ERROR(gpu_init_sa_index_dto(&index->sa));
  return (SUCCESS);
}

gpu_error_t gpu_init_index(gpu_index_buffer_t **index, const char* const indexRaw,
                           const uint64_t bwtSize, const uint32_t samplingRate, const gpu_index_coding_t indexCoding,
                           const uint32_t numSupportedDevices, const gpu_module_t activeModules)
{
  gpu_index_buffer_t* const iBuff = (gpu_index_buffer_t *) malloc(sizeof(gpu_index_buffer_t));

  if(activeModules & GPU_FMI) GPU_ERROR(gpu_init_fmi_index(&iBuff->fmi,indexRaw,bwtSize,indexCoding,numSupportedDevices));
  if(activeModules & GPU_SA)  GPU_ERROR(gpu_init_sa_index(&iBuff->sa,indexRaw,bwtSize,samplingRate,indexCoding,numSupportedDevices));
  iBuff->activeModules = activeModules;

  GPU_ERROR(gpu_transform_index(indexRaw, iBuff, indexCoding, activeModules));

  (* index) = iBuff;
  return (SUCCESS);
}


/************************************************************
 Functions to release the index data from the DEVICE & HOST
************************************************************/

gpu_error_t gpu_free_index_host(gpu_index_buffer_t* const index, const gpu_module_t activeModules)
{
  if(activeModules & GPU_FMI) GPU_ERROR(gpu_free_fmi_index_host(&index->fmi));
  if(activeModules & GPU_SA)  GPU_ERROR(gpu_free_sa_index_host(&index->sa));
  return(SUCCESS);
}

gpu_error_t gpu_free_unused_index_host(gpu_index_buffer_t* index, gpu_device_info_t** const devices, const gpu_module_t activeModules)
{
  if(activeModules & GPU_FMI) GPU_ERROR(gpu_free_unused_fmi_index_host(&index->fmi, devices));
  if(activeModules & GPU_SA)  GPU_ERROR(gpu_free_unused_sa_index_host(&index->sa, devices));
  return(SUCCESS);
}

gpu_error_t gpu_free_index_device(gpu_index_buffer_t* index, gpu_device_info_t** const devices, const gpu_module_t activeModules)
{
  if(activeModules & GPU_FMI) GPU_ERROR(gpu_free_fmi_index_device(&index->fmi, devices));
  if(activeModules & GPU_SA)  GPU_ERROR(gpu_free_sa_index_device(&index->sa, devices));
  return(SUCCESS);
}

gpu_error_t gpu_free_index_metainfo(gpu_index_buffer_t* index, const gpu_module_t activeModules)
{
  if(activeModules & GPU_FMI) GPU_ERROR(gpu_free_fmi_index_metainfo(&index->fmi));
  if(activeModules & GPU_SA)  GPU_ERROR(gpu_free_sa_index_metainfo(&index->sa));
  return(SUCCESS);
}

gpu_error_t gpu_free_index(gpu_index_buffer_t **index, gpu_device_info_t** const devices, const gpu_module_t activeModules)
{
  gpu_index_buffer_t* iBuff = (* index);

  GPU_ERROR(gpu_free_index_host(iBuff, activeModules));
  GPU_ERROR(gpu_free_index_device(iBuff, devices, activeModules));

  if(iBuff != NULL){
    GPU_ERROR(gpu_free_index_metainfo(iBuff, activeModules));
    free(iBuff);
    iBuff = NULL;
  }

  (* index) = iBuff;
  return(SUCCESS);
}

#endif /* GPU_INDEX_C_ */



