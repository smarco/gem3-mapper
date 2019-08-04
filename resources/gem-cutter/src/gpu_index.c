/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2011-2018 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_INDEX_C_
#define GPU_INDEX_C_

#include "../include/gpu_index.h"
#include "../include/gpu_io.h"


/************************************************************
Get information functions
************************************************************/

gpu_error_t gpu_index_get_size(const gpu_index_buffer_t* const index, size_t *bytesPerIndex, const gpu_module_t activeModules)
{
  (* bytesPerIndex) = 0;

  if(activeModules & GPU_FMI){
    size_t bytesPerFMI = 0, bytesPerFmiTable = 0;
    GPU_ERROR(gpu_fmi_index_get_size(&index->fmi, &bytesPerFMI));
    (* bytesPerIndex) += bytesPerFMI;
    GPU_ERROR(gpu_fmi_table_get_size(&index->fmi.table, &bytesPerFmiTable));
    (* bytesPerIndex) += bytesPerFmiTable;
  }

  if(activeModules & GPU_SA){
    size_t bytesPerSA = 0;
    GPU_ERROR(gpu_sa_index_get_size(&index->sa, &bytesPerSA));
    (* bytesPerIndex) += bytesPerSA;
  }

  return (SUCCESS);
}


/************************************************************
INPUT / OUPUT Functions
************************************************************/

gpu_error_t gpu_index_read_specs(int fp, gpu_index_buffer_t* const index, const gpu_module_t activeModules)
{
  if(activeModules & GPU_FMI){
    GPU_ERROR(gpu_fmi_index_read_specs(fp, &index->fmi));
    GPU_ERROR(gpu_fmi_table_read_specs(fp, &index->fmi.table));
  }
  if(activeModules & GPU_SA){
    GPU_ERROR(gpu_sa_index_read_specs(fp, &index->sa));
  }
  return (SUCCESS);
}

gpu_error_t gpu_index_read(int fp, gpu_index_buffer_t* const index, const gpu_module_t activeModules)
{
  if(activeModules & GPU_FMI){
    GPU_ERROR(gpu_fmi_index_read_specs(fp, &index->fmi));
    GPU_ERROR(gpu_fmi_table_read_specs(fp, &index->fmi.table));
    GPU_ERROR(gpu_fmi_index_read(fp, &index->fmi));
    GPU_ERROR(gpu_fmi_table_read(fp, &index->fmi.table));
  }
  if(activeModules & GPU_SA){
    GPU_ERROR(gpu_sa_index_read_specs(fp, &index->sa));
    GPU_ERROR(gpu_sa_index_read(fp, &index->sa));
  }
  return (SUCCESS);
}

gpu_error_t gpu_index_write_specs(int fp, const gpu_index_buffer_t* const index, const gpu_module_t activeModules)
{
  if(activeModules & GPU_FMI){
    GPU_ERROR(gpu_fmi_index_write_specs(fp, &index->fmi));
    GPU_ERROR(gpu_fmi_table_write_specs(fp, &index->fmi.table));
  }
  if(activeModules & GPU_SA){
    GPU_ERROR(gpu_sa_index_write_specs(fp, &index->sa));
  }
  return(SUCCESS);
}

gpu_error_t gpu_index_write(int fp, const gpu_index_buffer_t* const index, const gpu_module_t activeModules)
{
  if(activeModules & GPU_FMI){
    GPU_ERROR(gpu_fmi_index_write_specs(fp, &index->fmi));
    GPU_ERROR(gpu_fmi_table_write_specs(fp, &index->fmi.table));
    GPU_ERROR(gpu_fmi_index_write(fp, &index->fmi));
    GPU_ERROR(gpu_fmi_table_write(fp, &index->fmi.table));
  }
  if(activeModules & GPU_SA){
    GPU_ERROR(gpu_sa_index_write_specs(fp, &index->sa));
    GPU_ERROR(gpu_sa_index_write(fp, &index->sa));
  }
  return(SUCCESS);
}


/************************************************************
Functions to transference the index (HOST <-> DEVICES)
************************************************************/

gpu_error_t gpu_index_transfer_CPU_to_GPUs(gpu_index_buffer_t* const index, gpu_device_info_t** const devices, const gpu_module_t activeModules)
{
  if(activeModules & GPU_FMI){
    GPU_ERROR(gpu_fmi_index_transfer_CPU_to_GPUs(&index->fmi, devices));
    GPU_ERROR(gpu_fmi_table_transfer_CPU_to_GPUs(&index->fmi.table, devices));
  }
  if(activeModules & GPU_SA){
    GPU_ERROR(gpu_sa_index_transfer_CPU_to_GPUs(&index->sa, devices));
  }
  return (SUCCESS);
}


/************************************************************
Functions to transform the index
************************************************************/

gpu_error_t gpu_index_transform_ASCII(const gpu_index_dto_t* const textRaw, gpu_index_buffer_t* const index, const gpu_module_t activeModules)
{
  if(activeModules & GPU_FMI){
    GPU_ERROR(gpu_fmi_index_transform_ASCII(textRaw->fmi.h_plain, &index->fmi));
    GPU_ERROR(gpu_fmi_table_construction(&index->fmi.table, index->fmi.h_fmi, index->fmi.bwtSize));
  }
  if(activeModules & GPU_SA){
    GPU_ERROR(gpu_sa_index_transform_ASCII(textRaw->sa.h_plain, &index->sa));
  }
  return (SUCCESS);
}

gpu_error_t gpu_index_transform_GEM_FULL(const gpu_index_dto_t* const indexRaw, gpu_index_buffer_t* const index, const gpu_module_t activeModules)
{
  if(activeModules & GPU_FMI){
    GPU_ERROR(gpu_fmi_index_transform_GEM_FULL((gpu_gem_fmi_dto_t*)indexRaw, &index->fmi));
    GPU_ERROR(gpu_fmi_table_construction(&index->fmi.table, index->fmi.h_fmi, index->fmi.bwtSize));
  }
  if(activeModules & GPU_SA){
    GPU_ERROR(gpu_sa_index_transform_GEM_FULL((gpu_gem_sa_dto_t*)indexRaw, &index->sa));
  }
  return (SUCCESS);
}

gpu_error_t gpu_index_load_specs_MFASTA_FULL(const gpu_index_dto_t* const indexRaw, gpu_index_buffer_t* const index, const gpu_module_t activeModules)
{
  const char* const filename = indexRaw->filename;

  if(activeModules & GPU_FMI){
    GPU_ERROR(gpu_fmi_index_load_specs_MFASTA_FULL(filename, &index->fmi));
    GPU_ERROR(gpu_fmi_table_load_default_specs(&index->fmi.table));
  }
  if(activeModules & GPU_SA){
    GPU_ERROR(gpu_sa_index_load_specs_MFASTA_FULL(filename, &index->sa));
  }

  return (SUCCESS);
}

gpu_error_t gpu_index_transform_MFASTA_FULL(const gpu_index_dto_t* const indexRaw, gpu_index_buffer_t* const index, const gpu_module_t activeModules)
{
  const char* const filename = indexRaw->filename;

  if(activeModules & GPU_FMI){
    GPU_ERROR(gpu_fmi_index_transform_MFASTA_FULL(filename, &index->fmi));
    GPU_ERROR(gpu_fmi_table_construction(&index->fmi.table, index->fmi.h_fmi, index->fmi.bwtSize));
  }
  if(activeModules & GPU_SA){
    GPU_ERROR(gpu_sa_index_transform_MFASTA_FULL(filename, &index->sa));
  }

  return (SUCCESS);
}

gpu_error_t gpu_index_set_specs(gpu_index_buffer_t* const index, const gpu_index_dto_t* const indexRaw,
                                const gpu_index_coding_t indexCoding, const gpu_module_t activeModules)
{
  const char* const filename = indexRaw->filename;

  switch(indexCoding){
    case GPU_INDEX_ASCII:
      /* Not need special I/O initialization */
      break;
    case GPU_INDEX_GEM_FULL:
      /* Not need special I/O initialization */
      break;
    case GPU_INDEX_GEM_FILE:
      GPU_ERROR(gpu_io_load_index_specs_GEM_FULL(filename, index, activeModules));
      break;
    case GPU_INDEX_MFASTA_FILE:
      GPU_ERROR(gpu_index_load_specs_MFASTA_FULL(indexRaw, index, activeModules));
      break;
    case GPU_INDEX_PROFILE_FILE:
      GPU_ERROR(gpu_io_load_index_specs_PROFILE(filename, index, activeModules));
      break;
    default:
      GPU_ERROR(E_INDEX_CODING);
    break;
  }

  return (SUCCESS);
}

gpu_error_t gpu_index_transform(gpu_index_buffer_t* const index, const gpu_index_dto_t* const indexRaw,
                                const gpu_index_coding_t indexCoding, const gpu_module_t activeModules)
{
  const char* const filename = indexRaw->filename;

  switch(indexCoding){
    case GPU_INDEX_ASCII:
      GPU_ERROR(gpu_index_transform_ASCII(indexRaw, index, activeModules));
      break;
    case GPU_INDEX_GEM_FULL:
      GPU_ERROR(gpu_index_transform_GEM_FULL(indexRaw, index, activeModules));
      break;
    case GPU_INDEX_GEM_FILE:
      GPU_ERROR(gpu_io_load_index_GEM_FULL(filename, index, activeModules));
      //GPU_ERROR(gpu_save_index_PROFILE("internalIndexGEM", index, activeModules)); //DEBUG: backup the index
      break;
    case GPU_INDEX_MFASTA_FILE:
      GPU_ERROR(gpu_index_transform_MFASTA_FULL(indexRaw, index, activeModules));
      break;
    case GPU_INDEX_PROFILE_FILE:
      GPU_ERROR(gpu_io_load_index_PROFILE(filename, index, activeModules));
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

gpu_error_t gpu_index_init_dto(gpu_index_buffer_t *index, const gpu_module_t activeModules)
{
  if(activeModules & GPU_FMI){
    GPU_ERROR(gpu_fmi_index_init_dto(&index->fmi));
    GPU_ERROR(gpu_fmi_table_init_dto(&index->fmi.table));
  }
  if(activeModules & GPU_SA){
    GPU_ERROR(gpu_sa_index_init_dto(&index->sa));
  }

  return (SUCCESS);
}

gpu_error_t gpu_index_init(gpu_index_buffer_t** const index, const gpu_index_dto_t* const rawIndex,
                           const uint32_t numSupportedDevices,const gpu_module_t activeModules)
{
  gpu_index_buffer_t* const iBuff = (gpu_index_buffer_t *) malloc(sizeof(gpu_index_buffer_t));

  //Initialize the the active index modules
  iBuff->activeModules = activeModules & GPU_INDEX;

  GPU_ERROR(gpu_fmi_index_init_dto(&iBuff->fmi));
  GPU_ERROR(gpu_fmi_table_init_dto(&iBuff->fmi.table));
  GPU_ERROR(gpu_fmi_index_init(&iBuff->fmi, iBuff->fmi.bwtSize, numSupportedDevices));
  GPU_ERROR(gpu_fmi_table_init(&iBuff->fmi.table, iBuff->fmi.table.maxLevelsTableLUT, numSupportedDevices));
  GPU_ERROR(gpu_sa_index_init(&iBuff->sa, iBuff->sa.numEntries, iBuff->sa.sampligRate, numSupportedDevices));

  if(activeModules & GPU_FMI) GPU_ERROR(gpu_index_set_specs(iBuff, rawIndex, rawIndex->fmi.indexCoding, GPU_FMI));
  if(activeModules & GPU_SA)  GPU_ERROR(gpu_index_set_specs(iBuff, rawIndex, rawIndex->sa.indexCoding, GPU_SA));

  (* index) = iBuff;
  return (SUCCESS);
}

gpu_error_t gpu_index_allocate(gpu_index_buffer_t *index, const gpu_module_t activeModules)
{
  if(activeModules & GPU_FMI){
    GPU_ERROR(gpu_fmi_index_allocate(&index->fmi));
    GPU_ERROR(gpu_fmi_table_allocate(&index->fmi.table));
  }
  if(activeModules & GPU_SA){
    GPU_ERROR(gpu_sa_index_allocate(&index->sa));
  }

  return (SUCCESS);
}

gpu_error_t gpu_index_load(gpu_index_buffer_t *index, const gpu_index_dto_t * const rawIndex,
                           const gpu_module_t activeModules)
{
  if(activeModules & GPU_FMI){
    GPU_ERROR(gpu_fmi_index_allocate(&index->fmi));
    GPU_ERROR(gpu_fmi_table_allocate(&index->fmi.table));
    GPU_ERROR(gpu_index_transform(index, rawIndex, rawIndex->fmi.indexCoding, GPU_FMI));
  }

  if(activeModules & GPU_SA){
    GPU_ERROR(gpu_sa_index_allocate(&index->sa));
    GPU_ERROR(gpu_index_transform(index, rawIndex, rawIndex->sa.indexCoding, GPU_SA));
  }

  return (SUCCESS);
}


/************************************************************
 Functions to release the index data from the DEVICE & HOST
************************************************************/

gpu_error_t gpu_index_free_host(gpu_index_buffer_t* const index, const gpu_module_t activeModules)
{
  if(activeModules & GPU_FMI){
    GPU_ERROR(gpu_fmi_index_free_host(&index->fmi));
    GPU_ERROR(gpu_fmi_table_free_host(&index->fmi.table));
  }
  if(activeModules & GPU_SA){
    GPU_ERROR(gpu_sa_index_free_host(&index->sa));
  }
  return(SUCCESS);
}

gpu_error_t gpu_index_free_unused_host(gpu_index_buffer_t* index, gpu_device_info_t** const devices, const gpu_module_t activeModules)
{
  if(activeModules & GPU_FMI){
    GPU_ERROR(gpu_fmi_index_free_unused_host(&index->fmi, devices));
    GPU_ERROR(gpu_fmi_table_free_unused_host(&index->fmi.table, devices));
  }
  if(activeModules & GPU_SA){
    GPU_ERROR(gpu_sa_index_free_unused_host(&index->sa, devices));
  }
  return(SUCCESS);
}

gpu_error_t gpu_index_free_device(gpu_index_buffer_t* index, gpu_device_info_t** const devices, const gpu_module_t activeModules)
{
  if(activeModules & GPU_FMI){
    GPU_ERROR(gpu_fmi_index_free_device(&index->fmi, devices));
    GPU_ERROR(gpu_fmi_table_free_device(&index->fmi.table, devices));
  }
  if(activeModules & GPU_SA){
    GPU_ERROR(gpu_sa_index_free_device(&index->sa, devices));
  }
  return(SUCCESS);
}

gpu_error_t gpu_index_free_metainfo(gpu_index_buffer_t* index, const gpu_module_t activeModules)
{
  if(activeModules & GPU_FMI){
    GPU_ERROR(gpu_fmi_index_free_metainfo(&index->fmi));
    GPU_ERROR(gpu_fmi_table_free_metainfo(&index->fmi.table));
  }
  if(activeModules & GPU_SA){
    GPU_ERROR(gpu_sa_index_free_metainfo(&index->sa));
  }
  return(SUCCESS);
}

gpu_error_t gpu_index_free(gpu_index_buffer_t **index, gpu_device_info_t** const devices, const gpu_module_t activeModules)
{
  gpu_index_buffer_t* iBuff = (* index);

  GPU_ERROR(gpu_index_free_host(iBuff, activeModules));
  GPU_ERROR(gpu_index_free_device(iBuff, devices, activeModules));

  if(iBuff != NULL){
    GPU_ERROR(gpu_index_free_metainfo(iBuff, activeModules));
    free(iBuff);
    iBuff = NULL;
  }

  (* index) = iBuff;
  return(SUCCESS);
}

#endif /* GPU_INDEX_C_ */



