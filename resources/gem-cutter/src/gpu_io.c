/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2011-2018 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_IO_C_
#define GPU_IO_C_

#include "../include/gpu_io.h"

/************************************************************
Basic primitives for input/output
************************************************************/

gpu_error_t gpu_io_read_buffered(int fp, void* const buffer, const size_t bytesRequest)
{
  const uint32_t numRequests  = GPU_DIV_CEIL(bytesRequest, GPU_FILE_SIZE_BLOCK);
  size_t result, numBytesRequested = 0;
  uint32_t idRequest;
  for(idRequest = 0; idRequest < numRequests; ++idRequest){
    const size_t requestSize = GPU_MIN(GPU_FILE_SIZE_BLOCK, bytesRequest - numBytesRequested);
    result = read(fp, buffer + numBytesRequested, requestSize);
    if (result != requestSize) return (E_READING_FILE);
    numBytesRequested += requestSize;
  }
  // Succeed
  return (SUCCESS);
}

gpu_error_t gpu_io_write_buffered(int fp, void* const buffer, const size_t bytesRequest)
{
  const uint32_t numRequests  = GPU_DIV_CEIL(bytesRequest, GPU_FILE_SIZE_BLOCK);
  size_t result, numBytesRequested = 0;
  uint32_t idRequest;
  for(idRequest = 0; idRequest < numRequests; ++idRequest){
    const size_t requestSize = GPU_MIN(GPU_FILE_SIZE_BLOCK, bytesRequest - numBytesRequested);
    result = write(fp, buffer + numBytesRequested, requestSize);
    if (result != requestSize) return (E_WRITING_FILE);
    numBytesRequested += requestSize;
  }
  // Succeed
  return (SUCCESS);
}

/************************************************************
Primitives for input/output
************************************************************/

gpu_error_t gpu_io_load_specs_BWT_MFASTA(const char* const fn, gpu_index_buffer_t* const index,
										                     const gpu_module_t activeModules)
{
  FILE *fp = NULL;
  uint64_t sizeFile = 0;

  if (activeModules & GPU_FMI){
    fp = fopen(fn, "rb");
    if (fp == NULL) return (E_OPENING_FILE);

    fseek(fp, 0L, SEEK_END);
    sizeFile = ftell(fp);

    index->fmi.bwtSize    = sizeFile;
    index->fmi.numEntries = GPU_DIV_CEIL(index->fmi.bwtSize, GPU_FMI_ENTRY_SIZE) + 1;

    fclose(fp);
  }

  if (activeModules & GPU_SA){
    return(E_NOT_IMPLEMENTED);
  }

  return (SUCCESS);
}

gpu_error_t gpu_io_load_BWT_MFASTA(const char* const fn, gpu_index_buffer_t* const index, char **h_BWT)
{
  FILE *fp = NULL;
  char *h_ascii_BWT = NULL;
  char lineFile[GPU_FILE_SIZE_LINES];
  uint64_t sizeFile = 0, position = 0;
  int32_t charsRead = 0;

  fp = fopen(fn, "rb");
  if (fp == NULL) return (E_OPENING_FILE);

  fseek(fp, 0L, SEEK_END);
  sizeFile = ftell(fp);
  rewind(fp);

  h_ascii_BWT = (char*) malloc(sizeFile * sizeof(char));
  if (h_ascii_BWT == NULL) return (E_ALLOCATE_MEM);

  while((!feof(fp)) && (fgets(lineFile, GPU_FILE_SIZE_LINES, fp) != NULL)){
    if (lineFile[0] != '>'){
      charsRead = strlen(lineFile);
      if(charsRead) charsRead--;
      memcpy((h_ascii_BWT + position), lineFile, charsRead);
      position +=  charsRead;
    }
  }

  index->fmi.bwtSize    = position;
  index->fmi.numEntries = GPU_DIV_CEIL(index->fmi.bwtSize, GPU_FMI_ENTRY_SIZE) + 1;
  (* h_BWT)             = h_ascii_BWT;

  fclose(fp);
  return (SUCCESS);
}

gpu_error_t gpu_io_save_index_PROFILE(const char* const fn, const gpu_index_buffer_t* const index,
									                    const gpu_module_t activeModules)
{
  const uint32_t sizeFileName = 512;
  char fileName[sizeFileName];
  int fp = 0, openMode = GPU_FILE_BASIC_MODE | O_WRONLY;

  if((index->activeModules & activeModules) == 0)
    return(E_MODULE_NOT_FOUND);

  if (activeModules & GPU_SA){
    sprintf(fileName, "%s.%lu.%lu.sa", fn, index->sa.numEntries, index->sa.sampligRate);
    fp = open(fileName, openMode, GPU_FILE_PERMISIONS);
    if (fp < 0) return (E_WRITING_FILE);
    GPU_ERROR(gpu_index_write(fp, index, GPU_SA));
    close(fp);
  }

  if(activeModules & GPU_FMI){
    sprintf(fileName, "%s.%lu.%u.fmi", fn, index->fmi.bwtSize, GPU_FMI_ENTRY_SIZE);
    fp = open(fileName, openMode, GPU_FILE_PERMISIONS);
    if (fp < 0) return (E_WRITING_FILE);
    GPU_ERROR(gpu_index_write(fp, index, GPU_FMI));
    close(fp);
  }

  return (SUCCESS);
}

gpu_error_t gpu_io_load_index_specs_PROFILE(const char* const fn, gpu_index_buffer_t* const index,
											                      const gpu_module_t activeModules)
{
  int fp = 0, openMode = GPU_FILE_BASIC_MODE | O_RDONLY;

  if((activeModules & GPU_INDEX) == 0)
    return(E_MODULE_NOT_FOUND);

  fp = open(fn, openMode, GPU_FILE_PERMISIONS);
  if (fp < 0) return (E_OPENING_FILE);

  GPU_ERROR(gpu_index_read_specs(fp, index, activeModules));

  close(fp);
  return (SUCCESS);
}

gpu_error_t gpu_io_load_index_PROFILE(const char* const fn, gpu_index_buffer_t* const index,
									  const gpu_module_t activeModules)
{
  int fp = 0, openMode = GPU_FILE_BASIC_MODE | O_RDONLY;

  if((activeModules & GPU_INDEX) == 0)
    return(E_MODULE_NOT_FOUND);

  fp = open(fn, openMode, GPU_FILE_PERMISIONS);
  if (fp < 0) return (E_OPENING_FILE);

  GPU_ERROR(gpu_index_read(fp, index, activeModules));

  close(fp);
  return (SUCCESS);
}

gpu_error_t gpu_io_load_reference_specs_MFASTA(const char* const fn, gpu_reference_buffer_t* const reference,
											   const gpu_module_t activeModules)
{
  FILE *fp = NULL;
  uint64_t sizeFile = 0;

  fp = fopen(fn, "rb");
  if (fp == NULL) return (E_OPENING_FILE);

  fseek(fp, 0L, SEEK_END);
  sizeFile = ftell(fp);

  reference->size = sizeFile;
  reference->numEntriesPlain = GPU_DIV_CEIL(reference->size, GPU_REFERENCE_PLAIN__CHARS_PER_ENTRY) + GPU_REFERENCE_END_PADDING;

  fclose(fp);
  return (SUCCESS);
}

gpu_error_t gpu_io_load_reference_MFASTA(const char* const fn, gpu_reference_buffer_t* const reference,
									     const gpu_module_t activeModules)
{
  FILE *fp = NULL;
  char lineFile[GPU_FILE_SIZE_LINES], *tmp_reference;
  uint64_t sizeFile = 0, position = 0;
  int32_t charsRead = 0;

  fp = fopen(fn, "rb");
  if (fp == NULL) return (E_OPENING_FILE);

  fseek(fp, 0L, SEEK_END);
  sizeFile = ftell(fp);
  rewind(fp);

  tmp_reference = (char*) malloc(sizeFile * sizeof(char));
  if (tmp_reference == NULL) return (E_ALLOCATE_MEM);

  if ((fgets(lineFile, GPU_FILE_SIZE_LINES, fp) == NULL) || (lineFile[0] != '>'))
    return (E_READING_FILE);

  while((!feof(fp)) && (fgets(lineFile, GPU_FILE_SIZE_LINES, fp) != NULL)){
    if (lineFile[0] != '>'){
      charsRead = strlen(lineFile);
      if(charsRead) charsRead--;
      memcpy((tmp_reference + position), lineFile, charsRead);
      position +=  charsRead;
    }
  }

  reference->size = position;
  reference->numEntriesPlain = GPU_DIV_CEIL(reference->size, GPU_REFERENCE_PLAIN__CHARS_PER_ENTRY) + GPU_REFERENCE_END_PADDING;
  GPU_ERROR(gpu_reference_transform(reference, tmp_reference, GPU_REF_ASCII, activeModules));

  fclose(fp);
  free(tmp_reference);
  return (SUCCESS);
}

gpu_error_t gpu_io_load_reference_specs_PROFILE(const char* const fn, gpu_reference_buffer_t* const reference,
												const gpu_module_t activeModules)
{
  int fp = 0, openMode = GPU_FILE_BASIC_MODE | O_RDONLY;

  if((activeModules & GPU_REFERENCE) == 0)
    return(E_MODULE_NOT_FOUND);

  fp = open(fn, openMode, GPU_FILE_PERMISIONS);
  if (fp < 0) return (E_OPENING_FILE);

  GPU_ERROR(gpu_reference_read_specs(fp, reference, activeModules));

  close(fp);
  return (SUCCESS);
}

gpu_error_t gpu_io_load_reference_PROFILE(const char* const fn, gpu_reference_buffer_t* const reference,
										  const gpu_module_t activeModules)
{
  int fp = 0, openMode = GPU_FILE_BASIC_MODE | O_RDONLY;

  if((activeModules & GPU_REFERENCE) == 0)
    return(E_MODULE_NOT_FOUND);

  fp = open(fn, openMode, GPU_FILE_PERMISIONS);
  if (fp < 0) return (E_OPENING_FILE);

  GPU_ERROR(gpu_reference_read(fp, reference, activeModules));

  close(fp);
  return (SUCCESS);
}

gpu_error_t gpu_io_save_reference_PROFILE(const char* const fn, const gpu_reference_buffer_t* const reference,
										  const gpu_module_t activeModules)
{
  int fp = 0, openMode = GPU_FILE_BASIC_MODE | O_WRONLY;

  if((activeModules & GPU_REFERENCE) == 0)
    return(E_MODULE_NOT_FOUND);

  fp = open(fn, openMode, GPU_FILE_PERMISIONS);
  if (fp < 0) return (E_OPENING_FILE);

  GPU_ERROR(gpu_reference_write(fp, reference, activeModules));

  close(fp);
  return (SUCCESS);
}

gpu_error_t gpu_io_load_index_specs_GEM_FULL(const char* const fn, gpu_index_buffer_t* const index,
										     const gpu_module_t activeModules)
{
  int fp = 0, openMode = GPU_FILE_BASIC_MODE | O_RDONLY;
  off64_t currentOffset = 0, fileOffsetFMIndex = 0, fileOffsetSAIndex = 0, fileOffsetRef = 0;
  gpu_module_t fileActiveModules = GPU_NONE_MODULES;

  fp = open(fn, openMode, GPU_FILE_PERMISIONS);
  if (fp < 0) return (E_OPENING_FILE);

  GPU_ERROR(gpu_io_load_module_info_GEM_FULL(fp, &fileActiveModules));

  if((fileActiveModules & activeModules) == 0)
    return(E_MODULE_NOT_FOUND);

  GPU_ERROR(gpu_io_load_offsets_info_GEM_FULL(fp, &fileOffsetFMIndex, &fileOffsetSAIndex, &fileOffsetRef));

  if(activeModules & GPU_FMI){
    currentOffset = lseek64(fp, fileOffsetFMIndex, SEEK_SET);
    if (currentOffset < 0) return (E_READING_FILE);
    GPU_ERROR(gpu_index_read_specs(fp, index, GPU_FMI));
  }

  if(activeModules & GPU_SA){
    currentOffset = lseek64(fp, fileOffsetSAIndex, SEEK_SET);
    if (currentOffset < 0) return (E_READING_FILE);
    GPU_ERROR(gpu_index_read_specs(fp, index, GPU_SA));
  }

  // Sanity check, re-calculate the active modules
  index->activeModules |= fileActiveModules & activeModules & GPU_INDEX;

  close(fp);
  return (SUCCESS);
}

gpu_error_t gpu_io_load_index_GEM_FULL(const char* const fn, gpu_index_buffer_t* const index,
									   const gpu_module_t activeModules)
{
  int fp = 0, openMode = GPU_FILE_BASIC_MODE | O_RDONLY;
  off64_t currentOffset = 0, fileOffsetFMIndex = 0, fileOffsetSAIndex = 0, fileOffsetRef = 0;
  gpu_module_t fileActiveModules = GPU_NONE_MODULES;

  fp = open(fn, openMode, GPU_FILE_PERMISIONS);
  if (fp < 0) return (E_OPENING_FILE);

  GPU_ERROR(gpu_io_load_module_info_GEM_FULL(fp, &fileActiveModules));

  if((fileActiveModules & activeModules) == 0)
    return(E_MODULE_NOT_FOUND);

  GPU_ERROR(gpu_io_load_offsets_info_GEM_FULL(fp, &fileOffsetFMIndex, &fileOffsetSAIndex, &fileOffsetRef));

  if(activeModules & GPU_FMI){
    currentOffset = lseek64(fp, fileOffsetFMIndex, SEEK_SET);
    if (currentOffset < 0) return (E_READING_FILE);
    GPU_ERROR(gpu_index_read(fp, index, GPU_FMI));
  }

  if(activeModules & GPU_SA){
    currentOffset = lseek64(fp, fileOffsetSAIndex, SEEK_SET);
    if (currentOffset < 0) return (E_READING_FILE);
    GPU_ERROR(gpu_index_read(fp, index, GPU_SA));
  }

  // Sanity check, re-calculate the active modules
  index->activeModules |= fileActiveModules & activeModules & GPU_INDEX;

  close(fp);
  return (SUCCESS);
}

gpu_error_t gpu_io_save_index_GEM_FULL(const char* const fn, const gpu_index_buffer_t* const index,
									   const gpu_module_t activeModules)
{
  int fp = 0, openMode = GPU_FILE_BASIC_MODE | O_WRONLY;
  gpu_module_t storedModules = index->activeModules & activeModules;
  off64_t currentOffset = 0, fileOffsetFMIndex = 0, fileOffsetSAIndex = 0, fileOffsetRef = 0;

  if((index->activeModules & activeModules) == 0)
    return(E_MODULE_NOT_FOUND);

  fp = open(fn, openMode, GPU_FILE_PERMISIONS);
  if (fp < 0) return (E_OPENING_FILE);

  GPU_ERROR(gpu_io_save_module_info_GEM_FULL(fp, storedModules));
  GPU_ERROR(gpu_io_save_offsets_info_GEM_FULL(fp, fileOffsetFMIndex, fileOffsetSAIndex, fileOffsetRef));

  if(storedModules & GPU_FMI){
    fileOffsetFMIndex = lseek64(fp, 0, SEEK_CUR);
    GPU_ERROR(gpu_index_write(fp, index, GPU_FMI));
  }

  if(storedModules & GPU_SA){
    fileOffsetSAIndex = lseek64(fp, 0, SEEK_CUR);
    GPU_ERROR(gpu_index_write(fp, index, GPU_SA));
  }

  //Rewind the file
  currentOffset = lseek64(fp, 0, SEEK_SET);
  if(currentOffset < 0) return (E_WRITING_FILE);

  GPU_ERROR(gpu_io_save_module_info_GEM_FULL(fp, storedModules));
  GPU_ERROR(gpu_io_save_offsets_info_GEM_FULL(fp, fileOffsetFMIndex, fileOffsetSAIndex, fileOffsetRef));

  close(fp);
  return (SUCCESS);
}

gpu_error_t gpu_io_load_reference_specs_GEM_FULL(const char* const fn, gpu_reference_buffer_t* const reference,
												 const gpu_module_t activeModules)
{
  int fp = 0, openMode = GPU_FILE_BASIC_MODE | O_RDONLY;
  off64_t currentOffset = 0, fileOffsetFMIndex = 0, fileOffsetSAIndex = 0, fileOffsetRef = 0;
  gpu_module_t fileActiveModules = GPU_NONE_MODULES;

  fp = open(fn, openMode, GPU_FILE_PERMISIONS);
  if (fp < 0) return(E_OPENING_FILE);

  GPU_ERROR(gpu_io_load_module_info_GEM_FULL(fp, &fileActiveModules));

  if((fileActiveModules & activeModules) == 0)
    return (E_MODULE_NOT_FOUND);

  reference->activeModules = fileActiveModules & activeModules & GPU_REFERENCE;

  GPU_ERROR(gpu_io_load_offsets_info_GEM_FULL(fp, &fileOffsetFMIndex, &fileOffsetSAIndex, &fileOffsetRef));

  currentOffset = lseek64(fp, fileOffsetRef, SEEK_SET);
  if (currentOffset < 0) return(E_READING_FILE);

  GPU_ERROR(gpu_reference_read_specs(fp, reference, GPU_REFERENCE));

  close(fp);
  return (SUCCESS);
}

gpu_error_t gpu_io_load_reference_GEM_FULL(const char* const fn, gpu_reference_buffer_t* const reference,
										   const gpu_module_t activeModules)
{
  int fp = 0, openMode = GPU_FILE_BASIC_MODE | O_RDONLY;
  off64_t currentOffset = 0, fileOffsetFMIndex = 0, fileOffsetSAIndex = 0, fileOffsetRef = 0;
  gpu_module_t fileActiveModules = GPU_NONE_MODULES;

  fp = open(fn, openMode, GPU_FILE_PERMISIONS);
  if (fp < 0) return(E_OPENING_FILE);

  GPU_ERROR(gpu_io_load_module_info_GEM_FULL(fp, &fileActiveModules));

  if((fileActiveModules & activeModules) == 0)
    return(E_MODULE_NOT_FOUND);

  reference->activeModules = fileActiveModules & activeModules & GPU_REFERENCE;

  GPU_ERROR(gpu_io_load_offsets_info_GEM_FULL(fp, &fileOffsetFMIndex, &fileOffsetSAIndex, &fileOffsetRef));

  currentOffset = lseek64(fp, fileOffsetRef, SEEK_SET);
  if (currentOffset < 0) return(E_READING_FILE);

  GPU_ERROR(gpu_reference_read(fp, reference, GPU_REFERENCE));

  close(fp);
  return (SUCCESS);
}

gpu_error_t gpu_io_save_reference_GEM_FULL(const char* const fn, const gpu_reference_buffer_t* const reference,
										   const gpu_module_t activeModules)
{
  int fp = 0, openMode = GPU_FILE_BASIC_MODE | O_WRONLY;
  off64_t currentOffset = 0, fileOffsetFMIndex = 0, fileOffsetSAIndex = 0, fileOffsetRef = 0;

  if((activeModules & GPU_REFERENCE) == 0)
    return (E_MODULE_NOT_FOUND);

  fp = open(fn, openMode, GPU_FILE_PERMISIONS);
  if (fp < 0) return (E_OPENING_FILE);

  GPU_ERROR(gpu_io_save_module_info_GEM_FULL(fp, activeModules));
  GPU_ERROR(gpu_io_save_offsets_info_GEM_FULL(fp, fileOffsetFMIndex, fileOffsetSAIndex, fileOffsetRef));

  fileOffsetRef = lseek64(fp, 0, SEEK_CUR);
  if (fileOffsetRef < 0) return(E_READING_FILE);

  GPU_ERROR(gpu_reference_write(fp, reference, GPU_REFERENCE));

  //Rewind file
  currentOffset = lseek64(fp, 0, SEEK_SET);
  if (currentOffset < 0) return(E_READING_FILE);

  GPU_ERROR(gpu_io_save_module_info_GEM_FULL(fp, activeModules));
  GPU_ERROR(gpu_io_save_offsets_info_GEM_FULL(fp, fileOffsetFMIndex, fileOffsetSAIndex, fileOffsetRef));

  close(fp);
  return (SUCCESS);
}

gpu_error_t gpu_io_save_module_info_GEM_FULL(const int fp, const gpu_module_t fileActiveModules)
{
  size_t result, bytesRequest;
  bytesRequest = sizeof(gpu_module_t);
  result = write(fp, (void *)&fileActiveModules, bytesRequest);
  if (result != bytesRequest) return (E_WRITING_FILE);
  return(SUCCESS);
}

gpu_error_t gpu_io_save_offsets_info_GEM_FULL(const int fp, const off64_t fileOffsetFMIndex,
                                              const off64_t fileOffsetSAIndex, const off64_t fileOffsetRef)
{
  size_t result, bytesRequest;
  bytesRequest = sizeof(off64_t);
  result = write(fp, (void *)&fileOffsetFMIndex, bytesRequest);
  if (result != bytesRequest) return (E_WRITING_FILE);
  bytesRequest = sizeof(off64_t);
  result = write(fp, (void *)&fileOffsetSAIndex, bytesRequest);
  if (result != bytesRequest) return (E_WRITING_FILE);
  bytesRequest = sizeof(off64_t);
  result = write(fp, (void *)&fileOffsetRef, bytesRequest);
  if (result != bytesRequest) return (E_WRITING_FILE);
  return(SUCCESS);
}

gpu_error_t gpu_io_load_module_info_GEM_FULL(const int fp, gpu_module_t* const fileActiveModules)
{
  size_t result, bytesRequest;
  bytesRequest = sizeof(gpu_module_t);
  result = read(fp, (void*)fileActiveModules, bytesRequest);
  if (result != bytesRequest) return(E_READING_FILE);
  return(SUCCESS);
}

gpu_error_t gpu_io_load_offsets_info_GEM_FULL(const int fp, off64_t* const fileOffsetFMIndex,
                                              off64_t* const fileOffsetSAIndex, off64_t* const fileOffsetRef)
{
  size_t result, bytesRequest;
  bytesRequest = sizeof(off64_t);
  result = read(fp, (void *)fileOffsetFMIndex, bytesRequest);
  if (result != bytesRequest) return(E_READING_FILE);
  bytesRequest = sizeof(off64_t);
  result = read(fp, (void *)fileOffsetSAIndex, bytesRequest);
  if (result != bytesRequest) return(E_READING_FILE);
  bytesRequest = sizeof(off64_t);
  result = read(fp, (void *)fileOffsetRef, bytesRequest);
  if (result != bytesRequest) return(E_READING_FILE);
  return(SUCCESS);
}

void gpu_io_save_indexed_structures_GEM_(const char* const fn, const gpu_gem_fmi_dto_t* const gemFMindex,
                                         const gpu_gem_ref_dto_t* const gemRef, const gpu_gem_sa_dto_t* const gemSAindex,
                                         const gpu_module_t activeModules)
{
  int fp = 0, openMode = GPU_FILE_BASIC_MODE | O_WRONLY;
  off64_t currentOffset = 0, fileOffsetFMIndex = 0, fileOffsetSAIndex = 0, fileOffsetRef = 0;

  /* Objects to initialize */
  gpu_reference_buffer_t ref;
  gpu_index_buffer_t     index;

  // Initialize the reference structure
  GPU_ERROR(gpu_reference_init_dto(&ref));
  ref.activeModules 	= activeModules & GPU_REFERENCE;
  ref.size          	= gemRef->ref_length * 2; //Forward & reverse genomes
  ref.numEntriesPlain   = GPU_DIV_CEIL(ref.size, GPU_REFERENCE_PLAIN__CHARS_PER_ENTRY) + GPU_REFERENCE_END_PADDING;
  ref.numEntriesMasked  = GPU_DIV_CEIL(ref.size, GPU_REFERENCE_MASKED__CHARS_PER_ENTRY) + GPU_REFERENCE_END_PADDING;

  // Initialize the index (SA & FMI) structure
  GPU_ERROR(gpu_index_init_dto(&index, activeModules & GPU_INDEX));
  index.activeModules                  = activeModules & GPU_INDEX;
  // Initialize the FM-Index
  index.fmi.bwtSize                    = gemFMindex->bwt_length;
  index.fmi.numEntries                 = GPU_DIV_CEIL(index.fmi.bwtSize, GPU_FMI_ENTRY_SIZE) + 1;
  // Initialize the FMI Table
  index.fmi.table.maxLevelsTableLUT    = gemFMindex->num_levels_fmi_table;
  index.fmi.table.skipLevelsTableLUT   = gemFMindex->skip_levels_fmi_table;
  index.fmi.table.occThresholdTableLUT = gemFMindex->occ_threashold_fmi_table;
  index.fmi.table.formatTableLUT       = GPU_FMI_TABLE_MULTILEVEL_LINKED;
  // Initialize the Suffix-Array index
  index.sa.sampligRate                 = gemSAindex->sa_sampling;
  index.sa.numEntries                  = GPU_DIV_CEIL(gemSAindex->sa_length, gemSAindex->sa_sampling);

  fp = open(fn, openMode, GPU_FILE_PERMISIONS);
  if (fp < 0) GPU_ERROR(E_OPENING_FILE);

  //Rewind the file (overwrites the old existing content)
  currentOffset = lseek64(fp, 0, SEEK_SET);
  if(currentOffset < 0) GPU_ERROR(E_WRITING_FILE);

  GPU_ERROR(gpu_io_save_module_info_GEM_FULL(fp, activeModules));
  GPU_ERROR(gpu_io_save_offsets_info_GEM_FULL(fp, fileOffsetFMIndex, fileOffsetSAIndex, fileOffsetRef));

  fileOffsetFMIndex = lseek64(fp, 0, SEEK_CUR);
  fileOffsetSAIndex = fileOffsetFMIndex;
  fileOffsetRef     = fileOffsetSAIndex;

  if(index.activeModules & GPU_INDEX){
    if(index.activeModules & GPU_FMI){
      GPU_ERROR(gpu_index_allocate(&index, GPU_FMI));
      GPU_ERROR(gpu_index_transform(&index, (gpu_index_dto_t*)gemFMindex, gemFMindex->index_coding, GPU_FMI));
      GPU_ERROR(gpu_index_write(fp, &index, GPU_FMI));
      GPU_ERROR(gpu_index_free_host(&index,GPU_FMI));
      //GPU_ERROR(gpu_save_index_PROFILE("internalIndexGEM", fmi)); // Dumping the index
      fileOffsetSAIndex = lseek64(fp, 0, SEEK_CUR);
      fileOffsetRef     = fileOffsetSAIndex;
    }
    if(index.activeModules & GPU_SA){
      GPU_ERROR(gpu_index_allocate(&index, GPU_SA));
      GPU_ERROR(gpu_index_transform(&index, (gpu_index_dto_t*)gemSAindex, gemSAindex->index_coding, GPU_SA));
      GPU_ERROR(gpu_index_write(fp, &index, GPU_SA));
      GPU_ERROR(gpu_index_free_host(&index,GPU_SA));
      fileOffsetRef = lseek64(fp, 0, SEEK_CUR);
    }
  }

  if(ref.activeModules & GPU_REFERENCE){
    GPU_ERROR(gpu_reference_allocate(&ref, GPU_REFERENCE));
    GPU_ERROR(gpu_reference_transform(&ref, (char*)gemRef, gemRef->ref_coding, GPU_REFERENCE));
    GPU_ERROR(gpu_reference_write(fp, &ref, GPU_REFERENCE));
    GPU_ERROR(gpu_reference_free_host(&ref));
  }

  //Rewind the file
  currentOffset = lseek64(fp, 0, SEEK_SET);
  if(currentOffset < 0) GPU_ERROR(E_WRITING_FILE);

  GPU_ERROR(gpu_io_save_module_info_GEM_FULL(fp, activeModules));
  GPU_ERROR(gpu_io_save_offsets_info_GEM_FULL(fp, fileOffsetFMIndex, fileOffsetSAIndex, fileOffsetRef));

  close(fp);
}

#endif /* GPU_IO_C_ */

