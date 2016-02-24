/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_IO_C_
#define GPU_IO_C_

#include "../include/gpu_io.h"

/************************************************************
Primitives for input/output
************************************************************/

gpu_error_t gpu_io_load_specs_BWT_MFASTA(const char* const fn, gpu_index_buffer_t* const index, const gpu_module_t activeModules)
{
  FILE *fp = NULL;
  char *h_ascii_BWT = NULL;
  char lineFile[GPU_FILE_SIZE_LINES];
  uint64_t sizeFile = 0, position = 0;
  int32_t charsRead = 0;

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

gpu_error_t gpu_io_save_index_PROFILE(const char* const fn, const gpu_index_buffer_t* const index, const gpu_module_t activeModules)
{
  const uint32_t sizeFileName = 512;
  char fileName[sizeFileName];
  FILE *fp = NULL;

  if((index->activeModules & activeModules) == 0)
    return(E_MODULE_NOT_FOUND);

  if (activeModules & GPU_SA){
    sprintf(fileName, "%s.%lu.%u.sa", fn, index->sa.numEntries, index->sa.sampligRate);
    fp = fopen(fileName, "wb");
    if (fp == NULL) return (E_WRITING_FILE);
    GPU_ERROR(gpu_index_write(fp, index, GPU_SA));
    fclose(fp);
  }

  if(activeModules & GPU_FMI){
    sprintf(fileName, "%s.%lu.%u.fmi", fn, index->fmi.bwtSize, GPU_FMI_ENTRY_SIZE);
    fp = fopen(fileName, "wb");
    if (fp == NULL) return (E_WRITING_FILE);
    GPU_ERROR(gpu_index_write(fp, index, GPU_FMI));
    fclose(fp);
  }

  return (SUCCESS);
}

gpu_error_t gpu_io_load_index_specs_PROFILE(const char* const fn, gpu_index_buffer_t* const index, const gpu_module_t activeModules)
{
  FILE *fp = NULL;

  if((activeModules & GPU_INDEX) == 0)
    return(E_MODULE_NOT_FOUND);

  fp = fopen(fn, "rb");
  if (fp == NULL) return (E_OPENING_FILE);

  GPU_ERROR(gpu_index_read_specs(fp, index, activeModules));

  fclose(fp);
  return (SUCCESS);
}

gpu_error_t gpu_io_load_index_PROFILE(const char* const fn, gpu_index_buffer_t* const index, const gpu_module_t activeModules)
{
  FILE *fp = NULL;

  if((activeModules & GPU_INDEX) == 0)
    return(E_MODULE_NOT_FOUND);

  fp = fopen(fn, "rb");
  if (fp == NULL) return (E_OPENING_FILE);

  GPU_ERROR(gpu_index_read(fp, index, activeModules));

  fclose(fp);
  return (SUCCESS);
}

gpu_error_t gpu_io_load_reference_specs_MFASTA(const char *fn, gpu_reference_buffer_t* const reference)
{
  FILE *fp = NULL;
  uint64_t sizeFile = 0;

  fp = fopen(fn, "rb");
  if (fp == NULL) return (E_OPENING_FILE);

  fseek(fp, 0L, SEEK_END);
  sizeFile = ftell(fp);

  reference->size = sizeFile;
  reference->numEntries = GPU_DIV_CEIL(reference->size, GPU_REFERENCE_CHARS_PER_ENTRY) + GPU_REFERENCE_END_PADDING;

  fclose(fp);
  return (SUCCESS);
}

gpu_error_t gpu_io_load_reference_MFASTA(const char *fn, gpu_reference_buffer_t* const reference)
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
  reference->numEntries = GPU_DIV_CEIL(reference->size, GPU_REFERENCE_CHARS_PER_ENTRY) + GPU_REFERENCE_END_PADDING;
  GPU_ERROR(gpu_reference_transform(reference, tmp_reference, GPU_REF_ASCII, GPU_REFERENCE));

  fclose(fp);
  free(tmp_reference);
  return (SUCCESS);
}

gpu_error_t gpu_io_load_reference_specs_PROFILE(const char* const fn, gpu_reference_buffer_t* const reference)
{
  FILE *fp = NULL;
  size_t result;

  fp = fopen(fn, "rb");
  if (fp == NULL) return (E_OPENING_FILE);

  GPU_ERROR(gpu_reference_read_specs(fp, reference, GPU_REFERENCE));

  fclose(fp);
  return (SUCCESS);
}

gpu_error_t gpu_io_load_reference_PROFILE(const char* const fn, gpu_reference_buffer_t* const reference)
{
  FILE *fp = NULL;
  size_t result;

  fp = fopen(fn, "rb");
  if (fp == NULL) return (E_OPENING_FILE);

  GPU_ERROR(gpu_reference_read(fp, reference, GPU_REFERENCE));

  fclose(fp);
  return (SUCCESS);
}

gpu_error_t gpu_io_save_reference_PROFILE(const char* const fn, const gpu_reference_buffer_t* const reference)
{
  FILE *fp = NULL;
  size_t result;

  fp = fopen(fn, "wb");
  if (fp == NULL) return (E_OPENING_FILE);

  GPU_ERROR(gpu_reference_write(fp, reference, GPU_REFERENCE));

  fclose(fp);
  return (SUCCESS);
}

gpu_error_t gpu_io_load_index_specs_GEM_FULL(const char *fn, gpu_index_buffer_t* const index, const gpu_module_t activeModules)
{
  FILE *fp = NULL;
  size_t result;
  off64_t fileOffsetFMIndex = 0, fileOffsetSAIndex = 0, fileOffsetRef = 0;

  fp = fopen(fn, "rb");
  if (fp == NULL) return (E_OPENING_FILE);

  result = fread(&index->activeModules, sizeof(gpu_module_t), 1, fp);
  if (result != 1) return (E_READING_FILE);

  if((index->activeModules & activeModules) == 0)
    return(E_MODULE_NOT_FOUND);

  result = fread(&fileOffsetFMIndex, sizeof(off64_t), 1, fp);
  if (result != 1) return (E_READING_FILE);
  result = fread(&fileOffsetSAIndex, sizeof(off64_t), 1, fp);
  if (result != 1) GPU_ERROR(E_READING_FILE);
  result = fread(&fileOffsetRef, sizeof(off64_t), 1, fp);
  if (result != 1) return (E_READING_FILE);

  if(activeModules & GPU_FMI){
    result = fseeko64(fp, fileOffsetFMIndex, SEEK_SET);
    if (result != 0) return (E_READING_FILE);
    GPU_ERROR(gpu_index_read_specs(fp, index, GPU_FMI));
  }

  if(activeModules & GPU_SA){
    result = fseeko64(fp, fileOffsetSAIndex, SEEK_SET);
    if (result != 0) return (E_READING_FILE);
    GPU_ERROR(gpu_index_read_specs(fp, index, GPU_SA));
  }

  fclose(fp);
  return (SUCCESS);
}

gpu_error_t gpu_io_load_index_GEM_FULL(const char *fn, gpu_index_buffer_t* const index, const gpu_index_coding_t activeModules)
{
  FILE *fp = NULL;
  size_t result;
  off64_t fileOffsetFMIndex = 0, fileOffsetSAIndex = 0, fileOffsetRef = 0;

  fp = fopen(fn, "rb");
  if (fp == NULL) return (E_OPENING_FILE);

  result = fread(&index->activeModules, sizeof(gpu_module_t), 1, fp);
  if (result != 1) return (E_READING_FILE);

  if((index->activeModules & activeModules) == 0)
    return(E_MODULE_NOT_FOUND);

  result = fread(&fileOffsetFMIndex, sizeof(off64_t), 1, fp);
  if (result != 1) return (E_READING_FILE);
  result = fread(&fileOffsetSAIndex, sizeof(off64_t), 1, fp);
  if (result != 1) GPU_ERROR(E_READING_FILE);
  result = fread(&fileOffsetRef, sizeof(off64_t), 1, fp);
  if (result != 1) return (E_READING_FILE);

  if(activeModules & GPU_FMI){
    result = fseeko64(fp, fileOffsetFMIndex, SEEK_SET);
    if (result != 0) return (E_READING_FILE);
    GPU_ERROR(gpu_index_read(fp, index, GPU_FMI));
  }

  if(activeModules & GPU_SA){
    result = fseeko64(fp, fileOffsetSAIndex, SEEK_SET);
    if (result != 0) return (E_READING_FILE);
    GPU_ERROR(gpu_index_read(fp, index, GPU_SA));
  }

  fclose(fp);
  return (SUCCESS);
}

gpu_error_t gpu_io_save_index_GEM_FULL(const char* const fn, const gpu_index_buffer_t* const index, const gpu_index_coding_t activeModules)
{
  FILE *fp = NULL;
  size_t result;
  gpu_module_t storedModules = index->activeModules & activeModules;
  off64_t fileOffsetFMIndex = 0, fileOffsetSAIndex = 0, fileOffsetRef = 0;

  if((index->activeModules & activeModules) == 0)
    return(E_MODULE_NOT_FOUND);

  fp = fopen(fn, "wb");
  if (fp == NULL) return (E_OPENING_FILE);

  result = fwrite(&storedModules,     sizeof(gpu_module_t), 1, fp);
  if (result != 1) return (E_WRITING_FILE);
  result = fwrite(&fileOffsetFMIndex, sizeof(off64_t), 1, fp);
  if (result != 1) GPU_ERROR(E_WRITING_FILE);
  result = fwrite(&fileOffsetSAIndex, sizeof(off64_t), 1, fp);
  if (result != 1) GPU_ERROR(E_WRITING_FILE);
  result = fwrite(&fileOffsetRef,     sizeof(off64_t), 1, fp);
  if (result != 1) GPU_ERROR(E_WRITING_FILE);

  if(storedModules & GPU_FMI){
    fileOffsetFMIndex = ftello64(fp);
    GPU_ERROR(gpu_index_write(fp, index, GPU_FMI));
  }

  if(storedModules & GPU_SA){
    fileOffsetSAIndex = ftello64(fp);
    GPU_ERROR(gpu_index_write(fp, index, GPU_SA));
  }

  rewind(fp);
  result = fwrite(&storedModules,     sizeof(gpu_module_t), 1, fp);
  if (result != 1) GPU_ERROR(E_WRITING_FILE);
  result = fwrite(&fileOffsetFMIndex, sizeof(off64_t), 1, fp);
  if (result != 1) GPU_ERROR(E_WRITING_FILE);
  result = fwrite(&fileOffsetSAIndex, sizeof(off64_t), 1, fp);
  if (result != 1) GPU_ERROR(E_WRITING_FILE);
  result = fwrite(&fileOffsetRef,     sizeof(off64_t), 1, fp);
  if (result != 1) GPU_ERROR(E_WRITING_FILE);

  fclose(fp);
  return (SUCCESS);
}

gpu_error_t gpu_io_load_reference_specs_GEM_FULL(const char* const fn, gpu_reference_buffer_t* const reference)
{
  FILE *fp = NULL;
  size_t result;
  off64_t fileOffsetFMIndex = 0, fileOffsetSAIndex = 0, fileOffsetRef = 0;

  fp = fopen(fn, "rb");
  if (fp == NULL) return (E_OPENING_FILE);

  result = fread(&reference->activeModules, sizeof(gpu_module_t), 1, fp);
  if (result != 1) return (E_READING_FILE);

  if((reference->activeModules & GPU_REFERENCE) == 0)
    return (E_MODULE_NOT_FOUND);

  result = fread(&fileOffsetFMIndex, sizeof(off64_t), 1, fp);
  if (result != 1) return (E_READING_FILE);
  result = fread(&fileOffsetSAIndex, sizeof(off64_t), 1, fp);
  if (result != 1) GPU_ERROR(E_READING_FILE);
  result = fread(&fileOffsetRef, sizeof(off64_t), 1, fp);
  if (result != 1) return (E_READING_FILE);
  result = fseeko64(fp, fileOffsetRef, SEEK_SET);
  if (result != 0) return (E_READING_FILE);

  GPU_ERROR(gpu_reference_read_specs(fp, reference, GPU_REFERENCE));

  fclose(fp);
  return (SUCCESS);
}

gpu_error_t gpu_io_load_reference_GEM_FULL(const char* const fn, gpu_reference_buffer_t* const reference)
{
  FILE *fp = NULL;
  size_t result;
  off64_t fileOffsetFMIndex = 0, fileOffsetSAIndex = 0, fileOffsetRef = 0;

  fp = fopen(fn, "rb");
  if (fp == NULL) return (E_OPENING_FILE);

  result = fread(&reference->activeModules, sizeof(gpu_module_t), 1, fp);
  if (result != 1) return (E_READING_FILE);

  if((reference->activeModules & GPU_REFERENCE) == 0)
    return (E_MODULE_NOT_FOUND);

  result = fread(&fileOffsetFMIndex, sizeof(off64_t), 1, fp);
  if (result != 1) return (E_READING_FILE);
  result = fread(&fileOffsetSAIndex, sizeof(off64_t), 1, fp);
  if (result != 1) GPU_ERROR(E_READING_FILE);
  result = fread(&fileOffsetRef, sizeof(off64_t), 1, fp);
  if (result != 1) return (E_READING_FILE);
  result = fseeko64(fp, fileOffsetRef, SEEK_SET);
  if (result != 0) return (E_READING_FILE);

  GPU_ERROR(gpu_reference_read(fp, reference, GPU_REFERENCE));

  fclose(fp);
  return (SUCCESS);
}

gpu_error_t gpu_io_save_reference_GEM_FULL(const char* const fn, const gpu_reference_buffer_t* const reference)
{
  FILE *fp = NULL;
  size_t result;
  gpu_module_t activeModules = reference->activeModules;
  off64_t fileOffsetFMIndex = 0, fileOffsetSAIndex = 0, fileOffsetRef = 0;

  fp = fopen(fn, "wb");
  if (fp == NULL) return (E_OPENING_FILE);

  result = fwrite(&activeModules, sizeof(gpu_module_t), 1, fp);
  if (result != 1) return (E_READING_FILE);
  result = fwrite(&fileOffsetFMIndex, sizeof(off64_t), 1, fp);
  if (result != 1) return (E_READING_FILE);
  result = fwrite(&fileOffsetSAIndex, sizeof(off64_t), 1, fp);
  if (result != 1) return (E_READING_FILE);
  result = fwrite(&fileOffsetRef,   sizeof(off64_t), 1, fp);
  if (result != 1) return (E_READING_FILE);

  fileOffsetRef = ftello64(fp);
  GPU_ERROR(gpu_reference_write(fp, reference, GPU_REFERENCE));

  rewind(fp);
  result = fwrite(&activeModules,   sizeof(gpu_module_t), 1, fp);
  if (result != 1) GPU_ERROR(E_WRITING_FILE);
  result = fwrite(&fileOffsetFMIndex, sizeof(off64_t), 1, fp);
  if (result != 1) GPU_ERROR(E_WRITING_FILE);
  result = fwrite(&fileOffsetSAIndex, sizeof(off64_t), 1, fp);
  if (result != 1) GPU_ERROR(E_WRITING_FILE);
  result = fwrite(&fileOffsetRef,   sizeof(off64_t), 1, fp);
  if (result != 1) GPU_ERROR(E_WRITING_FILE);

  fclose(fp);
  return (SUCCESS);
}

void gpu_io_save_indexed_structures_GEM_(const char* const fileName, const gpu_gem_fmi_dto_t* const gemFMindex,
                                         const gpu_gem_ref_dto_t* const gemRef, const gpu_gem_sa_dto_t* const gemSAindex,
                                         const gpu_module_t activeModules)
{
  off64_t fileOffsetFMIndex = 0, fileOffsetSAIndex = 0, fileOffsetRef = 0;
  FILE *fp = NULL;
  size_t result;

  /* Objects to initialize */
  gpu_reference_buffer_t ref;
  gpu_index_buffer_t     index;

  // Initialize the reference structure
  GPU_ERROR(gpu_reference_init_dto(&ref));
  ref.activeModules = activeModules;
  ref.size          = gemRef->ref_length;
  ref.numEntries    = GPU_DIV_CEIL(ref.size, GPU_REFERENCE_CHARS_PER_ENTRY) + GPU_REFERENCE_END_PADDING;

  // Initialize the index (SA & FMI) structure
  GPU_ERROR(gpu_index_init_dto(&index, activeModules));
  index.activeModules  = activeModules;
  index.fmi.bwtSize    = gemFMindex->bwt_length;
  index.fmi.numEntries = GPU_DIV_CEIL(index.fmi.bwtSize, GPU_FMI_ENTRY_SIZE) + 1;
  index.sa.sampligRate = gemSAindex->sa_sampling;
  index.sa.numEntries  = GPU_DIV_CEIL(gemSAindex->sa_length, gemSAindex->sa_sampling);

  fp = fopen(fileName, "wb");
  if (fp == NULL) GPU_ERROR(E_OPENING_FILE);
  result = fwrite(&activeModules,   sizeof(gpu_module_t), 1, fp);
  if (result != 1) GPU_ERROR(E_WRITING_FILE);
  result = fwrite(&fileOffsetFMIndex, sizeof(off64_t), 1, fp);
  if (result != 1) GPU_ERROR(E_WRITING_FILE);
  result = fwrite(&fileOffsetSAIndex, sizeof(off64_t), 1, fp);
  if (result != 1) GPU_ERROR(E_WRITING_FILE);
  result = fwrite(&fileOffsetRef,   sizeof(off64_t), 1, fp);
  if (result != 1) GPU_ERROR(E_WRITING_FILE);

  fileOffsetFMIndex = ftello64(fp);
  fileOffsetSAIndex = fileOffsetFMIndex;
  fileOffsetRef     = fileOffsetSAIndex;

  if(index.activeModules & GPU_INDEX){
    if(index.activeModules & GPU_FMI){
      GPU_ERROR(gpu_index_allocate(&index, GPU_FMI));
      GPU_ERROR(gpu_index_transform(&index, (gpu_index_dto_t*)gemFMindex, gemFMindex->index_coding, GPU_FMI));
      GPU_ERROR(gpu_index_write(fp, &index, GPU_FMI));
      //GPU_ERROR(gpu_save_index_PROFILE("internalIndexGEM", fmi)); //DEBUG: backup the index
      fileOffsetSAIndex = ftello64(fp);
      fileOffsetRef     = fileOffsetSAIndex;
    }
    if(index.activeModules & GPU_SA){
      GPU_ERROR(gpu_index_allocate(&index, GPU_SA));
      GPU_ERROR(gpu_index_transform(&index, (gpu_index_dto_t*)gemSAindex, gemSAindex->index_coding, GPU_SA));
      GPU_ERROR(gpu_index_write(fp, &index, GPU_SA));
      fileOffsetRef = ftello64(fp);
    }
  }

  if(ref.activeModules & GPU_REFERENCE){
    GPU_ERROR(gpu_reference_allocate(&ref, GPU_REFERENCE));
    GPU_ERROR(gpu_reference_transform(&ref, (char*)gemRef, gemRef->ref_coding, GPU_REFERENCE));
    GPU_ERROR(gpu_reference_write(fp, &ref, GPU_REFERENCE));
  }

  rewind(fp);
  result = fwrite(&activeModules,    sizeof(gpu_module_t), 1, fp);
  if (result != 1) GPU_ERROR(E_WRITING_FILE);
  result = fwrite(&fileOffsetFMIndex, sizeof(off64_t), 1, fp);
  if (result != 1) GPU_ERROR(E_WRITING_FILE);
  result = fwrite(&fileOffsetSAIndex, sizeof(off64_t), 1, fp);
  if (result != 1) GPU_ERROR(E_WRITING_FILE);
  result = fwrite(&fileOffsetRef,     sizeof(off64_t), 1, fp);
  if (result != 1) GPU_ERROR(E_WRITING_FILE);

  fclose(fp);
}

#endif /* GPU_IO_C_ */

