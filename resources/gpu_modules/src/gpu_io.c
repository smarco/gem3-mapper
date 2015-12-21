#include "../include/gpu_io.h"


/************************************************************
Primitives for input/output
************************************************************/


gpu_error_t gpu_load_BWT_MFASTA(const char* const fn, gpu_index_buffer_t* const fmi, char **h_BWT)
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

  fmi->bwtSize    = position;
  fmi->numEntries = GPU_DIV_CEIL(fmi->bwtSize, GPU_FMI_ENTRY_SIZE) + 1;
  (* h_BWT)     = h_ascii_BWT;

  fclose(fp);
  return (SUCCESS);
}

gpu_error_t gpu_save_index_PROFILE(const char* const fn, const gpu_index_buffer_t* const index)
{
  const uint32_t sizeFileName = 512;
  char fileName[sizeFileName];
  FILE *fp = NULL;

  sprintf(fileName, "%s.%lu.%u.fmi", fn, index->bwtSize, GPU_FMI_ENTRY_SIZE);
  fp = fopen(fileName, "wb");
  if (fp == NULL) return (E_WRITING_FILE);

  GPU_ERROR(gpu_write_index(fp, index));

  fclose(fp);
  return (SUCCESS);
}

gpu_error_t gpu_load_index_PROFILE(const char* const fn, gpu_index_buffer_t* const index)
{
  FILE *fp = NULL;

  fp = fopen(fn, "rb");
  if (fp == NULL) return (E_OPENING_FILE);

  GPU_ERROR(gpu_read_index(fp, index));

  fclose(fp);
  return (SUCCESS);
}

gpu_error_t gpu_load_reference_MFASTA(const char *fn, gpu_reference_buffer_t* const reference)
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
  GPU_ERROR(gpu_transform_reference(tmp_reference, reference, GPU_REF_ASCII));

  fclose(fp);
  free(tmp_reference);
  return (SUCCESS);
}

gpu_error_t gpu_load_reference_PROFILE(const char* const fn, gpu_reference_buffer_t* const reference)
{
  FILE *fp = NULL;
  size_t result;

  fp = fopen(fn, "rb");
  if (fp == NULL) return (E_OPENING_FILE);

  GPU_ERROR(gpu_read_reference(fp, reference));

  fclose(fp);
  return (SUCCESS);
}

gpu_error_t gpu_save_reference_PROFILE(const char* const fn, const gpu_reference_buffer_t* const reference)
{
  FILE *fp = NULL;
  size_t result;

  fp = fopen(fn, "wb");
  if (fp == NULL) return (E_OPENING_FILE);

  GPU_ERROR(gpu_write_reference(fp, reference));

  fclose(fp);
  return (SUCCESS);
}

gpu_error_t gpu_load_index_GEM_FULL(const char *fn, gpu_index_buffer_t* const index)
{
  FILE *fp = NULL;
  size_t result;
  off64_t fileOffsetIndex = 0, fileOffsetRef = 0;

  fp = fopen(fn, "rb");
  if (fp == NULL) return (E_OPENING_FILE);

  result = fread(&index->activeModules, sizeof(gpu_module_t), 1, fp);
  if (result != 1) return (E_READING_FILE);
  if((index->activeModules & GPU_INDEX) == 0) return(E_MODULE_NOT_FOUND);

  result = fread(&fileOffsetIndex, sizeof(off64_t), 1, fp);
  if (result != 1) return (E_READING_FILE);
  result = fread(&fileOffsetRef, sizeof(off64_t), 1, fp);
  if (result != 1) return (E_READING_FILE);
  result = fseeko64(fp, fileOffsetIndex, SEEK_SET);
  if (result != 0) return (E_READING_FILE);

  GPU_ERROR(gpu_read_index(fp, index));

  fclose(fp);
  return (SUCCESS);
}

gpu_error_t gpu_load_reference_GEM_FULL(const char* const fn, gpu_reference_buffer_t* const reference)
{
  FILE *fp = NULL;
  size_t result;
  off64_t fileOffsetIndex = 0, fileOffsetRef = 0;

  fp = fopen(fn, "rb");
  if (fp == NULL) return (E_OPENING_FILE);

  result = fread(&reference->activeModules, sizeof(gpu_module_t), 1, fp);
  if (result != 1) return (E_READING_FILE);
  if((reference->activeModules & GPU_REFERENCE) == 0) return (E_MODULE_NOT_FOUND);

  result = fread(&fileOffsetIndex, sizeof(off64_t), 1, fp);
  if (result != 1) return (E_READING_FILE);
  result = fread(&fileOffsetRef, sizeof(off64_t), 1, fp);
  if (result != 1) return (E_READING_FILE);
  result = fseeko64(fp, fileOffsetRef, SEEK_SET);
  if (result != 0) return (E_READING_FILE);

  GPU_ERROR(gpu_read_reference(fp, reference));

  fclose(fp);
  return (SUCCESS);
}

gpu_error_t gpu_save_reference_GEM_FULL(const char* const fn, const gpu_reference_buffer_t* const reference)
{
  FILE *fp = NULL;
  size_t result;
  gpu_module_t activeModules = reference->activeModules & GPU_REFERENCE;
  off64_t fileOffsetIndex = 0, fileOffsetRef = 0;

  fp = fopen(fn, "wb");
  if (fp == NULL) return (E_OPENING_FILE);

  result = fwrite(&activeModules, sizeof(gpu_module_t), 1, fp);
  if (result != 1) return (E_READING_FILE);
  result = fwrite(&fileOffsetIndex, sizeof(off64_t), 1, fp);
  if (result != 1) return (E_READING_FILE);
  result = fwrite(&fileOffsetRef,   sizeof(off64_t), 1, fp);
  if (result != 1) return (E_READING_FILE);

  fileOffsetRef = ftello64(fp);
  GPU_ERROR(gpu_write_reference(fp, reference));

  rewind(fp);
  result = fwrite(&activeModules,   sizeof(gpu_module_t), 1, fp);
  if (result != 1) GPU_ERROR(E_WRITING_FILE);
  result = fwrite(&fileOffsetIndex, sizeof(off64_t), 1, fp);
  if (result != 1) GPU_ERROR(E_WRITING_FILE);
  result = fwrite(&fileOffsetRef,   sizeof(off64_t), 1, fp);
  if (result != 1) GPU_ERROR(E_WRITING_FILE);

  fclose(fp);
  return (SUCCESS);
}

gpu_error_t gpu_save_index_GEM_FULL(const char* const fn, const gpu_index_buffer_t* const index)
{
  FILE *fp = NULL;
  size_t result;
  gpu_module_t activeModules = index->activeModules & GPU_INDEX;
  off64_t fileOffsetIndex = 0, fileOffsetRef = 0;

  fp = fopen(fn, "wb");
  if (fp == NULL) return (E_OPENING_FILE);

  result = fwrite(&activeModules, sizeof(gpu_module_t), 1, fp);
  if (result != 1) return (E_WRITING_FILE);
  result = fwrite(&fileOffsetIndex, sizeof(off64_t), 1, fp);
  if (result != 1) GPU_ERROR(E_WRITING_FILE);
  result = fwrite(&fileOffsetRef,   sizeof(off64_t), 1, fp);
  if (result != 1) GPU_ERROR(E_WRITING_FILE);

  fileOffsetIndex = ftello64(fp);
  GPU_ERROR(gpu_write_index(fp, index));

  rewind(fp);
  result = fwrite(&activeModules,   sizeof(gpu_module_t), 1, fp);
  if (result != 1) GPU_ERROR(E_WRITING_FILE);
  result = fwrite(&fileOffsetIndex, sizeof(off64_t), 1, fp);
  if (result != 1) GPU_ERROR(E_WRITING_FILE);
  result = fwrite(&fileOffsetRef,   sizeof(off64_t), 1, fp);
  if (result != 1) GPU_ERROR(E_WRITING_FILE);

  fclose(fp);
  return (SUCCESS);
}

void gpu_save_indexed_structures_GEM_(const char* const fileName, const gpu_gem_fmi_dto_t* const gemIndex,
                                      const gpu_gem_ref_dto_t* const gemRef, const gpu_module_t activeModules)
{
  gpu_reference_buffer_t ref;
  gpu_index_buffer_t     index;
  off64_t fileOffsetIndex = 0, fileOffsetRef = 0;
  FILE *fp = NULL;
  size_t result;

  //Initialize the reference structure
  GPU_ERROR(gpu_init_reference_dto(&ref));
  ref.activeModules = activeModules & GPU_REFERENCE;

  //Initialize the index structure
  GPU_ERROR(gpu_init_index_dto(&index));
  index.activeModules = activeModules & GPU_INDEX;

  fp = fopen(fileName, "wb");
  if (fp == NULL) GPU_ERROR(E_OPENING_FILE);
  result = fwrite(&activeModules,   sizeof(gpu_module_t), 1, fp);
  if (result != 1) GPU_ERROR(E_WRITING_FILE);
  result = fwrite(&fileOffsetIndex, sizeof(off64_t), 1, fp);
  if (result != 1) GPU_ERROR(E_WRITING_FILE);
  result = fwrite(&fileOffsetRef,   sizeof(off64_t), 1, fp);
  if (result != 1) GPU_ERROR(E_WRITING_FILE);

  fileOffsetIndex = ftello64(fp);
  fileOffsetRef   = fileOffsetIndex;

  if(index.activeModules & GPU_INDEX){
    GPU_ERROR(gpu_transform_index((char*)gemIndex, &index, gemIndex->index_coding));
    GPU_ERROR(gpu_write_index(fp, &index));
    //GPU_ERROR(gpu_save_index_PROFILE("internalIndexGEM", fmi)); //DEBUG: backup the index
    fileOffsetRef = ftello64(fp);
  }

  if(ref.activeModules & GPU_REFERENCE){
    GPU_ERROR(gpu_transform_reference((char*)gemRef, &ref, gemRef->ref_coding));
    GPU_ERROR(gpu_write_reference(fp, &ref));
  }

  rewind(fp);
  result = fwrite(&activeModules,   sizeof(gpu_module_t), 1, fp);
  if (result != 1) GPU_ERROR(E_WRITING_FILE);
  result = fwrite(&fileOffsetIndex, sizeof(off64_t), 1, fp);
  if (result != 1) GPU_ERROR(E_WRITING_FILE);
  result = fwrite(&fileOffsetRef,   sizeof(off64_t), 1, fp);
  if (result != 1) GPU_ERROR(E_WRITING_FILE);

  fclose(fp);
}
