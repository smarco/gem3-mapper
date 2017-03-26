/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_FMI_INDEX_C_
#define GPU_FMI_INDEX_C_

#include "../include/gpu_fmi_index.h"


/************************************************************
Get information functions
************************************************************/

gpu_error_t gpu_fmi_index_get_size(const gpu_fmi_buffer_t* const fmi, size_t* const bytesPerFMI)
{
  (* bytesPerFMI) = fmi->numEntries * sizeof(gpu_fmi_entry_t);
  return (SUCCESS);
}


/************************************************************
 GLOBAL METHODS: INPUT / OUPUT Functions
************************************************************/

gpu_error_t gpu_fmi_index_read_specs(int fp, gpu_fmi_buffer_t* const fmi)
{
  size_t result, bytesRequest;

  bytesRequest = sizeof(uint64_t);
  result = read(fp, (void *)&fmi->numEntries, bytesRequest);
  if (result != bytesRequest) return (E_READING_FILE);

  bytesRequest = sizeof(uint64_t);
  result = read(fp, (void *)&fmi->bwtSize, bytesRequest);
  if (result != bytesRequest) return (E_READING_FILE);

  return (SUCCESS);
}

gpu_error_t gpu_fmi_index_read(int fp, gpu_fmi_buffer_t* const fmi)
{
  const size_t   bytesRequest = sizeof(gpu_fmi_entry_t) * fmi->numEntries;
  const uint64_t numRequests  = GPU_DIV_CEIL(bytesRequest, GPU_FILE_SIZE_BLOCK);
  size_t result, numBytesRequested = 0;
  uint64_t idRequest;

  for(idRequest = 0; idRequest < numRequests; ++idRequest){
    const size_t requestSize = GPU_MIN(GPU_FILE_SIZE_BLOCK, bytesRequest - numBytesRequested);
    result = read(fp, (void* )fmi->h_fmi + numBytesRequested, requestSize);
    if (result != requestSize) return (E_READING_FILE);
    numBytesRequested += requestSize;
  }

  return (SUCCESS);
}

gpu_error_t gpu_fmi_index_write_specs(int fp, const gpu_fmi_buffer_t* const fmi)
{
  size_t result, bytesRequest;

  bytesRequest = sizeof(uint64_t);
  result = write(fp, (void *)&fmi->numEntries, bytesRequest);
  if (result != bytesRequest) return (E_WRITING_FILE);

  bytesRequest = sizeof(uint64_t);
  result = write(fp, (void *)&fmi->bwtSize, bytesRequest);
  if (result != bytesRequest) return (E_WRITING_FILE);

  return (SUCCESS);
}

gpu_error_t gpu_fmi_index_write(int fp, const gpu_fmi_buffer_t* const fmi)
{
  const size_t   bytesRequest = sizeof(gpu_fmi_entry_t) * fmi->numEntries;
  const uint64_t numRequests  = GPU_DIV_CEIL(bytesRequest, GPU_FILE_SIZE_BLOCK);
  size_t result, numBytesRequested = 0;
  uint64_t idRequest;

  for(idRequest = 0; idRequest < numRequests; ++idRequest){
    const size_t requestSize = GPU_MIN(GPU_FILE_SIZE_BLOCK, bytesRequest - numBytesRequested);
    result = write(fp, (void* )fmi->h_fmi + numBytesRequested, requestSize);
    if (result != requestSize) return (E_WRITING_FILE);
    numBytesRequested += requestSize;
  }

  return (SUCCESS);
}

gpu_error_t gpu_fmi_index_load_specs_MFASTA_FULL(const char* const fn, gpu_fmi_buffer_t* const fmi)
{
  FILE *fp = NULL;
  uint64_t sizeFile = 0;

  fp = fopen(fn, "rb");
  if (fp == NULL) return (E_OPENING_FILE);

  fseek(fp, 0L, SEEK_END);
  sizeFile = ftell(fp);

  fmi->bwtSize    = sizeFile;
  fmi->numEntries = GPU_DIV_CEIL(fmi->bwtSize, GPU_FMI_ENTRY_SIZE) + 1;

  fclose(fp);
  return (SUCCESS);
}

gpu_error_t gpu_fmi_index_load_MFASTA_FULL(const char* const fn, gpu_fmi_buffer_t* const fmi, char **h_BWT)
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
  (* h_BWT)       = h_ascii_BWT;

  fclose(fp);
  return (SUCCESS);
}


/************************************************************
 GLOBAL METHODS: Transfer the index (HOST <-> DEVICES)
************************************************************/

gpu_error_t gpu_fmi_index_transfer_CPU_to_GPUs(gpu_fmi_buffer_t* const fmi, gpu_device_info_t** const devices)
{
  uint32_t deviceFreeMemory, idSupportedDevice;
  uint32_t numSupportedDevices = devices[0]->numSupportedDevices;

  for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
    if(fmi->memorySpace[idSupportedDevice] == GPU_DEVICE_MAPPED){
      const size_t cpySize = fmi->numEntries * sizeof(gpu_fmi_entry_t);
      deviceFreeMemory = gpu_device_get_free_memory(devices[idSupportedDevice]->idDevice);
      if ((GPU_CONVERT__B_TO_MB(cpySize)) > deviceFreeMemory) return(E_INSUFFICIENT_MEM_GPU);
        CUDA_ERROR(cudaSetDevice(devices[idSupportedDevice]->idDevice));
      //Synchronous allocate & transfer the FM-index to the GPU
      CUDA_ERROR(cudaMalloc((void**) &fmi->d_fmi[idSupportedDevice], cpySize));
      CUDA_ERROR(cudaMemcpy(fmi->d_fmi[idSupportedDevice], fmi->h_fmi, cpySize, cudaMemcpyHostToDevice));
    }else{
      fmi->d_fmi[idSupportedDevice] = fmi->h_fmi;
    }
  }

  return (SUCCESS);
}


/************************************************************
 GLOBAL METHODS: Functions to transform the index
************************************************************/

gpu_error_t gpu_fmi_index_transform_ASCII(const char* const textBWT, gpu_fmi_buffer_t* const fmi)
{
    gpu_index_bitmap_entry_t  *h_bitmaps_BWT  = NULL;
    gpu_index_counter_entry_t *h_counters_FMI = NULL;

    const uint32_t bwtNumEntries      = fmi->numEntries * (GPU_FMI_ENTRY_SIZE / GPU_UINT32_LENGTH);
    const uint32_t countersNumEntries = fmi->numEntries * (GPU_FMI_ENTRY_SIZE / GPU_UINT32_LENGTH);

    h_bitmaps_BWT = (gpu_index_bitmap_entry_t *) malloc(bwtNumEntries * sizeof(gpu_index_bitmap_entry_t));
    if (h_bitmaps_BWT == NULL) return (E_ALLOCATE_MEM);
    h_counters_FMI = (gpu_index_counter_entry_t *) malloc(countersNumEntries * sizeof(gpu_index_counter_entry_t));
    if (h_counters_FMI == NULL) return (E_ALLOCATE_MEM);

    GPU_ERROR(gpu_fmi_index_build_PEQ(fmi, textBWT, h_bitmaps_BWT));
    GPU_ERROR(gpu_fmi_index_build_COUNTERS(fmi, h_counters_FMI, textBWT));
    GPU_ERROR(gpu_fmi_index_build_FMI(fmi, h_bitmaps_BWT, h_counters_FMI));

    free(h_bitmaps_BWT);
    free(h_counters_FMI);
    return(SUCCESS);
}

gpu_error_t gpu_fmi_index_transform_GEM_FULL(const gpu_gem_fmi_dto_t* const gpu_gem_fmi_dto, gpu_fmi_buffer_t* const fmi)
{
  // BWT Parameters
  const uint64_t BWT_MINOR_BLOCKS_PER_MAYOR_BLOCK = (1<<10); /* 1024 */
  const uint64_t BWT_MINOR_BLOCK_LENGTH = 64;
  const uint64_t* c = gpu_gem_fmi_dto->c;
  const uint64_t* C = gpu_gem_fmi_dto->C;
  const uint64_t* mayor_counters = gpu_gem_fmi_dto->mayor_counters;
  const uint32_t* bwt_mem = (uint32_t*) gpu_gem_fmi_dto->bwt_mem;
  const uint64_t bwt_length = gpu_gem_fmi_dto->bwt_length;
  // Initialize FMI structures
  fmi->bwtSize    = bwt_length;
  fmi->numEntries = GPU_DIV_CEIL(fmi->bwtSize, GPU_FMI_ENTRY_SIZE) + 1;
  gpu_fmi_entry_t* h_fmi = fmi->h_fmi; // Host FMI
  // Iterate thought the BWT memory layout
  uint16_t* minor_counters;
  uint32_t x0, x1, x2;
  uint32_t y0, y1, y2;
  uint32_t z0, z1, z2;
  uint32_t w0, w1, w2;
  uint64_t minor_block = 0, h_fmi_entry = 0;
  uint64_t bwt_pos;
  for (bwt_pos=0;bwt_pos<bwt_length;) {
    // Get next block (Block0)
    minor_counters = (uint16_t*)bwt_mem; bwt_mem += 4; // Minor Counters
    x0 = *(bwt_mem); ++bwt_mem;
    y0 = *(bwt_mem); ++bwt_mem;
    x1 = *(bwt_mem); ++bwt_mem;
    y1 = *(bwt_mem); ++bwt_mem;
    x2 = *(bwt_mem); ++bwt_mem;
    y2 = *(bwt_mem); ++bwt_mem;
    bwt_mem += 2; // Sampling bitmap
    bwt_pos+=BWT_MINOR_BLOCK_LENGTH; // Next Block
    // Check end-of-bwt
    if (bwt_pos<bwt_length) {
      // Get next block (Block1)
      bwt_mem += 4; // Skip Minor Counters
      z0 = *(bwt_mem); ++bwt_mem;
      w0 = *(bwt_mem); ++bwt_mem;
      z1 = *(bwt_mem); ++bwt_mem;
      w1 = *(bwt_mem); ++bwt_mem;
      z2 = *(bwt_mem); ++bwt_mem;
      w2 = *(bwt_mem); ++bwt_mem;
      bwt_mem += 2; // Sampling bitmap
      bwt_pos+=BWT_MINOR_BLOCK_LENGTH; // Next Block
    } else {
      // Padding zeros
      z0 = 0; z1 = 0; z2 = ~0;
      w0 = 0; w1 = 0; w2 = ~0;
    }
    // Write Counters
    if (h_fmi_entry % 2 == 0) {
      h_fmi->counters[0] = minor_counters[0] + mayor_counters[0]; // 'A'
      h_fmi->counters[1] = minor_counters[1] + mayor_counters[1]; // 'C'
    } else {
      h_fmi->counters[0] = minor_counters[2] + mayor_counters[2]; // 'G'
      h_fmi->counters[1] = minor_counters[3] + mayor_counters[3]; // 'T'
    }
    // Write Bitmap
    h_fmi->bitmaps[0]  = gpu_bit_reverse(y0) & gpu_bit_reverse(~y2);
    h_fmi->bitmaps[1]  = gpu_bit_reverse(y1) & gpu_bit_reverse(~y2);
    h_fmi->bitmaps[2]  = gpu_bit_reverse(~y2);
    h_fmi->bitmaps[3]  = gpu_bit_reverse(x0) & gpu_bit_reverse(~x2);
    h_fmi->bitmaps[4]  = gpu_bit_reverse(z0) & gpu_bit_reverse(~z2);
    h_fmi->bitmaps[5]  = gpu_bit_reverse(z1) & gpu_bit_reverse(~z2);
    h_fmi->bitmaps[6]  = gpu_bit_reverse(~z2);
    h_fmi->bitmaps[7]  = gpu_bit_reverse(x1) & gpu_bit_reverse(~x2);
    h_fmi->bitmaps[8]  = gpu_bit_reverse(w0) & gpu_bit_reverse(~w2);
    h_fmi->bitmaps[9]  = gpu_bit_reverse(w1) & gpu_bit_reverse(~w2);
    h_fmi->bitmaps[10] = gpu_bit_reverse(~w2);
    h_fmi->bitmaps[11] = gpu_bit_reverse(~x2);
    // Next Entry
    ++h_fmi; ++h_fmi_entry;
    minor_block += 2;
    if (minor_block >= BWT_MINOR_BLOCKS_PER_MAYOR_BLOCK) {
      minor_block = 0;
      mayor_counters += 8;
    }
  }

  --h_fmi;
  int32_t padding_module = bwt_length % GPU_FMI_ENTRY_SIZE;
  uint32_t mask;
  mask = gpu_gen_mask(padding_module);
  padding_module -= GPU_UINT32_LENGTH;
  h_fmi->bitmaps[3]   &= mask;
  h_fmi->bitmaps[7]   &= mask;
  h_fmi->bitmaps[11]  &= mask;
  mask = gpu_gen_mask(padding_module);
  padding_module -= GPU_UINT32_LENGTH;
  h_fmi->bitmaps[0]  &= mask;
  h_fmi->bitmaps[1]  &= mask;
  h_fmi->bitmaps[2]  &= mask;
  mask = gpu_gen_mask(padding_module);
  padding_module -= GPU_UINT32_LENGTH;
  h_fmi->bitmaps[4]  &= mask;
  h_fmi->bitmaps[5]  &= mask;
  h_fmi->bitmaps[6]  &= mask;
  mask = gpu_gen_mask(padding_module);
  padding_module -= GPU_UINT32_LENGTH;
  h_fmi->bitmaps[8]  &= mask;
  h_fmi->bitmaps[9]  &= mask;
  h_fmi->bitmaps[10] &= mask;

  ++h_fmi;
  // Ghost entry for alternate counters
  if (h_fmi_entry % 2 == 0) {   // Write Counters
    h_fmi->counters[0] = c[0] + C[0]; // 'A'
    h_fmi->counters[1] = c[1] + C[1]; // 'C'
  } else {
    h_fmi->counters[0] = c[2] + C[2]; // 'G'
    h_fmi->counters[1] = c[3] + C[3]; // 'T'
  }
  // Write Bitmap
  h_fmi->bitmaps[0]  =  0; h_fmi->bitmaps[1]  =  0;
  h_fmi->bitmaps[2]  =  0; h_fmi->bitmaps[3]  =  0;
  h_fmi->bitmaps[4]  =  0; h_fmi->bitmaps[5]  =  0;
  h_fmi->bitmaps[6]  =  0; h_fmi->bitmaps[7]  =  0;
  h_fmi->bitmaps[8]  =  0; h_fmi->bitmaps[9]  =  0;
  h_fmi->bitmaps[10] =  0; h_fmi->bitmaps[11] =  0;
  // Return SUCCESS
  return (SUCCESS);
}


gpu_error_t gpu_fmi_index_transform_MFASTA_FULL(const char* const indexRaw, gpu_fmi_buffer_t* const fmi)
{
  char* h_BWT = NULL;
  GPU_ERROR(gpu_fmi_index_load_MFASTA_FULL(indexRaw, fmi, &h_BWT));
  GPU_ERROR(gpu_fmi_index_transform_ASCII(h_BWT, fmi));
  free(h_BWT);
  return (SUCCESS);
}


/************************************************************
 GLOBAL METHODS: Index initialization functions
************************************************************/

gpu_error_t gpu_fmi_index_init_dto(gpu_fmi_buffer_t* const fmi)
{
  //Initialize the FMI index structure
  fmi->d_fmi          = NULL;
  fmi->h_fmi          = NULL;
  fmi->hostAllocStats = GPU_PAGE_UNLOCKED;
  fmi->memorySpace    = NULL;
  fmi->bwtSize        = 0;
  fmi->numEntries     = 0;

  return (SUCCESS);
}

gpu_error_t gpu_fmi_index_init(gpu_fmi_buffer_t* const fmi, const uint64_t bwtSize, const uint32_t numSupportedDevices)
{
  uint32_t idSupDevice;
  GPU_ERROR(gpu_fmi_index_init_dto(fmi));

  fmi->bwtSize        = bwtSize;
  fmi->numEntries     = GPU_DIV_CEIL(fmi->bwtSize, GPU_FMI_ENTRY_SIZE) + 1;

  fmi->d_fmi = (gpu_fmi_entry_t **) malloc(numSupportedDevices * sizeof(gpu_fmi_entry_t *));
  if (fmi->d_fmi == NULL) GPU_ERROR(E_ALLOCATE_MEM);
  fmi->memorySpace = (memory_alloc_t *) malloc(numSupportedDevices * sizeof(memory_alloc_t));
  if (fmi->memorySpace == NULL) GPU_ERROR(E_ALLOCATE_MEM);

  for(idSupDevice = 0; idSupDevice < numSupportedDevices; ++idSupDevice){
    fmi->d_fmi[idSupDevice]       = NULL;
    fmi->memorySpace[idSupDevice] = GPU_NONE_MAPPED;
  }

  return (SUCCESS);
}

gpu_error_t gpu_fmi_index_allocate(gpu_fmi_buffer_t* const fmi)
{
  fmi->numEntries = GPU_DIV_CEIL(fmi->bwtSize, GPU_FMI_ENTRY_SIZE) + 1;
  if(fmi->hostAllocStats & GPU_PAGE_LOCKED){
    CUDA_ERROR(cudaHostAlloc((void**) &fmi->h_fmi, fmi->numEntries * sizeof(gpu_fmi_entry_t), cudaHostAllocMapped));
  }else{
    fmi->h_fmi = malloc (fmi->numEntries * sizeof(gpu_fmi_entry_t));
    if (fmi->h_fmi == NULL) return (E_ALLOCATE_MEM);
  }
  return(SUCCESS);
}

/************************************************************
 GLOBAL METHODS: Functions to release DEVICE & HOST indexes
************************************************************/

gpu_error_t gpu_fmi_index_free_host(gpu_fmi_buffer_t* const fmi)
{
  if(fmi->h_fmi != NULL){
      if(fmi->hostAllocStats == GPU_PAGE_LOCKED) CUDA_ERROR(cudaFreeHost(fmi->h_fmi));
      else free(fmi->h_fmi);
      fmi->h_fmi = NULL;
    }

    return(SUCCESS);
}

gpu_error_t gpu_fmi_index_free_unused_host(gpu_fmi_buffer_t* const fmi, gpu_device_info_t** const devices)
{
  uint32_t idSupportedDevice, numSupportedDevices;
  bool indexInHostSideUsed = false;

  numSupportedDevices = devices[0]->numSupportedDevices;
  //Free all the unused references in the host side
  for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
    if(fmi->memorySpace[idSupportedDevice] == GPU_HOST_MAPPED) indexInHostSideUsed = true;
  }

  if(!indexInHostSideUsed){
    GPU_ERROR(gpu_fmi_index_free_host(fmi));
  }

  return(SUCCESS);
}

gpu_error_t gpu_fmi_index_free_device(gpu_fmi_buffer_t* const fmi, gpu_device_info_t** const devices)
{
  const uint32_t numSupportedDevices = devices[0]->numSupportedDevices;
  uint32_t idSupportedDevice;

  //Free all the references in the devices
  for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
    CUDA_ERROR(cudaSetDevice(devices[idSupportedDevice]->idDevice));
    if(fmi->d_fmi[idSupportedDevice] != NULL){
      if(fmi->memorySpace[idSupportedDevice] == GPU_DEVICE_MAPPED)
        CUDA_ERROR(cudaFree(fmi->d_fmi[idSupportedDevice]));
      fmi->d_fmi[idSupportedDevice] = NULL;
    }
  }

  //Free the index list
  if(fmi->d_fmi != NULL){
    free(fmi->d_fmi);
    fmi->d_fmi = NULL;
  }

  return(SUCCESS);
}

gpu_error_t gpu_fmi_index_free_metainfo(gpu_fmi_buffer_t* const fmi)
{
  //Free the index list
  if(fmi->memorySpace != NULL){
    free(fmi->memorySpace);
    fmi->memorySpace = NULL;
  }
  return(SUCCESS);
}


/************************************************************
 LOCAL METHODS: Basic conversion primitives
************************************************************/

uint32_t gpu_char_to_bin(const char base)
{
  uint32_t indexBase = GPU_ENC_DNA_CHAR_X;
  indexBase = ((base =='A') || (base =='a')) ? GPU_ENC_DNA_CHAR_A : indexBase;
  indexBase = ((base =='C') || (base =='c')) ? GPU_ENC_DNA_CHAR_C : indexBase;
  indexBase = ((base =='G') || (base =='g')) ? GPU_ENC_DNA_CHAR_G : indexBase;
  indexBase = ((base =='T') || (base =='t')) ? GPU_ENC_DNA_CHAR_T : indexBase;
  return(indexBase);
}

char gpu_bin_to_char(const uint32_t indexBase)
{
  const char LUT[5] = {'A', 'C', 'G', 'T', 'N'};
  const char base   = (indexBase < 5) ? LUT[indexBase] : 'N';
  return (base);
}

void gpu_encode_entry_BWT_to_PEQ(gpu_index_bitmap_entry_t* const bwt_entry,
                                 const char base, const uint32_t size)
{
  const uint32_t binBase = gpu_char_to_bin(base);
  uint32_t idBase;
  for(idBase = 0; idBase < GPU_FMI_BWT_CHAR_LENGTH; ++idBase){
    bwt_entry->bitmaps[idBase] |= ((binBase >> idBase) & 0x1) << (GPU_UINT32_LENGTH - size - 1);
  }
}

void gpu_bining_bases(gpu_index_counter_entry_t* const counterEntry, const char base)
{
  const uint32_t indexBase = gpu_char_to_bin(base);
  if(indexBase < 4) counterEntry->counters[indexBase]++;
}

void gpu_index_set_layout_counter(const gpu_index_counter_entry_t* const h_counters_FMI,
                                  gpu_fmi_entry_t* const h_fmi, const uint64_t idEntry)
{
  uint32_t idCounter;
  for(idCounter = 0; idCounter < GPU_FMI_COUNTERS_PER_ENTRY; ++idCounter){
    const uint32_t FMI_SET_ALT_COUNTER = idEntry % GPU_FMI_COUNTERS_PER_ENTRY;
    h_fmi->counters[idCounter] = h_counters_FMI->counters[(FMI_SET_ALT_COUNTER * GPU_FMI_COUNTERS_PER_ENTRY) + idCounter];
  }
}

void gpu_index_set_layout_bitmap(gpu_index_bitmap_entry_t* const h_bitmap_BWT, gpu_fmi_entry_t* const h_fmi)
{
  const uint32_t LUT[12] = {3,7,11,0,1,2,4,5,6,8,9,10};                                                     // 1 1 (1) 0 2 2 (2) 0 3 3 (3) (0)
  uint32_t idPacket, idFMIBucket;
  const uint32_t FMI_NUM_BITMAPS        = GPU_FMI_ENTRY_SIZE * GPU_FMI_BWT_CHAR_LENGTH / GPU_UINT32_LENGTH; // 12 words             (4 FMI entries x 3 bits)
  const uint32_t NUM_PACKET_BMP_ENTRIES = GPU_FMI_ENTRY_SIZE  / GPU_UINT32_LENGTH;                          // 4 BMP entries        (128 bases / 32 bits)
  const uint32_t NUM_PACKET_FMI_ENTRIES = FMI_NUM_BITMAPS / NUM_PACKET_BMP_ENTRIES;                         // 3 FMI bitmap entries (12 bitmaps / 4 packets)

  for(idFMIBucket = 0; idFMIBucket < NUM_PACKET_BMP_ENTRIES; ++idFMIBucket)                                 // Iterate over BMP entries (4)
    h_bitmap_BWT[idFMIBucket].bitmaps[NUM_PACKET_FMI_ENTRIES - 1] = ~ h_bitmap_BWT[idFMIBucket].bitmaps[NUM_PACKET_FMI_ENTRIES - 1];

  for(idPacket = 0; idPacket < NUM_PACKET_BMP_ENTRIES; ++idPacket){                                       // Iterate over BMP entries (4)
    for(idFMIBucket = 0; idFMIBucket < NUM_PACKET_FMI_ENTRIES; ++idFMIBucket){                            // Iterate over FMI bitmaps (3)
      h_fmi->bitmaps[LUT[idPacket * NUM_PACKET_FMI_ENTRIES + idFMIBucket]] =  h_bitmap_BWT[idPacket].bitmaps[idFMIBucket];
    }
  }
}


/************************************************************
 LOCAL METHODS: Main conversion primitives
************************************************************/

gpu_error_t gpu_fmi_index_build_PEQ(const gpu_fmi_buffer_t* const fmi, const char* const h_ascii_BWT,
                                    gpu_index_bitmap_entry_t* const h_bitmap_BWT)
{
  uint64_t idEntry, i, bwtPosition;
  unsigned char bwtChar;

  //Padded to FMI_ENTRY_SIZE (128 bases)
  const uint32_t bwtNumEntries = fmi->numEntries * (GPU_FMI_ENTRY_SIZE / GPU_UINT32_LENGTH);

  for(idEntry = 0; idEntry < bwtNumEntries; ++idEntry){
    for(i = 0; i < GPU_UINT32_LENGTH; ++i){
      bwtPosition = (idEntry * GPU_UINT32_LENGTH) + i;
      if (bwtPosition < fmi->numEntries) bwtChar = h_ascii_BWT[bwtPosition];
        else bwtChar = 'N'; //filling BWT padding
      gpu_encode_entry_BWT_to_PEQ(&h_bitmap_BWT[idEntry], bwtChar, i);
    }
  }
  return(SUCCESS);
}

gpu_error_t gpu_fmi_index_build_COUNTERS(const gpu_fmi_buffer_t* const fmi, gpu_index_counter_entry_t* const h_counters_FMI,
                                         const char* const h_ascii_BWT)
{
  uint64_t idEntry, i, bwtPosition;
  uint32_t idBase;
  char bwtChar;

  gpu_index_counter_entry_t localCounters;
  const uint32_t            countersNumEntries = fmi->numEntries * (GPU_FMI_ENTRY_SIZE / GPU_UINT32_LENGTH);

  // Initialize first local BWT entry
  for(idEntry = 0; idEntry < countersNumEntries; ++idEntry){
    for(idBase = 0; idBase < GPU_FMI_NUM_COUNTERS; ++idBase){
      h_counters_FMI[idEntry].counters[idBase] = 0;
    }
  }

  // Accumulate values locally for each BWT entry
  for(idEntry = 0; idEntry < countersNumEntries - 1; ++idEntry){
    for(i = 0; i < GPU_UINT32_LENGTH; ++i){
      bwtPosition = (idEntry * GPU_UINT32_LENGTH) + i;
      if (bwtPosition < fmi->bwtSize) bwtChar = h_ascii_BWT[bwtPosition];
        else bwtChar = 'N'; //filling BWT padding (will be N)
      gpu_bining_bases(&h_counters_FMI[idEntry + 1], bwtChar);
    }
  }

  // Accumulate values globally for all the BWT
  for(idEntry = 1; idEntry < countersNumEntries; ++idEntry){
    for(idBase = 0; idBase < GPU_FMI_NUM_COUNTERS; ++idBase){
      h_counters_FMI[idEntry].counters[idBase] += h_counters_FMI[idEntry - 1].counters[idBase];
    }
  }

  // Prepare the local counters with the previous accumulative letters
  localCounters.counters[0] = 0;
  for(idBase = 1; idBase < GPU_FMI_NUM_COUNTERS; ++idBase){
    localCounters.counters[idBase] = localCounters.counters[idBase - 1] + h_counters_FMI[countersNumEntries - 4].counters[idBase - 1];
  }

  // Accumulate the previous alphabet letters to the global counters
  for(idEntry = 0; idEntry < countersNumEntries; ++idEntry){
    for(idBase = 1; idBase < GPU_FMI_NUM_COUNTERS; ++idBase){
      h_counters_FMI[idEntry].counters[idBase] += localCounters.counters[idBase];
    }
  }

  return (SUCCESS);
}

gpu_error_t gpu_fmi_index_build_FMI(gpu_fmi_buffer_t* const fmi, gpu_index_bitmap_entry_t* const h_bitmap_BWT,
                                    const gpu_index_counter_entry_t* const h_counters_FMI)
{
  const uint32_t BITMAPS_PER_FMI = GPU_FMI_ENTRY_SIZE / GPU_UINT32_LENGTH;  // 4 FMI bitmap entries
  uint64_t idEntry;

  for(idEntry = 0; idEntry < fmi->numEntries; ++idEntry){
    gpu_index_set_layout_counter(&h_counters_FMI[idEntry * BITMAPS_PER_FMI], &fmi->h_fmi[idEntry], idEntry);
    gpu_index_set_layout_bitmap(&h_bitmap_BWT[idEntry * BITMAPS_PER_FMI], &fmi->h_fmi[idEntry]);
  }

  return(SUCCESS);
}

#endif /* GPU_FMI_INDEX_C_ */

