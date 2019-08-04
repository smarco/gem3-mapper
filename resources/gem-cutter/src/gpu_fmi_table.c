/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2011-2018 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_FMI_TABLE_C_
#define GPU_FMI_TABLE_C_

#include "../include/gpu_fmi_table.h"

/* Local function for fmi-table rank queries */
uint32_t countBitmapCPU(const uint32_t bitmap, const int32_t shift, const uint32_t idxCounterGroup)
{
  // Definning the bitmap masking for the corresponding interval
  uint32_t mask = GPU_UINT32_ONES << (GPU_UINT32_LENGTH - shift);
  // Avoiding tha out-of-range interval bases
  mask = (shift > GPU_UINT32_LENGTH) ? GPU_UINT32_ONES : mask;
  mask = (shift > 0) ? mask : GPU_UINT32_ZEROS;
  // Counting the FM-index entry internal bases
  mask = (idxCounterGroup) ? ~mask : mask;
  return (__builtin_popcount(bitmap & mask));
}

/* Local function for fmi-table rank queries */
uint32_t computeBitmapsCPU(uint32_t* const vbitmap, const uint32_t bitmapPosition,
                           const uint32_t bit0, const uint32_t bit1,
                           const uint32_t missedEntry, const uint32_t idBitmap)
{
  // Compute the internal interval
  const int32_t  relativePosition = bitmapPosition - (idBitmap * GPU_UINT32_LENGTH);
  uint32_t bmpCollapsed, numCaracters;
  // Computing the wavelet bitmap
  vbitmap[0] = bit0 ? vbitmap[0] : ~vbitmap[0];
  vbitmap[1] = bit1 ? vbitmap[1] : ~vbitmap[1];
  // Counting the FM-index entry internal bases
  bmpCollapsed  = vbitmap[0] & vbitmap[1] & vbitmap[2];
  numCaracters  = countBitmapCPU(bmpCollapsed, relativePosition, missedEntry);
  // Returning the occurence caracters
  return (numCaracters);
}

/* Local function for fmi-table rank queries */
void LF_mapping_advance_step(const gpu_fmi_entry_t* const fmi, const uint64_t interval, uint64_t* const new_interval, const uint32_t base)
{
  // Defining bitmap layout
  const uint32_t NUM_BITMAPS = GPU_FMI_ENTRY_SIZE / GPU_UINT32_LENGTH;
  const uint32_t LUT[12] = {3,7,11,0,1,2,4,5,6,8,9,10};
  // Indexing the FMI entry
  const uint64_t entryIdx       = interval / GPU_FMI_ENTRY_SIZE;
  const uint32_t bitmapPosition = interval % GPU_FMI_ENTRY_SIZE;
  // Gathering the base of the seed
  const uint32_t bit0 =  base & 0x1L;
  const uint32_t bit1 = (base & 0x2L) >> 1;
  const uint32_t missedEntry = (entryIdx % GPU_FMI_ALTERNATE_COUNTERS == bit1) ? 0 : 1;
  const uint64_t bigCounter  = fmi[entryIdx + missedEntry].counters[bit0];
  // Reorder bitmaps layout for low flag
  uint32_t numCharacters = 0;
  for(uint32_t idBitmap = 0; idBitmap < NUM_BITMAPS; ++idBitmap){
    const uint32_t initBitmap  = idBitmap * GPU_FMI_BWT_CHAR_LENGTH;
    uint32_t vbitmap[GPU_FMI_BWT_CHAR_LENGTH] =
        {fmi[entryIdx].bitmaps[LUT[initBitmap]],fmi[entryIdx].bitmaps[LUT[initBitmap + 1]],fmi[entryIdx].bitmaps[LUT[initBitmap + 2]]};
    numCharacters += computeBitmapsCPU(vbitmap, bitmapPosition, bit0, bit1, missedEntry, idBitmap);
  }
  (* new_interval) = (missedEntry) ? bigCounter - numCharacters : bigCounter + numCharacters;
}

void gpu_fmi_table_process_entry(const gpu_fmi_entry_t* const fmi, gpu_sa_entry_t* const currentIntervals,
                                 const gpu_sa_entry_t L, const gpu_sa_entry_t R)
{
  const uint32_t numBases = GPU_FMI_TABLE_ALPHABET_SIZE;
  uint32_t idBase = 0;
  // Compute the next interval & save interval in LUT
  LF_mapping_advance_step(fmi, L, currentIntervals + idBase, idBase);
  // Corner case for the N base (top)
  for(idBase = 0; idBase < numBases; ++idBase)
    LF_mapping_advance_step(fmi, R, currentIntervals + idBase + 1, idBase);
}

gpu_error_t gpu_fmi_table_init_offsets(offset_table_t* const offsetsTableLUT, const uint32_t numLevels)
{
  // Init the fmi-table basic features
  const uint32_t numBases  = GPU_FMI_TABLE_ALPHABET_SIZE;
  uint32_t numElemTableLUT = GPU_FMI_TABLE_MIN_ELEMENTS;
  uint32_t idLevel = 0, numElemPerPartialLevel = 1, numElemPerLevel;
  // Init the offsets for the level 0
  offsetsTableLUT[0].init = 0;
  offsetsTableLUT[0].top  = 1;
  // Process the offsets for the resting levels
  for(idLevel = 1; idLevel < numLevels; ++idLevel){
    numElemPerPartialLevel *= numBases;
    numElemPerLevel = numElemPerPartialLevel + GPU_DIV_CEIL(numElemPerPartialLevel, numBases);
    offsetsTableLUT[idLevel].init = numElemTableLUT;
    offsetsTableLUT[idLevel].top  = numElemTableLUT + numElemPerPartialLevel;
    numElemTableLUT += numElemPerLevel;
  }
  // Succeed
  return(SUCCESS);
}

gpu_error_t gpu_fmi_table_process_backward_level(const gpu_fmi_entry_t* const h_fmi, const uint32_t idLevel,
                                                offset_table_t* const offsetsTableLUT, gpu_sa_entry_t* const fmiTableLUT)
{
  const uint32_t numBases            = GPU_FMI_TABLE_ALPHABET_SIZE;
  const uint32_t initOffsetPrevLevel = offsetsTableLUT[idLevel-1].init, topOffsetPrevLevel = offsetsTableLUT[idLevel-1].top;
  const uint32_t initOffsetCurrLevel = offsetsTableLUT[idLevel].init,   topOffsetCurrLevel = offsetsTableLUT[idLevel].top;
  const uint32_t numElements         = initOffsetCurrLevel - initOffsetPrevLevel;
  const uint32_t numLeftElements     = topOffsetPrevLevel  - initOffsetPrevLevel;
  const uint32_t numRightElements    = initOffsetCurrLevel - topOffsetPrevLevel;
  uint32_t idBase, idEntry;
  // Process all the level in a block way (each block is represented by a base)
  for(idBase = 0; idBase < numBases; ++idBase){
    // Calculates the internal block offsets for this level
    const uint32_t idBaseLeftOffset  = idBase * numLeftElements;
    const uint32_t idBaseRightOffset = idBase * numRightElements;
    // Process the next interval for all entries in the block
    for(idEntry = 0; idEntry < numElements; ++idEntry){
      gpu_sa_entry_t prevInterval, currInterval;
      const uint32_t inIdEntryLUT  = initOffsetPrevLevel + idEntry;
            uint32_t outIdEntryLUT = initOffsetCurrLevel + idBaseLeftOffset + idEntry;
      // Calculates the corner case offsets to keep a multiple of numBases layout
      if(idEntry >= numLeftElements) outIdEntryLUT = topOffsetCurrLevel + idBaseRightOffset + (idEntry - numLeftElements);
      // Process the specific entry for the new level
      prevInterval = fmiTableLUT[inIdEntryLUT];
      LF_mapping_advance_step(h_fmi, prevInterval, &currInterval, idBase);
      fmiTableLUT[outIdEntryLUT] = currInterval;
    }
  }
  // Succeed
  return(SUCCESS);
}

gpu_error_t gpu_fmi_table_get_positions(const uint32_t idLevel, const uint32_t idEntry, const offset_table_t* const offsetsTableLUT,
                                        uint32_t* const idGlobalL, uint32_t* const idGlobalR)
{
  const uint32_t initOffsetCurrLevel = offsetsTableLUT[idLevel].init, topOffsetCurrLevel = offsetsTableLUT[idLevel].top;
  // Creating the hash key to obtain the corresponding FMI interval
  const uint32_t idTableLeft  = idEntry;
  const uint32_t idTableRight = idTableLeft >> 2;
  // Gathering the table entries for L & R (specialized table layout)
  uint32_t idL = initOffsetCurrLevel + idTableLeft;
  uint32_t idR = initOffsetCurrLevel + idTableLeft + 1;
  if((idEntry & GPU_FMI_TABLE_KEY_MASK) == (GPU_FMI_TABLE_ALPHABET_SIZE - 1))
    idR = topOffsetCurrLevel + idTableRight;
  (* idGlobalL) = idL; (* idGlobalR) = idR;
  // Succeed
  return(SUCCESS);
}

gpu_error_t gpu_fmi_table_process_backward_links(const uint32_t idLevel, const uint64_t occThreshold,
                                                 const offset_table_t* const offsetsTableLUT, gpu_sa_entry_t* const fmiTableLUT)
{
  const uint32_t initOffsetCurrLevel = offsetsTableLUT[idLevel].init, topOffsetCurrLevel = offsetsTableLUT[idLevel].top;
  const uint32_t numEntries          = topOffsetCurrLevel - initOffsetCurrLevel;
  uint32_t idEntry;
  // Extracts the superior link mark
  for(idEntry = 0; idEntry < numEntries; ++idEntry){
    uint64_t L, R, occ, linkContent = ((uint64_t)idLevel) << GPU_FMI_TABLE_FIELD_LENGTH;
    uint32_t idL, idR;
    gpu_fmi_table_get_positions(idLevel, idEntry, offsetsTableLUT, &idL, &idR);
    L = fmiTableLUT[idL]; R = fmiTableLUT[idR];
    occ = R - L;
    if(occ <= occThreshold){
      // Extracts the superior link mark
      const uint32_t idParentLevel     = idLevel - 1;
      const uint32_t offsetParentEntry = offsetsTableLUT[idParentLevel].init;
      const uint32_t idParentEntry     = idEntry & (~(GPU_UINT32_ONES << (idParentLevel << 1)));
      linkContent = fmiTableLUT[offsetParentEntry + idParentEntry] & GPU_FMI_TABLE_LINK_MASK;
    }
    fmiTableLUT[idL] |= linkContent;
  }
  // Succeed
  return(SUCCESS);
}

gpu_error_t gpu_fmi_table_process_forward_links(const uint32_t idLevel, const uint64_t occThreshold,
                                                const offset_table_t* const offsetsTableLUT, gpu_sa_entry_t* const fmiTableLUT)
{
  const uint32_t numBases            = GPU_FMI_TABLE_ALPHABET_SIZE;
  const uint32_t initOffsetCurrLevel = offsetsTableLUT[idLevel].init, topOffsetCurrLevel = offsetsTableLUT[idLevel].top;
  const uint32_t numEntries          = topOffsetCurrLevel - initOffsetCurrLevel;

  uint32_t idEntry;
  //Extracts the superior link mark
  for(idEntry = 0; idEntry < numEntries; ++idEntry){
    uint64_t L, R, occ, linkContent = ((uint64_t)idLevel) << GPU_FMI_TABLE_FIELD_LENGTH;
    uint32_t idL, idR;
    gpu_fmi_table_get_positions(idLevel, idEntry, offsetsTableLUT, &idL, &idR);
    L = fmiTableLUT[idL]; R = fmiTableLUT[idR];
    occ = R - L;
    if(occ <= occThreshold){
      //Extracts the superior link mark
      const uint32_t idParentLevel     = idLevel - 1;
      const uint32_t idParentEntry     = idEntry >> (numBases >> 1);
      const uint32_t offsetParentEntry = offsetsTableLUT[idParentLevel].init;
      linkContent = fmiTableLUT[offsetParentEntry + idParentEntry] & GPU_FMI_TABLE_LINK_MASK;
    }
    fmiTableLUT[idL] |= linkContent;
  }
  // Succeed
  return(SUCCESS);
}

gpu_error_t gpu_fmi_table_process_forward_level(const gpu_fmi_entry_t* const fmi, const uint32_t idLevel,
                                                offset_table_t* const offsetsTableLUT, gpu_sa_entry_t* const fmiTableLUT)
{
  const uint32_t numBases            = GPU_FMI_TABLE_ALPHABET_SIZE;
  const uint32_t initOffsetCurrLevel = offsetsTableLUT[idLevel-1].init, topOffsetCurrLevel = offsetsTableLUT[idLevel-1].top;
  const uint32_t initOffsetNextLevel = offsetsTableLUT[idLevel].init,   topOffsetNextLevel = offsetsTableLUT[idLevel].top;
  const uint32_t numElementsLevel    = topOffsetCurrLevel - initOffsetCurrLevel;
  // Generates the next LUT level (GPU_FMI_TABLE_ALPHABET_SIZE+1 elements per entry)
  uint32_t idEntry, idBase;
  for(idEntry = 0; idEntry < numElementsLevel; ++idEntry){
    gpu_sa_entry_t currentIntervals[numBases + 1], L, R;
    uint32_t idL = initOffsetCurrLevel + idEntry, idR = initOffsetCurrLevel + idEntry + 1;
    if((idEntry % numBases) == (numBases - 1))
      idR = topOffsetCurrLevel + (idEntry / numBases);
    // Generates the next level intervals
    L = fmiTableLUT[idL]; R = fmiTableLUT[idR];
    gpu_fmi_table_process_entry(fmi, currentIntervals, L, R);
    // Store the next level intervals into the LUT
    for(idBase = 0; idBase < numBases; ++idBase)
      fmiTableLUT[initOffsetNextLevel + (idEntry * numBases + idBase)] = currentIntervals[idBase];
    fmiTableLUT[topOffsetNextLevel + idEntry] = currentIntervals[idBase];
  }
  // Succeed
  return(SUCCESS);
}

gpu_error_t gpu_fmi_table_get_num_elements(const uint32_t numLevels, uint32_t* const totalElemTableLUT)
{
  // Init the fmi-table basic features
  const uint32_t numBases  = GPU_FMI_TABLE_ALPHABET_SIZE;
  uint32_t numElemTableLUT = GPU_FMI_TABLE_MIN_ELEMENTS;
  uint32_t idLevel = 0, numElemPerPartialLevel = 1, numElemPerLevel;
  // Process the offsets for the resting levels
  for(idLevel = 1; idLevel < numLevels; ++idLevel){
    numElemPerPartialLevel *= numBases;
    numElemPerLevel = numElemPerPartialLevel + GPU_DIV_CEIL(numElemPerPartialLevel, numBases);
    numElemTableLUT += numElemPerLevel;
  }
  (* totalElemTableLUT) = numElemTableLUT;
  // Succeed
  return(SUCCESS);
}

gpu_error_t gpu_fmi_table_init_dto(gpu_fmi_table_t* const fmiTable)
{
  //Initialize the FMI index structure
  fmiTable->formatTableLUT       = GPU_FMI_TABLE_DISABLED;
  fmiTable->d_offsetsTableLUT    = NULL;
  fmiTable->h_offsetsTableLUT    = NULL;
  fmiTable->d_fmiTableLUT        = NULL;
  fmiTable->h_fmiTableLUT        = NULL;
  fmiTable->hostAllocStats       = GPU_PAGE_UNLOCKED;
  fmiTable->memorySpace          = NULL;
  fmiTable->maxLevelsTableLUT    = 0;
  fmiTable->skipLevelsTableLUT   = 0;
  fmiTable->totalElemTableLUT    = 0;
  fmiTable->occThresholdTableLUT = 0;
  // Succeed
  return (SUCCESS);
}

gpu_error_t gpu_fmi_table_init(gpu_fmi_table_t* const fmiTable, const uint32_t numLevels, const uint32_t numSupportedDevices)
{
  uint32_t idSupDevice, totalElemTableLUT;
  GPU_ERROR(gpu_fmi_table_init_dto(fmiTable));
  GPU_ERROR(gpu_fmi_table_get_num_elements(numLevels, &totalElemTableLUT));

  fmiTable->maxLevelsTableLUT = numLevels;
  fmiTable->totalElemTableLUT = totalElemTableLUT;
  fmiTable->formatTableLUT    = GPU_FMI_TABLE_MULTILEVEL_LINKED;

  fmiTable->d_offsetsTableLUT = (offset_table_t **) malloc(numSupportedDevices * sizeof(offset_table_t *));
  if (fmiTable->d_offsetsTableLUT == NULL) GPU_ERROR(E_ALLOCATE_MEM);
  fmiTable->d_fmiTableLUT = (gpu_sa_entry_t **) malloc(numSupportedDevices * sizeof(gpu_sa_entry_t *));
  if (fmiTable->d_fmiTableLUT == NULL) GPU_ERROR(E_ALLOCATE_MEM);
  fmiTable->memorySpace = (memory_alloc_t *) malloc(numSupportedDevices * sizeof(memory_alloc_t));
  if (fmiTable->memorySpace == NULL) GPU_ERROR(E_ALLOCATE_MEM);

  for(idSupDevice = 0; idSupDevice < numSupportedDevices; ++idSupDevice){
    fmiTable->d_offsetsTableLUT[idSupDevice] = NULL;
    fmiTable->d_fmiTableLUT[idSupDevice]     = NULL;
    fmiTable->memorySpace[idSupDevice]       = GPU_NONE_MAPPED;
  }
  // Succeed
  return (SUCCESS);
}

gpu_error_t gpu_fmi_table_allocate(gpu_fmi_table_t* const fmiTable)
{
  // Init the fmi-table specifications
  GPU_ERROR(gpu_fmi_table_get_num_elements(fmiTable->maxLevelsTableLUT, &fmiTable->totalElemTableLUT));
  if(fmiTable->hostAllocStats & GPU_PAGE_LOCKED){
    // Allocate the metadata used in the LUT fmi-table
    CUDA_ERROR(cudaHostAlloc((void**) &fmiTable->h_offsetsTableLUT, fmiTable->maxLevelsTableLUT * sizeof(offset_table_t), cudaHostAllocMapped));
    // Allocate the space used for the LUT fmi-table
    CUDA_ERROR(cudaHostAlloc((void**) &fmiTable->h_fmiTableLUT, fmiTable->totalElemTableLUT * sizeof(gpu_sa_entry_t), cudaHostAllocMapped));
  }else{
    // Allocate the metadata used in the LUT fmi-table
    fmiTable->h_offsetsTableLUT = (offset_table_t*) malloc(fmiTable->maxLevelsTableLUT * sizeof(offset_table_t));
    if (fmiTable->h_offsetsTableLUT == NULL) return (E_ALLOCATE_MEM);
    // Allocate the space used for the LUT fmi-table
    fmiTable->h_fmiTableLUT = (gpu_sa_entry_t*) malloc(fmiTable->totalElemTableLUT * sizeof(gpu_sa_entry_t));
    if (fmiTable->h_fmiTableLUT == NULL) return (E_ALLOCATE_MEM);
  }
  // Succeed
  return(SUCCESS);
}

gpu_error_t gpu_fmi_table_load_default_specs(gpu_fmi_table_t* const fmiTable)
{
  if (fmiTable->maxLevelsTableLUT  == 0) fmiTable->maxLevelsTableLUT  = GPU_FMI_TABLE_DEFAULT_LEVELS;
  if (fmiTable->skipLevelsTableLUT == 0) fmiTable->skipLevelsTableLUT = GPU_FMI_TABLE_DEFAULT_SKIP_LEVELS;
  GPU_ERROR(gpu_fmi_table_get_num_elements(fmiTable->maxLevelsTableLUT, &fmiTable->totalElemTableLUT));
  return (SUCCESS);
}

gpu_error_t gpu_fmi_table_read_specs(int fp, gpu_fmi_table_t* const fmiTable)
{
  size_t result, bytesRequest;
  // Read the specifications of the LUT fmi-table
  bytesRequest = sizeof(uint32_t);
  result = read(fp, (void *)&fmiTable->maxLevelsTableLUT, bytesRequest);
  if (result != bytesRequest) return (E_READING_FILE);
  bytesRequest = sizeof(uint32_t);
  result = read(fp, (void *)&fmiTable->skipLevelsTableLUT, bytesRequest);
  if (result != bytesRequest) return (E_READING_FILE);
  bytesRequest = sizeof(uint32_t);
  result = read(fp, (void *)&fmiTable->totalElemTableLUT, bytesRequest);
  if (result != bytesRequest) return (E_READING_FILE);
  // Succeed
  return (SUCCESS);
}

gpu_error_t gpu_fmi_table_read(int fp, gpu_fmi_table_t* const fmiTable)
{
  size_t result, bytesRequest, numBytesRequested = 0;
  uint64_t idRequest, numRequests;
  // Read the metadata used in the LUT fmi-table
  bytesRequest = sizeof(offset_table_t) * fmiTable->maxLevelsTableLUT;
  result = read(fp, (void* )fmiTable->h_offsetsTableLUT, bytesRequest);
  if (result != bytesRequest) return (E_READING_FILE);
  // Read the LUT fmi-table
  bytesRequest = sizeof(gpu_sa_entry_t) * fmiTable->totalElemTableLUT;
  numRequests  = GPU_DIV_CEIL(bytesRequest, GPU_FILE_SIZE_BLOCK);
  for(idRequest = 0; idRequest < numRequests; ++idRequest){
    const size_t requestSize = GPU_MIN(GPU_FILE_SIZE_BLOCK, bytesRequest - numBytesRequested);
    result = read(fp, (void* )fmiTable->h_fmiTableLUT + numBytesRequested, requestSize);
    if (result != requestSize) return (E_READING_FILE);
    numBytesRequested += requestSize;
  }
  // Succeed
  return (SUCCESS);
}

gpu_error_t gpu_fmi_table_write_specs(int fp, const gpu_fmi_table_t* const fmiTable)
{
  size_t result, bytesRequest;
  // Write the specifications of the LUT fmi-table
  bytesRequest = sizeof(uint32_t);
  result = write(fp, (void *)&fmiTable->maxLevelsTableLUT, bytesRequest);
  if (result != bytesRequest) return (E_WRITING_FILE);
  bytesRequest = sizeof(uint32_t);
  result = write(fp, (void *)&fmiTable->skipLevelsTableLUT, bytesRequest);
  if (result != bytesRequest) return (E_WRITING_FILE);
  bytesRequest = sizeof(uint32_t);
  result = write(fp, (void *)&fmiTable->totalElemTableLUT, bytesRequest);
  if (result != bytesRequest) return (E_WRITING_FILE);
  // Succeed
  return (SUCCESS);
}

gpu_error_t gpu_fmi_table_write(int fp, const gpu_fmi_table_t* const fmiTable)
{
  size_t result, bytesRequest, numBytesRequested = 0;
  uint64_t idRequest, numRequests;
  // Write the metadata used in the LUT fmi-table
  bytesRequest = sizeof(offset_table_t) * fmiTable->maxLevelsTableLUT;
  result = write(fp, (void* )fmiTable->h_offsetsTableLUT, bytesRequest);
  if (result != bytesRequest) return (E_WRITING_FILE);
  // Write the LUT fmi-table
  bytesRequest = sizeof(gpu_sa_entry_t) * fmiTable->totalElemTableLUT;
  numRequests  = GPU_DIV_CEIL(bytesRequest, GPU_FILE_SIZE_BLOCK);
  for(idRequest = 0; idRequest < numRequests; ++idRequest){
    const size_t requestSize = GPU_MIN(GPU_FILE_SIZE_BLOCK, bytesRequest - numBytesRequested);
    result = write(fp, (void* )fmiTable->h_fmiTableLUT + numBytesRequested, requestSize);
    if (result != requestSize) return (E_WRITING_FILE);
    numBytesRequested += requestSize;
  }
  // Succeed
  return (SUCCESS);
}

gpu_error_t gpu_fmi_table_transfer_CPU_to_GPUs(gpu_fmi_table_t* const fmiTable, gpu_device_info_t** const devices)
{
  uint32_t deviceFreeMemory, idSupportedDevice;
  uint32_t numSupportedDevices = devices[0]->numSupportedDevices;

  for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
    if(fmiTable->memorySpace[idSupportedDevice] == GPU_DEVICE_MAPPED){
      const size_t cpySizeMeta = fmiTable->maxLevelsTableLUT * sizeof(offset_table_t);
      const size_t cpySizeLUT  = fmiTable->totalElemTableLUT * sizeof(gpu_sa_entry_t);
      deviceFreeMemory = gpu_device_get_free_memory(devices[idSupportedDevice]->idDevice);
      if ((GPU_CONVERT__B_TO_MB(cpySizeMeta + cpySizeLUT)) > deviceFreeMemory)
        return(E_INSUFFICIENT_MEM_GPU);
      CUDA_ERROR(cudaSetDevice(devices[idSupportedDevice]->idDevice));
      //Synchronous allocate & transfer the FM-index to the GPU
      CUDA_ERROR(cudaMalloc((void**) &fmiTable->d_offsetsTableLUT[idSupportedDevice], cpySizeMeta));
      CUDA_ERROR(cudaMemcpy(fmiTable->d_offsetsTableLUT[idSupportedDevice], fmiTable->h_offsetsTableLUT, cpySizeMeta, cudaMemcpyHostToDevice));
      //Synchronous allocate & transfer the FM-index to the GPU
      CUDA_ERROR(cudaMalloc((void**) &fmiTable->d_fmiTableLUT[idSupportedDevice], cpySizeLUT));
      CUDA_ERROR(cudaMemcpy(fmiTable->d_fmiTableLUT[idSupportedDevice], fmiTable->h_fmiTableLUT, cpySizeLUT, cudaMemcpyHostToDevice));
    }else{
      fmiTable->d_offsetsTableLUT[idSupportedDevice] = fmiTable->h_offsetsTableLUT;
      fmiTable->d_fmiTableLUT[idSupportedDevice]     = fmiTable->h_fmiTableLUT;
    }
  }
  // Succeed
  return (SUCCESS);
}

/* (numLevels, Size) => (8, 0.8MB) (9, 3.2MB) (10, 13MB) (11, 52MB) (12, 210MB) (13, 840MB) */
gpu_error_t gpu_fmi_table_build(gpu_fmi_table_t* const fmiTable, const gpu_fmi_entry_t* const h_fmi, const uint64_t bwtSize)
{
  // Get FMI table specifications
  const uint32_t maxLevels = fmiTable->maxLevelsTableLUT;
  uint32_t idLevel = 0;
  // Pre-calculate internal offsets
  GPU_ERROR(gpu_fmi_table_init_offsets(fmiTable->h_offsetsTableLUT, maxLevels));
  // Generates and fills 0th LUT level
  fmiTable->h_fmiTableLUT[0] = 0;
  fmiTable->h_fmiTableLUT[1] = bwtSize;
  idLevel++;
  // Generates and fills 1st LUT level
  if(idLevel < maxLevels)
    GPU_ERROR(gpu_fmi_table_process_forward_level(h_fmi, idLevel, fmiTable->h_offsetsTableLUT, fmiTable->h_fmiTableLUT));
  idLevel++;
  // Generates and fills nth LUT levels
  for(; idLevel < maxLevels; ++idLevel)
    GPU_ERROR(gpu_fmi_table_process_backward_level(h_fmi, idLevel, fmiTable->h_offsetsTableLUT, fmiTable->h_fmiTableLUT));
  // Succeed
  return(SUCCESS);
}

gpu_error_t gpu_fmi_table_process_links(const uint32_t numLevels, const uint64_t occThreshold, const offset_table_t* const offsetsTableLUT, gpu_sa_entry_t* const fmiTableLUT)
{
  uint32_t idLevel = 1;
  // Attaching the 1st LUT level links
  GPU_ERROR(gpu_fmi_table_process_forward_links(idLevel, occThreshold, offsetsTableLUT, fmiTableLUT));
  // Attaching the nth LUT level links
  for(idLevel = 2; idLevel < numLevels; ++idLevel)
    GPU_ERROR(gpu_fmi_table_process_backward_links(idLevel, occThreshold, offsetsTableLUT, fmiTableLUT));
  // Succeed
  return(SUCCESS);
}

gpu_error_t gpu_fmi_table_print_links(const offset_table_t* const offsetsTableLUT, gpu_sa_entry_t* const fmiTableLUT)
{
  uint32_t idEntry, idLevel, maxLevel = 10;
  for(idLevel = 0; idLevel < maxLevel; ++idLevel){
    const uint32_t numEntries = offsetsTableLUT[idLevel + 1].init - offsetsTableLUT[idLevel].init;
    const uint32_t offEntries = offsetsTableLUT[idLevel].init;
    printf("====== idLevel=%u - numEntries=%u ====== \n", idLevel, numEntries);
    for(idEntry = 0; idEntry < numEntries; ++idEntry){
       const long long unsigned int idEntryRaw      = fmiTableLUT[offEntries + idEntry];
       const long long unsigned int idEntryLink     = (idEntryRaw & GPU_FMI_TABLE_LINK_MASK) >> GPU_FMI_TABLE_FIELD_LENGTH;
       const long long unsigned int idEntryInterval = idEntryRaw & GPU_FMI_TABLE_FIELD_MASK;
       printf("[%u] idEntry=%llu idEntryLink=%llu idEntryInterval=%llu \n", idEntry, idEntryRaw, idEntryLink, idEntryInterval);
    }
  }
  // Succeed
  return(SUCCESS);
}

gpu_error_t gpu_fmi_table_construction(gpu_fmi_table_t* const fmiTable, const gpu_fmi_entry_t* const h_fmi, const uint64_t bwtSize)
{
  // Get FMI table specifications
  const uint32_t maxLevels                    = fmiTable->maxLevelsTableLUT;
  const uint32_t occThreshold                 = fmiTable->occThresholdTableLUT;
  const offset_table_t* const offsetsTableLUT = fmiTable->h_offsetsTableLUT;
  gpu_sa_entry_t* const fmiTableLUT           = fmiTable->h_fmiTableLUT;
  // Build FMI Table & attach the multilevel links
  if(maxLevels >= GPU_FMI_TABLE_MIN_LEVELS){
    GPU_ERROR(gpu_fmi_table_build(fmiTable, h_fmi, bwtSize));
    GPU_ERROR(gpu_fmi_table_process_links(maxLevels, occThreshold, offsetsTableLUT, fmiTableLUT));
    //GPU_ERROR(gpu_fmi_table_print_links(offsetsTableLUT, fmiTableLUT));
    // Succeed
    return(SUCCESS);
  }else{
    return(E_FMI_TABLE_INCOMPATIBLE_SIZE);
  }
}

gpu_error_t gpu_fmi_table_get_size(const gpu_fmi_table_t* const fmiTable, size_t* const bytesPerFmiTable)
{
  const size_t fmIndexTableSize        = fmiTable->totalElemTableLUT * sizeof(gpu_sa_entry_t);
  const size_t fmIndexTableOffsetsSize = fmiTable->maxLevelsTableLUT * sizeof(offset_table_t);

  (* bytesPerFmiTable) = fmIndexTableSize + fmIndexTableOffsetsSize;
  return (SUCCESS);
}

gpu_error_t gpu_fmi_table_free_host(gpu_fmi_table_t* const fmiTable)
{
  if(fmiTable->h_offsetsTableLUT != NULL){
    if(fmiTable->hostAllocStats == GPU_PAGE_LOCKED)
      CUDA_ERROR(cudaFreeHost(fmiTable->h_offsetsTableLUT));
    else
      free(fmiTable->h_offsetsTableLUT);
    fmiTable->h_offsetsTableLUT = NULL;
  }

  if(fmiTable->h_fmiTableLUT != NULL){
    if(fmiTable->hostAllocStats == GPU_PAGE_LOCKED)
      CUDA_ERROR(cudaFreeHost(fmiTable->h_fmiTableLUT));
    else
      free(fmiTable->h_fmiTableLUT);
    fmiTable->h_fmiTableLUT = NULL;
  }
  // Succeed
  return(SUCCESS);
}

gpu_error_t gpu_fmi_table_free_unused_host(gpu_fmi_table_t* const fmiTable, gpu_device_info_t** const devices)
{
  uint32_t idSupportedDevice, numSupportedDevices;
  bool indexInHostSideUsed = false;

  numSupportedDevices = devices[0]->numSupportedDevices;
  //Free all the unused references in the host side
  for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
    if(fmiTable->memorySpace[idSupportedDevice] == GPU_HOST_MAPPED) indexInHostSideUsed = true;
  }

  if(!indexInHostSideUsed){
    GPU_ERROR(gpu_fmi_table_free_host(fmiTable));
  }

  return(SUCCESS);
}

gpu_error_t gpu_fmi_table_free_device(gpu_fmi_table_t* const fmiTable, gpu_device_info_t** const devices)
{
  const uint32_t numSupportedDevices = devices[0]->numSupportedDevices;
  uint32_t idSupportedDevice;

  //Free all the references in the devices
  for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
    CUDA_ERROR(cudaSetDevice(devices[idSupportedDevice]->idDevice));
    if(fmiTable->d_offsetsTableLUT[idSupportedDevice] != NULL){
      if(fmiTable->memorySpace[idSupportedDevice] == GPU_DEVICE_MAPPED)
        CUDA_ERROR(cudaFree(fmiTable->d_offsetsTableLUT[idSupportedDevice]));
      fmiTable->d_offsetsTableLUT[idSupportedDevice] = NULL;
    }
    if(fmiTable->d_fmiTableLUT[idSupportedDevice] != NULL){
      if(fmiTable->memorySpace[idSupportedDevice] == GPU_DEVICE_MAPPED)
        CUDA_ERROR(cudaFree(fmiTable->d_fmiTableLUT[idSupportedDevice]));
      fmiTable->d_fmiTableLUT[idSupportedDevice] = NULL;
    }
  }

  //Free the index list
  if(fmiTable->d_offsetsTableLUT != NULL){
    free(fmiTable->d_offsetsTableLUT);
    fmiTable->d_offsetsTableLUT = NULL;
  }

  //Free the index list
  if(fmiTable->d_fmiTableLUT != NULL){
    free(fmiTable->d_fmiTableLUT);
    fmiTable->d_fmiTableLUT = NULL;
  }

  return(SUCCESS);
}

gpu_error_t gpu_fmi_table_free_metainfo(gpu_fmi_table_t* const fmiTable)
{
  //Free the index list
  if(fmiTable->memorySpace != NULL){
    free(fmiTable->memorySpace);
    fmiTable->memorySpace = NULL;
  }
  return(SUCCESS);
}


#endif /* GPU_FMI_TABLE_C_ */

