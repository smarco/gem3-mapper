/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_TEXT_CORE_H_
#define GPU_TEXT_CORE_H_

extern "C" {
#include "gpu_commons.h"
}

__device__ const ulong2 GPU_TEXT_INIT = {GPU_UINT64_ZEROS, GPU_UINT64_ONES};


GPU_INLINE __device__ uint8_t gpu_text_lookup(const uint64_t* const text, const uint64_t textPosition, ulong2* const globalInfo, const uint32_t BASE_TEXT_LENGTH)
{
  uint8_t BASES_PER_TEXT_ENTRY = GPU_UINT64_LENGTH / BASE_TEXT_LENGTH;
  uint8_t BASE_TEXT_MASK       = GPU_UINT64_ONES >> (GPU_UINT64_LENGTH - BASE_TEXT_LENGTH);
  //Text meta information to reduce the main memory requests
  uint64_t infoBase = globalInfo->x; // Packet BASES_PER_TEXT_ENTRY bases (cached)
  uint64_t infoText = globalInfo->y; // Position cached
  //Address the cached block request (hash function)
  const uint64_t idEntryText = textPosition / BASES_PER_TEXT_ENTRY;
  const uint64_t idIntraText = textPosition % BASES_PER_TEXT_ENTRY;
  uint32_t shiftBits, base;
  //Make the query request if data was not cached
  if((infoText != idEntryText) || (infoText == GPU_UINT64_ONES)){
    infoBase = text[idEntryText];
    infoText = idEntryText;
  }
  //Logical base extraction from cached text requested
  shiftBits = idIntraText * BASE_TEXT_LENGTH;
  base      = (infoBase >> shiftBits) & BASE_TEXT_MASK;
  //Return the requested data
  globalInfo->x = infoBase; //Contain text cached data
  globalInfo->y = infoText; //Contain the text position cached data
  return(base);             //Requested base
}

/*
GPU_INLINE __device__ void gpu_text_update(uint64_t* const text, const uint8_t base, const uint64_t textPosition, ulong2* const globalInfo, const uint32_t BASE_TEXT_LENGTH)
{
  uint8_t BASES_PER_TEXT_ENTRY = GPU_UINT64_LENGTH / BASE_TEXT_LENGTH;
  uint8_t BASE_TEXT_MASK       = GPU_UINT64_ONES >> (GPU_UINT64_LENGTH - BASE_TEXT_LENGTH);
  //Text meta information to reduce the main memory requests
  uint64_t infoBase = globalInfo->x; // Packet BASES_PER_TEXT_ENTRY bases (cached)
  uint64_t infoText = globalInfo->y; // Position cached
  //Address the cached block request (hash function)
  const uint64_t idEntryText = refText / BASES_PER_TEXT_ENTRY;
  const uint64_t idIntraText = refText % BASES_PER_TEXT_ENTRY;
  uint32_t shiftBits, base;
  //Make the query request if data was not cached
  if((infoText != idEntry) || (infoText == GPU_UINT64_ONES)){
    infoBase = text[idEntryText];
    infoText = idEntryText;
  }
  //Logical base extraction from cached text requested
  shiftBits = idIntraText * BASE_TEXT_LENGTH;
  base      = (infoBase >> shiftBits) & BASE_TEXT_MASK;
  //Return the requested data
  globalInfo->x = infoBase; //Contain text cached data
  globalInfo->y = infoText; //Contain the text position cached data
  return(base);             //Requested base
}*/


#endif /* GPU_TEXT_CORE_H_ */

