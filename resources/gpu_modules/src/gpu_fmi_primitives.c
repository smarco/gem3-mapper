/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_FMI_PRIMITIVES_C_
#define GPU_FMI_PRIMITIVES_C_

#include "../include/gpu_fmi_primitives.h"
#include "../include/gpu_sa_primitives.h"

/************************************************************
Functions to get the GPU FMI buffers
************************************************************/

gpu_fmi_entry_t* gpu_fmi_buffer_get_index_(const void* const fmiBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) fmiBuffer;
  return(mBuff->index->fmi.h_fmi);
}

/************************************************************
Debug functions for the index
************************************************************/
/*
char gpu_fmi_search_bin_to_char(const uint32_t base)
{
  switch(base){
    case GPU_ENC_DNA_CHAR_A:
      return('A');
    case GPU_ENC_DNA_CHAR_C:
      return('C');
    case GPU_ENC_DNA_CHAR_G:
      return('G');
    case GPU_ENC_DNA_CHAR_T:
      return('T');
    default :
      return('X');
  }
}

uint32_t gpu_fmi_search_print_seed(const gpu_fmi_search_seed_t seed, const uint32_t seedSize)
{
  char plainSeed[GPU_FMI_SEED_MAX_CHARS] = {0};
  uint64_t bitmap = seed.hi;
  uint32_t idBase;

  for(idBase = 0; idBase < seedSize; ++idBase){
    uint32_t base = bitmap & 0x3;
    plainSeed[idBase] = gpu_fmi_search_bin_to_char(base);
    bitmap >>= GPU_FMI_SEED_CHAR_LENGTH;
    if(idBase == GPU_UINT32_LENGTH) bitmap = seed.low;
  }

  for(idBase = 0; idBase < seedSize; ++idBase)
    printf("%c", plainSeed[seedSize - idBase - 1]);

  return(SUCCESS);
}

uint32_t flag_print = 0;
uint32_t gpu_fmi_search_print_buffer(const void* const fmiBuffer)
{
  if(flag_print == 0){
    gpu_buffer_t* const mBuff  = (gpu_buffer_t *) fmiBuffer;
    const uint32_t maxSeeds    = 100; // Just check and print the first results
    const uint32_t numSeeds    = mBuff->data.search.seeds.numSeeds;
          uint32_t missMatches = 0;
          uint32_t idSeed;

    printf("Buffer: %d ------------------------------------\n", mBuff->idBuffer);
    for(idSeed = 0; idSeed < numSeeds; ++idSeed){
      const uint64_t hiSeedSection = mBuff->data.search.seeds.h_seeds[idSeed].low;
      const uint32_t seedSize = hiSeedSection >> (GPU_UINT64_LENGTH - GPU_FMI_SEED_FIELD_SIZE);
      printf("[%d] seed=", idSeed);
      gpu_fmi_search_print_seed(mBuff->data.search.seeds.h_seeds[idSeed], seedSize);
      printf("\t size=%d \t (GPU) lo=%lu \t hi=%lu \n",
             seedSize,
             mBuff->data.search.saIntervals.h_intervals[idSeed].low,
             mBuff->data.search.saIntervals.h_intervals[idSeed].hi);
      missMatches++;
    }
    printf("Buffer: %d ------------------------------------\n", mBuff->idBuffer);
    flag_print = 1;
  }
  return (SUCCESS);
}*/


#endif /* GPU_FMI_PRIMITIVES_C_ */

