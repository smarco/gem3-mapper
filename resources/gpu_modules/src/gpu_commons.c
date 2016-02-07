#ifndef GPU_COMMONS_C_
#define GPU_COMMONS_C_

#include "../include/gpu_commons.h"

uint32_t gpu_bit_reverse(uint32_t a)
{
    uint32_t t;
    a = (a << 15) | (a >> 17);
    t = (a ^ (a >> 10)) & 0x003f801f;
    a = (t + (t << 10)) ^ a;
    t = (a ^ (a >>  4)) & 0x0e038421;
    a = (t + (t <<  4)) ^ a;
    t = (a ^ (a >>  2)) & 0x22488842;
    a = (t + (t <<  2)) ^ a;
    return a;
}

uint32_t gpu_gen_mask(const int32_t shift)
{
  uint32_t mask = GPU_UINT32_ONES << (GPU_UINT32_LENGTH - shift);
  mask = (shift > GPU_UINT32_LENGTH) ? GPU_UINT32_ONES : mask;
  mask = (shift > 0) ? mask : GPU_UINT32_ZEROS;
  return(mask);
}

#endif /* GPU_COMMONS_C_ */
