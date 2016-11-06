/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

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

uint32_t gpu_count_active_bits(uint32_t a)
{
  a = a - ((a >> 1) & 0x55555555);
  a = (a & 0x33333333) + ((a >> 2) & 0x33333333);
  return (((a + (a >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

uint32_t gpu_gen_mask(const int32_t shift)
{
  uint32_t mask = GPU_UINT32_ONES << (GPU_UINT32_LENGTH - shift);
  mask = (shift > GPU_UINT32_LENGTH) ? GPU_UINT32_ONES : mask;
  mask = (shift > 0) ? mask : GPU_UINT32_ZEROS;
  return(mask);
}

uint8_t gpu_base2log(uint16_t value)
{
    uint8_t result = 0;
    while (value >>= 1)
        result++;
    return(result);
}

bool gpu_is_pow_two(uint32_t value)
{
    return ((value != 0) && ((value & (value - 1)) == 0));
}

#endif /* GPU_COMMONS_C_ */
