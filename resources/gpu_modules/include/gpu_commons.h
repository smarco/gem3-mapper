/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_COMMONS_H_
#define GPU_COMMONS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <cuda.h>
#include <cuda_runtime.h>

/********************************
Common constants for Device & Host
*********************************/

#define GPU_UINT8_LENGTH          8
#define GPU_UINT32_LENGTH         32
#define GPU_UINT64_LENGTH         64
#define GPU_UINT64_MAX_VALUE      ULONG_MAX
#define GPU_UINT32_ZEROS          0x00000000u
#define GPU_UINT32_ONES           0xFFFFFFFFu
#define GPU_UINT32_MASK_ONE_LOW   0x00000001u
#define GPU_UINT32_MASK_ONE_HIGH  0x80000000u
#define GPU_UINT8_ONES            0xFF
#define GPU_UINT64_ONES           0xFFFFFFFFFFFFFFFFu
#define GPU_UINT8_SIZE            1
#define GPU_UINT32_SIZE           4
#define GPU_UINT64_SIZE           8

#define GPU_INLINE                inline

/* Functions inline */
#define GPU_SELECT_OFFSET(NUM_A,NUM_B)      ((NUM_A/NUM_B)*NUM_B)
#define GPU_DIV_CEIL(NUMERATOR,DENOMINATOR) (((NUMERATOR)+((DENOMINATOR)-1))/(DENOMINATOR))
#define GPU_ROUND(NUM)                      ((int)((NUM) < 0 ? ((NUM) - 0.5) : ((NUM) + 0.5)))
#define GPU_MIN(NUM_A,NUM_B)                (((NUM_A) < (NUM_B)) ? (NUM_A) : (NUM_B))
#define GPU_MAX(NUM_A,NUM_B)                (((NUM_A) > (NUM_B)) ? (NUM_A) : (NUM_B))
#define GPU_ALIGN_TO(ADDR,BYTES)            ((void*)(GPU_DIV_CEIL((uint64_t)ADDR,BYTES) * BYTES))
#define BASE2LOG(NUM)                       (base2log(NUM))

/* Conversion utils */
#define GPU_CONVERT_BYTES_TO_BITS(number) ((number) * 8)
#define GPU_CONVERT_BITS_TO_BYTES(number) ((number) / 8)

#define GPU_CONVERT__B_TO_KB(number) ((number) / (1024))
#define GPU_CONVERT__B_TO_MB(number) ((number) / (1024 * 1024))
#define GPU_CONVERT__B_TO_GB(number) ((number) / (1024 * 1024 * 1024))

#define GPU_CONVERT_KB_TO__B(number) ((number) * (1024))
#define GPU_CONVERT_KB_TO_MB(number) ((number) / (1024))
#define GPU_CONVERT_KB_TO_GB(number) ((number) / (1024 * 1024))

#define GPU_CONVERT_MB_TO__B(number) ((number) * (1024 * 1024))
#define GPU_CONVERT_MB_TO_KB(number) ((number) * (1024))
#define GPU_CONVERT_MB_TO_GB(number) ((number) / (1024))

#define GPU_CONVERT_GB_TO__B(number) ((number) * (1024 * 1024 * 1024))
#define GPU_CONVERT_GB_TO_KB(number) ((number) * (1024 * 1024))
#define GPU_CONVERT_GB_TO_MB(number) ((number) * (1024))

#define GPU_CONVERT__HZ_TO_GHZ(number) ((number) / (1000 * 1000 * 1000))
#define GPU_CONVERT_KHZ_TO_GHZ(number) ((number) / (1000 * 1000))
#define GPU_CONVERT_MHZ_TO_GHZ(number) ((number) / (1000))

/* System */
#define GPU_FILE_SIZE_LINES   250
#define GPU_FILE_SIZE_BLOCK   (256 * 1024 * 1024)

/* Encoded DNA Nucleotides */
#define GPU_ENC_DNA_CHAR_A    0LL
#define GPU_ENC_DNA_CHAR_C    1LL
#define GPU_ENC_DNA_CHAR_G    2LL
#define GPU_ENC_DNA_CHAR_T    3LL

#define GPU_ENC_DNA_CHAR_N    4LL
#define GPU_ENC_DNA_CHAR_X    4LL
#define GPU_ENC_DNA_CHAR_SEP  5LL
#define GPU_ENC_DNA_CHAR_JUMP 6LL

#include "../gpu_interface.h"
#include "gpu_errors.h"
#include "gpu_sample.h"

uint32_t gpu_bit_reverse(uint32_t a);
uint32_t gpu_gen_mask(const int32_t shift);
uint32_t gpu_count_active_bits(uint32_t a);
uint8_t  base2log(const uint16_t value);

#endif /* GPU_COMMONS_H_ */
