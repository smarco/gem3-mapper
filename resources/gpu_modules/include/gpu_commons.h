/*
 * PROJECT: Bit-Parallel Myers on GPU
 * FILE: myers-interface.h
 * DATE: 4/7/2014
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 * DESCRIPTION: Common headers and data structures for BPM on GPU library
 */


#ifndef GPU_COMMONS_H_
#define GPU_COMMONS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include <cuda.h>
#include <cuda_runtime.h>

/********************************
Common constants for Device & Host
*********************************/

#define	GPU_UINT8_LENGTH			8
#define	GPU_UINT32_LENGTH			32
#define	GPU_UINT64_LENGTH			64
#define	GPU_UINT64_MAX_VALUE		ULONG_MAX
#define GPU_UINT32_ZEROS			0x00000000u
#define GPU_UINT32_ONES  			0xFFFFFFFFu
#define GPU_UINT32_MASK_ONE_LOW  	0x00000001u
#define GPU_UINT32_MASK_ONE_HIGH  	0x80000000u
#define GPU_UINT8_ONES  			0xFF
#define GPU_UINT64_ONES  			0xFFFFFFFFFFFFFFFFu
#define	GPU_UINT8_SIZE				1
#define	GPU_UINT32_SIZE				4
#define	GPU_UINT64_SIZE				8

#define	GPU_INLINE			inline

/* Functions inline */
#define GPU_DIV_CEIL(NUMERATOR,DENOMINATOR) (((NUMERATOR)+((DENOMINATOR)-1))/(DENOMINATOR))
#define GPU_ROUND(NUM) ((int)((NUM) < 0 ? ((NUM) - 0.5) : ((NUM) + 0.5)))
#define GPU_MIN(NUM_A,NUM_B) (((NUM_A) < (NUM_B)) ? (NUM_A) : (NUM_B))
#define GPU_MAX(NUM_A,NUM_B) (((NUM_A) > (NUM_B)) ? (NUM_A) : (NUM_B))
#define GPU_ALIGN_TO(ADDR,BYTES) ((void*)(GPU_DIV_CEIL((uint64_t)ADDR,BYTES) * BYTES))


/* Conversion utils */
#define GPU_CONVERT_B_TO_KB(number) ((number) / (1024))
#define GPU_CONVERT_B_TO_MB(number) ((number) / (1024 * 1024))
#define GPU_CONVERT_B_TO_GB(number) ((number) / (1024 * 1024 * 1024))
#define GPU_CONVERT_MB_TO_B(number) ((number) * (1024 * 1024))

/* System */
#define GPU_FILE_SIZE_LINES				250

/* Encoded DNA Nucleotides */
#define GPU_ENC_DNA_CHAR_A 0LL
#define GPU_ENC_DNA_CHAR_C 1LL
#define GPU_ENC_DNA_CHAR_G 2LL
#define GPU_ENC_DNA_CHAR_T 3LL

#define GPU_ENC_DNA_CHAR_N    4LL
#define GPU_ENC_DNA_CHAR_X    4LL
#define GPU_ENC_DNA_CHAR_SEP  5LL
#define GPU_ENC_DNA_CHAR_JUMP 6LL

#include "../gpu_interface.h"
#include "gpu_errors.h"
#include "gpu_sample.h"

#endif /* GPU_COMMONS_H_ */
