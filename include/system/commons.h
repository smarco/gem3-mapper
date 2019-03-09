/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Base module containing general purpose functions
 */

#ifndef COMMONS_H_
#define COMMONS_H_

/*
 * SETUP
 */
#define _GNU_SOURCE
//#define __USE_GNU 1

/*
 * GENERAL HEADERS
 */
#include <stdio.h>
#include <stdlib.h>

#include <stdbool.h>
#include <stdint.h>
#include <float.h>
#include <inttypes.h>
#include <ctype.h>
#include <sys/types.h>
#include <time.h>
#include <sys/time.h>

#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <limits.h>
#include <pwd.h>

#include <string.h>
#include <math.h>
#include <stdarg.h>

#include <err.h>
#include <errno.h>
#include <assert.h>
#include <signal.h>

//#ifndef __APPLE__
//#include <endian.h>
//#endif

/*
 * Debug
 */
#define GEM_DEEP_DEBUG false

/*
 * Common constants
 */
#define STREAM_FILE_NAME "stream"
#define STREAM_FILE_DESCRIPTOR "<<STREAM>>"
#define ALL UINT64_MAX
#define NO_STRATA ((int64_t)(-1))
// Special Characters
#define EOS '\0'
#define EOL '\n'
#define TAB '\t'
#define DOS_EOL '\r'
#define PLUS '+'
#define MINUS '-'
#define FORMAT '%'
#define SPACE ' '
#define SLASH '/'
#define STAR '*'
#define DOT '.'
#define EQUAL '='
#define COMA ','
#define SEMICOLON ';'
#define COLON ':'
#define HASH '#'
#define UNDERSCORE '_'
#define BANG '!'
#define TILDE '~'

// Buffer sizes
#define BUFFER_SIZE_1K   (1ul<<10)
#define BUFFER_SIZE_2K   (1ul<<11)
#define BUFFER_SIZE_4K   (1ul<<12)
#define BUFFER_SIZE_8K   (1ul<<13)
#define BUFFER_SIZE_16K  (1ul<<14)
#define BUFFER_SIZE_32K  (1ul<<15)
#define BUFFER_SIZE_64K  (1ul<<16)
#define BUFFER_SIZE_128K (1ul<<17)
#define BUFFER_SIZE_256K (1ul<<18)
#define BUFFER_SIZE_512K (1ul<<19)
#define BUFFER_SIZE_1M   (1ul<<20)
#define BUFFER_SIZE_2M   (1ul<<21)
#define BUFFER_SIZE_4M   (1ul<<22)
#define BUFFER_SIZE_8M   (1ul<<23)
#define BUFFER_SIZE_16M  (1ul<<24)
#define BUFFER_SIZE_32M  (1ul<<25)
#define BUFFER_SIZE_64M  (1ul<<26)
#define BUFFER_SIZE_128M (1ul<<27)
#define BUFFER_SIZE_256M (1ul<<28)
#define BUFFER_SIZE_512M (1ul<<29)
#define BUFFER_SIZE_1G   (1ul<<30)
#define BUFFER_SIZE_2G   (1ul<<31)
#define BUFFER_SIZE_4G   (1ul<<32)
#define BUFFER_SIZE_8G   (1ul<<33)
#define BUFFER_SIZE_16G  (1ul<<34)
#define BUFFER_SIZE_32G  (1ul<<35)
#define BUFFER_SIZE_64G  (1ul<<36)
#define BUFFER_SIZE_128G (1ul<<37)
#define BUFFER_SIZE_256G (1ul<<38)

// Metric Factors
#define METRIC_FACTOR_1K   (1000ul)
#define METRIC_FACTOR_1M   (1000000ul)
#define METRIC_FACTOR_1G   (1000000000ul)

// Number of lines
#define NUM_LINES_1K      (1000ul)
#define NUM_LINES_2K      (2000ul)
#define NUM_LINES_5K      (5000ul)
#define NUM_LINES_10K    (10000ul)
#define NUM_LINES_20K    (20000ul)
#define NUM_LINES_50K    (50000ul)
#define NUM_LINES_100K  (100000ul)
#define NUM_LINES_200K  (200000ul)
#define NUM_LINES_500K  (500000ul)
#define NUM_LINES_1M   (1000000ul)
#define NUM_LINES_2M   (2000000ul)
#define NUM_LINES_5M   (5000000ul)
#define NUM_LINES_10M (10000000ul)
#define NUM_LINES_20M (20000000ul)
#define NUM_LINES_50M (50000000ul)
// BM sizes
#define UINT512_LENGTH 512
#define UINT512_SIZE 64
#define UINT256_LENGTH 256
#define UINT256_SIZE 32
#define UINT128_LENGTH 128
#define UINT128_SIZE 16
#define UINT64_LENGTH 64
#define UINT64_SIZE 8
#define UINT32_LENGTH 32
#define UINT32_SIZE 4
#define UINT16_LENGTH 16
#define UINT16_SIZE 2
#define UINT8_LENGTH 8
#define UINT8_SIZE 1
// Types Length
#define INT_MAX_LENGTH 20

/*
 * Common Masks
 */
#define UINT64_ZEROS 0x0000000000000000ull
#define UINT64_ONES  0xFFFFFFFFFFFFFFFFull
#define UINT32_ZEROS 0x00000000ul
#define UINT32_ONES  0xFFFFFFFFul
// Extraction masks
#define UINT64_ONE_MASK        0x0000000000000001ull
#define UINT64_ZERO_MASK       0xFFFFFFFFFFFFFFFEull
#define UINT64_ONE_LAST_MASK   0x8000000000000000ull
#define UINT64_ZERO_LAST_MASK  0x7FFFFFFFFFFFFFFFull
#define UINT32_ONE_MASK        0x00000001ul
#define UINT32_ZERO_MASK       0xFFFFFFFEul
#define UINT32_ONE_LAST_MASK   0x80000000ul
#define UINT32_ZERO_LAST_MASK  0x7FFFFFFFul
// Counting Masks
extern const uint64_t uint64_mask_ones[];
extern const uint64_t uint64_mask_inverse_ones[];
#define uint64_erank_mask(position) uint64_mask_ones[(position)]
#define uint64_erank_inv_mask(position) uint64_mask_inverse_ones[(position)]

/*
 * Conversion utils
 */
#define CONVERT_B_TO_KB(number) ((number)/(1024))
#define CONVERT_B_TO_MB(number) ((number)/(1024*1024))
#define CONVERT_B_TO_GB(number) ((number)/(1024*1024*1024))

/*
 * Compilation/Code Miscellaneous
 */
// GemTools Inline
#ifdef __clang__
#define GEM_INLINE_C
#define GEM_INLINE_H
#else
//#define GEM_INLINE inline
#define GEM_INLINE_C inline
#define GEM_INLINE_H extern inline
#endif

// Conditional expect
#define gem_expect_true(condition)  __builtin_expect(condition,1)
#define gem_expect_false(condition) __builtin_expect(condition,0)
// Macro Stringify
#define QUOTE(value) #value

/*
 * Helper functions (OPERATIVE)
 */
#define SWAP(a,b) do {__typeof__(a) aux = a; a = b; b = aux;} while (0)
#define CFREE(handler) if (handler!=NULL) mm_free(handler);
#define DELEGATE_ERROR(funtion_call,error_code) if ((error_code=funtion_call)) { return error_code; }

/*
 * Random number generator
 */
#define gem_srand() srand(time(0))
// Pseudo-random numbers in [min, max]
#define gem_rand(min,max)   ( min + ( rand()%(max-min+1) ) )
#define gem_rand_f(min,max) ( min + ((double)rand()/(double)(RAND_MAX+1)) * (max-min+1) )
uint64_t gem_rand_IID(const uint64_t min,const uint64_t max);
// Pseudo-random numbers in [0, 1)
#define gem_drand() (double)rand()/(double)(RAND_MAX+1)

/*
 * Common numerical data processing/formating
 */
#define MIN(a,b) (((a)<=(b))?(a):(b))
#define MAX(a,b) (((a)>=(b))?(a):(b))
#define ABS(a) (((a)>=0)?(a):-(a))
#define IS_NUMBER(character) ('0' <= (character) && (character) <= '9')
#define IS_DIGIT(character) IS_NUMBER(character)
#define IS_LETTER(character) (('a' <= (character) && (character) <= 'z') || ('A' <= (character) && (character) <= 'Z'))
#define IS_ALPHANUMERIC(character) (IS_NUMBER(character) || IS_LETTER(character))
#define IS_BETWEEN(number,a,b) ((a)<=(number) && (number)<=(b))

#define BOUNDED_SUBTRACTION(minuend,subtrahend,limit) (((minuend)>((limit)+(subtrahend))) ? (minuend)-(subtrahend):(limit))
#define BOUNDED_ADDITION(summand_A,summand_B,limit) ((((summand_A)+(summand_B))<(limit)) ? (summand_A)+(summand_B):(limit))

#define PERCENTAGE(AMOUNT,TOTAL) ((TOTAL)?100.0*(float)(AMOUNT)/(float)(TOTAL):0.0)
#define DIV(NUMERATOR,DENOMINATOR) ((DENOMINATOR)?(NUMERATOR)/(DENOMINATOR):(0))
#define DIV_CEIL(NUMERATOR,DENOMINATOR) (((NUMERATOR)+((DENOMINATOR)-1))/(DENOMINATOR))
#define DIV_F(NUMERATOR,DENOMINATOR) ((DENOMINATOR)?(float)(NUMERATOR)/(float)(DENOMINATOR):(0))

#define POW4(number) (1<<((number)<<1))

#define DIV_POW2(NUMERATOR,DENOMINATOR_LOG2) ((NUMERATOR)>>(DENOMINATOR_LOG2))
extern const uint64_t uint64_mask_mod_pow2[];
#define MOD_POW2(NUMERATOR,DENOMINATOR_LOG2) ((NUMERATOR) & uint64_mask_mod_pow2[DENOMINATOR_LOG2])

uint64_t integer_num_ciphers(const uint64_t number);

uint64_t integer_proportion(const double proportion,const uint64_t length);
uint64_t integer_lower_power_of_two(uint64_t number);
uint64_t integer_upper_power_of_two(uint64_t number);

int integer_to_ascii(char* const buffer,uint64_t number);

float gem_log2(float number);
float gem_loge(float number);

/*
 * Statistical Utils
 */
double standard_normal_CDF(double x);

/*
 * CheckSum & BitDisplay
 */
uint64_t checksum_uint64(uint64_t* mem,const uint64_t num_words);
void checksum_incremental_uint64(uint64_t* const checksum,const uint64_t word);

#define CHAR_TO_PRINTABLE(character) ( ('!' <= (character) && (character) <= '~') ? (character) : '#')
void fprintf_uint64_binary(FILE* const stream,const uint64_t word);
void fprintf_uint64_footprint(FILE* const stream,const uint64_t word);

/*
 * Common parsing
 */
#define IS_EOL(character) gem_expect_false((character)==EOL)
#define IS_ANY_EOL(character) ((character)==EOL || (character)==DOS_EOL)
#define IS_HEX_DIGIT(character) (IS_NUMBER(character) || ('a' <= (character) && (character) <= 'f') || ('A' <= (character) && (character) <= 'F'))

#define GET_DIGIT(character) ((character) - '0')
#define GET_HEX_DIGIT(character) (IS_NUMBER(character) ? GET_DIGIT(character) : (toupper(character) - 'A' + 10))

#define IS_END_OF_RECORD(character) ( (character)==EOL || (character)==EOS )
#define IS_END_OF_FIELD(character) ( IS_END_OF_RECORD(character) || (character)==SPACE || (character)==TAB )

/*
 * Time (ms)
 */
#define TIME_DIFF_NS(start,end) ((end.tv_sec*1000000000 + end.tv_nsec) - (start.tv_sec*1000000000 + start.tv_nsec))
#define TIME_DIFF_S(start,end) ((end.tv_sec + end.tv_nsec/1E9) - (start.tv_sec + start.tv_nsec/1E9))

/*
 * Popcount macros
 */
#ifdef __SSE4_2__
  #include <nmmintrin.h>
  #define POPCOUNT_64(word64) _mm_popcnt_u64((word64))
  #define POPCOUNT_32(word32) _mm_popcnt_u32((word32))
#else
  #define POPCOUNT_64(word64) __builtin_popcountll((word64))
  #define POPCOUNT_32(word32) __builtin_popcount((word32))
#endif

/*
 * Prefetch macros
 */
#ifdef __SSE__
  #include <xmmintrin.h>
  #define PREFETCH(ADDR) _mm_prefetch(((const char*)ADDR),_MM_HINT_NTA)
#else
  #define PREFETCH(ADDR) __builtin_prefetch(((const char*)ADDR),0,0)
#endif

/*
 * System
 */
typedef enum {
  DEVICE_CPU,
  DEVICE_GPU,
} device_t;

uint64_t system_get_num_processors(void);
char* system_get_cwd(void);
char* system_get_hostname(void);
char* system_get_user_name(void);
void system_print_info(FILE* const stream);

void system_get_time(struct timespec *ts);

#endif /* COMMONS_H_ */
