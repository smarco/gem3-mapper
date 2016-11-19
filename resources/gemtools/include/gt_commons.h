/*
 * PROJECT: GEM-Tools library
 * FILE: gt_commons.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Base module containing general purpose functions
 */

#ifndef GT_COMMONS_H_
#define GT_COMMONS_H_

#define inline

/*
 * VERSION
 */
#define GT_VERSION "1.6"
#define GT_GIT_URL "https://github.com/gemtools/gemtools"

/*
 * SETUP
 */
#define _GNU_SOURCE
#define __USE_GNU

/*
 * GENERAL HEADERS
 */
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <float.h>
#include <inttypes.h>
#include <stdlib.h>

#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>

#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <getopt.h>

#include <ctype.h>
#include <sys/types.h>
#include <time.h>
#include <sys/time.h>

#include <errno.h>
#include <err.h>
#include <assert.h>
#include <signal.h>

#include <pthread.h>

/*
 * COMMONS
 */

// Data constants
#define UINT64_ZEROS 0x0000000000000000ul
#define UINT64_ONES  0xFFFFFFFFFFFFFFFFul
#define UINT64_ONE   0x0000000000000001ul

// Special characters
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

// Buffer sizes
#define GT_BUFFER_SIZE_1K   ((1<<10))
#define GT_BUFFER_SIZE_2K   ((1<<11))
#define GT_BUFFER_SIZE_4K   ((1<<12))
#define GT_BUFFER_SIZE_8K   ((1<<13))
#define GT_BUFFER_SIZE_16K  ((1<<14))
#define GT_BUFFER_SIZE_32K  ((1<<15))
#define GT_BUFFER_SIZE_64K  ((1<<16))
#define GT_BUFFER_SIZE_128K ((1<<17))
#define GT_BUFFER_SIZE_256K ((1<<18))
#define GT_BUFFER_SIZE_512K ((1<<19))
#define GT_BUFFER_SIZE_1M   ((1<<20))
#define GT_BUFFER_SIZE_2M   ((1<<21))
#define GT_BUFFER_SIZE_4M   ((1<<22))
#define GT_BUFFER_SIZE_8M   ((1<<23))
#define GT_BUFFER_SIZE_16M  ((1<<24))
#define GT_BUFFER_SIZE_32M  ((1<<25))
#define GT_BUFFER_SIZE_64M  ((1<<26))
#define GT_BUFFER_SIZE_128M ((1<<27))
#define GT_BUFFER_SIZE_256M ((1<<28))
#define GT_BUFFER_SIZE_512M ((1<<29))
#define GT_BUFFER_SIZE_1G   ((1<<30))
#define GT_BUFFER_SIZE_2G   ((1<<31))
#define GT_BUFFER_SIZE_4G   ((1<<32))
#define GT_BUFFER_SIZE_8G   ((1<<33))

// Number of lines
#define GT_NUM_LINES_1K      (1000)
#define GT_NUM_LINES_2K      (2000)
#define GT_NUM_LINES_5K      (5000)
#define GT_NUM_LINES_10K    (10000)
#define GT_NUM_LINES_20K    (20000)
#define GT_NUM_LINES_50K    (50000)
#define GT_NUM_LINES_100K  (100000)
#define GT_NUM_LINES_200K  (200000)
#define GT_NUM_LINES_500K  (500000)
#define GT_NUM_LINES_1M   (1000000)
#define GT_NUM_LINES_2M   (2000000)
#define GT_NUM_LINES_5M   (5000000)
#define GT_NUM_LINES_10M (10000000)
#define GT_NUM_LINES_20M (20000000)
#define GT_NUM_LINES_50M (50000000)

// Conditional expect
#define gt_expect_true(condition) __builtin_expect(condition,1)
#define gt_expect_false(condition) __builtin_expect(condition,0)

// GemTools Inline
#define GT_INLINE inline

// Macro Stringify
#define GT_QUOTE(value) #value

// Size of an static array
#define GT_CONST_ARRAY_SIZE(array,type) (sizeof(array)/sizeof(type))

/*
 * Is functions
 */
#define gt_is_number(character) ('0' <= (character) && (character) <= '9')
#define gt_is_hex_digit(character) (gt_is_number(character) || ('a' <= (character) && (character) <= 'f') || ('A' <= (character) && (character) <= 'F'))
#define gt_is_letter(character) (('a' <= (character) && (character) <= 'z') || ('A' <= (character) && (character) <= 'Z'))
#define gt_is_alphanumeric(character) (gt_is_number(character) || gt_is_letter(character))

#define gt_is_end_of_record(character) ( (character)==EOL || (character)==EOS )
#define gt_is_end_of_field(character) ( gt_is_end_of_record(character) || (character)==SPACE || (character)==TAB )

#define gt_get_cipher(character) ((character) - '0')
#define gt_get_hex_cipher(character) (gt_is_number(character)?gt_get_cipher(character):(toupper(character) - 'A' + 10))

/*
 * Helper functions (OPERATIVE)
 */
#define gt_cfree(handler) if (handler!=NULL) gt_free(handler);
#define GT_MIN(a,b) ((a)<=(b)?(a):(b))
#define GT_MAX(a,b) ((a)>=(b)?(a):(b))
#define GT_ABS(a) ((a)>=0?(a):-(a))
#define GT_SWAP(a,b) do {typeof(a) aux = a; a = b; b = aux; } while (0)
#define GT_BETWEEN(number,a,b) ((a)<=(number) && (number)<=(b))
GT_INLINE uint64_t gt_get_integer_proportion(const float proportion,const uint64_t length);

/*
 * Print's template helpers
 */
GT_INLINE uint64_t gt_calculate_memory_required_v(const char *template,va_list v_args);
GT_INLINE uint64_t gt_calculate_memory_required_va(const char *template,...);

/*
 * Error value return wrapper
 */
#define GT_DELEGATE_ERROR(funtion_call,error_code) if ((error_code=funtion_call)) { return error_code; }

/*
 * Mutex/Cond Helpers
 */

#define GT_MUTEX_INIT(mutex) \
  gt_cond_fatal_error__perror(pthread_mutex_init(&(mutex),NULL),SYS_MUTEX_INIT);
#define GT_BEGIN_MUTEX_SECTION(mutex) \
  gt_cond_fatal_error(pthread_mutex_lock(&(mutex)),SYS_MUTEX);
#define GT_END_MUTEX_SECTION(mutex) \
  gt_cond_fatal_error(pthread_mutex_unlock(&(mutex)),SYS_MUTEX);
#define GT_CV_SIGNAL(cv) \
  gt_cond_fatal_error(pthread_cond_signal(&(cv)),SYS_COND_VAR);
#define GT_CV_BROADCAST(cv) \
  gt_cond_fatal_error(pthread_cond_broadcast(&(cv)),SYS_COND_VAR);
#define GT_CV_WAIT(cv,mutex) \
  gt_cond_fatal_error(pthread_cond_wait(&(cv),&(mutex)),SYS_COND_VAR);

/*
 * Random number generator
 */
#define gt_srand() srand(time(0))
// pseudo-random numbers in [min, max]
#define gt_rand(min,max)   ( min + ( rand()%(max-min+1) ) )
#define gt_rand_f(min,max) ( min + ((double)rand()/(double)(RAND_MAX+1)) * (max-min+1) )
GT_INLINE uint64_t gt_rand_IID(const uint64_t min,const uint64_t max);
// pseudo-random numbers in [0, 1)
#define gt_drand() (double)rand()/(double)(RAND_MAX+1)

/*
 * Common constants
 */
#define GT_STREAM_FILE_NAME "<<STREAM>>"
#define GT_ALL UINT64_MAX
#define GT_NO_STRATA ((int64_t)(-1))

/*
 * Common data processing/formating
 */
#define GT_GET_PERCENTAGE(AMOUNT,TOTAL) ((TOTAL)?100.0*(float)(AMOUNT)/(float)(TOTAL):0.0)
#define GT_DIV(NUMERATOR,DENOMINATOR) ((DENOMINATOR)?(NUMERATOR)/(DENOMINATOR):(0))
#define GT_DIV_F(NUMERATOR,DENOMINATOR) ((DENOMINATOR)?(float)(NUMERATOR)/(float)(DENOMINATOR):(0))

/*
 * Target dependent functions
 */
// Popcount macros
#ifdef __SSE4_2__
  #include <nmmintrin.h>
  #define GT_POPCOUNT_64 _mm_popcnt_u64
  #define GT_POPCOUNT_32 _mm_popcnt_u32
#else
  #define GT_POPCOUNT_64 __builtin_popcountll
  #define GT_POPCOUNT_32 __builtin_popcount
#endif
// Prefetch macros
#ifdef __SSE__
  #include <xmmintrin.h>
  #define GT_PREFETCH(ADDR) _mm_prefetch(ADDR,_MM_HINT_NTA)
#else
  #define GT_PREFETCH(ADDR) __builtin_prefetch(ADDR,0,0)
#endif

#endif /* GT_COMMONS_H_ */
