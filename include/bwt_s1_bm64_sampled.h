/*
 * PROJECT: GEMMapper
 * FILE: bwt_s1_bm64_sampled.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Provides basic routines to encode a DNA text into a Burrows-Wheeler transform
 *              using a compact bit representation and counter buckets as to enhance Occ/rank queries
 */

#ifndef BWT_S1_BM64_SAMPLED_H_
#define BWT_S1_BM64_SAMPLED_H_

/*
 * Checkers
 */
#ifndef SAMPLING_SA_DIRECT
#error "BWT implementation only compatible with SAMPLING_SA_DIRECT ('sampling_sa.h')"
#endif

#define BWT_CHECK(bwt) GEM_CHECK_NULL(bwt)
#define BWT_BUILDER_CHECK(bwt) BWT_CHECK(bwt)

#define BWT_CHECK_INDEX(i) gem_cond_fatal_error(i > bwt->n,BWT_INDEX,i)
#define BWT_CHECK_ENC(c) gem_cond_fatal_error(!is_enc_dna(c),BWT_ENC_CHAR,C)
#define BWT_CHECK_INDEX__ENC(i,c) \
  DNA_PBWT_CHECK_INDEX(i); \
  DNA_PBWT_CHECK_ENCODED_CHAR(c)

/*
 * BWT Structure
 */
typedef struct {
  /* Meta-Data */
  uint64_t length;                    // Length of the BWT
  uint64_t* c;                        // Occurrences of each character
  uint64_t* C;                        // The cumulative occurrences ("ranks") of the symbols of the string
  /* BWT Buckets-Structure */
  uint64_t num_minor_blocks;          // Total number of minor blocks
  uint64_t num_mayor_blocks;          // Total number of mayor blocks
  uint64_t* mayor_counters;           // Pointer to the Mayor Counters (Rank)
  uint64_t* bwt_mem;                  // Pointer to the BWT structure in memory
  /* MM */
  mm_t* mm_counters;
  mm_t* mm_bwt_mem;
} bwt_t;
typedef struct {
  /* Meta-Data */
  bwt_t bwt;
  /* Auxiliary variables */
  uint64_t mayor_counter;             // Current mayorCounter position
  uint64_t* minor_block_mem;          // Pointer to current minorBlock position
  uint64_t* mayor_occ;                // Mayor accumulated counts
  uint64_t* minor_occ;                // Minor accumulated counts
  uint64_t sampling_mayor_count;      // Sampling-Mayor accumulated count
  uint64_t sampling_minor_count;      // Sampling-Minor accumulated count
} bwt_builder_t;
// BWT Elements
typedef struct {
  uint64_t block_pos;        // minorBlock position (@i / 64)
  uint64_t block_mod;        // Current position modulus letters per minorBlock (@i % 64)
  uint64_t* mayor_counters;  // Pointer to the proper MayorCounters (@i.MayorCounters)
  uint64_t* block_mem;       // Pointer to the proper MinorBlock (@i.MinorBlock)
} bwt_block_locator_t;
typedef struct {
  uint64_t bitmap_1__2[4]; /* 0 - x00 - {A,N}   => ~W[1] & ~W[2]
                            * 1 - x01 - {C,SEP} => ~W[1] &  W[2]
                            * 2 - x10 - {G,J}  =>  W[1] & ~W[2]
                            * 3 - x11 - {T,#8}  =>  W[1] &  W[2] */
  uint64_t bitmap_3[2];    /* 0 - 0xx - {A,C,G,T}     => ~W[3]
                            * 1 - 1xx - {N,SEP,J,#8} =>  W[3] */
  uint64_t gap_mask;       /* Masking bits until @lo position [-XXXXXXX-lo-xxxxxxx-hi-......] */
} bwt_block_elms_t;

/*
 * BWT Dimensions
 *   Alphabet = 3bits => 3 dense
 *   MinorBlock => 16B (minorCounters) + 24B (BM) + 8B (samplingBM) => 48B
 *   MayorBlock => 48KB (1024 MinorBlocks)
 *   |MinorCounters| = 8*UINT16
 *   |MayorCounters| = 8*UINT64  (Stored separately)
 *   MinorBlock.BM = 3*UINT64
 *     - Bitmap-LayOut Bitwise (3 x 64word)
 *     - [64bits-Bit0-{LSB}][64bits-Bit1][64bits-Bit2-{MSB}]
 *   SamplingBM = 1*UINT64
 *     - Bitmap of sampled positions (UINT64 within each MinorBlock.BM[3])
 *     - MayorCounters => UINT64 (Stored at MayorCounters[7]) [Piggybacking]
 *     - MinorCounters => UINT16 (Stored at MinorBlock.MinorCounters[7]) [Piggybacking]
 * Remarks::
 *   Number of MayorBlocks (3GB genome) => 3*(1024^3)/65536 = 49152 MBlocks
 *   Assuming SysPage = 4KB =>
 *       1 SysPage => 85 minorBlocks (1 unaligned, will cause 2 Pagefaults)
 *       3 SysPage => 256 minorBlocks (2 unaligned) => 1,011 PageFaults per access [OK]
 *       48KB/4KB = 12 SysPage per MayorBlock
 */
#define BWT_MINOR_COUNTER_RANGE         8
#define BWT_MINOR_BLOCK_COUNTER_LENGTH 16
#define BWT_MINOR_BLOCK_LENGTH         64
#define BWT_MINOR_BLOCK_BITS            3

#define BWT_MINOR_BLOCK_SIZE \
  (BWT_MINOR_BLOCK_COUNTER_LENGTH*BWT_MINOR_COUNTER_RANGE + \
   BWT_MINOR_BLOCK_LENGTH*BWT_MINOR_BLOCK_BITS + \
   BWT_MINOR_BLOCK_LENGTH)/8                              /* 48B */
#define BWT_MINOR_BLOCK_64WORDS (BWT_MINOR_BLOCK_SIZE/8)

// #define BWT_MAYOR_BLOCK_NUM_MINOR_BLOCKS ((((1<<(16))-1)/BWT_MINOR_BLOCK_LENGTH) + 1) /* (((2^16)−1)÷64)+1 = 1024 */
#define BWT_MAYOR_BLOCK_NUM_MINOR_BLOCKS (1<<10) /* 1024 */
#define BWT_MAYOR_BLOCK_LENGTH (BWT_MAYOR_BLOCK_NUM_MINOR_BLOCKS*BWT_MINOR_BLOCK_LENGTH) /* 1024*64 = 65536 */
#define BWT_MAYOR_BLOCK_SIZE   (BWT_MAYOR_BLOCK_NUM_MINOR_BLOCKS*BWT_MINOR_BLOCK_SIZE)   /* 1024*48B = 48KB */

#define BWT_MAYOR_COUNTER_SIZE  UINT64_SIZE
#define BWT_MAYOR_COUNTER_RANGE 8

#define SAMPLING_MAYOR_COUNTER_SIZE  UINT64_SIZE
#define SAMPLING_MINOR_COUNTERS_SIZE UINT16_SIZE
#define SAMPLING_ENC_CHAR 7

#endif /* BWT_S1_BM64_SAMPLED_H_ */
