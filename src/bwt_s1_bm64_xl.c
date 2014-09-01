/*
 * PROJECT: GEMMapper
 * FILE: bwt.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Provides basic routines to encode a DNA text into a Burrows-Wheeler transform
 *              using a compact bit representation and counter buckets as to enhance Occ/rank queries
 */

/*
 * TODO
 *   - Marked bitmap shortcut into the BLOCK ¿how?
 *   - (char_enc-DNA_RANGE) to (char_enc & 011)
 */

#include "bwt.h"
#include "sparse_array_locator.h"

/*
 * Checkers
 */
#define BWT_CHECK(bwt) GEM_CHECK_NULL(bwt)
#define BWT_BUILDER_CHECK(bwt) BWT_CHECK(bwt)

#define BWT_CHECK_INDEX(i) gem_cond_fatal_error(i > bwt->n,BWT_INDEX,i)
#define BWT_CHECK_ENC(c) gem_cond_fatal_error(!is_enc_dna(c),BWT_ENC_CHAR,C)
#define BWT_CHECK_INDEX__ENC(i,c) \
  DNA_PBWT_CHECK_INDEX(i); \
  DNA_PBWT_CHECK_ENCODED_CHAR(c)

/*
 * Profiler/Stats
 */
#define BWT_CHAR_TICK()
#define BWT_ERANK_TICK()
#define BWT_ERANK_INTERVAL_TICK()
#define BWT_LF_TICK()

#define BWT_PREFETCH_TICK()
#define BWT_PRECOMPUTE_TICK()
#define BWT_PREFETCH__PRECOMPUTE_TICK()


// BWT Structure
typedef struct {
  uint16_t counter[4];                // Extra counters {N,SEP,...}
  uint64_t layer_2;                   // Third bitmap layer (MSB)
} xl_block_t;
struct _bwt_t {
  /* Meta-Data */
  uint64_t length;                    // Length of the BWT
  uint64_t* c;                        // Occurrences of each character
  uint64_t* C;                        // The cumulative occurrences ("ranks") of the symbols of the string
  /* BWT Buckets-Structure */
  uint64_t num_minor_blocks;          // Total number of minor blocks
  uint64_t num_mayor_blocks;          // Total number of mayor blocks
  uint64_t* mayor_counters;           // Pointer to the Mayor Counters
  uint64_t* bwt_mem;                  // Pointer to the BWT structure in memory
  /* (MSB) Extra-Layer {'|','N',...} */
  sparse_array_locator_t* xl_locator; // 3-Bit layer locator
  uint64_t num_xl_bitmaps;            // Total 3-Bit layer Bitmaps
  xl_block_t* xl_mem;                 // 3-Bit layer structure [16Counter]x4 + [64BitMap]
  /* MM */
  mm_t* mm_mayor_counters;
  mm_t* mm_bwt_mem;
  mm_t* mm_xl_mem;
};
struct _bwt_builder_t {
  /* Meta-Data */
  struct _bwt_t bwt;
  /* (MSB) Extra-Layer {'|','N',...} */
  sparse_array_locator_t* xl_locator; // 3-Bit layer locator
  svector_t* xl_memory;               // (xl_block_t) 3-Bit layer structure [16Counter]x2 + [64BitMap]?
  svector_iterator_t xl_iterator;     // @xl_memory Iterator
  /* Auxiliary variables */
  uint64_t mayor_counter;             // Current mayorCounter position
  uint64_t minor_block;               // Current minorBlock position
  uint64_t* minor_block_mem;          // Pointer to current minorBlock position
  uint64_t* mayor_occ;                // Mayor accumulated counts
  uint64_t* minor_occ;                // Minor accumulated counts
  uint64_t* minor_occ_ext;            // Minor accumulated counts (for extended characters)
};
// BWT Elements
struct _bwt_block_locator_t {
  uint64_t block_pos;        // minorBlock position (@i / 64)
  uint64_t block_mod;        // Current position modulus letters per minorBlock (@i % 64)
  uint64_t* mayor_counters;  // Pointer to the proper MayorCounters (@i.MayorCounters)
  uint64_t* block_mem;       // Pointer to the proper MinorBlock (@i.MinorBlock)
  bool has_xl_block;         // Has the proper MinorBlock a xl_block?
  xl_block_t* xl_block;      // Pointer to the proper XL-Block (@i.XL)
};
struct _bwt_block_elms_t {
  uint64_t bitmap_1__2[4]; /* 0 - x00 - {A,N}   => ~W[1] & ~W[2]
                            * 1 - x01 - {C,SEP} => ~W[1] &  W[2]
                            * 2 - x10 - {G,#7}  =>  W[1] & ~W[2]
                            * 3 - x11 - {T,#8}  =>  W[1] &  W[2] */
  uint64_t bitmap_3[2];    /* 0 - 0xx - {A,C,G,T}     => ~W[3]
                            * 1 - 1xx - {N,SEP,#7,#8} =>  W[3] */
  uint64_t gap_mask;       /* Masking bits until @lo position [-XXXXXXX-lo-xxxxxxx-hi-......] */
};

/*
 * BWT Dimensions
 *   Alphabet = 3bits => (2 dense + 1 space) [Annotated in SparseBitmap]
 *   MinorBlock => 24B
 *   MayorBlock => 24KB (1024 MinorBlocks)
 *   |MinorCounters| = 4*UINT16  (Sparse 4*UINT16)
 *   |MayorCounters| = 8*UINT64  (Stored separately)
 *   MinorBlock.BM = 2*UINT64 (Sparse-Bit2 1*UINT64)
 *     - Bitmap-LayOut Bitwise (3 x 64word)
 *     - [64bits-Bit0-{LSB}][64bits-Bit1][64bits-Bit2-{MSB}]
 * Remarks::
 *   Number of MayorBlocks (3GB genome) => 3*1024^3/65536 = 49152 MBlocks => 6144 Bytes of SparseBitmap-Annotation
 *   Assuming SysPage = 4KB =>
 *       1 SysPage => 170 minorBlocks (1 unaligned, will cause 2 Pagefaults)
 *       3 SysPage => 3*170 minorBlocks (2 unaligned) => (3*170+2/170) => 1,0039 PageFaults per access [OK]
 *       24KB/4KB = 6 SysPage per MayorBlock
 */
#define BWT_MINOR_COUNTER_RANGE        DNA_RANGE /* 4 */
#define BWT_MINOR_BLOCK_COUNTER_LENGTH 16
#define BWT_MINOR_BLOCK_LENGTH         64
#define BWT_MINOR_BLOCK_BITS           DNA_RANGE_BITS /* 2 */

#define BWT_MINOR_BLOCK_SIZE \
  (BWT_MINOR_BLOCK_COUNTER_LENGTH*BWT_MINOR_COUNTER_RANGE + \
   BWT_MINOR_BLOCK_LENGTH*BWT_MINOR_BLOCK_BITS)/8  /* 24B */
#define BWT_MINOR_BLOCK_64WORDS (BWT_MINOR_BLOCK_SIZE/8)

// #define BWT_MAYOR_BLOCK_NUM_MINOR_BLOCKS ((((1<<(16))-1)/BWT_MINOR_BLOCK_LENGTH) + 1) /* (((2^16)−1)÷64)+1 = 1024 */
#define BWT_MAYOR_BLOCK_NUM_MINOR_BLOCKS (1<<10) /* 1024 */
#define BWT_MAYOR_BLOCK_LENGTH (BWT_MAYOR_BLOCK_NUM_MINOR_BLOCKS*BWT_MINOR_BLOCK_LENGTH) /* 1024*64 = 65536 */
#define BWT_MAYOR_BLOCK_SIZE   (BWT_MAYOR_BLOCK_NUM_MINOR_BLOCKS*BWT_MINOR_BLOCK_SIZE)   /* 1024*24B = 24576B = 24KB */

#define BWT_MAYOR_COUNTER_SIZE  UINT64_SIZE
#define BWT_MAYOR_COUNTER_RANGE 8

/*
 * BWT Builder Auxiliary functions
 */
GEM_INLINE void bwt_builder_initialize(
    bwt_builder_t* const bwt_builder,
    const uint64_t bwt_text_length,const uint64_t* const character_occurrences) {
  struct _bwt_t* const bwt = &bwt_builder->bwt;
  /*
   * Meta-Data
   */
  bwt->length = bwt_text_length;
  // Set characters count
  uint8_t i;
  bwt->c = mm_calloc(BWT_MAYOR_COUNTER_RANGE,uint64_t,true);
  for (i=0;i<BWT_MAYOR_COUNTER_RANGE;++i) bwt->c[i] = character_occurrences[i];
  // Set the cumulative occurrences (the relative ordering of the symbols is encoded)
  bwt->C = mm_calloc(BWT_MAYOR_COUNTER_RANGE,uint64_t,true);
  bwt->C[ENC_DNA_CHAR_A]   = 0;
  bwt->C[ENC_DNA_CHAR_C]   = bwt->C[ENC_DNA_CHAR_A] + bwt->c[ENC_DNA_CHAR_A];
  bwt->C[ENC_DNA_CHAR_G]   = bwt->C[ENC_DNA_CHAR_C] + bwt->c[ENC_DNA_CHAR_C];
  bwt->C[ENC_DNA_CHAR_T]   = bwt->C[ENC_DNA_CHAR_G] + bwt->c[ENC_DNA_CHAR_G];
  bwt->C[ENC_DNA_CHAR_N]   = bwt->C[ENC_DNA_CHAR_T] + bwt->c[ENC_DNA_CHAR_T];
  bwt->C[ENC_DNA_CHAR_SEP] = bwt->C[ENC_DNA_CHAR_N] + bwt->c[ENC_DNA_CHAR_N];
  bwt->num_minor_blocks = DIV_CEIL(bwt->length,BWT_MINOR_BLOCK_LENGTH);
  bwt->num_mayor_blocks = DIV_CEIL(bwt->num_minor_blocks,BWT_MAYOR_BLOCK_NUM_MINOR_BLOCKS);
  /*
   * Allocate BWT Structures
   */
  bwt->mayor_counters = mm_calloc(bwt->num_mayor_blocks*BWT_MAYOR_COUNTER_RANGE,uint64_t,true);
  bwt->bwt_mem = (uint64_t*) mm_calloc(bwt->num_minor_blocks*BWT_MINOR_BLOCK_SIZE,uint8_t,true);
  // Allocate eXtraLayer
  bwt_builder->xl_locator = sparse_array_locator_new(0,bwt->num_minor_blocks);
  bwt_builder->xl_memory = svector_new(mm_pool_get_slab(mm_pool_32MB),xl_block_t);
  svector_iterator_new(&bwt_builder->xl_iterator,bwt_builder->xl_memory,SVECTOR_WRITE_ITERATOR,0);
  // Auxiliary variables
  bwt_builder->mayor_counter = 0;
  bwt_builder->minor_block = 0;
  bwt_builder->minor_block_mem = bwt->bwt_mem;
  bwt_builder->mayor_occ = mm_calloc(BWT_MAYOR_COUNTER_RANGE,uint64_t,true);
  bwt_builder->minor_occ = mm_calloc(BWT_MAYOR_COUNTER_RANGE,uint64_t,true);
  bwt_builder->minor_occ_ext = mm_calloc(BWT_MAYOR_COUNTER_RANGE,uint64_t,true);
}
#define bwt_builder_write_enc(enc,bit_mask,layer_0,layer_1,layer_2) \
  switch (enc) { \
    case ENC_DNA_CHAR_A:   /* 000 */ \
      /* layer_0 |= bit_mask; */ \
      /* layer_1 |= bit_mask; */ \
      /* layer_2 |= bit_mask; */ \
      break; \
    case ENC_DNA_CHAR_C:   /* 001 */ \
      layer_0 |= bit_mask; \
      /* layer_1 |= bit_mask; */ \
      /* layer_2 |= bit_mask; */ \
      break; \
    case ENC_DNA_CHAR_G:   /* 010 */ \
      /* layer_0 |= bit_mask; */ \
      layer_1 |= bit_mask; \
      /* layer_2 |= bit_mask; */ \
      break; \
    case ENC_DNA_CHAR_T:   /* 011 */ \
      layer_0 |= bit_mask; \
      layer_1 |= bit_mask; \
      /* layer_2 |= bit_mask; */ \
      break; \
    case ENC_DNA_CHAR_N:   /* 100 */ \
      /* layer_0 |= bit_mask; */ \
      /* layer_1 |= bit_mask; */ \
      layer_2 |= bit_mask; \
      break; \
    case ENC_DNA_CHAR_SEP: /* 101 */ \
      layer_0 |= bit_mask; \
      /* layer_1 |= bit_mask; */ \
      layer_2 |= bit_mask; \
      break; \
    case ENC_DNA_CHAR_JUMP: /* 110 */ \
      /* layer_0 |= bit_mask; */ \
      layer_1 |= bit_mask; \
      layer_2 |= bit_mask; \
      break; \
    default: \
      GEM_INVALID_CASE(); \
      break; \
  } \
  bit_mask <<= 1

GEM_INLINE void bwt_builder_write_mayor_counters(bwt_builder_t* const bwt_builder) {
  uint64_t* const mayor_counters =
      bwt_builder->bwt.mayor_counters + bwt_builder->mayor_counter;
  uint8_t i;
  for (i=0;i<BWT_MAYOR_COUNTER_RANGE;++i) {
    // Update sentinel counters
    bwt_builder->mayor_occ[i] += bwt_builder->minor_occ[i];
    bwt_builder->minor_occ[i] = 0;
    // Dump counters
    mayor_counters[i] = bwt_builder->mayor_occ[i] + bwt_builder->bwt.C[i];
  }
  // Update MayorCounter position
  bwt_builder->mayor_counter += BWT_MAYOR_COUNTER_RANGE;
}
GEM_INLINE void bwt_builder_write_minor_counters(bwt_builder_t* const bwt_builder) {
  uint16_t* const minor_counters_mem = (uint16_t*) bwt_builder->minor_block_mem;
  uint8_t i;
  // Dump main-4 counters
  for (i=0;i<BWT_MINOR_COUNTER_RANGE;++i) {
    minor_counters_mem[i] = bwt_builder->minor_occ[i]; // Write minorCounter
    bwt_builder->minor_occ[i] = 0; // Reset
  }
  // Store ext-4 counters
  for (;i<BWT_MAYOR_COUNTER_RANGE;++i) {
    bwt_builder->minor_occ_ext[i] = bwt_builder->minor_occ[i];
    bwt_builder->minor_occ[i] = 0;
  }
  // Update minorBlock pointer
  ++(bwt_builder->minor_block_mem);
}
GEM_INLINE void bwt_builder_write_minor_block(
    bwt_builder_t* const bwt_builder,
    const uint64_t layer_0,const uint64_t layer_1,const uint64_t layer_2) {
  *bwt_builder->minor_block_mem = layer_0; ++(bwt_builder->minor_block_mem);
  *bwt_builder->minor_block_mem = layer_1; ++(bwt_builder->minor_block_mem);
  // Store eXtraLayer (if non empty)
  if (layer_2 != 0 || bwt_builder->minor_block==0) { // First minorBlock always written (XL not empty)
    sparse_array_locator_mark(bwt_builder->xl_locator,bwt_builder->minor_block);
    xl_block_t* const xl_block = svector_iterator_get_elm(&bwt_builder->xl_iterator);
    // Store ext-4 counters
    uint8_t i;
    for (i=BWT_MINOR_COUNTER_RANGE;i<BWT_MAYOR_COUNTER_RANGE;++i) {
      xl_block->counter[i-BWT_MINOR_COUNTER_RANGE] = bwt_builder->minor_occ_ext[i];
    }
    // Store layer_2
    xl_block->layer_2 = layer_2;
  }
  // Next
  ++(bwt_builder->minor_block);
  svector_write_iterator_next(&bwt_builder->xl_iterator);
}
/*
 * BWT Builder
 */
GEM_INLINE bwt_builder_t* bwt_builder_new(
    dna_text_t* const bwt_text,const uint64_t* const character_occurrences,
    const bool check,const bool verbose) {
  // TODO Checks
  /*
   * Allocate & initialize builder
   */
  bwt_builder_t* const bwt_builder = mm_alloc(bwt_builder_t);
  const uint64_t bwt_length = bwt_text->length;
  bwt_builder_initialize(bwt_builder,bwt_length,character_occurrences);
  /*
   * Compute BWT & write
   */
  // Cursors/Locators
  uint64_t mayorBlock_pos = BWT_MAYOR_BLOCK_LENGTH;
  uint64_t minorBlock_pos = BWT_MINOR_BLOCK_LENGTH;
  // BM Layers
  uint64_t layer_0=0, layer_1=0, layer_2=0, bit_mask=UINT64_ONE_MASK;
  // Iterate over the BWT  // TODO Implement Ticker
  const uint8_t* const bwt = bwt_text->text;
  uint64_t bwt_pos = 0;
  while (bwt_pos < bwt_length) {
    // Get BWT character
    const uint8_t enc = bwt[bwt_pos++];
    // Print MayorCounters
    if (mayorBlock_pos==BWT_MAYOR_BLOCK_LENGTH) {
      bwt_builder_write_mayor_counters(bwt_builder);
      mayorBlock_pos = 0;
    }
    // Print MinorCounters
    if (minorBlock_pos==BWT_MINOR_BLOCK_LENGTH) {
      bwt_builder_write_minor_counters(bwt_builder);
      minorBlock_pos = 0;
    }
    // Write character to BitMap
    bwt_builder_write_enc(enc,bit_mask,layer_0,layer_1,layer_2);
    // Account current char & Next
    ++(bwt_builder->minor_occ[enc]);
    ++minorBlock_pos; ++mayorBlock_pos;
    // Close MinorBlock (dump)
    if (minorBlock_pos==BWT_MINOR_BLOCK_LENGTH) {
      bwt_builder_write_minor_block(bwt_builder,layer_0,layer_1,layer_2);
      // Reset
      layer_0=0; layer_1=0; layer_2=0;
      bit_mask = UINT64_ONE_MASK;
    }
  }
  // Close last MinorBlock (dump)
  if (minorBlock_pos!=BWT_MINOR_BLOCK_LENGTH) {
    bwt_builder_write_minor_block(bwt_builder,layer_0,layer_1,layer_2);
  }
//
//
//  // Possible consistency checks
//  if (do_tests) {
//    if (be_verbose) {
//      tprintf(Time, "Beginning BWT test (%ld calls)...\n", (uint64_t) bwt->n * 8);
//    }
//    /* We check that everything be consistent */
//    idx_t dist = 0;
//    memset(occs, 0, DNA_T_RANGE * sizeof(uint64_t));
//    pos = 0;
//    ticker_reset(ticker, bwt->n - 1);
//    for (i = 0; i < bwt->n; ++i) {
//      const ch_t c = init(i);
//      ch_t j;
//      for (j = 0; j < 7; ++j) {
//        assert((bwt_LF_erank_encoded(bwt,j,i)-bwt->C[j])==occs[j]);
//      }
//      assert(bwt_char_encoded(bwt,i)==bwt_dna_p_encode[c]);
//      assert(bwt_LF(bwt,i)==bwt_LF_erank(bwt,c,i));
//      if (!GEM_MOD_BY_POWER_OF_TWO_64(i,sampling_rate_log)) {
//        assert(bwt_sampled_LF_bit(bwt,pos));
//        assert( (!i) || (dist < (1<<sampling_rate_log)) );
//        dist = 0;
//      } else {
//        dist++;
//      }
//      pos = bwt_LF(bwt, pos);
//      ++occs[bwt_dna_p_encode[c]];
//      if (be_verbose) ticker_update(ticker, i);
//    }
//  }
//  ticker_delete(ticker);
  return bwt_builder;
}
GEM_INLINE bwt_t* bwt_builder_get_bwt(bwt_builder_t* const bwt_builder) {
  return &bwt_builder->bwt;
}
GEM_INLINE void bwt_builder_write(fm_t* const file_manager,bwt_builder_t* const bwt_builder) {
  FM_CHECK(file_manager);
  BWT_BUILDER_CHECK(bwt_builder);
  /* Meta-Data */
  fm_write_uint64(file_manager,bwt_builder->bwt.length); // Length of the BWT
  fm_write_mem(file_manager,bwt_builder->bwt.c,BWT_MAYOR_COUNTER_RANGE*UINT64_SIZE); // Occurrences of each character
  fm_write_mem(file_manager,bwt_builder->bwt.C,BWT_MAYOR_COUNTER_RANGE*UINT64_SIZE); // The cumulative occurrences
  /* BWT Buckets-Structure */
  fm_write_uint64(file_manager,bwt_builder->bwt.num_minor_blocks); // Total number of minor blocks
  fm_write_uint64(file_manager,bwt_builder->bwt.num_mayor_blocks); // Total number of mayor blocks
  /* Mayor Counters [Aligned to 4KB] */
  fm_skip_align_4KB(file_manager);
  fm_write_mem(file_manager,bwt_builder->bwt.mayor_counters,
      bwt_builder->bwt.num_mayor_blocks*BWT_MAYOR_COUNTER_RANGE*BWT_MAYOR_COUNTER_SIZE);
  /* BWT Structure [Aligned to 4KB] */
  fm_skip_align_4KB(file_manager);
  fm_write_mem(file_manager,bwt_builder->bwt.bwt_mem,bwt_builder->bwt.num_minor_blocks*BWT_MINOR_BLOCK_SIZE);
  /* Extra-Layer {'|','N'} */
  sparse_array_locator_write(file_manager,bwt_builder->xl_locator); // 3-Bit layer locator
  fm_write_uint64(file_manager,bwt_builder->xl_memory->elements_used); // Total 3-Bit layer Bitmaps
  svector_write(file_manager,bwt_builder->xl_memory);
}
GEM_INLINE void bwt_builder_delete(bwt_builder_t* const bwt_builder) {
  BWT_BUILDER_CHECK(bwt_builder);
  mm_free(bwt_builder->bwt.c);
  mm_free(bwt_builder->bwt.C);
  mm_free(bwt_builder->bwt.mayor_counters);
  mm_free(bwt_builder->bwt.bwt_mem);
  sparse_array_locator_delete(bwt_builder->xl_locator);
  svector_delete(bwt_builder->xl_memory);
  mm_free(bwt_builder->mayor_occ);
  mm_free(bwt_builder->minor_occ);
  mm_free(bwt_builder->minor_occ_ext);
  mm_free(bwt_builder);
}
/*
 * BWT Loader
 */
GEM_INLINE bwt_t* bwt_read(fm_t* const file_manager,const bool check) {
  // Allocate handler
  bwt_t* const bwt = mm_alloc(struct _bwt_t);
  /* Meta-Data */
  bwt->length = fm_read_uint64(file_manager); // Length of the BWT
  bwt->c = mm_calloc(BWT_MAYOR_COUNTER_RANGE,uint64_t,true);
  fm_read_mem(file_manager,bwt->c,BWT_MAYOR_COUNTER_RANGE*UINT64_SIZE); // Occurrences of each character
  bwt->C = mm_calloc(BWT_MAYOR_COUNTER_RANGE,uint64_t,true);
  fm_read_mem(file_manager,bwt->C,BWT_MAYOR_COUNTER_RANGE*UINT64_SIZE); // The cumulative occurrences
  /* BWT Buckets-Structure */
  bwt->num_minor_blocks = fm_read_uint64(file_manager); // Total number of minor blocks
  bwt->num_mayor_blocks = fm_read_uint64(file_manager); // Total number of mayor blocks
  /* Mayor Counters [Aligned to 4KB] */
  fm_skip_align_4KB(file_manager);
  bwt->mm_mayor_counters = fm_load_mem(file_manager,
      bwt->num_mayor_blocks*BWT_MAYOR_COUNTER_RANGE*BWT_MAYOR_COUNTER_SIZE);
  bwt->mayor_counters = mm_get_base_mem(bwt->mm_mayor_counters);
  /* BWT Structure [Aligned to 4KB] */
  fm_skip_align_4KB(file_manager);
  bwt->mm_bwt_mem = fm_load_mem(file_manager,bwt->num_minor_blocks*BWT_MINOR_BLOCK_SIZE);
  bwt->bwt_mem = mm_get_base_mem(bwt->mm_mayor_counters);
  /* Extra-Layer {'|','N'} */
  bwt->xl_locator = sparse_array_locator_read(file_manager); // 3-Bit layer locator
  bwt->num_xl_bitmaps = fm_read_uint64(file_manager); // Total 3-Bit layer Bitmaps
  bwt->mm_xl_mem = fm_load_mem(file_manager,bwt->num_xl_bitmaps*sizeof(xl_block_t));
  bwt->xl_mem = mm_get_base_mem(bwt->mm_xl_mem);
  // Return
  return bwt;
}
GEM_INLINE bwt_t* bwt_read_mem(mm_t* const memory_manager,const bool check) {
  // Allocate handler
  bwt_t* const bwt = mm_alloc(struct _bwt_t);
  /* Meta-Data */
  bwt->length = mm_read_uint64(memory_manager); // Length of the BWT
  bwt->c = mm_read_mem(memory_manager,BWT_MAYOR_COUNTER_RANGE*UINT64_SIZE); // Occurrences of each character
  bwt->C = mm_read_mem(memory_manager,BWT_MAYOR_COUNTER_RANGE*UINT64_SIZE); // The cumulative occurrences
  /* BWT Buckets-Structure */
  bwt->num_minor_blocks = mm_read_uint64(memory_manager); // Total number of minor blocks
  bwt->num_mayor_blocks = mm_read_uint64(memory_manager); // Total number of mayor blocks
  /* Mayor Counters [Aligned to 4KB] */
  mm_skip_align_4KB(memory_manager);
  bwt->mm_mayor_counters = NULL;
  bwt->mayor_counters = mm_read_mem(memory_manager,
      bwt->num_mayor_blocks*BWT_MAYOR_COUNTER_RANGE*BWT_MAYOR_COUNTER_SIZE);
  /* BWT Structure [Aligned to 4KB] */
  mm_skip_align_4KB(memory_manager);
  bwt->mm_bwt_mem = NULL;
  bwt->bwt_mem = mm_read_mem(memory_manager,bwt->num_minor_blocks*BWT_MINOR_BLOCK_SIZE);
  /* Extra-Layer {'|','N'} */
  bwt->xl_locator = sparse_array_locator_read_mem(memory_manager); // 3-Bit layer locator
  bwt->num_xl_bitmaps = mm_read_uint64(memory_manager); // Total 3-Bit layer Bitmaps
  bwt->mm_xl_mem = NULL;
  bwt->xl_mem = mm_read_mem(memory_manager,bwt->num_xl_bitmaps*sizeof(xl_block_t));
  // Return
  return bwt;
}
GEM_INLINE void bwt_delete(bwt_t* const bwt) {
  BWT_CHECK(bwt);
  // Free counters memory
  if (bwt->mm_mayor_counters) {
    mm_free(bwt->c);
    mm_free(bwt->c);
    mm_bulk_free(bwt->mm_mayor_counters);
  }
  // Free bwt memory
  if (bwt->mm_bwt_mem) mm_bulk_free(bwt->mm_bwt_mem);
  // Free Extra-Layer
  sparse_array_locator_delete(bwt->xl_locator);
  if (bwt->mm_xl_mem) mm_bulk_free(bwt->mm_xl_mem);
}
/*
 * BWT General Accessors
 */
GEM_INLINE uint64_t bwt_builder_get_length(const bwt_builder_t* const bwt_builder) {
  BWT_BUILDER_CHECK(bwt_builder);
  return bwt_get_length(&bwt_builder->bwt);
}
GEM_INLINE uint64_t bwt_get_size_(bwt_t* const bwt,
    const uint64_t xl_locator_size,const uint64_t num_xl_bitmaps) {
  // Compute sizes
  const uint64_t minor_blocks_size = bwt->num_minor_blocks*BWT_MINOR_BLOCK_SIZE;
  const uint64_t mayor_counters_size = bwt->num_mayor_blocks*BWT_MAYOR_COUNTER_RANGE*UINT64_SIZE;
  const uint64_t xl_bitmap_size = num_xl_bitmaps*sizeof(xl_block_t);
  return (BWT_MAYOR_COUNTER_RANGE*UINT64_SIZE) +   /* bwt->c */
         (BWT_MAYOR_COUNTER_RANGE*UINT64_SIZE) +   /* bwt->C */
         minor_blocks_size +             /* bwt->bwt_mem */
         mayor_counters_size +           /* bwt->mayor_counters */
         xl_locator_size+xl_bitmap_size; /* bwt->xl */
}
GEM_INLINE uint64_t bwt_builder_get_total_size(bwt_builder_t* const bwt_builder) {
  return bwt_get_size_(&bwt_builder->bwt,
      bwt_builder->xl_locator->total_size,bwt_builder->xl_memory->elements_used);
}
GEM_INLINE uint64_t bwt_get_length(const bwt_t* const bwt) {
  BWT_CHECK(bwt);
  return bwt->length;
}
GEM_INLINE uint64_t bwt_get_size(bwt_t* const bwt) {
  return bwt_get_size_(bwt,
      bwt->xl_locator->total_size,bwt->num_xl_bitmaps);
}
/*
 * XOR table to mask bitmap depending based on the character (enc)
 */
const int64_t xor_table_1[] = {-1ll, -1ll, -1ll, -1ll, -0ll, -0ll, -0ll, -0ll};
const int64_t xor_table_2[] = {-1ll, -1ll, -0ll, -0ll, -1ll, -1ll, -0ll, -0ll};
const int64_t xor_table_3[] = {-1ll, -0ll, -1ll, -0ll, -1ll, -0ll, -1ll, -0ll};
/* Gets the mayor_counters & minor_block corresponding to position @i */
#define BWT_LOCATE_BLOCK_MAIN(bwt,position,block_pos,block_mod,mayor_counters,block_mem) \
  const uint64_t block_pos = position / BWT_MINOR_BLOCK_LENGTH; \
  const uint64_t block_mod = position % BWT_MINOR_BLOCK_LENGTH; \
  const uint64_t* const mayor_counters = bwt->mayor_counters + (position/BWT_MAYOR_BLOCK_LENGTH)*BWT_MAYOR_COUNTER_RANGE; \
  const uint64_t* const block_mem = bwt->bwt_mem + block_pos*BWT_MINOR_BLOCK_64WORDS
#define BWT_LOCATE_BLOCK_DNA(bwt,position,char_enc,block_pos,block_mod,mayor_counters,block_mem,has_xl_block,xl_block) \
  BWT_LOCATE_BLOCK_MAIN(bwt,position,block_pos,block_mod,mayor_counters,block_mem); \
  /* Locate XL */ \
  const bool has_xl_block = sparse_array_locator_is_marked(bwt->xl_locator,block_pos); \
  const xl_block_t* const xl_block = (gem_expect_false(has_xl_block)) ? \
    bwt->xl_mem + sparse_array_locator_get_erank(bwt->xl_locator,block_pos) : NULL
#define BWT_LOCATE_BLOCK_EXT(bwt,position,char_enc,block_pos,block_mod,mayor_counters,block_mem,has_xl_block,xl_block) \
  BWT_LOCATE_BLOCK_MAIN(bwt,position,block_pos,block_mod,mayor_counters,block_mem); \
  /* Locate XL */ \
  uint64_t erank; \
  const bool has_xl_block = sparse_array_locator_get_erank__marked(bwt->xl_locator,block_pos,&erank); \
  const xl_block_t* const xl_block = bwt->xl_mem + erank
/* Gets the mayor_counters & minor_block corresponding to position @i */
GEM_INLINE void bwt_get_block_location_dna(
    const bwt_t* const bwt,const uint64_t position,
    bwt_block_locator_t* const block_loc) {
  block_loc->block_pos = position / BWT_MINOR_BLOCK_LENGTH;
  block_loc->block_mod = position % BWT_MINOR_BLOCK_LENGTH;
  block_loc->mayor_counters = bwt->mayor_counters + (position/BWT_MAYOR_BLOCK_LENGTH)*BWT_MAYOR_COUNTER_RANGE;
  block_loc->block_mem = bwt->bwt_mem + block_loc->block_pos*BWT_MINOR_BLOCK_64WORDS;
  // Locate XL
  block_loc->has_xl_block = sparse_array_locator_is_marked(bwt->xl_locator,block_loc->block_pos);
  if (gem_expect_false(block_loc->has_xl_block)) {
    block_loc->xl_block = bwt->xl_mem+sparse_array_locator_get_erank(bwt->xl_locator,block_loc->block_pos);
  }
}
GEM_INLINE void bwt_get_block_location_ext(
    const bwt_t* const bwt,const uint64_t position,
    bwt_block_locator_t* const block_loc) {
  block_loc->block_pos = position / BWT_MINOR_BLOCK_LENGTH;
  block_loc->block_mod = position % BWT_MINOR_BLOCK_LENGTH;
  block_loc->mayor_counters = bwt->mayor_counters + (position/BWT_MAYOR_BLOCK_LENGTH)*BWT_MAYOR_COUNTER_RANGE;
  block_loc->block_mem = bwt->bwt_mem + block_loc->block_pos*BWT_MINOR_BLOCK_64WORDS;
  // Locate XL
  uint64_t erank;
  block_loc->has_xl_block = sparse_array_locator_get_erank__marked(bwt->xl_locator,block_loc->block_pos,&erank);
  block_loc->xl_block = bwt->xl_mem+erank;
}
/* Prefetches the bwt-block */
#define BWT_PREFETCH_BLOCK(block_loc) \
  PREFETCH(block_loc->mayor_counters); \
  PREFETCH(block_loc->block_mem); \
  PREFETCH(block_loc->block_mem+1); \
  PREFETCH(block_loc->block_mem+2)
/* Computes and returns the encoded letter */
GEM_INLINE uint8_t bwt_char_(
    const bwt_t* const bwt,
    const uint64_t block_pos,const uint64_t block_mod,const uint64_t* const block_mem,
    const bool has_xl_block,const xl_block_t* const xl_block) {
  const uint64_t letter_mask = 1ull << block_mod;
  const uint8_t bit_1 = ((block_mem[1] & letter_mask) != 0);
  const uint8_t bit_2 = ((block_mem[2] & letter_mask) != 0);
  if (gem_expect_false(has_xl_block)) {
    const uint8_t bit_3 = ((xl_block->layer_2 & letter_mask) != 0);
    return (bit_3 << 2) | (bit_2 << 1) | bit_1;
  } else {
    return (bit_2 << 1) | bit_1;
  }
}
/* Computes and returns the rank */
GEM_INLINE uint64_t bwt_erank_dna(
    const bwt_t* const bwt,const uint8_t char_enc,
    const uint64_t block_pos,const uint64_t block_mod,
    const uint64_t* const mayor_counters,const uint64_t* const block_mem,
    const bool has_xl_block,const xl_block_t* const xl_block) {
  // Fetching Regular DNA Characters
  const uint64_t sum_counters = mayor_counters[char_enc] + ((uint16_t*)block_mem)[char_enc];
  uint64_t bitmap = (block_mem[1]^xor_table_3[char_enc]) & (block_mem[2]^xor_table_2[char_enc]);
  if (gem_expect_false(has_xl_block)) bitmap &= (xl_block->layer_2^xor_table_1[char_enc]);
  // Return rank
  return sum_counters + POPCOUNT_64(bitmap & uint64_erank_mask(block_mod));
}
GEM_INLINE uint64_t bwt_erank_ext(
    const bwt_t* const bwt,const uint8_t char_enc,
    const uint64_t block_pos,const uint64_t block_mod,
    const uint64_t* const mayor_counters,const uint64_t* const block_mem,
    const bool has_xl_block,const xl_block_t* const xl_block) {
  // Fetching XL (Extra Layer)
  const uint64_t sum_counters = mayor_counters[char_enc] + xl_block->counter[char_enc-DNA_RANGE];
  if (gem_expect_false(has_xl_block)) {
    return sum_counters +
           POPCOUNT_64((block_mem[1]^xor_table_3[char_enc]) &
                       (block_mem[2]^xor_table_2[char_enc]) &
                       (xl_block->layer_2^xor_table_1[char_enc]) &
                        uint64_erank_mask(block_mod));
  } else {
    return sum_counters;
  }
}
/* Computes and returns the rank of an interval located in the same block */
GEM_INLINE void bwt_erank_interval_dna(
    const bwt_t* const bwt,const uint8_t char_enc,const uint64_t lo_value,
    const uint64_t block_pos,const uint64_t block_mod,
    const uint64_t* const mayor_counters,const uint64_t* const block_mem,
    const bool has_xl_block,const xl_block_t* const xl_block,
    uint64_t* const lo,uint64_t* const hi) {
  // Fetching Regular DNA Characters
  const uint64_t sum_counters = mayor_counters[char_enc] + ((uint16_t*)block_mem)[char_enc];
  uint64_t bitmap = (block_mem[1]^xor_table_3[char_enc]) & (block_mem[2]^xor_table_2[char_enc]);
  if (gem_expect_false(has_xl_block)) bitmap &= (xl_block->layer_2^xor_table_1[char_enc]);
  *hi = sum_counters +
        ( POPCOUNT_64(bitmap & uint64_erank_mask(block_mod)) );
  *lo = sum_counters +
        ( POPCOUNT_64(bitmap & uint64_erank_mask(lo_value % BWT_MINOR_BLOCK_LENGTH)) );
}
GEM_INLINE void bwt_erank_interval_ext(
    const bwt_t* const bwt,const uint8_t char_enc,const uint64_t lo_value,
    const uint64_t block_pos,const uint64_t block_mod,
    const uint64_t* const mayor_counters,const uint64_t* const block_mem,
    const bool has_xl_block,const xl_block_t* const xl_block,
    uint64_t* const lo,uint64_t* const hi) {
  // Fetching XL (Extra Layer)
  const uint64_t sum_counters = mayor_counters[char_enc] + xl_block->counter[char_enc-DNA_RANGE];
  if (gem_expect_false(has_xl_block)) {
    const uint64_t bitmap = (block_mem[1]^xor_table_3[char_enc]) &
                            (block_mem[2]^xor_table_2[char_enc]) &
                            (xl_block->layer_2^xor_table_1[char_enc]);
    *hi = sum_counters +
          ( POPCOUNT_64(bitmap & uint64_erank_mask(block_mod)) );
    *lo = sum_counters +
          ( POPCOUNT_64(bitmap & uint64_erank_mask(lo_value % BWT_MINOR_BLOCK_LENGTH)) );
  } else {
    *hi = sum_counters;
    *lo = sum_counters;
  }
}
/* Pre-computes the block elements (faster computation of all possible ranks) */
GEM_INLINE void bwt_precompute_(
    const bwt_t* const bwt,
    const uint64_t block_pos,const uint64_t block_mod,
    const uint64_t* const mayor_counters,const uint64_t* const block_mem,
    const bool has_xl_block,const xl_block_t* const xl_block,
    bwt_block_elms_t* const bwt_block_elms) {
  const uint64_t erase_mask = uint64_erank_mask(block_mod);
  const uint64_t bitmap_1 = block_mem[1] & erase_mask;
  const uint64_t bitmap_N1 = (block_mem[1]^(-1ull)) & erase_mask;
  const uint64_t bitmap_2 = block_mem[2];
  const uint64_t bitmap_N2 = bitmap_2^(-1ull);
  bwt_block_elms->bitmap_1__2[0] = bitmap_N1 & bitmap_N2;
  bwt_block_elms->bitmap_1__2[1] = bitmap_N1 & bitmap_2;
  bwt_block_elms->bitmap_1__2[2] = bitmap_1  & bitmap_N2;
  bwt_block_elms->bitmap_1__2[3] = bitmap_1  & bitmap_2;
  // Fetching Extra Characters
  if (gem_expect_false(has_xl_block)) {
    const uint64_t bitmap_3 = xl_block->layer_2;
    bwt_block_elms->bitmap_3[0] = (bitmap_3^(-1ull));
    bwt_block_elms->bitmap_3[1] = bitmap_3;
  }
}
/*
 * BWT Character Accessors
 */
GEM_INLINE uint8_t bwt_char(const bwt_t* const bwt,const uint64_t position) {
  BWT_CHAR_TICK();
  /* Locate Block */
  const uint64_t block_pos = position / BWT_MINOR_BLOCK_LENGTH;
  const uint64_t block_mod = position % BWT_MINOR_BLOCK_LENGTH;
  const uint64_t* const block_mem = bwt->bwt_mem + block_pos*BWT_MINOR_BLOCK_64WORDS;
  /* Locate XL */
  const bool has_xl_block = sparse_array_locator_is_marked(bwt->xl_locator,block_pos);
  const xl_block_t* const xl_block = (gem_expect_false(has_xl_block)) ?
    bwt->xl_mem + sparse_array_locator_get_erank(bwt->xl_locator,block_pos) : NULL;
  return bwt_char_(bwt,block_pos,block_mod,block_mem,has_xl_block,xl_block);
}
GEM_INLINE char bwt_char_character(const bwt_t* const bwt,const uint64_t position) {
  return dna_decode(bwt_char(bwt,position));
}
/*
 * BWT ERank (Exclusive Rank Function)
 */
GEM_INLINE uint64_t bwt_builder_erank(const bwt_builder_t* const bwt_builder,const uint8_t char_enc,const uint64_t position) {
  BWT_ERANK_TICK();
  BWT_LOCATE_BLOCK_MAIN((&bwt_builder->bwt),position,block_pos,block_mod,mayor_counters,block_mem);
  /* Locate XL */
  const bool has_xl_block = sparse_array_locator_is_marked(bwt_builder->xl_locator,block_pos);
  const xl_block_t* const xl_block = (gem_expect_false(has_xl_block)) ?
    svector_get_element(bwt_builder->xl_memory,
        sparse_array_locator_get_erank(bwt_builder->xl_locator,block_pos),xl_block_t) : NULL;
  return bwt_erank_dna(&bwt_builder->bwt,char_enc,block_pos,block_mod,mayor_counters,block_mem,has_xl_block,xl_block);
}
GEM_INLINE uint64_t bwt_erank(const bwt_t* const bwt,const uint8_t char_enc,const uint64_t position) {
  BWT_ERANK_TICK();
  if (gem_expect_true(char_enc < DNA_RANGE)) {
    BWT_LOCATE_BLOCK_DNA(bwt,position,char_enc,block_pos,block_mod,mayor_counters,block_mem,has_xl_block,xl_block);
    return bwt_erank_dna(bwt,char_enc,block_pos,block_mod,mayor_counters,block_mem,has_xl_block,xl_block);
  } else {
    BWT_LOCATE_BLOCK_EXT(bwt,position,char_enc,block_pos,block_mod,mayor_counters,block_mem,has_xl_block,xl_block);
    return bwt_erank_ext(bwt,char_enc,block_pos,block_mod,mayor_counters,block_mem,has_xl_block,xl_block);
  }
}
GEM_INLINE uint64_t bwt_erank_character(const bwt_t* const bwt,const char character,const uint64_t position) {
  return bwt_erank(bwt,dna_encode(character),position);
}
GEM_INLINE void bwt_erank_interval( // FIXME: Should return bool?
    const bwt_t* const bwt,const uint8_t char_enc,
    const uint64_t lo_in,const uint64_t hi_in,uint64_t* const lo_out,uint64_t* const hi_out) {
  BWT_ERANK_INTERVAL_TICK();
  if (gem_expect_true(char_enc < DNA_RANGE)) {
    BWT_LOCATE_BLOCK_DNA(bwt,hi_in,char_enc,block_pos,block_mod,mayor_counters,block_mem,has_xl_block,xl_block);
    bwt_erank_interval_dna(bwt,char_enc,lo_in,block_pos,block_mod,
        mayor_counters,block_mem,has_xl_block,xl_block,lo_out,hi_out);
  } else {
    BWT_LOCATE_BLOCK_EXT(bwt,hi_in,char_enc,block_pos,block_mod,mayor_counters,block_mem,has_xl_block,xl_block);
    bwt_erank_interval_ext(bwt,char_enc,lo_in,block_pos,block_mod,
        mayor_counters,block_mem,has_xl_block,xl_block,lo_out,hi_out);
  }
}
/*
 * BWT Prefetched ERank
 */
GEM_INLINE void bwt_prefetch(const bwt_t* const bwt,const uint64_t position,bwt_block_locator_t* const block_loc) {
  BWT_PREFETCH_TICK();
  bwt_get_block_location_ext(bwt,position,block_loc);
  BWT_PREFETCH_BLOCK(block_loc);
}
GEM_INLINE uint64_t bwt_prefetched_erank(
    const bwt_t* const bwt,const uint8_t char_enc,
    const uint64_t position,const bwt_block_locator_t* const block_loc) {
  BWT_ERANK_TICK();
  return (gem_expect_true(char_enc < DNA_RANGE)) ?
    bwt_erank_dna(bwt,char_enc,
        block_loc->block_pos,block_loc->block_mod,
        block_loc->mayor_counters,block_loc->block_mem,
        block_loc->has_xl_block,block_loc->xl_block) :
    bwt_erank_ext(bwt,char_enc,
        block_loc->block_pos,block_loc->block_mod,
        block_loc->mayor_counters,block_loc->block_mem,
        block_loc->has_xl_block,block_loc->xl_block);
}
GEM_INLINE void bwt_prefetched_erank_interval(  // FIXME: Should be returning bool
    const bwt_t* const bwt,const uint8_t char_enc,
    const uint64_t lo_in,const uint64_t hi_in,uint64_t* const lo_out,uint64_t* const hi_out,
    const bwt_block_locator_t* const block_loc) {
  BWT_ERANK_INTERVAL_TICK();
  return (gem_expect_true(char_enc < DNA_RANGE)) ?
    bwt_erank_interval_dna(bwt,char_enc,lo_in,
        block_loc->block_pos,block_loc->block_mod,
        block_loc->mayor_counters,block_loc->block_mem,
        block_loc->has_xl_block,block_loc->xl_block,
        lo_out,hi_out) :
    bwt_erank_interval_ext(bwt,char_enc,lo_in,
        block_loc->block_pos,block_loc->block_mod,
        block_loc->mayor_counters,block_loc->block_mem,
        block_loc->has_xl_block,block_loc->xl_block,
        lo_out,hi_out);
}
/*
 *  BWT Precomputed ERank (Precomputation of the block's elements)
 */
GEM_INLINE void bwt_precompute(
    const bwt_t* const bwt,const uint64_t position,
    bwt_block_locator_t* const block_loc,bwt_block_elms_t* const block_elms) {
  BWT_PRECOMPUTE_TICK();
  bwt_get_block_location_ext(bwt,position,block_loc);
  bwt_precompute_(bwt,
      block_loc->block_pos,block_loc->block_mod,
      block_loc->mayor_counters,block_loc->block_mem,
      block_loc->has_xl_block,block_loc->xl_block,
      block_elms);
}
GEM_INLINE void bwt_precompute_interval(
    const bwt_t* const bwt,const uint64_t lo,const uint64_t hi,
    bwt_block_locator_t* const block_loc,bwt_block_elms_t* const block_elms) {
  BWT_PRECOMPUTE_TICK();
  bwt_precompute(bwt,hi,block_loc,block_elms);
  block_elms->gap_mask = uint64_erank_inv_mask(lo % BWT_MINOR_BLOCK_LENGTH);
}
GEM_INLINE void bwt_prefetched_precompute(
    const bwt_t* const bwt,
    const bwt_block_locator_t* const block_loc,bwt_block_elms_t* const block_elms) {
  BWT_PREFETCH__PRECOMPUTE_TICK();
  bwt_precompute_(bwt,
      block_loc->block_pos,block_loc->block_mod,
      block_loc->mayor_counters,block_loc->block_mem,
      block_loc->has_xl_block,block_loc->xl_block,
      block_elms);
}
GEM_INLINE void bwt_prefetched_precompute_interval(
    const bwt_t* const bwt,const uint64_t lo,
    const bwt_block_locator_t* const block_loc,bwt_block_elms_t* const block_elms) {
  BWT_PREFETCH__PRECOMPUTE_TICK();
  bwt_prefetched_precompute(bwt,block_loc,block_elms);
  block_elms->gap_mask = uint64_erank_inv_mask(lo % BWT_MINOR_BLOCK_LENGTH);
}
GEM_INLINE uint64_t bwt_precomputed_erank(
    const bwt_t* const bwt,const uint8_t char_enc,
    const bwt_block_locator_t* const block_loc,const bwt_block_elms_t* const block_elms) {
  BWT_ERANK_TICK();
  // Check character
  if (gem_expect_true(char_enc < DNA_RANGE)) {
    // Fetching Regular DNA Characters
    return block_loc->mayor_counters[char_enc] + ((uint16_t*)block_loc->block_mem)[char_enc] +
           POPCOUNT_64(block_elms->bitmap_1__2[char_enc & 011]);
  } else {
    // Fetching Extra Characters
    if (gem_expect_false(block_loc->has_xl_block)) {
      return block_loc->mayor_counters[char_enc] + block_loc->xl_block->counter[char_enc-DNA_RANGE] +
             POPCOUNT_64(block_elms->bitmap_1__2[char_enc & 011] & block_elms->bitmap_3[char_enc>>2]);
    } else {
      return block_loc->mayor_counters[char_enc] + block_loc->xl_block->counter[char_enc-DNA_RANGE];
    }
  }
}
GEM_INLINE bool bwt_precomputed_erank_interval(
    const bwt_t* const bwt,const uint8_t char_enc,
    uint64_t* const lo_out,uint64_t* const hi_out,
    const bwt_block_locator_t* const block_loc,const bwt_block_elms_t* const block_elms) {
  BWT_ERANK_INTERVAL_TICK();
  // Check character
  if (gem_expect_true(char_enc < DNA_RANGE)) {
    // Fetching Regular DNA Characters
    const uint64_t bitmap = block_elms->bitmap_1__2[char_enc & 011];
    const uint64_t bitmap_gap = bitmap & block_elms->gap_mask;
    if (gem_expect_false(bitmap_gap==0)) return false; // Empty interval
    *hi_out = block_loc->mayor_counters[char_enc] +
              ((uint16_t*)block_loc->block_mem)[char_enc] + POPCOUNT_64(bitmap);
    *lo_out = *hi_out - POPCOUNT_64(bitmap_gap);
    return true;
  } else {
    // Fetching Extra Characters
    if (gem_expect_false(!block_loc->has_xl_block)) return false; // Empty interval
    const uint64_t bitmap = block_elms->bitmap_1__2[char_enc & 011] & block_elms->bitmap_3[char_enc>>2];
    const uint64_t bitmap_gap = bitmap & block_elms->gap_mask;
    if (gem_expect_false(bitmap_gap==0)) return false;
    *hi_out = block_loc->mayor_counters[char_enc] +
              block_loc->xl_block->counter[char_enc-DNA_RANGE] +
              POPCOUNT_64(bitmap);
    *lo_out = *hi_out - POPCOUNT_64(bitmap_gap);
    return true;
  }
}
/*
 * BWT LF (Last to first)
 */
GEM_INLINE uint64_t bwt_LF(
    const bwt_t* const bwt,const uint8_t position) {
  uint8_t char_enc;
  return bwt_LF__enc(bwt,position,&char_enc);
}
GEM_INLINE uint64_t bwt_prefetched_LF(
    const bwt_t* const bwt,const uint8_t position,
    const bwt_block_locator_t* const block_loc) {
  uint8_t char_enc;
  return bwt_prefetched_LF__enc(bwt,position,&char_enc,block_loc);
}
GEM_INLINE uint64_t bwt_LF__enc(
    const bwt_t* const bwt,const uint8_t position,uint8_t* const char_enc) {
  BWT_LF_TICK();
  BWT_LOCATE_BLOCK_EXT(bwt,position,char_enc,block_pos,block_mod,mayor_counters,block_mem,has_xl_block,xl_block);
  *char_enc = bwt_char_(bwt,block_pos,block_mod,block_mem,has_xl_block,xl_block);
  return (gem_expect_true(*char_enc < DNA_RANGE)) ?
      bwt_erank_dna(bwt,*char_enc,block_pos,block_mod,mayor_counters,block_mem,has_xl_block,xl_block) :
      bwt_erank_ext(bwt,*char_enc,block_pos,block_mod,mayor_counters,block_mem,has_xl_block,xl_block);
}
GEM_INLINE uint64_t bwt_LF__character(
    const bwt_t* const bwt,const uint64_t position,char* const character) {
  BWT_LF_TICK();
  uint8_t char_enc;
  const uint64_t rank_LF = bwt_LF__enc(bwt,position,&char_enc);
  *character = dna_decode(char_enc);
  return rank_LF;
}
GEM_INLINE uint64_t bwt_prefetched_LF__enc(
    const bwt_t* const bwt,const uint8_t position,uint8_t* const char_enc,
    const bwt_block_locator_t* const block_loc) {
  BWT_LF_TICK();
  *char_enc = bwt_char_(bwt,
      block_loc->block_pos,block_loc->block_mod,block_loc->block_mem,
      block_loc->has_xl_block,block_loc->xl_block);
  return (gem_expect_true(*char_enc < DNA_RANGE)) ?
      bwt_erank_dna(bwt,*char_enc,
          block_loc->block_pos,block_loc->block_mod,block_loc->mayor_counters,
          block_loc->block_mem,block_loc->has_xl_block,block_loc->xl_block) :
      bwt_erank_ext(bwt,*char_enc,
          block_loc->block_pos,block_loc->block_mod,block_loc->mayor_counters,
          block_loc->block_mem,block_loc->has_xl_block,block_loc->xl_block);
}
/*
 * Display
 */
GEM_INLINE void bwt_print_(FILE* const stream,bwt_t* const bwt,
    const uint64_t xl_locator_size,const uint64_t num_xl_bitmaps) {
  // Compute sizes
  const uint64_t minor_blocks_size = bwt->num_minor_blocks*BWT_MINOR_BLOCK_SIZE;
  const uint64_t mayor_counters_size = bwt->num_mayor_blocks*BWT_MAYOR_COUNTER_RANGE*UINT64_SIZE;
  const uint64_t xl_bitmap_size = num_xl_bitmaps*sizeof(xl_block_t);
  const uint64_t bwt_total_size =
      (BWT_MAYOR_COUNTER_RANGE*UINT64_SIZE)+(BWT_MAYOR_COUNTER_RANGE*UINT64_SIZE)+
      minor_blocks_size+mayor_counters_size+xl_locator_size+xl_bitmap_size;
  const uint64_t minor_counters_size = bwt->num_minor_blocks*BWT_MINOR_COUNTER_RANGE*BWT_MINOR_BLOCK_COUNTER_LENGTH/8;
  const uint64_t minor_bitmap_size = bwt->num_minor_blocks*BWT_MINOR_BLOCK_LENGTH*BWT_MINOR_BLOCK_BITS/8;
  // Display BWT
  tab_fprintf(stream,"[GEM]>BWT\n");
  tab_fprintf(stream,"  => Architecture\tGBWT.2l.1s.C64.c16.2bm64.xl\n");
  tab_fprintf(stream,"    => Buckets 2-levels\n");
  tab_fprintf(stream,"    => Step    1-step\n");
  tab_fprintf(stream,"    => MayorCounters.length 64bits\n");
  tab_fprintf(stream,"    => MinorCounters.length 16bits\n");
  tab_fprintf(stream,"    => Bitmap.Bitwise 2x64bits\n");
  tab_fprintf(stream,"      => Extra.Layer  Sparse.1x64bits (3rd-bit + xl.counters)\n");
  tab_fprintf(stream,"  => Total.length %lu\n",bwt->length);
  tab_fprintf(stream,"  => Total.Size %lu\n",bwt_total_size);
  tab_fprintf(stream,"    => Mayor.counters %lu (%lu MB) [%2.3f%%]\n",
      bwt->num_mayor_blocks,CONVERT_B_TO_MB(mayor_counters_size),
      PERCENTAGE(mayor_counters_size,bwt_total_size));
  tab_fprintf(stream,"    => Minor.Blocks %lu (%lu MB) [%2.3f%%]\n",
      bwt->num_minor_blocks,CONVERT_B_TO_MB(minor_blocks_size),
      PERCENTAGE(minor_blocks_size,bwt_total_size));
  tab_fprintf(stream,"      => Minor.Counters (%lu MB) [%2.3f%%]\n",
      CONVERT_B_TO_MB(minor_counters_size),
      PERCENTAGE(minor_counters_size,bwt_total_size));
  tab_fprintf(stream,"      => Minor.Bitmap (%lu MB) [%2.3f%%]\n",
      CONVERT_B_TO_MB(minor_bitmap_size),
      PERCENTAGE(minor_bitmap_size,bwt_total_size));
  tab_fprintf(stream,"    => Sparse.Locator (%lu MB) [%2.3f%%]\n",
      CONVERT_B_TO_MB(xl_locator_size),
      PERCENTAGE(xl_locator_size,bwt_total_size));
  tab_fprintf(stream,"    => Sparse.Bitmap (%lu MB) [%2.3f%%]\n",
      CONVERT_B_TO_MB(xl_bitmap_size),
      PERCENTAGE(xl_bitmap_size,bwt_total_size));
  tab_fprintf(stream,"  => Occurrences\tA=[%lu] C=[%lu] G=[%lu] T=[%lu] N=[%lu] |=[%lu]\n",
      bwt->c[ENC_DNA_CHAR_A],bwt->c[ENC_DNA_CHAR_C],bwt->c[ENC_DNA_CHAR_G],
      bwt->c[ENC_DNA_CHAR_T],bwt->c[ENC_DNA_CHAR_N],bwt->c[ENC_DNA_CHAR_SEP]);
  tab_fprintf(stream,"  => Cumulative.Occ\tA=[%lu] C=[%lu] G=[%lu] T=[%lu] N=[%lu] |=[%lu]\n",
      bwt->C[ENC_DNA_CHAR_A],bwt->C[ENC_DNA_CHAR_C],bwt->C[ENC_DNA_CHAR_G],
      bwt->C[ENC_DNA_CHAR_T],bwt->C[ENC_DNA_CHAR_N],bwt->C[ENC_DNA_CHAR_SEP]);
  /*
   * Stats // todo
   */

  // Flush
  fflush(stream);
}
GEM_INLINE void bwt_builder_print(FILE* const stream,bwt_builder_t* const bwt_builder) {
  bwt_print_(stream,&bwt_builder->bwt,bwt_builder->xl_locator->total_size,bwt_builder->xl_memory->elements_used);
  tab_global_inc();
  sparse_array_locator_print(stream,bwt_builder->xl_locator,false);
  tab_global_dec();
}
GEM_INLINE void bwt_print(FILE* const stream,bwt_t* const bwt) {
  bwt_print_(stream,bwt,bwt->xl_locator->total_size,bwt->num_xl_bitmaps);
  tab_global_inc();
  sparse_array_locator_print(stream,bwt->xl_locator,false);
  tab_global_dec();
}

