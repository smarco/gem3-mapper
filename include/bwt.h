/*
 * PROJECT: GEMMapper
 * FILE: bwt.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Provides basic routines to encode a DNA text into a Burrows-Wheeler transform
 *              using a compact bit representation and counter buckets as to enhance Occ/rank queries
 */

#ifndef BWT_H_
#define BWT_H_

#include "essentials.h"

#include "dna_string.h"
#include "dna_text.h"

/*
 * Check
 */
typedef struct _bwt_t bwt_t;
typedef struct _bwt_builder_t bwt_builder_t;
typedef struct _bwt_block_locator_t bwt_block_locator_t;
typedef struct _bwt_block_elms_t bwt_block_elms_t;
typedef struct {
  uint64_t lo;
  uint64_t hi;
} bwt_interval_t;

/*
 * BWT Builder
 */
GEM_INLINE bwt_builder_t* bwt_builder_new(
    dna_text_t* const bwt_text,const uint64_t* const character_occurrences,
    const bool check,const bool verbose,mm_slab_t* const mm_slab);
GEM_INLINE bwt_t* bwt_builder_get_bwt(bwt_builder_t* const bwt_builder);
GEM_INLINE void bwt_builder_write(fm_t* const file_manager,bwt_builder_t* const bwt_builder);
GEM_INLINE void bwt_builder_delete(bwt_builder_t* const bwt_builder);

/*
 * BWT Loader
 */
GEM_INLINE bwt_t* bwt_read(fm_t* const file_manager,const bool check);
GEM_INLINE bwt_t* bwt_read_mem(mm_t* const memory_manager,const bool check);
GEM_INLINE void bwt_write(bwt_t* const bwt,fm_t* const file_manager);
GEM_INLINE void bwt_delete(bwt_t* const bwt);

/*
 * BWT General Accessors
 */
GEM_INLINE uint64_t bwt_builder_get_length(const bwt_builder_t* const bwt_builder);
GEM_INLINE uint64_t bwt_builder_get_size(bwt_builder_t* const bwt_builder);

GEM_INLINE uint64_t bwt_get_length(const bwt_t* const bwt);
GEM_INLINE uint64_t bwt_get_size(bwt_t* const bwt);

/*
 * BWT Character Accessors
 */
GEM_INLINE uint8_t bwt_char(const bwt_t* const bwt,const uint64_t position);
GEM_INLINE char bwt_char_character(const bwt_t* const bwt,const uint64_t position);

/*
 * BWT ERank (Exclusive Rank Function)
 */
GEM_INLINE uint64_t bwt_builder_erank(const bwt_builder_t* const bwt_builder,const uint8_t char_enc,const uint64_t position);
GEM_INLINE uint64_t bwt_erank(const bwt_t* const bwt,const uint8_t char_enc,const uint64_t position);
GEM_INLINE uint64_t bwt_erank_character(const bwt_t* const bwt,const char character,const uint64_t position);
GEM_INLINE void bwt_erank_interval(
    const bwt_t* const bwt,const uint8_t char_enc,
    const uint64_t lo_in,const uint64_t hi_in,uint64_t* const lo_out,uint64_t* const hi_out);

/*
 * BWT Prefetched ERank
 */
GEM_INLINE void bwt_prefetch(const bwt_t* const bwt,const uint64_t position,bwt_block_locator_t* const block_loc);
GEM_INLINE uint64_t bwt_prefetched_erank(
    const bwt_t* const bwt,const uint8_t char_enc,
    const uint64_t position,const bwt_block_locator_t* const block_loc);
GEM_INLINE void bwt_prefetched_erank_interval(
    const bwt_t* const bwt,const uint8_t char_enc,
    const uint64_t lo_in,const uint64_t hi_in,uint64_t* const lo_out,uint64_t* const hi_out,
    const bwt_block_locator_t* const block_loc);

/*
 *  BWT Precomputed ERank (Precomputation of the block's elements)
 */
GEM_INLINE void bwt_precompute(
    const bwt_t* const bwt,const uint64_t position,
    bwt_block_locator_t* const block_loc,bwt_block_elms_t* const block_elms);
GEM_INLINE void bwt_precompute_interval(
    const bwt_t* const bwt,const uint64_t lo,const uint64_t hi,
    bwt_block_locator_t* const block_loc,bwt_block_elms_t* const block_elms);
GEM_INLINE void bwt_prefetched_precompute(
    const bwt_t* const bwt,
    const bwt_block_locator_t* const block_loc,bwt_block_elms_t* const block_elms);
GEM_INLINE void bwt_prefetched_precompute_interval(
    const bwt_t* const bwt,const uint64_t lo,
    const bwt_block_locator_t* const block_loc,bwt_block_elms_t* const block_elms);

GEM_INLINE uint64_t bwt_precomputed_erank(
    const bwt_t* const bwt,const uint8_t char_enc,
    const bwt_block_locator_t* const block_loc,const bwt_block_elms_t* const block_elms);
GEM_INLINE bool bwt_precomputed_erank_interval(
    const bwt_t* const bwt,const uint8_t char_enc,
    uint64_t* const lo_out,uint64_t* const hi_out,
    const bwt_block_locator_t* const block_loc,const bwt_block_elms_t* const block_elms);

/*
 * BWT LF (Last to first)
 */
GEM_INLINE uint64_t bwt_LF(
    const bwt_t* const bwt,const uint8_t position);
GEM_INLINE uint64_t bwt_prefetched_LF(
    const bwt_t* const bwt,const uint8_t position,
    const bwt_block_locator_t* const block_loc);

GEM_INLINE uint64_t bwt_LF__enc(
    const bwt_t* const bwt,const uint8_t position,uint8_t* const char_enc);
GEM_INLINE uint64_t bwt_LF__character(
    const bwt_t* const bwt,const uint64_t position,char* const character);
GEM_INLINE uint64_t bwt_prefetched_LF__enc(
    const bwt_t* const bwt,const uint8_t position,uint8_t* const char_enc,
    const bwt_block_locator_t* const block_loc);

/*
 * Display
 */
GEM_INLINE void bwt_builder_print(FILE* const stream,bwt_builder_t* const bwt_builder);
GEM_INLINE void bwt_print(FILE* const stream,bwt_t* const bwt);

#endif /* BWT_H_ */
