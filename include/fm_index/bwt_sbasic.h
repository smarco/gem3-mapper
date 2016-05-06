/*
 * PROJECT: GEMMapper
 * FILE: bwt_sbasic.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Provides basic routines to encode a DNA text into a
 *              Burrows-Wheeler transform using a compact bit representation
 *              and counter buckets as to enhance Occ/rank queries
 * MODEL:
 *   NAME:      bwt_sbasic
 *   ALPHABET:  3bits (8 chars)
 *   LEVELS:    2-levels
 *   STEP:      1s
 *   COUNTERS:  (8x64b,8x16b)
 *   BITMAP:    64bits
 *   SAMPLED:   Yes
 *   OTHERS:    -
 */

#ifndef BWT_SBASIC_H_
#define BWT_SBASIC_H_

#include "fm_index/bwt_commons.h"
#include "data_structures/dna_text.h"

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
} bwt_sbasic_t;
typedef struct {
  /* Meta-Data */
  bwt_sbasic_t bwt;
  /* Auxiliary variables */
  uint64_t mayor_counter;             // Current mayorCounter position
  uint64_t* minor_block_mem;          // Pointer to current minorBlock position
  uint64_t* mayor_occ;                // Mayor accumulated counts
  uint64_t* minor_occ;                // Minor accumulated counts
  uint64_t sampling_mayor_count;      // Sampling-Mayor accumulated count
  uint64_t sampling_minor_count;      // Sampling-Minor accumulated count
} bwt_sbasic_builder_t;

/*
 * BWT Builder
 */
bwt_sbasic_builder_t* bwt_sbasic_builder_new(
    dna_text_t* const bwt_text,
    const uint64_t* const character_occurrences,
    sampled_sa_builder_t* const sampled_sa,
    const bool check,
    const bool verbose);
void bwt_sbasic_builder_write(
    fm_t* const file_manager,
    bwt_sbasic_builder_t* const bwt_builder);
void bwt_sbasic_builder_delete(bwt_sbasic_builder_t* const bwt_builder);

/*
 * BWT Loader
 */
bwt_sbasic_t* bwt_sbasic_read_mem(mm_t* const memory_manager,const bool check);
void bwt_sbasic_delete(bwt_sbasic_t* const bwt);

/*
 * BWT General Accessors
 */
uint64_t bwt_sbasic_builder_get_length(const bwt_sbasic_builder_t* const bwt_builder);
uint64_t bwt_sbasic_builder_get_size(bwt_sbasic_builder_t* const bwt_builder);

uint64_t bwt_sbasic_get_length(const bwt_sbasic_t* const bwt);
uint64_t bwt_sbasic_get_size(bwt_sbasic_t* const bwt);

bool bwt_sbasic_is_same_bucket(const uint64_t lo,const uint64_t hi);

/*
 * BWT Character Accessors
 */
uint8_t bwt_sbasic_char(
    const bwt_sbasic_t* const bwt,
    const uint64_t position);
char bwt_sbasic_char_character(
    const bwt_sbasic_t* const bwt,
    const uint64_t position);

/*
 * BWT ERank (Exclusive Rank Function)
 */
uint64_t bwt_sbasic_builder_erank(
    const bwt_sbasic_builder_t* const bwt_builder,
    const uint8_t char_enc,
    const uint64_t position);
uint64_t bwt_sbasic_erank(
    const bwt_sbasic_t* const bwt,
    const uint8_t char_enc,
    const uint64_t position);
uint64_t bwt_sbasic_erank_character(
    const bwt_sbasic_t* const bwt,
    const char character,
    const uint64_t position);
void bwt_sbasic_erank_interval(
    const bwt_sbasic_t* const bwt,
    const uint8_t char_enc,
    const uint64_t lo_in,
    const uint64_t hi_in,
    uint64_t* const lo_out,
    uint64_t* const hi_out);
uint64_t bwt_sbasic_sampling_erank(
    const bwt_sbasic_t* const bwt,
    const uint64_t position);

/*
 * BWT Prefetched ERank
 */
void bwt_sbasic_prefetch(
    const bwt_sbasic_t* const bwt,
    const uint64_t position,
    bwt_block_locator_t* const block_loc);
uint64_t bwt_sbasic_prefetched_erank(
    const bwt_sbasic_t* const bwt,
    const uint8_t char_enc,
    const uint64_t position,
    const bwt_block_locator_t* const block_loc);
void bwt_sbasic_prefetched_erank_interval(
    const bwt_sbasic_t* const bwt,
    const uint8_t char_enc,
    const uint64_t lo_in,
    const uint64_t hi_in,
    uint64_t* const lo_out,
    uint64_t* const hi_out,
    const bwt_block_locator_t* const block_loc);

/*
 *  BWT Precomputed ERank (Precomputation of the block's elements)
 */
void bwt_sbasic_precompute(
    const bwt_sbasic_t* const bwt,
    const uint64_t position,
    bwt_block_locator_t* const block_loc,
    bwt_block_elms_t* const block_elms);
void bwt_sbasic_precompute_interval(
    const bwt_sbasic_t* const bwt,
    const uint64_t lo,
    const uint64_t hi,
    bwt_block_locator_t* const block_loc,
    bwt_block_elms_t* const block_elms);
void bwt_sbasic_prefetched_precompute(
    const bwt_sbasic_t* const bwt,
    const bwt_block_locator_t* const block_loc,
    bwt_block_elms_t* const block_elms);
void bwt_sbasic_prefetched_precompute_interval(
    const bwt_sbasic_t* const bwt,
    const uint64_t lo,
    const bwt_block_locator_t* const block_loc,
    bwt_block_elms_t* const block_elms);

uint64_t bwt_sbasic_precomputed_erank(
    const bwt_sbasic_t* const bwt,
    const uint8_t char_enc,
    const bwt_block_locator_t* const block_loc,
    const bwt_block_elms_t* const block_elms);
void bwt_sbasic_precomputed_erank_interval(
    const bwt_sbasic_t* const bwt,
    const uint8_t char_enc,
    uint64_t* const lo_out,
    uint64_t* const hi_out,
    const bwt_block_locator_t* const block_loc,
    const bwt_block_elms_t* const block_elms);

/*
 * BWT LF (Last to first)
 */
uint64_t bwt_sbasic_LF_(
    const bwt_sbasic_t* const bwt,
    const uint64_t position,
    uint8_t* const char_enc);
uint64_t bwt_sbasic_LF(
    const bwt_sbasic_t* const bwt,
    const uint64_t position,
    bool* const is_sampled);
uint64_t bwt_sbasic_prefetched_LF(
    const bwt_sbasic_t* const bwt,
    const uint64_t position,
    bool* const is_sampled,
    const bwt_block_locator_t* const block_loc);

uint64_t bwt_sbasic_LF__enc(
    const bwt_sbasic_t* const bwt,
    const uint64_t position,
    uint8_t* const char_enc,
    bool* const is_sampled);
uint64_t bwt_sbasic_LF__character(
    const bwt_sbasic_t* const bwt,
    const uint64_t position,
    char* const character,
    bool* const is_sampled);
uint64_t bwt_sbasic_prefetched_LF__enc(
    const bwt_sbasic_t* const bwt,
    const uint64_t position,
    uint8_t* const char_enc,
    bool* const is_sampled,
    const bwt_block_locator_t* const block_loc);

/*
 * Display
 */
void bwt_sbasic_builder_print(
    FILE* const stream,
    bwt_sbasic_builder_t* const bwt_builder);
void bwt_sbasic_print(
    FILE* const stream,
    bwt_sbasic_t* const bwt);

#endif /* BWT_SBASIC_H_ */
