/*
 * PROJECT: GEMMapper
 * FILE: bwt_commons.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef BWT_COMMONS_H_
#define BWT_COMMONS_H_

#include "utils/essentials.h"
#include "fm_index/sampled_sa.h"

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
extern uint64_t _bwt_ranks;

#define BWT_CHAR_TICK()                 PROF_BLOCK() { ++_bwt_ranks; }
#define BWT_ERANK_TICK()                PROF_BLOCK() { ++_bwt_ranks; }
#define BWT_ERANK_INTERVAL_TICK()       PROF_BLOCK() { ++_bwt_ranks; }
#define BWT_LF_TICK()                   PROF_BLOCK() { ++_bwt_ranks; }
#define BWT_SAMPLED_TICK()              PROF_BLOCK() { ++_bwt_ranks; }

#define BWT_PREFETCH_TICK()             PROF_BLOCK() { /* NOP */ }
#define BWT_PRECOMPUTE_TICK()           PROF_BLOCK() { /* NOP */ }

/*
 * Shared Constants/Tables
 */
extern const int64_t xor_table_1[];
extern const int64_t xor_table_2[];
extern const int64_t xor_table_3[];

/*
 * BWT Data-structures
 */
typedef struct {
  uint64_t lo;
  uint64_t hi;
} bwt_interval_t;
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
 * Shared Functions
 */
// Prefetch the BWT-block
#define BWT_PREFETCH_BLOCK(block_loc) \
  PREFETCH(block_loc->mayor_counters); \
  PREFETCH(block_loc->block_mem)

/*
 * Errors
 */
#define GEM_ERROR_BWT_WRONG_MODEL_NO "BWT error. Wrong BWT-Index Model %"PRIu64" (Expected model %"PRIu64")"

#endif /* BWT_COMMONS_H_ */
