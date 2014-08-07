/*
 * PROJECT: GEMMapper
 * FILE: rank_mtable.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Rank memoizated queries
 */

#ifndef RANK_MTABLE_H_
#define RANK_MTABLE_H_

#include "essentials.h"
#include "bwt.h"

/*
 * Constants
 */
#define RANK_MTABLE_SEARCH_DEPTH  11                          // Number of character that can be searched up in the table
#define RANK_MTABLE_LEVELS       (RANK_MTABLE_SEARCH_DEPTH+1) // One fake level (zero-HI)

/*
 * Check
 */
#define RANK_MTABLE_CHECK(rank_mtable) GEM_CHECK_NULL(rank_mtable)
#define RANK_MQUERY_CHECK(rank_query)  GEM_CHECK_NULL(rank_query)

typedef struct {
  // Meta-info
  uint64_t table_size;         // Total number of ranks stored
  uint64_t num_levels;         // Total depth of the table
  // Table
  uint64_t* level_skip;        // Skip from levels-to-level (Pre-computed)
  uint64_t** sa_ranks_levels;  // Pointers to the levels
  /* MM */
  mm_t* mm_sa_ranks;
} rank_mtable_t;

typedef struct {
  uint64_t hi_position; // Effective HI-position in the table
  uint64_t level;       // Level on the table
} rank_mquery_t;

/*
 * Loader/Setup
 */
GEM_INLINE rank_mtable_t* rank_mtable_read(fm_t* const file_manager);
GEM_INLINE rank_mtable_t* rank_mtable_read_mem(mm_t* const memory_manager);
GEM_INLINE void rank_mtable_write(fm_t* const file_manager,rank_mtable_t* const rank_mtable);
GEM_INLINE void rank_mtable_delete(rank_mtable_t* const rank_mtable);

/*
 * Builder
 */
GEM_INLINE rank_mtable_t* rank_mtable_builder_new(const bwt_builder_t* const bwt_builder,const bool verbose);
GEM_INLINE void rank_mtable_builder_delete(rank_mtable_t* const rank_mtable);

/*
 * Accessors
 */
GEM_INLINE uint64_t rank_mtable_get_size(const rank_mtable_t* const rank_mtable);

/*
 * Query
 */
GEM_INLINE void rank_mquery_new(rank_mquery_t* const query);
GEM_INLINE void rank_mquery_add_char(const rank_mtable_t* const rank_mtable,rank_mquery_t* const query,uint8_t const enc_char);
GEM_INLINE uint64_t rank_mquery_get_level(const rank_mquery_t* const query);
GEM_INLINE uint64_t rank_mquery_is_exhausted(const rank_mquery_t* const query);

/*
 * Fetch rank value
 */
GEM_INLINE void rank_mtable_fetch(
    const rank_mtable_t* const rank_mtable,const rank_mquery_t* const query,
    uint64_t* const lo,uint64_t* const hi);

/*
 * Display
 */
GEM_INLINE void rank_mtable_print(FILE* const stream,rank_mtable_t* const rank_mtable);

#endif /* RANK_MTABLE_H_ */
