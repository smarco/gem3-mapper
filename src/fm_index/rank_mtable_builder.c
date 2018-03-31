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
 * DESCRIPTION:
 *   FM-Index module provides a builder for a rank-memoization table
 *   that stores pre-computed bwt-intevals for all the k-mer (k=11)
 *   words from {ACGT}
 */

#include "text/dna_text.h"
#include "fm_index/rank_mtable_builder.h"

/*
 * Write mtable
 */
void rank_mtable_builder_write(
    fm_t* const file_manager,
    rank_mtable_t* const rank_mtable) {
  // Write Meta-info
  fm_write_uint64(file_manager,rank_mtable->num_levels);
  fm_write_uint64(file_manager,rank_mtable->table_size);
  fm_write_uint64(file_manager,rank_mtable->min_matching_depth);
  // Write table
  fm_skip_align_4KB(file_manager);
  fm_write_mem(file_manager,rank_mtable->sa_ranks_levels[0],rank_mtable->table_size*UINT64_SIZE);
}
/*
 * Find minimum matching depth
 */
void rank_mtable_builder_find_mmd(
    rank_mtable_t* const rank_mtable,
    const uint64_t level,
    rank_mquery_t* const query,
    uint64_t* const min_matching_depth) {
  // Check number of matches
  uint64_t lo, hi;
  rank_mtable_fetch(rank_mtable,query,&lo,&hi);
  if (hi-lo <= RANK_MTABLE_MMD_THRESHOLD) {
    uint64_t matching_level = level;
    if (hi-lo==0 && matching_level>0) --matching_level;
    *min_matching_depth = MIN(*min_matching_depth,matching_level);
  } else if (gem_expect_false(level+1 < rank_mtable->num_levels)) { // Control recursion level
    rank_mquery_t next_query;
    // Update 'A'
    next_query = *query;
    rank_mquery_add_char(rank_mtable,&next_query,ENC_DNA_CHAR_A);
    rank_mtable_builder_find_mmd(rank_mtable,level+1,&next_query,min_matching_depth);
    // Update 'C'
    next_query = *query;
    rank_mquery_add_char(rank_mtable,&next_query,ENC_DNA_CHAR_C);
    rank_mtable_builder_find_mmd(rank_mtable,level+1,&next_query,min_matching_depth);
    // Update 'G'
    next_query = *query;
    rank_mquery_add_char(rank_mtable,&next_query,ENC_DNA_CHAR_G);
    rank_mtable_builder_find_mmd(rank_mtable,level+1,&next_query,min_matching_depth);
    // Update 'T'
    next_query = *query;
    rank_mquery_add_char(rank_mtable,&next_query,ENC_DNA_CHAR_T);
    rank_mtable_builder_find_mmd(rank_mtable,level+1,&next_query,min_matching_depth);
  }
}
/*
 * Fill Rank-MTable
 */
void rank_mtable_builder_fill_ranks_lo(
    const bwt_builder_t* const bwt_builder,
    rank_mtable_t* const rank_mtable,
    uint64_t offset_lo,
    const uint64_t lo,
    const uint64_t level,
    ticker_t* const ticker) {
  // Control recursion level
  if (gem_expect_false(level < rank_mtable->num_levels)) {
    // Get level
    uint64_t* const rank_level = rank_mtable->sa_ranks_levels[level];
    // Calculate level offset
    const uint64_t next_level = level+1;
    const uint64_t level_skip = rank_mtable->level_skip[level];
    // Update 'A'
    rank_level[offset_lo] = bwt_builder_erank(bwt_builder,ENC_DNA_CHAR_A,lo);
    rank_mtable_builder_fill_ranks_lo(bwt_builder,rank_mtable,offset_lo,rank_level[offset_lo],next_level,ticker);
    // Update 'C'
    offset_lo+=level_skip;
    rank_level[offset_lo] = bwt_builder_erank(bwt_builder,ENC_DNA_CHAR_C,lo);
    rank_mtable_builder_fill_ranks_lo(bwt_builder,rank_mtable,offset_lo,rank_level[offset_lo],next_level,ticker);
    // Update 'G'
    offset_lo+=level_skip;
    rank_level[offset_lo] = bwt_builder_erank(bwt_builder,ENC_DNA_CHAR_G,lo);
    rank_mtable_builder_fill_ranks_lo(bwt_builder,rank_mtable,offset_lo,rank_level[offset_lo],next_level,ticker);
    // Update 'T'
    offset_lo+=level_skip;
    rank_level[offset_lo] = bwt_builder_erank(bwt_builder,ENC_DNA_CHAR_T,lo);
    rank_mtable_builder_fill_ranks_lo(bwt_builder,rank_mtable,offset_lo,rank_level[offset_lo],next_level,ticker);
    // Update ticker
    ticker_update(ticker,4);
  }
}
void rank_mtable_builder_fill_ranks_hi(
    const bwt_builder_t* const bwt_builder,
    rank_mtable_t* const rank_mtable,
    uint64_t offset_hi,
    const uint64_t hi,
    const uint64_t level,
    ticker_t* const ticker) {
  // Control recursion level
  if (gem_expect_false(level < rank_mtable->num_levels)) {
    // Get level
    uint64_t* const rank_level = rank_mtable->sa_ranks_levels[level];
    // Calculate level offset
    const uint64_t next_level = level+1;
    const uint64_t level_skip = rank_mtable->level_skip[level];
    // Update 'A'
    rank_level[offset_hi] = bwt_builder_erank(bwt_builder,ENC_DNA_CHAR_A,hi);
    rank_mtable_builder_fill_ranks_hi(bwt_builder,rank_mtable,offset_hi,rank_level[offset_hi],next_level,ticker);
    // Update 'C'
    offset_hi+=level_skip;
    rank_level[offset_hi] = bwt_builder_erank(bwt_builder,ENC_DNA_CHAR_C,hi);
    rank_mtable_builder_fill_ranks_hi(bwt_builder,rank_mtable,offset_hi,rank_level[offset_hi],next_level,ticker);
    // Update 'G'
    offset_hi+=level_skip;
    rank_level[offset_hi] = bwt_builder_erank(bwt_builder,ENC_DNA_CHAR_G,hi);
    rank_mtable_builder_fill_ranks_hi(bwt_builder,rank_mtable,offset_hi,rank_level[offset_hi],next_level,ticker);
    // Update 'T'
    offset_hi+=level_skip;
    rank_level[offset_hi] = bwt_builder_erank(bwt_builder,ENC_DNA_CHAR_T,hi);
    rank_mtable_builder_fill_ranks_hi(bwt_builder,rank_mtable,offset_hi,rank_level[offset_hi],next_level,ticker);
    // Update ticker
    ticker_update(ticker,4);
  }
}
void rank_mtable_builder_fill_ranks(
    const bwt_builder_t* const bwt_builder,
    rank_mtable_t* const rank_mtable,
    ticker_t* const ticker) {
  rank_mtable->sa_ranks_levels[0][1] = bwt_builder_get_length(bwt_builder); // Init_hi
  rank_mtable_builder_fill_ranks_hi(bwt_builder,rank_mtable,1,bwt_builder_get_length(bwt_builder),1,ticker);
  rank_mtable->sa_ranks_levels[0][0] = 0; // Init_lo
  rank_mtable->sa_ranks_levels[1][0] = 0; // 'A'
  rank_mtable_builder_fill_ranks_lo(bwt_builder,rank_mtable,0,0,2,ticker);
  // Find Minimum Matching Depth
  rank_mquery_t query;
  rank_mquery_new(&query);
  uint64_t min_matching_depth = RANK_MTABLE_SEARCH_DEPTH;
  rank_mtable_builder_find_mmd(rank_mtable,0,&query,&min_matching_depth);
  rank_mtable->min_matching_depth = min_matching_depth;
}
/*
 * Generate Rank-MTable
 */
rank_mtable_t* rank_mtable_builder_generate(
    const bwt_builder_t* const bwt_builder,
    const bool verbose) {
  // Alloc
  rank_mtable_t* const rank_mtable = mm_alloc(rank_mtable_t);
  // Set Meta-Inf
  rank_mtable->num_levels = RANK_MTABLE_LEVELS;
  rank_mtable->table_size = RANK_MTABLE_SIZE(RANK_MTABLE_SEARCH_DEPTH);
  // Allocate table
  rank_mtable->mm_sa_ranks = NULL;
  rank_mtable->sa_ranks_levels = mm_calloc(RANK_MTABLE_LEVELS,uint64_t*,true);
  rank_mtable->sa_ranks_levels[0] = mm_calloc(rank_mtable->table_size,uint64_t,true);
  rank_mtable->level_skip = mm_calloc(RANK_MTABLE_LEVELS,uint64_t,true);
  // Assign levels
  rank_mtable_init_levels(rank_mtable);
  // Initialize
  ticker_t ticker;
  ticker_percentage_reset(&ticker,verbose,"Building rank_mtable",rank_mtable->table_size,10,true);
  rank_mtable_builder_fill_ranks(bwt_builder,rank_mtable,&ticker);
  ticker_finish(&ticker);
  // Return
  return rank_mtable;
}
rank_mtable_t* rank_mtable_builder_new(
    const bwt_builder_t* const bwt_builder,
    const bool verbose) {
  return rank_mtable_builder_generate(bwt_builder,verbose);
}
/*
 * Delete
 */
void rank_mtable_builder_delete(rank_mtable_t* const rank_mtable) {
  mm_free(rank_mtable->sa_ranks_levels[0]);
  mm_free(rank_mtable->sa_ranks_levels);
  mm_free(rank_mtable->level_skip);
  mm_free(rank_mtable);
}
