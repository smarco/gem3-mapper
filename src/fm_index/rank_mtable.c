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
 *   FM-Index module provides a rank-memoization table that stores
 *   pre-computed bwt-intevals for all the k-mer (k=11) words from {ACGT}
 */

#include "fm_index/rank_mtable.h"
#include "fm_index/bwt/bwt.h"

/*
 * Setup
 */
void rank_mtable_init_levels(rank_mtable_t* const rank_mtable) {
  uint64_t* const level_skip = rank_mtable->level_skip;
  uint64_t i;
  level_skip[0] = 0; // Null
  level_skip[1] = 1; // Base skip
  rank_mtable->sa_ranks_levels[1] = rank_mtable->sa_ranks_levels[0] + 2;
  for (i=2;i<rank_mtable->num_levels;++i) {
    level_skip[i] = RANK_MTABLE_LEVEL_SIZE(i-1);
    rank_mtable->sa_ranks_levels[i] = rank_mtable->sa_ranks_levels[i-1] + level_skip[i];
  }
}
rank_mtable_t* rank_mtable_read(fm_t* const file_manager) {
  // Alloc
  rank_mtable_t* const rank_mtable = mm_alloc(rank_mtable_t);
  // Read Meta-info
  rank_mtable->num_levels = fm_read_uint64(file_manager);
  rank_mtable->table_size = fm_read_uint64(file_manager);
  rank_mtable->min_matching_depth = fm_read_uint64(file_manager);
  // Read ranks
  fm_skip_align_4KB(file_manager);
  rank_mtable->mm_sa_ranks = fm_load_mem(file_manager,rank_mtable->table_size*UINT64_SIZE);
  rank_mtable->sa_ranks_levels = mm_calloc(rank_mtable->num_levels,uint64_t*,true);
  rank_mtable->sa_ranks_levels[0] = mm_get_base_mem(rank_mtable->mm_sa_ranks);
  rank_mtable->level_skip = mm_calloc(rank_mtable->num_levels,uint64_t,true);
  // Assign levels
  rank_mtable_init_levels(rank_mtable);
  // Return
  return rank_mtable;
}
rank_mtable_t* rank_mtable_read_mem(mm_t* const memory_manager) {
  // Alloc
  rank_mtable_t* const rank_mtable = mm_alloc(rank_mtable_t);
  // Read Meta-info
  rank_mtable->num_levels = mm_read_uint64(memory_manager);
  rank_mtable->table_size = mm_read_uint64(memory_manager);
  rank_mtable->min_matching_depth = mm_read_uint64(memory_manager);
  // Read intervals
  mm_skip_align_4KB(memory_manager);
  rank_mtable->mm_sa_ranks = NULL;
  rank_mtable->sa_ranks_levels = mm_calloc(rank_mtable->num_levels,uint64_t*,true);
  rank_mtable->sa_ranks_levels[0] = mm_read_mem(memory_manager,rank_mtable->table_size*UINT64_SIZE);
  rank_mtable->level_skip = mm_calloc(rank_mtable->num_levels,uint64_t,true);
  // Assign levels
  rank_mtable_init_levels(rank_mtable);
  // Return
  return rank_mtable;
}
void rank_mtable_delete(rank_mtable_t* const rank_mtable) {
  mm_free(rank_mtable->sa_ranks_levels);
  mm_free(rank_mtable->level_skip);
  if (rank_mtable->mm_sa_ranks!=NULL) mm_bulk_free(rank_mtable->mm_sa_ranks);
  mm_free(rank_mtable);
}
/*
 * Accessors
 */
uint64_t rank_mtable_get_size(const rank_mtable_t* const rank_mtable) {
  return rank_mtable->table_size * UINT64_SIZE;
}
/*
 * Query
 */
void rank_mquery_new(rank_mquery_t* const query) {
  query->hi_position = 1;
  query->level = 0;
}
void rank_mquery_add_char(
    const rank_mtable_t* const rank_mtable,
    rank_mquery_t* const query,
    uint8_t const enc_char) {
  // Update HI => hi(n+1) = hi(n) + c*4^(level)
  ++(query->level);
  query->hi_position = query->hi_position + enc_char*rank_mtable->level_skip[query->level];
}
uint64_t rank_mquery_get_level(const rank_mquery_t* const query) {
  return query->level;
}
uint64_t rank_mquery_is_exhausted(const rank_mquery_t* const query) {
  return query->level >= RANK_MTABLE_SEARCH_DEPTH;
}
/*
 * Fetch rank value
 */
void rank_mtable_prefetch(
    const rank_mtable_t* const rank_mtable,
    const rank_mquery_t* const query) {
  PREFETCH((rank_mtable->sa_ranks_levels[query->level] + (query->hi_position-1)));
}
void rank_mtable_fetch(
    const rank_mtable_t* const rank_mtable,
    const rank_mquery_t* const query,
    uint64_t* const lo,
    uint64_t* const hi) {
  *hi = rank_mtable->sa_ranks_levels[query->level][query->hi_position];
  *lo = rank_mtable->sa_ranks_levels[query->level][query->hi_position-1];
}
/*
 * Display
 */
void rank_mtable_print_content_rec(
    FILE* const stream,
    rank_mtable_t* const rank_mtable,
    const uint64_t current_offset,
    const uint64_t current_level,
    const uint64_t max_level,
    char* const string) {
  // Control recursion level
  if (gem_expect_false(current_level < max_level)) {
    // Get level & Calculate level offset
    uint64_t* const rank_level = rank_mtable->sa_ranks_levels[current_level];
    const uint64_t next_level = current_level+1;
    const uint64_t level_skip = rank_mtable->level_skip[current_level];
    char* const string_level = string + (RANK_MTABLE_SEARCH_DEPTH-current_level);
    // Update 'A'
    const uint64_t offset_A = current_offset;
    *string_level = 'A';
    tab_fprintf(stream,"[%s] %lu %lu\n",string_level,rank_level[offset_A],rank_level[offset_A+1]);
    rank_mtable_print_content_rec(stream,rank_mtable,offset_A,next_level,max_level,string);
    // Update 'C'
    const uint64_t offset_C = offset_A + level_skip;
    *string_level = 'C';
    tab_fprintf(stream,"[%s] %lu %lu\n",string_level,rank_level[offset_C],rank_level[offset_C+1]);
    rank_mtable_print_content_rec(stream,rank_mtable,offset_C,next_level,max_level,string);
    // Update 'G'
    const uint64_t offset_G = offset_C + level_skip;
    *string_level = 'G';
    tab_fprintf(stream,"[%s] %lu %lu\n",string_level,rank_level[offset_G],rank_level[offset_G+1]);
    rank_mtable_print_content_rec(stream,rank_mtable,offset_G,next_level,max_level,string);
    // Update 'T'
    const uint64_t offset_T = offset_G + level_skip;
    *string_level = 'T';
    tab_fprintf(stream,"[%s] %lu %lu\n",string_level,rank_level[offset_T],rank_level[offset_T+1]);
    rank_mtable_print_content_rec(stream,rank_mtable,offset_T,next_level,max_level,string);
  }
}
void rank_mtable_print(
    FILE* const stream,
    rank_mtable_t* const rank_mtable,
    const bool print_content) {
  tab_fprintf(stream,"[GEM]>Rank.Table\n");
  tab_fprintf(stream,"  => Total.Cells        %"PRIu64"\n",rank_mtable->table_size);
  tab_fprintf(stream,"  => Total.Size         %"PRIu64" MB\n",CONVERT_B_TO_MB(rank_mtable->table_size*UINT64_SIZE));
  tab_fprintf(stream,"  => Num.Levels         %"PRIu64"\n",rank_mtable->num_levels);
  tab_fprintf(stream,"  => Min.matching.depth %"PRIu64"\n",rank_mtable->min_matching_depth);
  // Print content
  if (print_content) {
    tab_fprintf(stream,"  => Content\n");
    char string[RANK_MTABLE_SEARCH_DEPTH+1];
    const uint64_t max_level = 4; // MAX==rank_mtable->num_levels;
    tab_global_inc(); tab_global_inc();
    string[RANK_MTABLE_SEARCH_DEPTH] = '\0';
    rank_mtable_print_content_rec(stream,rank_mtable,0,1,max_level,string);
    tab_global_dec(); tab_global_dec();
  }
  // Flush
  fflush(stream);
}

