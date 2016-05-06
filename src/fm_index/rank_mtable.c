/*
 * PROJECT: GEMMapper
 * FILE: rank_mtable.c
 * DATE: 01/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "fm_index/rank_mtable.h"
#include "fm_index/bwt.h"

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
void rank_mtable_print(
    FILE* const stream,
    rank_mtable_t* const rank_mtable) {
  tab_fprintf(stream,"[GEM]>Rank.Table\n");
  tab_fprintf(stream,"  => Total.Cells        %"PRIu64"\n",rank_mtable->table_size);
  tab_fprintf(stream,"  => Total.Size         %"PRIu64" MB\n",CONVERT_B_TO_MB(rank_mtable->table_size*UINT64_SIZE));
  tab_fprintf(stream,"  => Num.Levels         %"PRIu64"\n",rank_mtable->num_levels);
  tab_fprintf(stream,"  => Min.matching.depth %"PRIu64"\n",rank_mtable->min_matching_depth);
  // Flush
  fflush(stream);
}
void rank_mtable_print_content(
    FILE* const stream,
    rank_mtable_t* const rank_mtable,
    const uint64_t text_length) {
  uint64_t i;
  for (i=0;i<rank_mtable->table_size;++i) {
    if (rank_mtable->sa_ranks_levels[0][i] > text_length) {
      printf("Here\n");
    }
    fprintf(stream,"%"PRIu64"\n",rank_mtable->sa_ranks_levels[0][i]);
  }
}
