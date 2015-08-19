/*
 * PROJECT: GEMMapper
 * FILE: rank_mtable.c
 * DATE: 01/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "rank_mtable.h"
#include "bwt.h"

#define RANK_MTABLE_SIZE(search_depth) ( (5*(POW4(search_depth)-1))/3 + 2 )  /* S := (5*(4^(n)-1))/3 + 2 */
#define RANK_MTABLE_LEVEL_SIZE(level)  (5*POW4(level-1))

/*
 * Loader/Setup
 */
GEM_INLINE void rank_mtable_init_levels(rank_mtable_t* const rank_mtable) {
  RANK_MTABLE_CHECK(rank_mtable);
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
GEM_INLINE rank_mtable_t* rank_mtable_read(fm_t* const file_manager) {
  FM_CHECK(file_manager);
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
GEM_INLINE rank_mtable_t* rank_mtable_read_mem(mm_t* const memory_manager) {
  MM_CHECK(memory_manager);
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
GEM_INLINE void rank_mtable_write(fm_t* const file_manager,rank_mtable_t* const rank_mtable) {
  FM_CHECK(file_manager);
  // Write Meta-info
  fm_write_uint64(file_manager,rank_mtable->num_levels);
  fm_write_uint64(file_manager,rank_mtable->table_size);
  fm_write_uint64(file_manager,rank_mtable->min_matching_depth);
  // Write table
  fm_skip_align_4KB(file_manager);
  fm_write_mem(file_manager,rank_mtable->sa_ranks_levels[0],rank_mtable->table_size*UINT64_SIZE);
}
GEM_INLINE void rank_mtable_delete(rank_mtable_t* const rank_mtable) {
  RANK_MTABLE_CHECK(rank_mtable);
  mm_free(rank_mtable->sa_ranks_levels);
  mm_free(rank_mtable->level_skip);
  if (rank_mtable->mm_sa_ranks!=NULL) mm_bulk_free(rank_mtable->mm_sa_ranks);
  mm_free(rank_mtable);
}
/*
 * Builder
 */
GEM_INLINE void rank_mtable_builder_fill_ranks_lo(
    const bwt_builder_t* const bwt_builder,rank_mtable_t* const rank_mtable,
    uint64_t offset_lo,const uint64_t lo,const uint64_t level,ticker_t* const ticker) {
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
GEM_INLINE void rank_mtable_builder_fill_ranks_hi(
    const bwt_builder_t* const bwt_builder,rank_mtable_t* const rank_mtable,
    uint64_t offset_hi,const uint64_t hi,const uint64_t level,ticker_t* const ticker) {
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
GEM_INLINE void rank_mtable_builder_find_mmd(
    const bwt_builder_t* const bwt_builder,rank_mtable_t* const rank_mtable,
    const uint64_t level,rank_mquery_t* const query,uint64_t* const min_matching_depth) {
  // Check number of matches
  uint64_t lo, hi;
  rank_mtable_fetch(rank_mtable,query,&lo,&hi);
  if (hi-lo <= RANK_MTABLE_MMD_THRESHOLD) {
    *min_matching_depth = MIN(*min_matching_depth,level);
  } else if (gem_expect_false(level+1 < rank_mtable->num_levels)) { // Control recursion level
    rank_mquery_t next_query;
    // Update 'A'
    next_query = *query;
    rank_mquery_add_char(rank_mtable,&next_query,ENC_DNA_CHAR_A);
    rank_mtable_builder_find_mmd(bwt_builder,rank_mtable,level+1,&next_query,min_matching_depth);
    // Update 'C'
    next_query = *query;
    rank_mquery_add_char(rank_mtable,&next_query,ENC_DNA_CHAR_C);
    rank_mtable_builder_find_mmd(bwt_builder,rank_mtable,level+1,&next_query,min_matching_depth);
    // Update 'G'
    next_query = *query;
    rank_mquery_add_char(rank_mtable,&next_query,ENC_DNA_CHAR_G);
    rank_mtable_builder_find_mmd(bwt_builder,rank_mtable,level+1,&next_query,min_matching_depth);
    // Update 'T'
    next_query = *query;
    rank_mquery_add_char(rank_mtable,&next_query,ENC_DNA_CHAR_T);
    rank_mtable_builder_find_mmd(bwt_builder,rank_mtable,level+1,&next_query,min_matching_depth);
  }
}
GEM_INLINE void rank_mtable_builder_fill_ranks(
    const bwt_builder_t* const bwt_builder,rank_mtable_t* const rank_mtable,ticker_t* const ticker) {
  rank_mtable->sa_ranks_levels[0][1] = bwt_builder_get_length(bwt_builder); // Init_hi
  rank_mtable_builder_fill_ranks_hi(bwt_builder,rank_mtable,1,bwt_builder_get_length(bwt_builder),1,ticker);
  rank_mtable->sa_ranks_levels[0][0] = 0; // Init_lo
  rank_mtable->sa_ranks_levels[1][0] = 0; // 'A'
  rank_mtable_builder_fill_ranks_lo(bwt_builder,rank_mtable,0,0,2,ticker);
  // Find Minimum Matching Depth
  rank_mquery_t query;
  rank_mquery_new(&query);
  uint64_t min_matching_depth = UINT64_MAX;
  rank_mtable_builder_find_mmd(bwt_builder,rank_mtable,0,&query,&min_matching_depth);
  rank_mtable->min_matching_depth = min_matching_depth;
}
GEM_INLINE rank_mtable_t* rank_mtable_builder_new(const bwt_builder_t* const bwt_builder,const bool verbose) {
  // Alloc
  rank_mtable_t* const rank_mtable = mm_alloc(rank_mtable_t);
  // Set Meta-Inf
  rank_mtable->num_levels=RANK_MTABLE_LEVELS;
  rank_mtable->table_size=RANK_MTABLE_SIZE(RANK_MTABLE_SEARCH_DEPTH);
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
GEM_INLINE void rank_mtable_builder_delete(rank_mtable_t* const rank_mtable) {
  RANK_MTABLE_CHECK(rank_mtable);
  mm_free(rank_mtable->sa_ranks_levels[0]);
  mm_free(rank_mtable->sa_ranks_levels);
  mm_free(rank_mtable->level_skip);
  mm_free(rank_mtable);
}
/*
 * Accessors
 */
GEM_INLINE uint64_t rank_mtable_get_size(const rank_mtable_t* const rank_mtable) {
  RANK_MTABLE_CHECK(rank_mtable);
  return rank_mtable->table_size * UINT64_SIZE;
}
/*
 * Query
 */
GEM_INLINE void rank_mquery_new(rank_mquery_t* const query) {
  RANK_MQUERY_CHECK(query);
  query->hi_position = 1;
  query->level = 0;
}
GEM_INLINE void rank_mquery_add_char(const rank_mtable_t* const rank_mtable,rank_mquery_t* const query,uint8_t const enc_char) {
  RANK_MQUERY_CHECK(query);
  // Update HI => hi(n+1) = hi(n) + c*4^(level)
  ++(query->level);
  query->hi_position = query->hi_position + enc_char*rank_mtable->level_skip[query->level];
}
GEM_INLINE uint64_t rank_mquery_get_level(const rank_mquery_t* const query) {
  RANK_MQUERY_CHECK(query);
  return query->level;
}
GEM_INLINE uint64_t rank_mquery_is_exhausted(const rank_mquery_t* const query) {
  RANK_MQUERY_CHECK(query);
  return query->level >= RANK_MTABLE_SEARCH_DEPTH;
}
/*
 * Fetch rank value
 */
GEM_INLINE void rank_mtable_fetch(
    const rank_mtable_t* const rank_mtable,const rank_mquery_t* const query,
    uint64_t* const lo,uint64_t* const hi) {
  RANK_MTABLE_CHECK(rank_mtable);
  RANK_MQUERY_CHECK(query);
  *hi = rank_mtable->sa_ranks_levels[query->level][query->hi_position];
  *lo = rank_mtable->sa_ranks_levels[query->level][query->hi_position-1];
}
/*
 * Display
 */
GEM_INLINE void rank_mtable_print(FILE* const stream,rank_mtable_t* const rank_mtable) {
  tab_fprintf(stream,"[GEM]>Rank.Table\n");
  tab_fprintf(stream,"  => Total.Cells        %"PRIu64"\n",rank_mtable->table_size);
  tab_fprintf(stream,"  => Total.Size         %"PRIu64" MB\n",CONVERT_B_TO_MB(rank_mtable->table_size*UINT64_SIZE));
  tab_fprintf(stream,"  => Num.Levels         %"PRIu64"\n",rank_mtable->num_levels);
  tab_fprintf(stream,"  => Min.matching.depth %"PRIu64"\n",rank_mtable->min_matching_depth);
  // Flush
  fflush(stream);
}
GEM_INLINE void rank_mtable_print_content(FILE* const stream,rank_mtable_t* const rank_mtable,const uint64_t text_length) {
  uint64_t i;
  for (i=0;i<rank_mtable->table_size;++i) {
    if (rank_mtable->sa_ranks_levels[0][i] > text_length) {
      printf("Here\n");
    }
    fprintf(stream,"%"PRIu64"\n",rank_mtable->sa_ranks_levels[0][i]);
  }
}
