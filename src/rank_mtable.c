/*
 * PROJECT: GEMMapper
 * FILE: rank_mtable.c
 * DATE: 01/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "rank_mtable.h"
#include "bwt.h"
#include "dna_string.h"

#define RANK_MTABLE_SIZE(search_depth) ( (POW4(search_depth+1)-1)/3 )  /* S := (4^(n+1)-1)/3 */
#define RANK_MTABLE_LEVEL_SIZE(level) POW4(level)

/*
 * Loader/Setup
 */
GEM_INLINE void rank_mtable_init_levels(rank_mtable_t* const rank_mtable) {
  RANK_MTABLE_CHECK(rank_mtable);
  uint64_t i, level_size = 1;
  for (i=1;i<rank_mtable->num_levels;++i) {
    rank_mtable->sa_ranks_levels[i] = rank_mtable->sa_ranks_levels[i-1] + level_size;
    level_size = level_size*DNA_RANGE;
  }
}
GEM_INLINE rank_mtable_t* rank_mtable_read(fm_t* const file_manager) {
  FM_CHECK(file_manager);
  // Alloc
  rank_mtable_t* const rank_mtable = mm_alloc(rank_mtable_t);
  // Read Meta-info
  rank_mtable->num_levels = fm_read_uint64(file_manager);
  rank_mtable->table_size = fm_read_uint64(file_manager);
  // Read intervals
  fm_skip_align_4KB(file_manager);
  rank_mtable->mm_sa_ranks = fm_load_mem(file_manager,rank_mtable->table_size*UINT64_SIZE);
  rank_mtable->sa_ranks_levels = mm_calloc(rank_mtable->num_levels,uint64_t*,true);
  rank_mtable->sa_ranks_levels[0] = mm_get_base_mem(rank_mtable->mm_sa_ranks);
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
  // Read intervals
  mm_skip_align_4KB(memory_manager);
  rank_mtable->mm_sa_ranks = NULL;
  rank_mtable->sa_ranks_levels = mm_calloc(rank_mtable->num_levels,uint64_t*,true);
  rank_mtable->sa_ranks_levels[0] = mm_read_mem(memory_manager,rank_mtable->table_size*UINT64_SIZE);
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
  // Write table
  fm_skip_align_4KB(file_manager);
  fm_write_mem(file_manager,rank_mtable->sa_ranks_levels[0],rank_mtable->table_size*UINT64_SIZE);
}
GEM_INLINE void rank_mtable_delete(rank_mtable_t* const rank_mtable) {
  RANK_MTABLE_CHECK(rank_mtable);
  mm_free(rank_mtable->sa_ranks_levels);
  if (rank_mtable->mm_sa_ranks!=NULL) mm_bulk_free(rank_mtable->mm_sa_ranks);
}
/*
 * Builder
 */
GEM_INLINE void rank_mtable_builder_fill_ranks(
    const bwt_builder_t* const bwt_builder,rank_mtable_t* const rank_mtable,
  uint64_t offset,const uint64_t hi,const uint64_t level,ticker_t* const ticker) {
  // Control recursion level
  if (gem_expect_false(level < rank_mtable->num_levels)) {
    // Get level
    uint64_t* const rank_level = rank_mtable->sa_ranks_levels[level];
    // Calculate level offset
    const uint64_t next_level = level+1;
    const uint64_t level_skip = POW4(level-1); /* 4^(level) */
    // Update 'A'
    rank_level[offset] = bwt_builder_erank(bwt_builder,ENC_DNA_CHAR_A,hi);
    rank_mtable_builder_fill_ranks(bwt_builder,rank_mtable,offset,rank_level[offset],next_level,ticker);
    // Update 'C'
    offset+=level_skip;
    rank_level[offset] = bwt_builder_erank(bwt_builder,ENC_DNA_CHAR_C,hi);
    rank_mtable_builder_fill_ranks(bwt_builder,rank_mtable,offset,rank_level[offset],next_level,ticker);
    // Update 'G'
    offset+=level_skip;
    rank_level[offset] = bwt_builder_erank(bwt_builder,ENC_DNA_CHAR_G,hi);
    rank_mtable_builder_fill_ranks(bwt_builder,rank_mtable,offset,rank_level[offset],next_level,ticker);
    // Update 'T'
    offset+=level_skip;
    rank_level[offset] = bwt_builder_erank(bwt_builder,ENC_DNA_CHAR_T,hi);
    rank_mtable_builder_fill_ranks(bwt_builder,rank_mtable,offset,rank_level[offset],next_level,ticker);
    // Update ticker
    ticker_update(ticker,4);
  }
}
GEM_INLINE rank_mtable_t* rank_mtable_builder_new(const bwt_builder_t* const bwt_builder,const bool verbose) {
  // Alloc
  rank_mtable_t* const rank_mtable = mm_alloc(rank_mtable_t);
  // Set Meta-Inf
  rank_mtable->num_levels=RANK_MTABLE_LEVELS__;
  rank_mtable->table_size=RANK_MTABLE_SIZE(RANK_MTABLE_SEARCH_DEPTH);
  // Allocate table
  rank_mtable->mm_sa_ranks = NULL;
  rank_mtable->sa_ranks_levels = mm_calloc(RANK_MTABLE_LEVELS__,uint64_t*,true);
  rank_mtable->sa_ranks_levels[0] = mm_calloc(rank_mtable->table_size,uint64_t,true);
  // Assign levels
  rank_mtable_init_levels(rank_mtable);
  // Initialize
  ticker_t ticker;
  ticker_percentage_reset(&ticker,verbose,"Building rank_mtable",rank_mtable->table_size-1,10,true);
  // rank_mtable->sa_ranks_levels[0][0] = 0; [Implicit]
  rank_mtable_builder_fill_ranks(bwt_builder,rank_mtable,0,bwt_builder_get_length(bwt_builder),1,&ticker);
  // Return
  return rank_mtable;
}
GEM_INLINE void rank_mtable_builder_delete(rank_mtable_t* const rank_mtable) {
  RANK_MTABLE_CHECK(rank_mtable);
  mm_free(rank_mtable->sa_ranks_levels[0]);
  mm_free(rank_mtable->sa_ranks_levels);
  mm_free(rank_mtable);
}
/*
 * Accessors
 */
GEM_INLINE uint64_t rank_mtable_get_size(rank_mtable_t* const rank_mtable) {
  RANK_MTABLE_CHECK(rank_mtable);
  return rank_mtable->table_size * UINT64_SIZE;
}
/*
 * Query
 */
GEM_INLINE void rank_mquery_new(rank_mquery_t* const query) {
  RANK_MQUERY_CHECK(query);
  query->hi_position = 0;
  query->lo_position = 0;
  query->hi_level = 0;
  query->lo_level = 0;
}
GEM_INLINE void rank_mquery_add_char(rank_mquery_t* const query,uint8_t const enc_char) {
  RANK_MQUERY_CHECK(query);
  // Update HI => hi(n+1) = hi(n) + c*4^(level)
  query->hi_position = query->hi_position + enc_char*RANK_MTABLE_LEVEL_SIZE(query->hi_level);
  ++(query->hi_level);
  // Update LO
  if (enc_char != ENC_DNA_CHAR_A) {
    query->lo_position = query->hi_position-1;
    query->lo_level = query->hi_level;
  }
}
GEM_INLINE uint64_t rank_mquery_get_level(rank_mquery_t* const query) {
  RANK_MQUERY_CHECK(query);
  return query->hi_level;
}
GEM_INLINE uint64_t rank_mquery_is_exhausted(rank_mquery_t* const query) {
  RANK_MQUERY_CHECK(query);
  return query->hi_level >= RANK_MTABLE_SEARCH_DEPTH;
}
/*
 * Fetch rank value
 */
GEM_INLINE void rank_mtable_fetch(
    rank_mtable_t* const rank_mtable,rank_mquery_t* const query,
    uint64_t* const lo,uint64_t* const hi) {
  RANK_MTABLE_CHECK(rank_mtable);
  RANK_MQUERY_CHECK(query);
  // RANK_MTABLE_QUERY_CHECK(rank_table,); // TODO
  *hi = rank_mtable->sa_ranks_levels[query->hi_level][query->hi_position];
  *lo = rank_mtable->sa_ranks_levels[query->lo_level][query->lo_position];
}
/*
 * Display
 */
GEM_INLINE void rank_mtable_print_cell_label(FILE* const stream,const uint64_t level,uint64_t position) {
  GEM_CHECK_NULL(stream);
  char label[level];
  int64_t i;
  // Decode position <-> label (base 4)
  for (i=0;i<level;++i) {
    label[i] = dna_decode(position % DNA_RANGE);
    position = position/4;
  }
  // Print label
  for (i=level-1;i>=0;--i) {
    fprintf(stream,"%c",label[i]);
  }
}
GEM_INLINE void rank_mtable_print_level(FILE* const stream,rank_mtable_t* const rank_mtable,const uint64_t level) {
  GEM_CHECK_NULL(stream);
  RANK_MTABLE_CHECK(rank_mtable);
  const uint64_t level_size = RANK_MTABLE_LEVEL_SIZE(level);
  tab_fprintf(stream,"  => Rank.Table.Level.%lu (%lu cells)\n",level,level_size);
  uint64_t i;
  tab_global_inc(); tab_global_inc();
  tab_fprintf(stream,"");
  for (i=0;i<level_size;++i) {
    if (i && i%4==0) { fprintf(stream,"\n"); tab_fprintf(stream,""); }
    fprintf(stream,"[");
    rank_mtable_print_cell_label(stream,level,i);
    fprintf(stream,"]=%lu",rank_mtable->sa_ranks_levels[level][i]);
  }
  fprintf(stream,"\n");
  tab_global_dec(); tab_global_dec();
}
GEM_INLINE void rank_mtable_print(FILE* const stream,rank_mtable_t* const rank_mtable) {
  tab_fprintf(stream,"[GEM]>Rank.Table\n");
  tab_fprintf(stream,"  => Total.Cells %lu\n",rank_mtable->table_size);
  tab_fprintf(stream,"  => Total.Size %lu MB\n",CONVERT_B_TO_MB(rank_mtable->table_size*UINT64_SIZE));
  tab_fprintf(stream,"  => Num.Levels %lu\n",rank_mtable->num_levels);
  rank_mtable_print_level(stream,rank_mtable,1);
  rank_mtable_print_level(stream,rank_mtable,2);
  // rank_mtable_print_level(stream,rank_mtable,3);
  // Flush
  fflush(stream);
}

