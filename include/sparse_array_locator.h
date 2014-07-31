/*
 * PROJECT: GEMMapper
 * FILE: sparse_array_locator.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef SPARSE_ARRAY_LOCATOR_H_
#define SPARSE_ARRAY_LOCATOR_H_

#include "essentials.h"
#include "stats_vector.h"
#include "segmented_vector.h"

/*
 * Checkers
 */
#define SPARSE_ARRAY_LOCATOR_CHECK(sparse_array_locator) \
  GEM_CHECK_NULL(sparse_array_locator)

#define SPARSE_ARRAY_LOCATOR_CHECK_POSITION(sparse_array_locator,position) \
  gem_check(position>=locator->total_length,SPARSE_ARRAY_LOCATOR_INDEX,position,locator->total_length)

#define SPARSE_ARRAY_LOCATOR_BUILDER_CHECK(locator_stats) \
  GEM_CHECK_NULL(locator_stats)

#define SPARSE_ARRAY_LOCATOR_STATS_CHECK(locator_stats) \
  GEM_CHECK_NULL(locator_stats)

/*
 * Sparse Array Locator
 */
typedef struct {
  /* Meta-Data */
  uint64_t total_length;     // Total length of the bitmap
  uint64_t total_size;       // Total Size (Bytes)
  uint64_t num_mayor_blocks; // Total Mayor Blocks
  uint64_t num_minor_blocks; // Total Minor Blocks
  /* Data Stored */
  uint64_t* mayor_counters;  // MAYOR_COUNTERS := [64MayorCounter]...[64MayorCounter]
  uint64_t* bitmap;          // MINOR_BLOCKS   := [48Bitmap][16MinorCounter]...[48Bitmap][16MinorCounter]
  /* MM */
  mm_t* mm;
  /* Space Range [Idx_start,Idx_end) */
  uint64_t first_block_offset; // First block offset
  uint64_t idx_offset;         // Global offset
} sparse_array_locator_t;
/*
 * Sparse Array Locator Builder
 */
typedef struct {
  /* Position & current Block */
  uint64_t position;
  uint64_t current_block;
  uint64_t current_count;
  /* Counting */
  uint64_t minor_counter_accumulated;
  uint64_t mayor_counter_accumulated;
  /* Memory */
  svector_t* minor_blocks;
  svector_iterator_t minor_blocks_iterator;
  vector_t* mayor_counters;
} sparse_array_locator_builder_t;
/*
 * Sparse Array Locator Stats
 */
typedef struct {
  /* Locator Stats */
  uint64_t locator_length;
  uint64_t marked_positions;
  /* Bitmaps */
  stats_vector_t* ones_density;
} sparse_array_locator_stats_t;

/*
 * Loader/Setup
 */
GEM_INLINE sparse_array_locator_t* sparse_array_locator_read(fm_t* const file_manager);
GEM_INLINE sparse_array_locator_t* sparse_array_locator_read_mem(mm_t* const memory_manager);
GEM_INLINE void sparse_array_locator_delete(sparse_array_locator_t* const locator);

/*
 * Accessors
 */
GEM_INLINE uint64_t sparse_array_locator_get_size(sparse_array_locator_t* const locator);

GEM_INLINE bool sparse_array_locator_is_marked(sparse_array_locator_t* const locator,const uint64_t position);
GEM_INLINE uint64_t sparse_array_locator_get_erank(sparse_array_locator_t* const locator,const uint64_t position);
GEM_INLINE bool sparse_array_locator_get_erank_if_marked(
    sparse_array_locator_t* const locator,const uint64_t position,uint64_t* const erank);
GEM_INLINE bool sparse_array_locator_get_erank__marked(
    sparse_array_locator_t* const locator,const uint64_t position,uint64_t* const erank);

/*
 * Builder
 */
// Static builder
GEM_INLINE sparse_array_locator_t* sparse_array_locator_new(const uint64_t idx_begin,const uint64_t idx_end);
GEM_INLINE void sparse_array_locator_mark(sparse_array_locator_t* const locator,const uint64_t position);
GEM_INLINE void sparse_array_locator_write(fm_t* const file_manager,sparse_array_locator_t* const locator);
GEM_INLINE void sparse_array_locator_merge__write(
    fm_t* const file_manager,sparse_array_locator_t** const locator,const uint64_t num_locators);
// Dynamic Builder
GEM_INLINE sparse_array_locator_builder_t* sparse_array_locator_builder_new(mm_slab_t* const mm_slab);
GEM_INLINE void sparse_array_locator_builder_delete(sparse_array_locator_builder_t* const locator_builder);
GEM_INLINE void sparse_array_locator_builder_next(
    sparse_array_locator_builder_t* const locator_builder,const bool mark_position);
GEM_INLINE void sparse_array_locator_builder_write(
    fm_t* const file_manager,sparse_array_locator_builder_t* const locator_builder);

/*
 * Display
 */
GEM_INLINE void sparse_array_locator_print(
    FILE* const stream,const sparse_array_locator_t* const locator,const bool display_content);

/*
 * Stats // TODO Enhance
 */
GEM_INLINE sparse_array_locator_stats_t* sparse_array_locator_stats_new();
GEM_INLINE void sparse_array_locator_stats_delete(sparse_array_locator_stats_t* const sparse_array_locator_stats);
GEM_INLINE void sparse_array_locator_stats_calculate(
    sparse_array_locator_stats_t* const sparse_array_locator_stats,sparse_array_locator_t* const locator);
GEM_INLINE void sparse_array_locator_stats_print(
    FILE* const stream,const char* const sparse_array_locator_stats_tag,
    sparse_array_locator_stats_t* const locator_stats);

/*
 * Error Messages
 */
#define GEM_ERROR_SPARSE_ARRAY_LOCATOR_INDEX "SparseArrayLocator::Index (%"PRIu64") out of bounds [0,%"PRIu64")"

#endif /* SPARSE_ARRAY_LOCATOR_H_ */
