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
 */

#include "utils/sparse_array_locator.h"

/*
 * Error Messages
 */
#define GEM_ERROR_SPARSE_ARRAY_LOCATOR_INDEX "SparseArrayLocator::Index (%"PRIu64") out of bounds [0,%"PRIu64")"

/*
 * Constants
 */
#define SAL_MINOR_BLOCK_SIZE          UINT64_SIZE   /* 64 */
#define SAL_MINOR_COUNTER_LENGTH      UINT16_LENGTH /* 16 */
#define SAL_MINOR_BLOCK_BITMAP_MASK   0x0000FFFFFFFFFFFFull
#define SAL_MINOR_BLOCK_ONE_LAST_MASK 0x0000800000000000ull
#define SAL_MINOR_BLOCK_LENGTH        (UINT64_LENGTH-SAL_MINOR_COUNTER_LENGTH) /* 64-16 = 48*/

#define SAL_MAYOR_COUNTER_SIZE UINT64_SIZE
#define SAL_MAYOR_BLOCK_LENGTH (1<<SAL_MINOR_COUNTER_LENGTH) /* ((2^16)-1)+1 = 65536 */
#define SAL_MAYOR_BLOCK_NUM_MINOR_BLOCKS (SAL_MAYOR_BLOCK_LENGTH/SAL_MINOR_BLOCK_LENGTH) /* 65536/48 = 1365 */

/*
 * Calculate SparseArrayLocator Dimensions (Number of Minor/Mayor Blocks)
 */
#define SAL_CALCULATE_DIMENSIONS(total_length,num_minor_blocks,num_mayor_blocks) \
  const uint64_t num_minor_blocks = (total_length+(SAL_MINOR_BLOCK_LENGTH-1))/SAL_MINOR_BLOCK_LENGTH; \
  const uint64_t num_mayor_blocks = (total_length+(SAL_MAYOR_BLOCK_LENGTH-1))/SAL_MAYOR_BLOCK_LENGTH

/*
 * Loader/Setup
 */
sparse_array_locator_t* sparse_array_locator_read(fm_t* const file_manager) {
  // Allocate handler
  sparse_array_locator_t* const locator = mm_alloc(sparse_array_locator_t);
  // Meta-Data
  locator->total_length = fm_read_uint64(file_manager);
  locator->total_size = fm_read_uint64(file_manager);
  locator->num_mayor_blocks = fm_read_uint64(file_manager);
  locator->num_minor_blocks = fm_read_uint64(file_manager);
  // Read locator
  locator->mm = fm_load_mem(file_manager,locator->total_size);
  locator->bitmap = mm_read_mem(locator->mm,locator->num_minor_blocks*SAL_MINOR_BLOCK_SIZE);
  locator->mayor_counters = mm_read_mem(locator->mm,locator->num_mayor_blocks*SAL_MAYOR_COUNTER_SIZE);
  // Return
  return locator;
}
sparse_array_locator_t* sparse_array_locator_read_mem(mm_t* const memory_manager) {
  // Allocate handler
  sparse_array_locator_t* const locator = mm_alloc(sparse_array_locator_t);
  // Meta-Data
  locator->total_length = mm_read_uint64(memory_manager);
  locator->total_size = mm_read_uint64(memory_manager);
  locator->num_mayor_blocks = mm_read_uint64(memory_manager);
  locator->num_minor_blocks = mm_read_uint64(memory_manager);
  // Read locator
  locator->mm = NULL;
  locator->bitmap = mm_read_mem(memory_manager,locator->num_minor_blocks*SAL_MINOR_BLOCK_SIZE);
  locator->mayor_counters = mm_read_mem(memory_manager,locator->num_mayor_blocks*SAL_MAYOR_COUNTER_SIZE);
  // Return
  return locator;
}
void sparse_array_locator_delete(sparse_array_locator_t* const locator) {
  if (locator->mm!=NULL) {
    if (locator->mm==(mm_t*)(-1)) {
      mm_free(locator->mayor_counters);
      mm_free(locator->bitmap);
    } else {
      mm_bulk_free(locator->mm);
    }
  }
  mm_free(locator);
}
/*
 * Accessors
 */
uint64_t sparse_array_locator_get_size(sparse_array_locator_t* const locator) {
  return locator->total_size;
}
#define SPARSE_ARRAY_LOCATOR_GET_MINOR_BLOCK_LOCATION(position,block_pos,block_mod) \
  /* Locate Minor Block */ \
  const uint64_t block_pos = position / SAL_MINOR_BLOCK_LENGTH; \
  const uint64_t block_mod = position % SAL_MINOR_BLOCK_LENGTH
#define SPARSE_ARRAY_LOCATOR_GET_MAYOR_BLOCK_LOCATION(position,mayor_block_pos) \
  /* Locate Mayor Block */ \
  const uint64_t mayor_block_pos = position / SAL_MAYOR_BLOCK_LENGTH
bool sparse_array_locator_is_marked(
    sparse_array_locator_t* const locator,
    const uint64_t position) {
  // Locate Block
  SPARSE_ARRAY_LOCATOR_GET_MINOR_BLOCK_LOCATION(position,block_pos,block_mod);
  return (locator->bitmap[block_pos] & (UINT64_ONE_MASK << block_mod));
}
uint64_t sparse_array_locator_get_erank(
    sparse_array_locator_t* const locator,
    const uint64_t position) {
  SPARSE_ARRAY_LOCATOR_GET_MINOR_BLOCK_LOCATION(position,block_pos,block_mod);
  const uint64_t* const block = locator->bitmap + block_pos;
  const uint16_t* const block_counter = ((const uint16_t* const)block) + 3;
  const uint64_t block_masked = *block & uint64_mask_ones[block_mod];
  SPARSE_ARRAY_LOCATOR_GET_MAYOR_BLOCK_LOCATION(position,mayor_block);
  return locator->mayor_counters[mayor_block] +
      ((gem_expect_true(block_masked)) ? POPCOUNT_64(block_masked) + *block_counter : *block_counter);
}
bool sparse_array_locator_get_erank_if_marked(
    sparse_array_locator_t* const locator,
    const uint64_t position,
    uint64_t* const erank) {
  SPARSE_ARRAY_LOCATOR_GET_MINOR_BLOCK_LOCATION(position,block_pos,block_mod);
  const uint64_t* const block_addr = locator->bitmap + block_pos;
  const uint64_t block = *block_addr;
  // Check if marked
  if ((block & (UINT64_ONE_MASK << block_mod))) {
    const uint16_t* const block_counter = ((const uint16_t* const)block_addr) + 3;
    const uint64_t block_masked = block & uint64_erank_mask(block_mod);
    SPARSE_ARRAY_LOCATOR_GET_MAYOR_BLOCK_LOCATION(position,mayor_block);
    *erank = locator->mayor_counters[mayor_block] +
        ((gem_expect_true(block_masked)) ? POPCOUNT_64(block_masked) + *block_counter : *block_counter);
    return true;
  } else {
    return false;
  }
}
bool sparse_array_locator_get_erank__marked(
    sparse_array_locator_t* const locator,
    const uint64_t position,
    uint64_t* const erank) {
  SPARSE_ARRAY_LOCATOR_GET_MINOR_BLOCK_LOCATION(position,block_pos,block_mod);
  const uint64_t* const block_addr = locator->bitmap + block_pos;
  const uint64_t block = *block_addr;
  // Check if marked
  const uint16_t* const block_counter = ((const uint16_t* const)block_addr) + 3;
  const uint64_t block_masked = block & uint64_erank_mask(block_mod);
  SPARSE_ARRAY_LOCATOR_GET_MAYOR_BLOCK_LOCATION(position,mayor_block);
  *erank = locator->mayor_counters[mayor_block] +
      ((gem_expect_true(block_masked)) ? POPCOUNT_64(block_masked) + *block_counter : *block_counter);
  return block & (UINT64_ONE_MASK << block_mod);
}
/*
 * Static Builder
 */
sparse_array_locator_t* sparse_array_locator_new(
    const uint64_t idx_begin,
    const uint64_t idx_end) {
  // Allocate handler
  sparse_array_locator_t* const locator = mm_alloc(sparse_array_locator_t);
  // Meta-Data
  const uint64_t first_block_offset = idx_begin % SAL_MINOR_BLOCK_LENGTH;
  const uint64_t total_length = idx_end-idx_begin + first_block_offset;
  SAL_CALCULATE_DIMENSIONS(total_length,num_minor_blocks,num_mayor_blocks);
  locator->total_length = total_length;
  locator->total_size = num_minor_blocks*UINT64_SIZE+num_mayor_blocks*UINT64_SIZE;
  locator->num_mayor_blocks = num_mayor_blocks;
  locator->num_minor_blocks = num_minor_blocks;
  locator->first_block_offset = first_block_offset;
  locator->idx_offset = idx_begin - first_block_offset;
  // Allocate Locator
  locator->mm = (mm_t*) (-1);
  locator->mayor_counters = mm_calloc(num_mayor_blocks,uint64_t,true);
  locator->bitmap = mm_calloc(num_minor_blocks,uint64_t,true);
  // Return
  return locator;
}
void sparse_array_locator_mark(
    sparse_array_locator_t* const locator,
    const uint64_t position) {
  const uint64_t effective_position = position - locator->idx_offset;
  // Locate Block
  SPARSE_ARRAY_LOCATOR_GET_MINOR_BLOCK_LOCATION(effective_position,block_pos,block_mod);
  // Mark position
  const uint64_t block_mask = UINT64_ONE_MASK << block_mod;
  locator->bitmap[block_pos] |= block_mask;
}
/*
 * Writer
 */
typedef struct {
  /* Current Locator Chunk */
  uint64_t* mayor_counter;
  uint64_t* block;
  uint64_t num_blocks;
  /* Global writing information */
  uint64_t global_block_position;
  uint64_t minor_block_position;
  uint64_t mayor_counter_accumulated;
  uint64_t minor_counter_accumulated;
} sparse_array_write_state_t;
void sparse_array_locator_write_chunk(
    fm_t* const file_manager,
    sparse_array_write_state_t* const write_state) {
  // Iterate over all Minor Blocks & write
  uint64_t i;
  for (i=0;i<write_state->num_blocks;++i) {
    // Update MayorCounter
    if (write_state->minor_block_position==SAL_MAYOR_BLOCK_NUM_MINOR_BLOCKS) {
      // Dump Mayor Counter
      write_state->mayor_counter_accumulated += write_state->minor_counter_accumulated;
      *(write_state->mayor_counter++) = write_state->mayor_counter_accumulated;
      // Reset
      write_state->minor_block_position=0;
      write_state->minor_counter_accumulated=0;
    }
    // Update MinorCounter & Write
    uint64_t* const block = write_state->block;
    uint16_t* const block_counter = ((uint16_t*)block) + 3;
    *block_counter = write_state->minor_counter_accumulated;
    fm_write_uint64(file_manager,*block);
    // Update Count
    *block_counter = 0;
    write_state->minor_counter_accumulated += POPCOUNT_64(*block);
    // Next
    ++write_state->block;
    ++write_state->minor_block_position;
  }
}
void sparse_array_locator_write_metadata(
    fm_t* const file_manager,
    const uint64_t total_length,
    const uint64_t total_size,
    const uint64_t num_mayor_blocks,
    const uint64_t num_minor_blocks) {
  fm_write_uint64(file_manager,total_length);
  fm_write_uint64(file_manager,total_size);
  fm_write_uint64(file_manager,num_mayor_blocks);
  fm_write_uint64(file_manager,num_minor_blocks);
}
void sparse_array_locator_write(
    fm_t* const file_manager,
    sparse_array_locator_t* const locator) {
  // Write Meta-Data
  sparse_array_locator_write_metadata(
      file_manager,locator->total_length,locator->total_size,
      locator->num_mayor_blocks,locator->num_minor_blocks);
  // Setup write chunk state
  sparse_array_write_state_t write_state = {
      /* Current Locator Chunk */
      .mayor_counter=locator->mayor_counters,
      .block=locator->bitmap,
      .num_blocks=locator->num_minor_blocks,
      /* Global writing information */
      .global_block_position=0,
      .minor_block_position=SAL_MAYOR_BLOCK_NUM_MINOR_BLOCKS,
      .mayor_counter_accumulated=0,
      .minor_counter_accumulated=0
  };
  // Write Chunk
  sparse_array_locator_write_chunk(file_manager,&write_state);
  // Write Mayor Counters
  fm_write_mem(file_manager,locator->mayor_counters,locator->num_mayor_blocks*SAL_MAYOR_COUNTER_SIZE);
}
void sparse_array_locator_merge__write(
    fm_t* const file_manager,
    sparse_array_locator_t** const locator,
    const uint64_t num_locators) {
  // Calculate Meta-Data
  uint64_t i, total_length;
  for (i=0,total_length=0;i<num_locators;++i) {
    total_length += locator[i]->total_length - locator[i]->first_block_offset;
  }
  SAL_CALCULATE_DIMENSIONS(total_length,num_minor_blocks,num_mayor_blocks);
  const uint64_t total_size = num_minor_blocks*UINT64_SIZE+num_mayor_blocks*UINT64_SIZE;
  // Write Meta-Data
  sparse_array_locator_write_metadata(
      file_manager,total_length,total_size,num_mayor_blocks,num_minor_blocks);
  // Setup write chunk state
  uint64_t* const mayor_counters = mm_calloc(num_mayor_blocks,uint64_t,true);
  sparse_array_write_state_t write_state = {
      /* Global mayor counters */
      .mayor_counter=mayor_counters,
      /* Global writing information */
      .global_block_position=0,
      .minor_block_position=SAL_MAYOR_BLOCK_NUM_MINOR_BLOCKS,
      .mayor_counter_accumulated=0,
      .minor_counter_accumulated=0
  };
  // Traverse all locators
  for (i=0;i<num_locators;++i) {
    const bool overlapping_blocks = (i < num_locators-1) && (locator[i+1]->first_block_offset > 0);
    const uint64_t num_blocks = overlapping_blocks ?
        locator[i]->num_minor_blocks-1 : locator[i]->num_minor_blocks;
    // Setup chunk
    write_state.block=locator[i]->bitmap;
    write_state.num_blocks=num_blocks;
    // Write Chunk
    sparse_array_locator_write_chunk(file_manager,&write_state);
    // OR overlapping block // MEM-HINT: After writing to have the block cached
    if (overlapping_blocks) {
      locator[i+1]->bitmap[0] |= locator[i]->bitmap[num_blocks];
    }
  }
  // Write Mayor Counters
  fm_write_mem(file_manager,mayor_counters,num_mayor_blocks*SAL_MAYOR_COUNTER_SIZE);
  // Free
  mm_free(mayor_counters);
}
/*
 * Dynamic Builder
 */
sparse_array_locator_builder_t* sparse_array_locator_builder_new(mm_slab_t* const mm_slab) {
  // Allocate
  sparse_array_locator_builder_t* const locator_builder = mm_alloc(sparse_array_locator_builder_t);
  // Initialize locator info
  locator_builder->position = 0;
  locator_builder->current_block = 0;
  locator_builder->current_count = 0;
  locator_builder->minor_counter_accumulated = 0;
  locator_builder->mayor_counter_accumulated = 0;
  // Initialize locator memory
  locator_builder->minor_blocks = svector_new(mm_slab,uint64_t);
  svector_iterator_new(&(locator_builder->minor_blocks_iterator),locator_builder->minor_blocks,SVECTOR_WRITE_ITERATOR,0);
  locator_builder->mayor_counters = vector_new(1000,uint64_t);
  // Return
  return locator_builder;
}
void sparse_array_locator_builder_delete(sparse_array_locator_builder_t* const locator_builder) {
  svector_delete(locator_builder->minor_blocks);
  vector_delete(locator_builder->mayor_counters);
  mm_free(locator_builder);
}
void sparse_array_locator_builder_next_word(sparse_array_locator_builder_t* const locator_builder) {
  // Counter OR Bitmap (16b+48b) = UINT64_T
  locator_builder->current_block |= locator_builder->minor_counter_accumulated << SAL_MINOR_BLOCK_LENGTH;
  locator_builder->minor_counter_accumulated = locator_builder->current_count;
  // Store
  *svector_iterator_get_element(&(locator_builder->minor_blocks_iterator),uint64_t) = locator_builder->current_block;
  svector_write_iterator_next(&(locator_builder->minor_blocks_iterator));
  // Reset
  locator_builder->current_block = 0;
}
void sparse_array_locator_builder_next(
    sparse_array_locator_builder_t* const locator_builder,
    const bool mark_position) {
  // Write Mayor Counter
  if (locator_builder->position%SAL_MAYOR_BLOCK_NUM_MINOR_BLOCKS==0) {
    locator_builder->mayor_counter_accumulated += locator_builder->current_count;
    vector_insert(locator_builder->mayor_counters,locator_builder->mayor_counter_accumulated,uint64_t);
    locator_builder->current_count = 0;
    locator_builder->minor_counter_accumulated = 0;
  }
  // Mark as present in the locator
  if (mark_position) {
    locator_builder->current_block = (locator_builder->current_block | SAL_MINOR_BLOCK_ONE_LAST_MASK);
    ++(locator_builder->current_count);
  }
  // Store & reset locator (if needed)
  ++(locator_builder->position);
  if (gem_expect_false(locator_builder->position%SAL_MINOR_BLOCK_LENGTH==0)) {
    sparse_array_locator_builder_next_word(locator_builder);
  } else {
    locator_builder->current_block >>= 1;
  }
}
void sparse_array_locator_builder_write(
    fm_t* const file_manager,
    sparse_array_locator_builder_t* const locator_builder) {
  // Flush
  const uint64_t block_mod = locator_builder->position % SAL_MINOR_BLOCK_LENGTH;
  if (block_mod!=0) {
    locator_builder->current_block >>= SAL_MINOR_BLOCK_LENGTH-(block_mod+1);
    sparse_array_locator_builder_next_word(locator_builder);
  }
  // Write Meta-Data
  const uint64_t total_length = svector_get_used(locator_builder->minor_blocks) * SAL_MINOR_BLOCK_LENGTH;
  SAL_CALCULATE_DIMENSIONS(total_length,num_minor_blocks,num_mayor_blocks);
  const uint64_t total_size = num_minor_blocks*UINT64_SIZE+num_mayor_blocks*UINT64_SIZE;
  sparse_array_locator_write_metadata(
      file_manager,total_length,total_size,num_mayor_blocks,num_minor_blocks);
  // Write Blocks (Bitmaps)
  svector_write(file_manager,locator_builder->minor_blocks);
  // Write Mayor Counters
  fm_write_mem(file_manager,vector_get_mem(locator_builder->mayor_counters,uint64_t),
      num_mayor_blocks*SAL_MAYOR_COUNTER_SIZE);
}
/*
 * Display
 */
void sparse_array_locator_print(
    FILE* const stream,
    const sparse_array_locator_t* const locator,
    const bool display_content) {
  // Display Meta-Data information
  const uint64_t total_size = locator->total_size;
  tab_fprintf(stream,"[GEM]>SparseArrayLocator\n");
  tab_fprintf(stream,"  => Length %"PRIu64"\n",locator->total_length);
  tab_fprintf(stream,"  => Blocks.Mayor %"PRIu64"\n",locator->num_mayor_blocks);
  tab_fprintf(stream,"  => Blocks.Minor %"PRIu64"\n",locator->num_minor_blocks);
  tab_fprintf(stream,"  => Total.Size %"PRIu64" MB (%2.3f%%)\n",CONVERT_B_TO_MB(total_size),PERCENTAGE(total_size,total_size));
  const uint64_t mayor_counters_size = locator->num_mayor_blocks*SAL_MAYOR_COUNTER_SIZE;
  tab_fprintf(stream,"    => MayorCounter.Size %"PRIu64" MB (%2.3f%%)\n",
      CONVERT_B_TO_MB(mayor_counters_size),PERCENTAGE(mayor_counters_size,total_size));
  const uint64_t minor_counters_size = locator->num_minor_blocks*SAL_MINOR_COUNTER_LENGTH/8;
  tab_fprintf(stream,"    => MinorCounter.Size %"PRIu64" MB (%2.3f%%)\n",
      CONVERT_B_TO_MB(minor_counters_size),PERCENTAGE(minor_counters_size,total_size));
  const uint64_t bitmap_size = locator->num_minor_blocks*SAL_MINOR_BLOCK_LENGTH/8;
  tab_fprintf(stream,"    => Bitmap.Size %"PRIu64" MB (%2.3f%%)\n",
      CONVERT_B_TO_MB(bitmap_size),PERCENTAGE(bitmap_size,total_size));
  // Display checksums
  uint64_t i, bitmap_checksum=0, mayor_counters_checksum=0;
  for (i=0;i<locator->num_minor_blocks;++i) checksum_incremental_uint64(&bitmap_checksum,locator->bitmap[i]);
  for (i=0;i<locator->num_mayor_blocks;++i) checksum_incremental_uint64(&mayor_counters_checksum,locator->mayor_counters[i]);
  tab_fprintf(stream,"  => Checksums\n");
  tab_fprintf(stream,"  =>   CS.Bitmaps %"PRIu64"\n",bitmap_checksum);
  tab_fprintf(stream,"  =>   CS.MayorCounters %"PRIu64"\n",mayor_counters_checksum);
  // Display content
  if (display_content) {
    for (i=0;i<locator->num_minor_blocks;++i) {
      const uint64_t block = locator->bitmap[i] & SAL_MINOR_BLOCK_BITMAP_MASK;
      const uint16_t block_counter = *(((uint16_t*)(locator->bitmap+i)) + 3);
      fprintf(stream,"(%u)[",block_counter);
      fprintf_uint64_binary(stream,block);fprintf(stream,"]\n");
    }
  }
  // Flush
  fflush(stream);
}
/*
 * Stats
 */
sparse_array_locator_stats_t* sparse_array_locator_stats_new(void) {
  // Allocate
  sparse_array_locator_stats_t* const locator_stats = mm_alloc(sparse_array_locator_stats_t);
  // Locator Stats
  locator_stats->marked_positions = 0;
  // Bitmaps
  locator_stats->ones_density = stats_vector_raw_new(SAL_MINOR_BLOCK_LENGTH,1);
  // Return
  return locator_stats;
}
void sparse_array_locator_stats_delete(sparse_array_locator_stats_t* const locator_stats) {
  // Free
  stats_vector_delete(locator_stats->ones_density);
  mm_free(locator_stats);
}
void sparse_array_locator_stats_calculate(
    sparse_array_locator_stats_t* const locator_stats,
    sparse_array_locator_t* const locator) {
  uint64_t i;
  // Calculate locator stats
  locator_stats->locator_length = locator->total_length;
  for (i=0;i<locator->num_minor_blocks;++i) {
    const uint64_t block = locator->bitmap[i] & SAL_MINOR_BLOCK_BITMAP_MASK;
    const uint64_t ones_count = POPCOUNT_64(block);
    locator_stats->marked_positions += ones_count;
    stats_vector_inc(locator_stats->ones_density,ones_count);
  }
}
void sparse_array_locator_stats_print(
    FILE* const stream,
    const char* const sparse_array_locator_stats_tag,
    sparse_array_locator_stats_t* const locator_stats) {
  // Display Statistics
  fprintf(stream,"SparseArrayLocator.Stats %s\n",sparse_array_locator_stats_tag);
  fprintf(stream,"  => Locator.Marked %"PRIu64" (%2.3f%%)\n",
      locator_stats->marked_positions,
      PERCENTAGE(locator_stats->marked_positions,locator_stats->locator_length));
  fprintf(stream,"  => Bitmaps.48b.OnesDensity\n");
  tab_global_inc(); tab_global_inc();
  stats_vector_display(stream,locator_stats->ones_density,false,false,NULL);
  tab_global_dec(); tab_global_dec();
  // Flush
  fflush(stream);
}

