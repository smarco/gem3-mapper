/*
 * PROJECT: GEMMapper
 * FILE: compact_dna_text.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   Provides functionality to store 64bit-words. It doesn't explicitly store zero bitmaps (sparse)
 *   Used to store the 3rd-bit layer of a 3bits alphabet. Due to the nature of the problem,
 *   this bit is seldom set. So its representation in a sparse-fashion results
 *   in a mayor saving of memory space.
 */

#include "sparse_bitmap.h"

#define SPARSE_BITMAP_BLOCK_LENGTH UINT128_LENGTH
#define SPARSE_BITMAP_BLOCK_MOD UINT64_LENGTH

#define SPARSE_BITMAP_EXTRACT_MASK UINT64_ONE_MASK

/*
 * Loader/Setup
 */
GEM_INLINE sparse_bitmap_t* sparse_bitmap_read(fm_t* const file_manager) {
  FM_CHECK(file_manager);
  // Allocate
  sparse_bitmap_t* const sparse_bitmap = mm_alloc(sparse_bitmap_t);
  // Read Bitmaps
  sparse_bitmap->num_bitmaps  = fm_read_uint64(file_manager);
  sparse_bitmap->mm_bitmaps = fm_load_mem(file_manager,sparse_bitmap->num_bitmaps*UINT64_SIZE);
  sparse_bitmap->bitmaps = mm_get_base_mem(sparse_bitmap->mm_bitmaps);
  // Read Locator
  sparse_bitmap->locator = sparse_array_locator_read(file_manager);
  // Return
  return sparse_bitmap;
}
GEM_INLINE sparse_bitmap_t* sparse_bitmap_read_mem(mm_t* const memory_manager) {
  MM_CHECK(memory_manager);
  // Allocate
  sparse_bitmap_t* const sparse_bitmap = mm_alloc(sparse_bitmap_t);
  // Read Bitmaps
  sparse_bitmap->num_bitmaps  = mm_read_uint64(memory_manager);
  sparse_bitmap->mm_bitmaps = NULL;
  sparse_bitmap->bitmaps = mm_read_mem(memory_manager,sparse_bitmap->num_bitmaps*UINT64_SIZE);
  // Read Locator
  sparse_bitmap->locator = sparse_array_locator_read_mem(memory_manager);
  // Return
  return sparse_bitmap;
}
GEM_INLINE void sparse_bitmap_delete(sparse_bitmap_t* const sparse_bitmap) {
  SPARSE_BITMAP_CHECK(sparse_bitmap);
  if (sparse_bitmap->mm_bitmaps!=NULL) mm_bulk_free(sparse_bitmap->mm_bitmaps);
  mm_free(sparse_bitmap);
}
/*
 * Accessors
 */
GEM_INLINE uint64_t sparse_bitmap_get_size(sparse_bitmap_t* const sparse_bitmap) {
  const uint64_t bitmaps_size = sparse_bitmap->num_bitmaps*UINT64_SIZE;
  return sparse_array_locator_get_size(sparse_bitmap->locator) + bitmaps_size;
}
GEM_INLINE bool sparse_bitmap_is_contained(sparse_bitmap_t* const sparse_bitmap,const uint64_t position) {
  SPARSE_BITMAP_CHECK(sparse_bitmap);
  return sparse_array_locator_is_marked(sparse_bitmap->locator,position);
}
GEM_INLINE uint64_t sparse_bitmap_get_bitmap(sparse_bitmap_t* const sparse_bitmap,const uint64_t position) {
  SPARSE_BITMAP_CHECK(sparse_bitmap);
  uint64_t bimap_position;
  if (sparse_array_locator_get_erank_if_marked(sparse_bitmap->locator,position,&bimap_position)) {
    return sparse_bitmap->bitmaps[bimap_position];
  } else {
    return 0;
  }
}
/*
 * Builder
 */
GEM_INLINE sparse_bitmap_builder_t* sparse_bitmap_builder_new(mm_slab_t* const mm_slab) {
  // Allocate
  sparse_bitmap_builder_t* const sparse_bitmap_builder = mm_alloc(sparse_bitmap_builder_t);
  // Initialize bitmaps
  sparse_bitmap_builder->bitmaps = svector_new(mm_slab,uint64_t);
  svector_iterator_new(&(sparse_bitmap_builder->bitmaps_iterator),sparse_bitmap_builder->bitmaps,SVECTOR_WRITE_ITERATOR,0);
  // Initialize locator
  sparse_bitmap_builder->locator_builder = sparse_array_locator_builder_new(mm_slab);
  // Return
  return sparse_bitmap_builder;
}
GEM_INLINE void sparse_bitmap_builder_delete(sparse_bitmap_builder_t* const sparse_bitmap_builder) {
  SPARSE_BITMAP_BUILDER_CHECK(sparse_bitmap_builder);
  sparse_array_locator_builder_delete(sparse_bitmap_builder->locator_builder);
  svector_delete(sparse_bitmap_builder->bitmaps);
  mm_free(sparse_bitmap_builder);
}
GEM_INLINE void sparse_bitmap_builder_add_bitmap(sparse_bitmap_builder_t* const sparse_bitmap_builder,const uint64_t bitmap) {
  SPARSE_BITMAP_BUILDER_CHECK(sparse_bitmap_builder);
  // Mark as present in the locator
  sparse_array_locator_builder_next(sparse_bitmap_builder->locator_builder,true);
  // Store bitmap
  *svector_iterator_get_element(&(sparse_bitmap_builder->bitmaps_iterator),uint64_t) = bitmap;
  svector_write_iterator_next(&(sparse_bitmap_builder->bitmaps_iterator));
}
GEM_INLINE void sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder_t* const sparse_bitmap_builder) {
  SPARSE_BITMAP_BUILDER_CHECK(sparse_bitmap_builder);
  // Skip position in the locator
  sparse_array_locator_builder_next(sparse_bitmap_builder->locator_builder,false);
}
GEM_INLINE void sparse_bitmap_builder_write(fm_t* const file_manager,sparse_bitmap_builder_t* const sparse_bitmap_builder) {
  SPARSE_BITMAP_BUILDER_CHECK(sparse_bitmap_builder);
  FM_CHECK(file_manager);
  // Write header
  fm_write_uint64(file_manager,svector_get_used(sparse_bitmap_builder->bitmaps)); // num_bitmaps
  // Write bitmaps
  svector_write(file_manager,sparse_bitmap_builder->bitmaps);
  // Write locator
  sparse_array_locator_builder_write(file_manager,sparse_bitmap_builder->locator_builder);
}
/*
 * Display
 */
GEM_INLINE void sparse_bitmap_print(FILE* const stream,sparse_bitmap_t* const sparse_bitmap,const bool display_content) {
  GEM_CHECK_NULL(stream);
  SPARSE_BITMAP_CHECK(sparse_bitmap);
  const uint64_t bitmaps_size = sparse_bitmap->num_bitmaps*UINT64_SIZE;
  const uint64_t total_size = sparse_array_locator_get_size(sparse_bitmap->locator) + bitmaps_size;
  tab_fprintf(stream,"[GEM]>SparseBitmap\n");
  tab_fprintf(stream,"  => Total.Size %lu (%2.3f%%)\n",total_size,PERCENTAGE(total_size,total_size));
  tab_fprintf(stream,"  => Bitmaps.Num %lu\n",sparse_bitmap->num_bitmaps);
  tab_fprintf(stream,"  => Bitmaps.Size %lu (%2.3f%%)\n",bitmaps_size,PERCENTAGE(bitmaps_size,total_size));
  sparse_array_locator_print(stream,sparse_bitmap->locator,display_content);
  if (display_content) {
    tab_fprintf(stream,"  => SparseBitmap.Content\n");
    /* Display full bitmaps */
    uint64_t i;
    const uint64_t num_bitmaps = sparse_bitmap->num_bitmaps;
    for (i=0;i<num_bitmaps;++i) {
      fprintf(stream," %8lu",sparse_bitmap->bitmaps[i]); // Print bitmap
      PRINT_COLUMN(stream,i,6," "); // Display in columns
    }
    PRINT_COLUMN_CLOSE(stream,i,6);
  }
  // Flush
  fflush(stream);
}

