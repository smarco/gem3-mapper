/*
 * PROJECT: GEMMapper
 * FILE: sparse_bitmap.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   Provides functionality to store 64bit-words. It doesn't explicitly store zero bitmaps (sparse)
 *   Used to store the 3rd-bit layer of a 3bits alphabet. Due to the nature of the problem,
 *   this bit is seldom set. So its representation in a sparse-fashion results
 *   in a mayor saving of memory space.
 */

#ifndef SPARSE_BITMAP_H_
#define SPARSE_BITMAP_H_

#include "utils/essentials.h"
#include "utils/segmented_vector.h"
#include "utils/sparse_array_locator.h"

/*
 * Sparse Bitmap Storage for 64bit bitmaps
 *   STRUCTURE := LOCATOR + BITMAPS
 *   BITMAPS := [64Bitmap]...[64Bitmap]
 */
typedef struct {
  /* Bitmaps */
  uint64_t* bitmaps;       // Bitmaps stored contiguously
  uint64_t num_bitmaps;    // Total number of bitmaps
  /* Locator */
  sparse_array_locator_t* locator;
  /* MM */
  mm_t* mm_bitmaps;
} sparse_bitmap_t;
/*
 * Sparse Bitmap Builder
 */
typedef struct {
  /* Bitmaps */
  svector_t* bitmaps;
  svector_iterator_t bitmaps_iterator;
  /* Locator */
  sparse_array_locator_builder_t* locator_builder;
} sparse_bitmap_builder_t;

/*
 * Loader/Setup
 */
sparse_bitmap_t* sparse_bitmap_read(fm_t* const file_manager);
sparse_bitmap_t* sparse_bitmap_read_mem(mm_t* const memory_manager);
void sparse_bitmap_delete(sparse_bitmap_t* const sparse_bitmap);

/*
 * Accessors
 */
uint64_t sparse_bitmap_get_size(sparse_bitmap_t* const sparse_bitmap);

bool sparse_bitmap_is_contained(
    sparse_bitmap_t* const sparse_bitmap,
    const uint64_t position);
uint64_t sparse_bitmap_get_bitmap(
    sparse_bitmap_t* const sparse_bitmap,
    const uint64_t position);
bool sparse_bitmap_get_bitmap_if_contained(
    sparse_bitmap_t* const sparse_bitmap,
    const uint64_t position);


/*
 * Builder
 */
sparse_bitmap_builder_t* sparse_bitmap_builder_new(mm_slab_t* const mm_slab);
void sparse_bitmap_builder_delete(
    sparse_bitmap_builder_t* const sparse_bitmap_builder);
void sparse_bitmap_builder_add_bitmap(
    sparse_bitmap_builder_t* const sparse_bitmap_builder,
    const uint64_t bitmap);
void sparse_bitmap_builder_skip_bitmap(
    sparse_bitmap_builder_t* const sparse_bitmap_builder);
void sparse_bitmap_builder_write(
    fm_t* const file_manager,
    sparse_bitmap_builder_t* const sparse_bitmap_builder);

/*
 * Display
 */
void sparse_bitmap_print(
    FILE* const stream,
    sparse_bitmap_t* const sparse_bitmap,
    const bool display_content);

#endif /* SPARSE_BITMAP_H_ */
