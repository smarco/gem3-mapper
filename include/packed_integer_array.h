/*
 * PROJECT: GEMMapper
 * FILE: packed_integer_array.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef PACKED_INTEGER_ARRAY_H_
#define PACKED_INTEGER_ARRAY_H_

#include "essentials.h"

/*
 * Checkers
 */
#define PACKED_INTEGER_ARRAY_CHECK(array) \
  GEM_CHECK_NULL(array); \
  GEM_CHECK_ZERO(array->integer_length); \
  GEM_CHECK_NULL(array->bitmap)

#define PACKED_INTEGER_ARRAY_CHECK_NO_NULL(array) \
  GEM_CHECK_NULL(array); \
  GEM_CHECK_ZERO(array->array_size); \
  GEM_CHECK_ZERO(array->num_elements); \
  GEM_CHECK_ZERO(array->integer_length); \
  GEM_CHECK_NULL(array->bitmap)

/*
 * Packed Integer Array
 */
typedef struct {
  /* Meta-Data */
  uint64_t num_elements;       // Number of elements
  uint64_t integer_length;     // Bits per integer
  uint64_t array_size;         // Total Size (Bytes)
  /* Array */
  uint64_t* bitmap;            // Packed Integer Array
  /* MM */
  mm_t* mm_bitmap;
  /* Space Range [Idx_start,Idx_end) */
  uint64_t first_block_offset; // First block offset
  uint64_t idx_offset;         // Global offset
} packed_integer_array_t;

/*
 * Loader/Setup
 */
GEM_INLINE packed_integer_array_t* packed_integer_array_new(
    const uint64_t idx_begin,const uint64_t idx_end,const uint64_t integer_length_bits);
GEM_INLINE packed_integer_array_t* packed_integer_array_read(fm_t* const file_manager);
GEM_INLINE packed_integer_array_t* packed_integer_array_read_mem(mm_t* const memory_manager);
GEM_INLINE void packed_integer_array_write(
    fm_t* const file_manager,packed_integer_array_t* const array);
GEM_INLINE void packed_integer_array_write_adaptor(
    fm_t* const file_manager,const uint64_t max_integer,const uint64_t num_integers,
    const uint64_t* const raw_integer_array);
GEM_INLINE void packed_integer_array_sliced_write(
    fm_t* const file_manager,packed_integer_array_t** const array,const uint64_t num_array);
GEM_INLINE void packed_integer_array_delete(packed_integer_array_t* const array);

/*
 * Accessors
 */
GEM_INLINE uint64_t packed_integer_array_load(packed_integer_array_t* const array,const uint64_t position);
GEM_INLINE void packed_integer_array_store(packed_integer_array_t* const array,const uint64_t position,const uint64_t integer);

GEM_INLINE uint64_t packed_integer_array_get_size(packed_integer_array_t* const array);
GEM_INLINE uint64_t packed_integer_array_get_length(packed_integer_array_t* const array);
GEM_INLINE uint64_t packed_integer_array_sliced_get_size(
    packed_integer_array_t** const array,const uint64_t num_array);

/*
 * Display
 */
GEM_INLINE void packed_integer_array_print(FILE* const stream,packed_integer_array_t* const array,const bool display_data);

#endif /* PACKED_INTEGER_ARRAY_H_ */
