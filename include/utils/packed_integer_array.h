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

#ifndef PACKED_INTEGER_ARRAY_H_
#define PACKED_INTEGER_ARRAY_H_

#include "utils/essentials.h"

/*
 * Sampling Rate
 */
typedef enum {
  SAMPLING_RATE_1   = 0,
  SAMPLING_RATE_2   = 1,
  SAMPLING_RATE_4   = 2,
  SAMPLING_RATE_8   = 3,
  SAMPLING_RATE_16  = 4,
  SAMPLING_RATE_32  = 5,
  SAMPLING_RATE_64  = 6,
  SAMPLING_RATE_128 = 7,
  SAMPLING_RATE_256 = 8,
  SAMPLING_RATE_NONE = UINT64_MAX-1,
  SAMPLING_RATE_RANGE = UINT64_MAX,
} sampling_rate_t;

/*
 * Packed Integer Array
 */
typedef struct {
  /* Meta-Data */
  uint64_t num_elements;       // Number of elements
  uint64_t integer_length;     // Bits per integer
  uint64_t array_size;         // Total Size (Bytes)
  /* Array */
  uint64_t* bitmap;            // Static Packed-Integer Array
  /* MM */
  mm_t* mm_bitmap;
} packed_integer_array_t;
typedef struct {
  /* Meta-Data */
  uint64_t num_elements;            // Number of elements
  uint64_t integer_length;          // Bits per integer
  /* Array */
  svector_t* bitmap;                // Dynamic Packed-Integer Array
  /* Writer */
  uint64_t lo64_offset;             // lo-word offset
  svector_iterator_t bitmap_writer; // svector write iterator (Bitmap)
} packed_integer_array_builder_t;

/*
 * Loader/Setup
 */
packed_integer_array_t* packed_integer_array_new(
    const uint64_t num_elements,
    const uint64_t integer_length_bits);
void packed_integer_array_delete(packed_integer_array_t* const array);
packed_integer_array_t* packed_integer_array_read(fm_t* const file_manager);
packed_integer_array_t* packed_integer_array_read_mem(mm_t* const memory_manager);
void packed_integer_array_write(
    fm_t* const file_manager,
    packed_integer_array_t* const array);

/*
 * Builder
 */
packed_integer_array_builder_t* packed_integer_array_builder_new(
    const uint64_t integer_length_bits,
    mm_slab_t* const mm_slab);
void packed_integer_array_builder_delete(
    packed_integer_array_builder_t* const array);
void packed_integer_array_builder_store(
    packed_integer_array_builder_t* const array,
    const uint64_t integer);
void packed_integer_array_builder_write(
    fm_t* const file_manager,
    packed_integer_array_builder_t** array_builders,
    const uint64_t num_array_builders);

/*
 * Accessors
 */
uint64_t packed_integer_array_get_size(const packed_integer_array_t* const array);
uint64_t packed_integer_array_get_length(const packed_integer_array_t* const array);

void packed_integer_array_prefetch(
    const packed_integer_array_t* const array,
    const uint64_t position);
uint64_t packed_integer_array_load(
    const packed_integer_array_t* const array,
    const uint64_t position);

void packed_integer_array_store(
    packed_integer_array_t* const array,
    const uint64_t position,
    const uint64_t integer);

/*
 * Display
 */
void packed_integer_array_print(
    FILE* const stream,
    const packed_integer_array_t* const array,
    const bool display_data);

#endif /* PACKED_INTEGER_ARRAY_H_ */
