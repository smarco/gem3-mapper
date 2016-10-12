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

#include "utils/packed_integer_array.h"

// Packed Array Size ((+8) to simplify the access to the last integer)
#define PACKED_INT_ARRAY_SIZE(num_elements,integer_length_bits) \
  (DIV_CEIL(((num_elements)*(integer_length_bits)),UINT64_LENGTH)*UINT64_SIZE + UINT64_SIZE)

/*
 * Loader/Setup
 */
packed_integer_array_t* packed_integer_array_new(
    const uint64_t num_elements,
    const uint64_t integer_length_bits) {
  // Allocate handler
  packed_integer_array_t* const array = mm_alloc(packed_integer_array_t);
  // Set Dimensions
  array->num_elements = num_elements;
  array->integer_length = integer_length_bits;
  // Array
  array->array_size = PACKED_INT_ARRAY_SIZE(array->num_elements,integer_length_bits);
  array->bitmap = (uint64_t*) mm_calloc(array->array_size,uint8_t,true);
  array->mm_bitmap = (mm_t*)(-1); // Trick to guide deallocation
  // Return
  return array;
}
void packed_integer_array_delete(packed_integer_array_t* const array) {
  if (array->mm_bitmap!=NULL) {
    if (array->mm_bitmap==(mm_t*)(-1)) {
      mm_free(array->bitmap);
    } else {
      mm_bulk_free(array->mm_bitmap);
    }
  }
  mm_free(array);
}
packed_integer_array_t* packed_integer_array_read(fm_t* const file_manager) {
  // Allocate handler
  packed_integer_array_t* const array = mm_alloc(packed_integer_array_t);
  // Read Meta-Data
  array->num_elements = fm_read_uint64(file_manager);
  array->integer_length = fm_read_uint64(file_manager);
  array->array_size = fm_read_uint64(file_manager);
  // Read array
  array->mm_bitmap = fm_load_mem(file_manager,array->array_size);
  array->bitmap = mm_get_base_mem(array->mm_bitmap);
  // Return
  return array;
}
packed_integer_array_t* packed_integer_array_read_mem(mm_t* const memory_manager) {
  // Allocate handler
  packed_integer_array_t* const array = mm_alloc(packed_integer_array_t);
  // Read Meta-Data
  array->num_elements = mm_read_uint64(memory_manager);
  array->integer_length = mm_read_uint64(memory_manager);
  array->array_size = mm_read_uint64(memory_manager);
  // Read array
  array->mm_bitmap = NULL;
  array->bitmap = mm_read_mem(memory_manager,array->array_size);
  // Return
  return array;
}
void packed_integer_array_write_metadata(
    fm_t* const file_manager,
    const uint64_t num_elements,
    const uint64_t integer_length,
    const uint64_t array_size) {
  // Write Meta-Data
  fm_write_uint64(file_manager,num_elements);
  fm_write_uint64(file_manager,integer_length);
  fm_write_uint64(file_manager,array_size);
}
void packed_integer_array_write(
    fm_t* const file_manager,
    packed_integer_array_t* const array) {
  // Write Meta-Data
  packed_integer_array_write_metadata(file_manager,
      array->num_elements,array->integer_length,array->array_size);
  // Write array
  fm_write_mem(file_manager,array->bitmap,array->array_size);
}
/*
 * Builder
 */
packed_integer_array_builder_t* packed_integer_array_builder_new(
    const uint64_t integer_length_bits,
    mm_slab_t* const mm_slab) {
  // Allocate handler
  packed_integer_array_builder_t* const array = mm_alloc(packed_integer_array_builder_t);
  // Set Dimensions
  array->num_elements = 0;
  array->integer_length = integer_length_bits;
  // Array
  array->bitmap = svector_new(mm_slab,uint64_t);
  // Writer
  array->lo64_offset = 0;
  svector_iterator_new(&array->bitmap_writer,array->bitmap,SVECTOR_WRITE_ITERATOR,0);
  // Return
  return array;
}
void packed_integer_array_builder_delete(packed_integer_array_builder_t* const array) {
  svector_delete(array->bitmap);
  mm_free(array);
}
void packed_integer_array_builder_store(
    packed_integer_array_builder_t* const array,
    const uint64_t integer) {
  uint64_t* next_word = svector_iterator_get_element(&array->bitmap_writer,uint64_t);
  // Store LO part
  *next_word |= (integer << array->lo64_offset);
  // Store HI part
  const uint64_t hi64_offset = UINT64_LENGTH-array->lo64_offset;
  if (hi64_offset < array->integer_length) {
    svector_write_iterator_next(&array->bitmap_writer);
    next_word = svector_iterator_get_element(&array->bitmap_writer,uint64_t);
    *next_word |= (integer >> hi64_offset);
  } else if (hi64_offset == array->integer_length) {
    svector_write_iterator_next(&array->bitmap_writer);
  }
  // Next
  ++array->num_elements;
  array->lo64_offset = (array->lo64_offset+array->integer_length) % UINT64_LENGTH;
}
typedef struct {
  /* Packed Integer Arrays Builders */
  uint64_t num_array_builders;
  packed_integer_array_builder_t** array_builders;
  uint64_t total_elements;
  /* Iterator */
  packed_integer_array_builder_t* current_array_builder;
  uint64_t array_builder_position;
  uint64_t array_builder_pending;
  uint64_t lo64_offset;
  svector_iterator_t bitmap_iterator;
} array_builder_hub_t;
void array_builder_hub_load_builder(
    array_builder_hub_t* const array_builder_hub,
    const uint64_t array_builder_position) {
  array_builder_hub->current_array_builder = array_builder_hub->array_builders[array_builder_position];
  array_builder_hub->array_builder_pending = array_builder_hub->current_array_builder->num_elements;
  if (array_builder_hub->array_builder_pending > 0) {
    array_builder_hub->lo64_offset = 0;
    svector_iterator_new(&array_builder_hub->bitmap_iterator,
        array_builder_hub->current_array_builder->bitmap,SVECTOR_READ_ITERATOR,0);
  } else {
    ++(array_builder_hub->array_builder_position);
    if (array_builder_hub->array_builder_position < array_builder_hub->num_array_builders) {
      array_builder_hub_load_builder(array_builder_hub,array_builder_hub->array_builder_position);
    } else {
      array_builder_hub->array_builder_pending = 0;
    }
  }
}
void array_builder_hub_new(
    array_builder_hub_t* const array_builder_hub,
    packed_integer_array_builder_t** array_builders,
    const uint64_t num_array_builders) {
  // Packed Integer Arrays Builders
  array_builder_hub->num_array_builders = num_array_builders;
  array_builder_hub->array_builders = array_builders;
  // Calculate dimensions
  array_builder_hub->total_elements = 0;
  uint64_t i;
  for (i=0;i<num_array_builders;++i) {
    array_builder_hub->total_elements += array_builders[i]->num_elements;
    if (array_builders[i]->lo64_offset > 0) svector_write_iterator_next(&(array_builders[i]->bitmap_writer));
  }
  // Iterator
  array_builder_hub->array_builder_position = 0;
  array_builder_hub_load_builder(array_builder_hub,0);
}
uint64_t array_builder_hub_get_num_elements(array_builder_hub_t* const array_builder_hub) {
  return array_builder_hub->total_elements;
}
bool array_builder_hub_eoi(array_builder_hub_t* const array_builder_hub) {
  return (array_builder_hub->array_builder_pending==0);
}
uint64_t array_builder_hub_next(array_builder_hub_t* const array_builder_hub) {
  // Locate next integer
  const uint64_t lo64_offset = array_builder_hub->lo64_offset;
  const uint64_t hi64_offset = UINT64_LENGTH - lo64_offset;
  const uint64_t current_word = *svector_iterator_get_element(&array_builder_hub->bitmap_iterator,uint64_t);
  const uint64_t integer_length = array_builder_hub->current_array_builder->integer_length;
  // Calculate the next integer
  uint64_t next_integer;
  if (hi64_offset >= integer_length) {
    next_integer = uint64_mask_ones[integer_length] & (current_word >> lo64_offset);
    if (hi64_offset==integer_length) {
      svector_read_iterator_next(&array_builder_hub->bitmap_iterator);
    }
  } else {
    svector_read_iterator_next(&array_builder_hub->bitmap_iterator);
    const uint64_t next_word = *svector_iterator_get_element(&array_builder_hub->bitmap_iterator,uint64_t);
    next_integer = uint64_mask_ones[integer_length] & ((current_word >> lo64_offset) | (next_word << hi64_offset));
  }
  // Update iterator
  --(array_builder_hub->array_builder_pending);
  if (array_builder_hub->array_builder_pending==0) {
    ++(array_builder_hub->array_builder_position);
    if (array_builder_hub->array_builder_position < array_builder_hub->num_array_builders) {
      array_builder_hub_load_builder(array_builder_hub,array_builder_hub->array_builder_position);
    }
  } else {
    array_builder_hub->lo64_offset = (array_builder_hub->lo64_offset+integer_length) % UINT64_LENGTH;
  }
  // Return
  return next_integer;
}
void packed_integer_array_builder_write(
    fm_t* const file_manager,
    packed_integer_array_builder_t** array_builders,
    const uint64_t num_array_builders) {
  // Create Hub
  array_builder_hub_t array_builder_hub;
  array_builder_hub_new(&array_builder_hub,array_builders,num_array_builders);
  // Write Meta-Data
  const uint64_t integer_length = array_builders[0]->integer_length;
  const uint64_t total_elements = array_builder_hub_get_num_elements(&array_builder_hub);
  const uint64_t array_size = PACKED_INT_ARRAY_SIZE(total_elements,integer_length);
  packed_integer_array_write_metadata(file_manager,total_elements,integer_length,array_size);
  // Write array
  uint64_t i, offset=0, bitmap_word=0;
  for (i=0;i<total_elements;++i) {
    // Get next array element
    const uint64_t integer = array_builder_hub_next(&array_builder_hub);
    // Store LO part
    bitmap_word |= (integer << offset);
    // Store HI part
    offset += integer_length;
    if (offset >= UINT64_LENGTH) {
      fm_write_uint64(file_manager,bitmap_word);
      bitmap_word = 0;
      offset -= UINT64_LENGTH;
      if (offset > 0) {
        bitmap_word |= (integer >> (integer_length-offset));
      }
    }
  }
  // Write
  if (offset > 0) {
    fm_write_uint64(file_manager,bitmap_word);
  }
  // Write ((+8) padding)
  fm_write_uint64(file_manager,0);
}
/*
 * Accessors
 */
uint64_t packed_integer_array_get_size(const packed_integer_array_t* const array) {
  return array->array_size;
}
uint64_t packed_integer_array_get_length(const packed_integer_array_t* const array) {
  return array->num_elements;
}
#define PACKED_INT_ARRAY_GET_LOCATION(lo64_bit_position,lo64_word_position,lo64_offset,hi64_offset) \
  const uint64_t lo64_bit_position = position*array->integer_length; \
  const uint64_t lo64_word_position = lo64_bit_position/UINT64_LENGTH; \
  const uint64_t lo64_offset = lo64_bit_position % UINT64_LENGTH; \
  const uint64_t hi64_offset = UINT64_LENGTH-lo64_offset
void packed_integer_array_prefetch(
    const packed_integer_array_t* const array,
    const uint64_t position) {
  // Prefetch
  const uint64_t lo64_bit_position = position*array->integer_length;
  const uint64_t lo64_word_position = lo64_bit_position/UINT64_LENGTH;
  PREFETCH(array->bitmap+lo64_word_position);
}
uint64_t packed_integer_array_load(
    const packed_integer_array_t* const array,
    const uint64_t position) {
  // Load LO & HI part
  PACKED_INT_ARRAY_GET_LOCATION(lo64_bit_position,lo64_word_position,lo64_offset,hi64_offset);
  if (hi64_offset >= array->integer_length) {
    return uint64_mask_ones[array->integer_length] & (array->bitmap[lo64_word_position] >> lo64_offset);
  } else {
    return uint64_mask_ones[array->integer_length] &
        ((array->bitmap[lo64_word_position] >> lo64_offset) | (array->bitmap[lo64_word_position+1] << hi64_offset));
  }
}
void packed_integer_array_store(
    packed_integer_array_t* const array,
    const uint64_t position,
    const uint64_t integer) {
  // Locate
  PACKED_INT_ARRAY_GET_LOCATION(lo64_bit_position,lo64_word_position,lo64_offset,hi64_offset);
  // Store LO part
  array->bitmap[lo64_word_position] |= (integer << lo64_offset);
  // Store HI part
  if (hi64_offset < array->integer_length) {
    array->bitmap[lo64_word_position+1] |= (integer >> hi64_offset);
  }
}
/*
 * Display
 */
void packed_integer_array_print(
    FILE* const stream,
    const packed_integer_array_t* const array,
    const bool display_data) {
  // Print meta-info
  tab_fprintf(stream,"[GEM]>PackedArray.info\n");
  tab_fprintf(stream,"  => PackedArray.Size %"PRIu64" MB\n",CONVERT_B_TO_MB(array->array_size));
  tab_fprintf(stream,"  => PackedArray.Length %"PRIu64"\n",array->num_elements);
  tab_fprintf(stream,"  => PackedArray.Integer %"PRIu64" bits\n",array->integer_length);
  // Print data
  if (display_data) {
    tab_fprintf(stream,"  => PackedArray.Data\n");
    uint64_t i;
    for (i=0;i<array->num_elements;++i) {
      fprintf(stream,"%8"PRIu64,packed_integer_array_load(array,i));
      PRINT_COLUMN(stream,i,10,"\t");
    }
    PRINT_COLUMN_CLOSE(stream,i,10);
  }
  // Flush
  fflush(stream);
}

