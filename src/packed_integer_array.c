/*
 * PROJECT: GEMMapper
 * FILE: packed_integer_array.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "packed_integer_array.h"

// Packed Array Size ((+8) to simplify the access to the last integer)
#define PACKED_INT_ARRAY_SIZE(num_elements,integer_length_bits) \
  (DIV_CEIL(((num_elements)*(integer_length_bits)),UINT64_LENGTH)*UINT64_SIZE + UINT64_SIZE)

/*
 * Loader/Setup
 */
GEM_INLINE packed_integer_array_t* packed_integer_array_new(
    const uint64_t idx_begin,const uint64_t idx_end,const uint64_t integer_length_bits) {
  // Allocate handler
  packed_integer_array_t* const array = mm_alloc(packed_integer_array_t);
  // Set Dimensions
  array->num_elements = idx_end-idx_begin;
  array->integer_length = integer_length_bits;
  // Array
  array->array_size = PACKED_INT_ARRAY_SIZE(array->num_elements,integer_length_bits);
  array->bitmap = (uint64_t*) mm_calloc(array->array_size,uint8_t,true);
  array->mm_bitmap = (mm_t*)(-1); // Trick to guide deallocation
  // Return
  return array;
}
GEM_INLINE packed_integer_array_t* packed_integer_array_read(fm_t* const file_manager) {
  FM_CHECK(file_manager);
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
GEM_INLINE packed_integer_array_t* packed_integer_array_read_mem(mm_t* const memory_manager) {
  MM_CHECK(memory_manager);
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
GEM_INLINE void packed_integer_array_write_metadata(
    fm_t* const file_manager,const uint64_t num_elements,
    const uint64_t integer_length,const uint64_t array_size) {
  // Write Meta-Data
  fm_write_uint64(file_manager,num_elements);
  fm_write_uint64(file_manager,integer_length);
  fm_write_uint64(file_manager,array_size);
}
GEM_INLINE void packed_integer_array_write(fm_t* const file_manager,packed_integer_array_t* const array) {
  FM_CHECK(file_manager);
  PACKED_INTEGER_ARRAY_CHECK_NO_NULL(array);
  // Write Meta-Data
  packed_integer_array_write_metadata(file_manager,
      array->num_elements,array->integer_length,array->array_size);
  // Write array
  fm_write_mem(file_manager,array->bitmap,array->array_size);
}
GEM_INLINE void packed_integer_array_write_adaptor(
    fm_t* const file_manager,const uint64_t max_integer,const uint64_t num_integers,
    const uint64_t* const raw_integer_array) {
  // Write Meta-Data
  const uint64_t sa_bit_length = integer_upper_power_of_two(max_integer);
  const uint64_t array_size = PACKED_INT_ARRAY_SIZE(num_integers,sa_bit_length);
  packed_integer_array_write_metadata(file_manager,num_integers,sa_bit_length,array_size);
  // Write array
  uint64_t i, offset, bitmap_word = 0;
  for (i=0,offset=0;i<num_integers;++i) {
    const uint64_t integer = raw_integer_array[i];
    // Store LO part
    bitmap_word |= (integer << offset);
    // Store HI part
    offset += sa_bit_length;
    if (offset >= UINT64_LENGTH) {
      fm_write_uint64(file_manager,bitmap_word);
      bitmap_word = 0;
      offset -= UINT64_LENGTH;
      if (offset > 0) {
        bitmap_word |= (integer >> (sa_bit_length-offset));
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
GEM_INLINE void packed_integer_array_delete(packed_integer_array_t* const array) {
  PACKED_INTEGER_ARRAY_CHECK(array);
  if (array->mm_bitmap!=NULL) {
    if (array->mm_bitmap==(mm_t*)(-1)) {
      mm_free(array->bitmap);
    } else {
      mm_bulk_free(array->mm_bitmap);
    }
  }
  mm_free(array);
}
/*
 * Accessors
 */
#define PACKED_INT_ARRAY_GET_LOCATION(lo64_bit_position,lo64_word_position,lo64_offset,hi64_offset) \
  const uint64_t lo64_bit_position = position*array->integer_length; \
  const uint64_t lo64_word_position = lo64_bit_position/UINT64_LENGTH; \
  const uint64_t lo64_offset = lo64_bit_position % UINT64_LENGTH; \
  const uint64_t hi64_offset = UINT64_LENGTH-lo64_offset
GEM_INLINE void packed_integer_array_prefetch(const packed_integer_array_t* const array,const uint64_t position) {
  PACKED_INTEGER_ARRAY_CHECK_NO_NULL(array);
  // Prefetch
  const uint64_t lo64_bit_position = position*array->integer_length;
  const uint64_t lo64_word_position = lo64_bit_position/UINT64_LENGTH;
  PREFETCH(array->bitmap+lo64_word_position);
}
GEM_INLINE uint64_t packed_integer_array_load(const packed_integer_array_t* const array,const uint64_t position) {
  PACKED_INTEGER_ARRAY_CHECK_NO_NULL(array);
  // Load LO & HI part
  PACKED_INT_ARRAY_GET_LOCATION(lo64_bit_position,lo64_word_position,lo64_offset,hi64_offset);
  if (hi64_offset >= array->integer_length) {
    return uint64_mask_ones[array->integer_length] & (array->bitmap[lo64_word_position] >> lo64_offset);
  } else {
    return uint64_mask_ones[array->integer_length] &
        ((array->bitmap[lo64_word_position] >> lo64_offset) | (array->bitmap[lo64_word_position+1] << hi64_offset));
  }
}
GEM_INLINE void packed_integer_array_store(packed_integer_array_t* const array,const uint64_t position,const uint64_t integer) {
  PACKED_INTEGER_ARRAY_CHECK_NO_NULL(array);
  // Locate
  PACKED_INT_ARRAY_GET_LOCATION(lo64_bit_position,lo64_word_position,lo64_offset,hi64_offset);
  // Store LO part
  array->bitmap[lo64_word_position] |= (integer << lo64_offset);
  // Store HI part
  if (hi64_offset < array->integer_length) {
    array->bitmap[lo64_word_position+1] |= (integer >> hi64_offset);
  }
}
GEM_INLINE uint64_t packed_integer_array_get_size(const packed_integer_array_t* const array) {
  PACKED_INTEGER_ARRAY_CHECK(array);
  return array->array_size;
}
GEM_INLINE uint64_t packed_integer_array_get_length(const packed_integer_array_t* const array) {
  PACKED_INTEGER_ARRAY_CHECK(array);
  return array->num_elements;
}
/*
 * Display
 */
GEM_INLINE void packed_integer_array_print(FILE* const stream,const packed_integer_array_t* const array,const bool display_data) {
  PACKED_INTEGER_ARRAY_CHECK(array);
  // Print meta-info
  tab_fprintf(stream,"[GEM]>PackedArray.info\n");
  tab_fprintf(stream,"  => PackedArray.Size %lu MB\n",CONVERT_B_TO_MB(array->array_size));
  tab_fprintf(stream,"  => PackedArray.Length %lu\n",array->num_elements);
  tab_fprintf(stream,"  => PackedArray.Integer %lu bits\n",array->integer_length);
  // Print data
  if (display_data) {
    tab_fprintf(stream,"  => PackedArray.Data\n");
    uint64_t i;
    for (i=0;i<array->num_elements;++i) {
      fprintf(stream,"%8lu",packed_integer_array_load(array,i));
      PRINT_COLUMN(stream,i,10,"\t");
    }
    PRINT_COLUMN_CLOSE(stream,i,10);
  }
  // Flush
  fflush(stream);
}

