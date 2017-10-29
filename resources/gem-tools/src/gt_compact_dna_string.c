/*
 * PROJECT: GEM-Tools library
 * FILE: gt_compact_dna_string.c
 * DATE: 20/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Bitmap (compact) representation of DNA-strings (using 8 characters alphabet)
 */

#include "gt_compact_dna_string.h"

/*
 * CDNA Bitmaps Encoding
 */
#define GT_CDNA_ENCODED_CHAR_BM0 ((uint8_t)1)
#define GT_CDNA_ENCODED_CHAR_BM1 ((uint8_t)2)
#define GT_CDNA_ENCODED_CHAR_BM2 ((uint8_t)4)

const uint8_t gt_cdna_encoded_char_bm[3] = {
  GT_CDNA_ENCODED_CHAR_BM0, GT_CDNA_ENCODED_CHAR_BM1, GT_CDNA_ENCODED_CHAR_BM2
};

/*
 * CDNA translation tables
 */
const char gt_cdna_decode[8] = {
  GT_DNA_CHAR_A, GT_DNA_CHAR_C, GT_DNA_CHAR_G, GT_DNA_CHAR_T,
  GT_DNA_CHAR_N, GT_DNA_CHAR_N, GT_DNA_CHAR_N, GT_DNA_CHAR_N
};
const uint8_t gt_cdna_encode[256] = {
    [0 ... 255] = GT_CDNA_ENC_CHAR_N,
    ['A'] = GT_CDNA_ENC_CHAR_A,['C'] = GT_CDNA_ENC_CHAR_C,['G'] = GT_CDNA_ENC_CHAR_G,['T'] = GT_CDNA_ENC_CHAR_T,
    ['a'] = GT_CDNA_ENC_CHAR_A,['c'] = GT_CDNA_ENC_CHAR_C,['g'] = GT_CDNA_ENC_CHAR_G,['t'] = GT_CDNA_ENC_CHAR_T,
};

/*
 * Block locator functions
 */
#define GT_CDNA_GET_NUM_BLOCKS(num_chars) ((num_chars+(GT_CDNA_BLOCK_CHARS-1))/GT_CDNA_BLOCK_CHARS)
#define GT_CDNA_GET_NUM_CHARS(num_blocks) (GT_CDNA_BLOCK_CHARS*num_blocks)
#define GT_CDNA_GET_BLOCKS_MEM(num_blocks) (GT_CDNA_BLOCK_SIZE*num_blocks)

/*
 * CDNA Internal Bitmap Handlers
 */
#define GT_CDNA_GET_BLOCK_POS(global_pos,block_num,block_pos) \
  block_num = global_pos/GT_CDNA_BLOCK_CHARS; \
  block_pos = global_pos % GT_CDNA_BLOCK_CHARS

#define GT_CDNA_INIT_BLOCK(block_mem) \
  ((uint64_t*)block_mem)[0]=UINT64_ZEROS; \
  ((uint64_t*)block_mem)[1]=UINT64_ZEROS; \
  ((uint64_t*)block_mem)[2]=UINT64_ONES

#define GT_CDNA_GET_MEM_BLOCK(bitmaps_mem,block_num) ((uint64_t*)bitmaps_mem)+((block_num)*GT_CDNA_BLOCK_BITMAPS)

#define GT_CDNA_LOAD_BLOCKS(block_mem,bm_0,bm_1,bm_2) { \
  bm_0 = *(block_mem); \
  bm_1 = *(block_mem+1); \
  bm_2 = *(block_mem+2); \
}

#define GT_CDNA_GET_BLOCKS(bitmaps_mem,block_num,bm_0,bm_1,bm_2) { \
  uint64_t* const block_mem = GT_CDNA_GET_MEM_BLOCK(bitmaps_mem,block_num); \
  GT_CDNA_LOAD_BLOCKS(block_mem,bm_0,bm_1,bm_2); \
}

#define GT_CDNA_SHIFT_FORWARD_CHARS(block_pos,bm_0,bm_1,bm_2) \
  bm_0 >>= block_pos; bm_1 >>= block_pos; bm_2 >>= block_pos

#define GT_CDNA_SHIFT_BACKWARD_CHARS(block_pos,bm_0,bm_1,bm_2) \
  bm_0 <<= block_pos; bm_1 <<= block_pos; bm_2 <<= block_pos

#define GT_CDNA_EXTRACT_CHAR(bm_0,bm_1,bm_2) \
  (((bm_0&GT_CDNA_EXTRACT_MASK))    |        \
   ((bm_1&GT_CDNA_EXTRACT_MASK)<<1) |        \
   ((bm_2&GT_CDNA_EXTRACT_MASK)<<2))

#define GT_CDNA_EXTRACT_LAST_CHAR(bm_0,bm_1,bm_2) \
  (((bm_0&GT_CDNA_EXTRACT_LAST_MASK))    |        \
   ((bm_1&GT_CDNA_EXTRACT_LAST_MASK)<<1) |        \
   ((bm_2&GT_CDNA_EXTRACT_LAST_MASK)<<2))

#define GT_CDNA_PROYECT_CHAR(block_mem,enc_char,bm_pos,bm_one_mask,bm_zero_mask) \
  if ((enc_char) & gt_cdna_encoded_char_bm[bm_pos]) { \
    block_mem[bm_pos] |= bm_one_mask; \
  } else { \
    block_mem[bm_pos] &= bm_zero_mask; \
  }

#define GT_CDNA_SET_CHAR(block_mem,block_pos,enc_char)  \
  const uint64_t bm_one_mask = GT_CDNA_ONE_MASK<<(block_pos); \
  const uint64_t bm_zero_mask = ~(bm_one_mask); \
  GT_CDNA_PROYECT_CHAR(block_mem,enc_char,0,bm_one_mask,bm_zero_mask); \
  GT_CDNA_PROYECT_CHAR(block_mem,enc_char,1,bm_one_mask,bm_zero_mask); \
  GT_CDNA_PROYECT_CHAR(block_mem,enc_char,2,bm_one_mask,bm_zero_mask)

/*
 * Constructor
 */
GT_INLINE gt_compact_dna_string* gt_cdna_string_new(const uint64_t initial_chars) {
  gt_compact_dna_string* cdna_string = gt_alloc(gt_compact_dna_string);
  const uint64_t initial_blocks = initial_chars ? GT_CDNA_GET_NUM_BLOCKS(initial_chars) : 1;
  cdna_string->bitmaps = gt_malloc(GT_CDNA_GET_BLOCKS_MEM(initial_blocks));
  cdna_string->allocated = GT_CDNA_GET_NUM_CHARS(initial_blocks);
  cdna_string->length = 0;
  GT_CDNA_INIT_BLOCK(cdna_string->bitmaps); // Init 0-block
  return cdna_string;
}
GT_INLINE void gt_cdna_string_resize(gt_compact_dna_string* const cdna_string,const uint64_t num_chars) {
  GT_COMPACT_DNA_STRING_CHECK(cdna_string);
  if (num_chars > cdna_string->allocated) {
    const uint64_t num_blocks = GT_CDNA_GET_NUM_BLOCKS(num_chars);
    cdna_string->bitmaps=realloc(cdna_string->bitmaps,GT_CDNA_GET_BLOCKS_MEM(num_blocks));
    cdna_string->allocated = GT_CDNA_GET_NUM_CHARS(num_blocks);
  }
}
GT_INLINE void gt_cdna_string_clear(gt_compact_dna_string* const cdna_string) {
  GT_COMPACT_DNA_STRING_CHECK(cdna_string);
  cdna_string->length = 0;
  GT_CDNA_INIT_BLOCK(cdna_string->bitmaps); // Init 0-block
}
GT_INLINE void gt_cdna_string_delete(gt_compact_dna_string* const cdna_string) {
  GT_COMPACT_DNA_STRING_CHECK(cdna_string);
  gt_free(cdna_string->bitmaps);
  gt_free(cdna_string);
}
/*
 * Handlers
 */
GT_INLINE void gt_cdna_allocate__init_blocks(gt_compact_dna_string* const cdna_string,const uint64_t pos) {
  GT_COMPACT_DNA_STRING_CHECK(cdna_string);
  if (pos >= cdna_string->length) {
    // Check allocated blocks
    if (pos >= cdna_string->allocated) gt_cdna_string_resize(cdna_string,pos+(GT_CDNA_BLOCK_CHARS*10));
    // Initialize new accessed blocks
    const uint64_t next_block_num = gt_expect_true(cdna_string->length>0) ? ((cdna_string->length-1)/GT_CDNA_BLOCK_CHARS)+1 : 1;
    const uint64_t top_block_num = pos/GT_CDNA_BLOCK_CHARS;
    if (next_block_num <= top_block_num) {
      uint64_t i;
      uint64_t* block_mem = GT_CDNA_GET_MEM_BLOCK(cdna_string->bitmaps,next_block_num);
      for (i=next_block_num; i<=top_block_num; ++i) {
        GT_CDNA_INIT_BLOCK(block_mem);
        block_mem+=GT_CDNA_BLOCK_BITMAPS;
      }
    }
    // Update total length
    cdna_string->length = pos+1;
  }
}
GT_INLINE char gt_cdna_string_get_char_at(gt_compact_dna_string* const cdna_string,const uint64_t position) {
  GT_COMPACT_DNA_STRING_CHECK(cdna_string);
  GT_COMPACT_DNA_STRING_POSITION_CHECK(cdna_string,position);
  uint64_t block_num, block_pos;
  uint64_t bm_0, bm_1, bm_2;
  GT_CDNA_GET_BLOCK_POS(position,block_num,block_pos);
  GT_CDNA_GET_BLOCKS(cdna_string->bitmaps,block_num,bm_0,bm_1,bm_2);
  GT_CDNA_SHIFT_FORWARD_CHARS(block_pos,bm_0,bm_1,bm_2);
  return gt_cdna_decode[GT_CDNA_EXTRACT_CHAR(bm_0,bm_1,bm_2)];
}
GT_INLINE void gt_cdna_string_set_char_at(gt_compact_dna_string* const cdna_string,const uint64_t position,const char character) {
  GT_COMPACT_DNA_STRING_CHECK(cdna_string);
  // Check allocated bitmaps
  gt_cdna_allocate__init_blocks(cdna_string,position);
  // Encode char
  uint64_t block_num, block_pos;
  GT_CDNA_GET_BLOCK_POS(position,block_num,block_pos);
  uint64_t* const block_mem = GT_CDNA_GET_MEM_BLOCK(cdna_string->bitmaps,block_num);
  const uint8_t enc_char = gt_cdna_encode(character);
  GT_CDNA_SET_CHAR(block_mem,block_pos,enc_char);
}
GT_INLINE uint64_t gt_cdna_string_get_length(gt_compact_dna_string* const cdna_string) {
  GT_COMPACT_DNA_STRING_CHECK(cdna_string);
  return cdna_string->length;
}
GT_INLINE void gt_cdna_string_append_string(gt_compact_dna_string* const cdna_string,const char* const string,const uint64_t length) {
  GT_COMPACT_DNA_STRING_CHECK(cdna_string);
  // Check allocated bitmaps
  const uint64_t total_chars = cdna_string->length+length-1;
  if (total_chars >= cdna_string->allocated) {
    gt_cdna_string_resize(cdna_string,total_chars);
  }
  // Copy string
  uint64_t block_num, block_pos, i;
  GT_CDNA_GET_BLOCK_POS(cdna_string->length,block_num,block_pos);
  uint64_t* block_mem = GT_CDNA_GET_MEM_BLOCK(cdna_string->bitmaps,block_num);
  for (i=0; i<length; ++i,++block_pos) {
    if (gt_expect_false(block_pos==GT_CDNA_BLOCK_CHARS)) {
      block_pos=0;
      block_mem+=GT_CDNA_BLOCK_BITMAPS;
    }
    const uint8_t enc_char = gt_cdna_encode(string[i]);
    GT_CDNA_SET_CHAR(block_mem,block_pos,enc_char);
  }
  // Update total length
  cdna_string->length = total_chars+1;
}
/*
 * Compact DNA String Sequence Iterator
 */
GT_INLINE void gt_cdna_string_new_iterator(
    gt_compact_dna_string* const cdna_string,const uint64_t position,gt_string_traversal const direction,
    gt_compact_dna_string_iterator* const cdna_string_iterator) {
  GT_COMPACT_DNA_STRING_CHECK(cdna_string);
  GT_NULL_CHECK(cdna_string_iterator);
  // Initialize the iterator
  cdna_string_iterator->cdna_string = cdna_string;
  // Seek to the current position
  if (position<cdna_string->length) {
    gt_cdna_string_iterator_seek(cdna_string_iterator,position,direction);
  } else {
    cdna_string_iterator->current_pos=cdna_string->length;
  }
}
GT_INLINE void gt_cdna_string_iterator_seek(
    gt_compact_dna_string_iterator* const cdna_string_iterator,
    const uint64_t position,gt_string_traversal const direction) {
  GT_NULL_CHECK(cdna_string_iterator);
  GT_NULL_CHECK(cdna_string_iterator->cdna_string);
  gt_compact_dna_string* const cdna_string = cdna_string_iterator->cdna_string;
  GT_COMPACT_DNA_STRING_POSITION_CHECK(cdna_string,position);
  // Set the iterator
  cdna_string_iterator->current_pos = position;
  cdna_string_iterator->current_pos_mod = position%GT_CDNA_BLOCK_CHARS;
  cdna_string_iterator->direction = direction;
  // Set the current bitmap and seek to current position
  uint64_t block_num, block_pos;
  GT_CDNA_GET_BLOCK_POS(position,block_num,block_pos);
  cdna_string_iterator->current_bitmap = GT_CDNA_GET_MEM_BLOCK(cdna_string->bitmaps,block_num);
  GT_CDNA_GET_BLOCKS(cdna_string->bitmaps,block_num,cdna_string_iterator->bm_0,cdna_string_iterator->bm_1,cdna_string_iterator->bm_2);
  if (cdna_string_iterator->direction==GT_ST_FORWARD) {
    GT_CDNA_SHIFT_FORWARD_CHARS(block_pos,
        cdna_string_iterator->bm_0,cdna_string_iterator->bm_1,cdna_string_iterator->bm_2);
  } else if (block_pos < GT_CDNA_BLOCK_CHARS-1) { // GT_ST_BACKWARD
    GT_CDNA_SHIFT_BACKWARD_CHARS((GT_CDNA_BLOCK_CHARS-1-block_pos),
        cdna_string_iterator->bm_0,cdna_string_iterator->bm_1,cdna_string_iterator->bm_2);
  }
}
GT_INLINE bool gt_cdna_string_iterator_eos(gt_compact_dna_string_iterator* const cdna_string_iterator) {
  GT_COMPACT_DNA_STRING_ITERATOR_CHECK(cdna_string_iterator);
  return cdna_string_iterator->current_pos>=cdna_string_iterator->cdna_string->length;
}
GT_INLINE char gt_cdna_string_iterator_following(gt_compact_dna_string_iterator* const cdna_string_iterator) {
  GT_COMPACT_DNA_STRING_ITERATOR_CHECK(cdna_string_iterator);
  GT_COMPACT_DNA_STRING_POSITION_CHECK(cdna_string_iterator->cdna_string,cdna_string_iterator->current_pos);
  // Extract character
  const char character = gt_cdna_decode[GT_CDNA_EXTRACT_CHAR(cdna_string_iterator->bm_0,cdna_string_iterator->bm_1,cdna_string_iterator->bm_2)];
  // Update position
  ++cdna_string_iterator->current_pos;
  // Seek to proper position (load block if necessary)
  if (cdna_string_iterator->current_pos_mod==GT_CDNA_BLOCK_CHARS-1) {
    cdna_string_iterator->current_pos_mod = 0;
    cdna_string_iterator->current_bitmap += GT_CDNA_BLOCK_BITMAPS;
    GT_CDNA_LOAD_BLOCKS(cdna_string_iterator->current_bitmap,cdna_string_iterator->bm_0,cdna_string_iterator->bm_1,cdna_string_iterator->bm_2);
  } else {
    ++cdna_string_iterator->current_pos_mod;
    GT_CDNA_SHIFT_FORWARD_CHARS(1,cdna_string_iterator->bm_0,cdna_string_iterator->bm_1,cdna_string_iterator->bm_2);
  }
  // Return the character decoded
  return character;
}
GT_INLINE char gt_cdna_string_iterator_previous(gt_compact_dna_string_iterator* const cdna_string_iterator) {
  GT_COMPACT_DNA_STRING_ITERATOR_CHECK(cdna_string_iterator);
  GT_COMPACT_DNA_STRING_POSITION_CHECK(cdna_string_iterator->cdna_string,cdna_string_iterator->current_pos);
  // Extract character
  const char character = gt_cdna_decode[GT_CDNA_EXTRACT_LAST_CHAR(cdna_string_iterator->bm_0,cdna_string_iterator->bm_1,cdna_string_iterator->bm_2)];
  // Update position
  --cdna_string_iterator->current_pos;
  // Seek to proper position (load block if necessary)
  if (cdna_string_iterator->current_pos_mod==0) {
    cdna_string_iterator->current_pos_mod = GT_CDNA_BLOCK_CHARS-1;
    cdna_string_iterator->current_bitmap -= GT_CDNA_BLOCK_BITMAPS;
    GT_CDNA_LOAD_BLOCKS(cdna_string_iterator->current_bitmap,cdna_string_iterator->bm_0,cdna_string_iterator->bm_1,cdna_string_iterator->bm_2);
  } else {
    --cdna_string_iterator->current_pos_mod;
    GT_CDNA_SHIFT_BACKWARD_CHARS(1,cdna_string_iterator->bm_0,cdna_string_iterator->bm_1,cdna_string_iterator->bm_2);
  }
  // Return the character decoded
  return character;
}
GT_INLINE char gt_cdna_string_iterator_next(gt_compact_dna_string_iterator* const cdna_string_iterator) {
  GT_COMPACT_DNA_STRING_ITERATOR_CHECK(cdna_string_iterator);
  return (cdna_string_iterator->direction==GT_ST_FORWARD) ?
      gt_cdna_string_iterator_following(cdna_string_iterator) :
      gt_cdna_string_iterator_previous(cdna_string_iterator);
}

