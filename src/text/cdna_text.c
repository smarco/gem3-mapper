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

/*
 * Pragmas
 */
#ifdef __clang__
#pragma GCC diagnostic ignored "-Winitializer-overrides"
#endif

/*
 * Include
 */
#include "text/cdna_text.h"
#include "text/dna_text.h"

/*
 * Errors
 */
#define GEM_ERROR_CDNA_INDEX_OUT_OF_RANGE "Compact DNA-Text. Requested index (%"PRIu64") out of range [0,%"PRIu64")"
#define GEM_ERROR_CDNA_NOT_VALID_CHARACTER "Compact DNA-Text Builder. Invalid character provided ASCII='%d'"

/*
 * Blocks encoding dimensions
 */
#define CDNA_BLOCK_CHARS 21 /* UINT64_LENGTH/DNA_EXT_RANGE_BITS */

/*
 * CDNA Masks
 */
#define CDNA_TEXT_READ_MASK         0x7000000000000000ull
#define CDNA_TEXT_INVERSE_READ_MASK 0x0000000000000007ull
#define CDNA_TEXT_WRITE_PADDING     0x0000000000000007ull

/*
 * CDNA Block Masks
 */
const uint64_t cdna_text_block_mask_left[256] =
{
  [0 ... 255] = 0xFFFFFFFFFFFFFFFFull,
  [ 0] = 0x7FFFFFFFFFFFFFFFull,
  [ 1] = 0x0FFFFFFFFFFFFFFFull,
  [ 2] = 0x01FFFFFFFFFFFFFFull,
  [ 3] = 0x003FFFFFFFFFFFFFull,
  [ 4] = 0x0007FFFFFFFFFFFFull,
  [ 5] = 0x0000FFFFFFFFFFFFull,
  [ 6] = 0x00001FFFFFFFFFFFull,
  [ 7] = 0x000003FFFFFFFFFFull,
  [ 8] = 0x0000007FFFFFFFFFull,
  [ 9] = 0x0000000FFFFFFFFFull,
  [10] = 0x00000001FFFFFFFFull,
  [11] = 0x000000003FFFFFFFull,
  [12] = 0x0000000007FFFFFFull,
  [13] = 0x0000000000FFFFFFull,
  [14] = 0x00000000001FFFFFull,
  [15] = 0x000000000003FFFFull,
  [16] = 0x0000000000007FFFull,
  [17] = 0x0000000000000FFFull,
  [18] = 0x00000000000001FFull,
  [19] = 0x000000000000003Full,
  [20] = 0x0000000000000007ull,
};
const uint64_t cdna_text_block_mask_right[256] =
{
  [0 ... 255] = 0xFFFFFFFFFFFFFFFFull,
  [ 0] = 0x0000000000000000ull,
  [ 1] = 0x7000000000000000ull,
  [ 2] = 0x7E00000000000000ull,
  [ 3] = 0x7FC0000000000000ull,
  [ 4] = 0x7FF8000000000000ull,
  [ 5] = 0x7FFF000000000000ull,
  [ 6] = 0x7FFFE00000000000ull,
  [ 7] = 0x7FFFFC0000000000ull,
  [ 8] = 0x7FFFFF8000000000ull,
  [ 9] = 0x7FFFFFF000000000ull,
  [10] = 0x7FFFFFFE00000000ull,
  [11] = 0x7FFFFFFFC0000000ull,
  [12] = 0x7FFFFFFFF8000000ull,
  [13] = 0x7FFFFFFFFF000000ull,
  [14] = 0x7FFFFFFFFFE00000ull,
  [15] = 0x7FFFFFFFFFFC0000ull,
  [16] = 0x7FFFFFFFFFFF8000ull,
  [17] = 0x7FFFFFFFFFFFF000ull,
  [18] = 0x7FFFFFFFFFFFFE00ull,
  [19] = 0x7FFFFFFFFFFFFFC0ull,
  [20] = 0x7FFFFFFFFFFFFFF8ull,
};
/*
 * CDNA Text (Builder incorporated)
 */
cdna_text_t* cdna_text_new(mm_slab_t* const mm_slab) {
  // Allocate
  cdna_text_t* const cdna_text = mm_alloc(cdna_text_t);
  // Text */
  cdna_text->text = svector_new(mm_slab,uint64_t);
  cdna_text->text_length = 0;
  // Iterator */
  svector_iterator_new(&(cdna_text->text_iterator),cdna_text->text,SVECTOR_WRITE_ITERATOR,0);
  // Internals (Builder) */
  cdna_text->position_mod21 = 0;
  cdna_text->current_word = 0;
  cdna_text->mm_slab = mm_slab;
  // Return
  return cdna_text;
}
void cdna_text_delete(cdna_text_t* const cdna_text) {
  svector_delete(cdna_text->text); // Free text
  mm_free(cdna_text); // Free handler
}
void cdna_text_flush(cdna_text_t* const cdna_text) {
  // Write Bitmap
  *svector_iterator_get_element(&(cdna_text->text_iterator),uint64_t) = cdna_text->current_word;
  svector_write_iterator_next(&(cdna_text->text_iterator));
  // Clear
  cdna_text->position_mod21 = 0;
  cdna_text->current_word = 0;
}
void cdna_text_add_char(cdna_text_t* const cdna_text,const uint8_t enc_char) {
  gem_fatal_check(!is_extended_dna_encoded(enc_char),CDNA_NOT_VALID_CHARACTER,enc_char);
  // Write new character
  cdna_text->current_word |= (uint64_t) enc_char;
  // Inc position
  ++(cdna_text->position_mod21);
  ++(cdna_text->text_length);
  // Dump (if needed)
  if (gem_expect_false(cdna_text->position_mod21==CDNA_BLOCK_CHARS)) {
    cdna_text_flush(cdna_text);
  } else {
    // Shift to allocate next character
    cdna_text->current_word <<= DNA_EXT_RANGE_BITS;
  }
}
void cdna_text_close(cdna_text_t* const cdna_text) {
  if (cdna_text->position_mod21 > 0) {
    // Pad the rest of the word
    cdna_text->current_word |= CDNA_TEXT_WRITE_PADDING;
    ++(cdna_text->position_mod21);
    while (cdna_text->position_mod21!=CDNA_BLOCK_CHARS) {
      cdna_text->current_word <<= DNA_EXT_RANGE_BITS;
      cdna_text->current_word |= CDNA_TEXT_WRITE_PADDING;
      ++(cdna_text->position_mod21);
    }
    // Flush & close
    cdna_text_flush(cdna_text);
  }
}
void cdna_text_write_as_bitwise_text(fm_t* const file_manager,cdna_text_t* const cdna_text) {
  cdna_bitwise_text_builder_t* const cdna_bitwise_text =
      cdna_bitwise_text_builder_new(file_manager,cdna_text->text_length,cdna_text->mm_slab);
  // Index Text Iterator
  cdna_text_iterator_t iterator;
  cdna_text_iterator_init(&iterator,cdna_text,0);
  while (!cdna_text_iterator_eoi(&iterator)) {
    cdna_bitwise_text_builder_add_char(cdna_bitwise_text,
        cdna_text_iterator_get_char_encoded(&iterator)); // Store Indexed-Text
    cdna_text_iterator_next_char(&iterator);
  }
  cdna_bitwise_text_builder_close(cdna_bitwise_text); // Finish writing Indexed-Text
  cdna_bitwise_text_builder_delete(cdna_bitwise_text); // Free
}
/*
 * CDNA Accessors
 */
uint64_t cdna_text_get_length(cdna_text_t* const cdna_text) {
  return cdna_text->text_length;
}
/*
 * CDNA Text Iterator
 */
void cdna_text_iterator_seek(cdna_text_iterator_t* const iterator,const uint64_t position) {
  gem_fatal_check(position>=iterator->cdna_text->text_length,CDNA_INDEX_OUT_OF_RANGE,position,iterator->cdna_text->text_length);
  // Locate
  iterator->position = position;
  iterator->position_mod21 = position%CDNA_BLOCK_CHARS;
  const uint64_t word_position = position/CDNA_BLOCK_CHARS;
  svector_read_iterator_seek(&(iterator->text_iterator),word_position);
  iterator->current_word = *svector_iterator_get_element(&(iterator->text_iterator),uint64_t);
  // Seek to the proper character
  if (gem_expect_true(iterator->position_mod21 > 0)) {
    iterator->current_word <<= (DNA_EXT_RANGE_BITS*iterator->position_mod21);
  }
}
void cdna_text_iterator_init(
    cdna_text_iterator_t* const iterator,
    cdna_text_t* const cdna_text,
    const uint64_t position) {
  // Text
  iterator->cdna_text = cdna_text;
  if (gem_expect_false(position>=iterator->cdna_text->text_length)) {
    iterator->eoi = true;
  } else {
    iterator->eoi = false;
    // Locate
    iterator->position = position;
    iterator->position_mod21 = position%CDNA_BLOCK_CHARS;
    const uint64_t word_position = position/CDNA_BLOCK_CHARS;
    svector_iterator_new(&(iterator->text_iterator),cdna_text->text,SVECTOR_READ_ITERATOR,word_position);
    iterator->current_word = *svector_iterator_get_element(&(iterator->text_iterator),uint64_t);
    // Seek to the proper character
    if (gem_expect_true(iterator->position_mod21 > 0)) {
      iterator->current_word <<= (DNA_EXT_RANGE_BITS*iterator->position_mod21);
    }
  }
}
bool cdna_text_iterator_eoi(cdna_text_iterator_t* const iterator) {
  return iterator->eoi;
}
uint8_t cdna_text_iterator_get_char_encoded(cdna_text_iterator_t* const iterator) {
  return (uint8_t) ((iterator->current_word & CDNA_TEXT_READ_MASK) >> 60);
}
void cdna_text_iterator_next_char(cdna_text_iterator_t* const iterator) {
  ++(iterator->position);
  ++(iterator->position_mod21);
  if (gem_expect_false(iterator->position>=iterator->cdna_text->text_length)) {
    iterator->eoi = true;
  } else {
    if (iterator->position_mod21==CDNA_BLOCK_CHARS) {
      iterator->position_mod21 = 0;
      svector_read_iterator_next(&(iterator->text_iterator));
      if (gem_expect_false(svector_read_iterator_eoi(&(iterator->text_iterator)))) {
        iterator->eoi = true;
      } else {
        iterator->current_word = *svector_iterator_get_element(&(iterator->text_iterator),uint64_t);
      }
    } else {
      iterator->current_word <<= DNA_EXT_RANGE_BITS;
    }
  }
}
/*
 * CDNA Text Reverse Iterator
 */
void cdna_text_reverse_iterator_init(
    cdna_text_iterator_t* const iterator,
    cdna_text_t* const cdna_text,
    const uint64_t position) {
  iterator->cdna_text = cdna_text;
  const uint64_t text_length = cdna_text->text_length;
  if (gem_expect_false(position >= text_length)) {
    iterator->eoi = true;
  } else {
    iterator->eoi = false;
    // Locate
    const uint64_t position_mod21 = position % CDNA_BLOCK_CHARS;
    const uint64_t position_block = position / CDNA_BLOCK_CHARS;
    const uint64_t text_length_block = (text_length-1) / CDNA_BLOCK_CHARS;
    iterator->position = position_block; // NOTE: @position denotes number of 64b-bitmaps left
    // Fetch bitmap
    if (text_length_block == position_block && cdna_text->position_mod21 > 0) {
      // Load cached bitmap (Adjusted)
      iterator->current_word = cdna_text->current_word >>
                               (DNA_EXT_RANGE_BITS*(cdna_text->position_mod21-position_mod21));
    } else {
      // Load stored bitmap (Adjusted)
      iterator->current_word = *svector_get_element(cdna_text->text,iterator->position,uint64_t) >>
                               (DNA_EXT_RANGE_BITS*(CDNA_BLOCK_CHARS-position_mod21-1));
    }
    // Number of characters left buffer
    iterator->position_mod21 = position_mod21;
  }
}
bool cdna_text_reverse_iterator_eoi(cdna_text_iterator_t* const iterator) {
  return iterator->eoi;
}
uint8_t cdna_text_reverse_iterator_get_char_encoded(cdna_text_iterator_t* const iterator) {
  return (uint8_t) (iterator->current_word & CDNA_TEXT_INVERSE_READ_MASK);
}
void cdna_text_reverse_iterator_next_char(cdna_text_iterator_t* const iterator) {
  if (gem_expect_false(iterator->position_mod21==0)) {
    if (iterator->position==0) {
      iterator->eoi = true;
    } else {
      --(iterator->position);
      iterator->current_word = *svector_get_element(iterator->cdna_text->text,iterator->position,uint64_t);
      iterator->position_mod21 = CDNA_BLOCK_CHARS-1;
    }
  } else {
    --(iterator->position_mod21);
    iterator->current_word >>= DNA_EXT_RANGE_BITS;
  }
}
/* CDNA Text Compare
 *   return (int) -1 if (a < b)
 *   return (int)  0 if (a == b)
 *   return (int)  1 if (a > b)
 */
void cdna_text_block_print(FILE* const stream,uint64_t word_block) {
  uint64_t i;
  for (i=0;i<CDNA_BLOCK_CHARS;++i) {
    fprintf(stream,"%c",dna_decode((uint8_t) ((word_block & CDNA_TEXT_READ_MASK) >> 60)));
    word_block <<= DNA_EXT_RANGE_BITS;
  }
}
void cdna_text_block_iterator_init(
    cdna_text_iterator_t* const iterator,
    cdna_text_t* const cdna_text,
    const uint64_t position) {
  // Text
  iterator->cdna_text = cdna_text;
  // Locate
  iterator->position = position;
  iterator->position_mod21 = position%CDNA_BLOCK_CHARS;
  const uint64_t word_position = position/CDNA_BLOCK_CHARS;
  svector_iterator_new(&(iterator->text_iterator),cdna_text->text,SVECTOR_READ_ITERATOR,word_position);
}
void cdna_text_block_iterator_next_block(cdna_text_iterator_t* const iterator) {
  svector_read_iterator_next(&(iterator->text_iterator));
}
uint64_t cdna_text_block_iterator_get_block(cdna_text_iterator_t* const iterator) {
  if (gem_expect_true(!svector_read_iterator_eoi(&(iterator->text_iterator)))) {
    return *svector_iterator_get_element(&(iterator->text_iterator),uint64_t);
  } else {
    cdna_text_iterator_seek(iterator,0);
    return *svector_iterator_get_element(&(iterator->text_iterator),uint64_t);
  }
}
#define cdna_text_cmp_load_block(iterator,mod,left,right,word,next_word,block) { \
  word = next_word; \
  cdna_text_block_iterator_next_block(iterator); \
  next_word = cdna_text_block_iterator_get_block(iterator); \
  const uint64_t mask = cdna_text_block_mask_left[mod]; \
  block = ((word & mask) << (left)) | ((next_word & (~mask)) >> (right)) ; \
}
int cdna_text_block_cmp(cdna_text_iterator_t* const iterator_a,cdna_text_iterator_t* const iterator_b) {
  // Compare successive blocks
  const uint64_t mod_a = iterator_a->position_mod21;
  const uint64_t mod_b = iterator_b->position_mod21;
  const uint64_t left_a = DNA_EXT_RANGE_BITS*mod_a;
  const uint64_t left_b = DNA_EXT_RANGE_BITS*mod_b;
  const uint64_t right_a = UINT64_LENGTH - 1 - left_a;
  const uint64_t right_b = UINT64_LENGTH - 1 - left_b;
  uint64_t word_a, word_b, block_a, block_b;
  uint64_t next_word_a = cdna_text_block_iterator_get_block(iterator_a);
  uint64_t next_word_b = cdna_text_block_iterator_get_block(iterator_b);
  do {
    cdna_text_cmp_load_block(iterator_a,mod_a,left_a,right_a,word_a,next_word_a,block_a);
    cdna_text_cmp_load_block(iterator_b,mod_b,left_b,right_b,word_b,next_word_b,block_b);
  } while (block_a==block_b);
  // Return
  return (block_a < block_b) ? -1 : 1;
}
/*
 * CDNA Text Character Iterator-HUB
 */
cdna_text_iterator_hub_t* cdna_text_iterator_hub_new(
    cdna_text_t** cdna_texts,
    const uint64_t num_texts,
    const uint64_t starting_position) {
  // Alloc
  cdna_text_iterator_hub_t* const iterator_hub = mm_alloc(cdna_text_iterator_hub_t);
  // Text(s)
  iterator_hub->num_texts = num_texts;
  // Set limits & iterators
  iterator_hub->limits = mm_calloc(num_texts,uint64_t,true);
  iterator_hub->iterators = mm_calloc(num_texts,cdna_text_iterator_t,true);
  uint64_t i, acc_offset;
  iterator_hub->iterator_num = UINT64_MAX;
  for (i=0,acc_offset=0;i<num_texts;++i) {
    const uint64_t new_acc_offset = acc_offset + cdna_texts[i]->text_length;
    if (acc_offset <= starting_position && starting_position < new_acc_offset) {
      iterator_hub->iterator_num = i;
      cdna_text_iterator_init(iterator_hub->iterators+i,cdna_texts[i],starting_position);
    } else {
      cdna_text_iterator_init(iterator_hub->iterators+i,cdna_texts[i],0);
    }
    acc_offset = new_acc_offset;
    iterator_hub->limits[i] = acc_offset;
  }
  // Current State
  iterator_hub->text_position = starting_position;
  iterator_hub->eoi = (iterator_hub->iterator_num==UINT64_MAX);
  // Return
  return iterator_hub;
}
void cdna_text_iterator_hub_delete(cdna_text_iterator_hub_t* const iterator_hub) {
  mm_free(iterator_hub->limits);
  mm_free(iterator_hub->iterators);
  mm_free(iterator_hub);
}
bool cdna_text_iterator_hub_eoi(cdna_text_iterator_hub_t* const iterator_hub) {
  return iterator_hub->eoi;
}
uint8_t cdna_text_iterator_hub_get_char_encoded(cdna_text_iterator_hub_t* const iterator_hub) {
  return cdna_text_iterator_get_char_encoded(iterator_hub->iterators+iterator_hub->iterator_num);
}
void cdna_text_iterator_hub_next_char(cdna_text_iterator_hub_t* const iterator_hub) {
  ++(iterator_hub->text_position);
  if (gem_expect_false(iterator_hub->text_position >= iterator_hub->limits[iterator_hub->iterator_num])) {
    ++(iterator_hub->iterator_num);
    while (gem_expect_false(iterator_hub->iterator_num < iterator_hub->num_texts) &&
        cdna_text_iterator_eoi(iterator_hub->iterators+iterator_hub->iterator_num)) {
        ++(iterator_hub->iterator_num);
    }
    if (gem_expect_false(iterator_hub->iterator_num >= iterator_hub->num_texts)) {
      iterator_hub->eoi = true;
    }
  } else {
    cdna_text_iterator_next_char(iterator_hub->iterators+iterator_hub->iterator_num);
  }
}
/*
 * Display
 */
void cdna_text_print(FILE* const stream,cdna_text_t* const cdna_text) {
  fprintf(stream,"[GEM]>Compacted DNA-text\n");
  fprintf(stream,"  => Architecture 3b.1bm64.21c\n");
  fprintf(stream,"  => Text.Length %"PRIu64"\n",cdna_text->text_length);
  const uint64_t cdna_text_size = DIV_CEIL(cdna_text->text_length,CDNA_BLOCK_CHARS);
  fprintf(stream,"  => Text.Size %"PRIu64" MB\n",CONVERT_B_TO_MB(cdna_text_size*UINT64_SIZE));
  // Flush
  fflush(stream);
}
