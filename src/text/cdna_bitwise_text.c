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

#include "text/cdna_bitwise_text.h"
#include "text/dna_text.h"

/*
 * Blocks encoding dimensions
 */
#define CDNA_BITWISE_BLOCK_CHARS UINT64_LENGTH

/*
 * CDNA Bitwise Masks
 */
#define CDNA_BITWISE_WRITE_LAST_MASK UINT64_ONE_LAST_MASK

/*
 * CDNA Bitwise Text (Loader/Setup)
 */
cdna_bitwise_text_t* cdna_bitwise_text_read(fm_t* const file_manager) {
  // Allocate
  cdna_bitwise_text_t* const cdna_text = mm_alloc(cdna_bitwise_text_t);
  // Read header
  cdna_text->text_length = fm_read_uint64(file_manager);
  cdna_text->text_size = fm_read_uint64(file_manager);
  // Read Text
  fm_skip_align_4KB(file_manager); // Align to 4KB-MP
  cdna_text->mm_text = fm_load_mem(file_manager,cdna_text->text_size);
  cdna_text->text = mm_get_base_mem(cdna_text->mm_text);
  // Read Sparse 3rd-layer
  cdna_text->sparse_bitmap = sparse_bitmap_read(file_manager);
  // Return
  return cdna_text;
}
cdna_bitwise_text_t* cdna_bitwise_text_read_mem(mm_t* const memory_manager) {
  // Allocate
  cdna_bitwise_text_t* const cdna_text = mm_alloc(cdna_bitwise_text_t);
  // Read header
  cdna_text->text_length = mm_read_uint64(memory_manager);
  cdna_text->text_size = mm_read_uint64(memory_manager);
  // Read Text
  mm_skip_align_4KB(memory_manager); // Align to 4KB-MP
  cdna_text->mm_text = NULL;
  cdna_text->text = mm_read_mem(memory_manager,cdna_text->text_size);
  // Read Sparse 3rd-layer
  cdna_text->sparse_bitmap = sparse_bitmap_read_mem(memory_manager);
  // Return
  return cdna_text;
}
void cdna_bitwise_text_delete(cdna_bitwise_text_t* const cdna_text) {
  // Free Text
  if (cdna_text->mm_text!=NULL) mm_bulk_free(cdna_text->mm_text);
  // Free sparse text
  sparse_bitmap_delete(cdna_text->sparse_bitmap);
  // Free handler
  mm_free(cdna_text);
}
/*
 * CDNA Bitwise Text (Builder)
 */
cdna_bitwise_text_builder_t* cdna_bitwise_text_builder_new(
    fm_t* const file_manager,
    const uint64_t text_length,
    mm_slab_t* const mm_slab) {
  // Allocate handler
  cdna_bitwise_text_builder_t* const cdna_text_fm = mm_alloc(cdna_bitwise_text_builder_t);
  // Initialize file manager
  cdna_text_fm->file_manager = file_manager;
  // Initialize sparseBM-builder
  cdna_text_fm->sparse_bitmap_builder = sparse_bitmap_builder_new(mm_slab);
  // Initialize
  cdna_text_fm->layer_0 = 0;
  cdna_text_fm->layer_1 = 0;
  cdna_text_fm->layer_2 = 0;
  cdna_text_fm->position_mod64 = 0;
  // Write Header
  fm_write_uint64(cdna_text_fm->file_manager,text_length); // Text Length
  const uint64_t text_size = DIV_CEIL(text_length,UINT64_LENGTH)*2*UINT64_SIZE;
  fm_write_uint64(cdna_text_fm->file_manager,text_size); // Text Size
  // Align to 4KB-MemPage
  fm_skip_align_4KB(cdna_text_fm->file_manager);
  // Return
  return cdna_text_fm;
}
void cdna_bitwise_text_builder_flush(cdna_bitwise_text_builder_t* const cdna_text) {
  // Write Bitmaps
  fm_write_uint64(cdna_text->file_manager,cdna_text->layer_0);
  fm_write_uint64(cdna_text->file_manager,cdna_text->layer_1);
  // Conditionally add the third-layer (Sparse bitmap storage)
  if (gem_expect_true(cdna_text->layer_2==0)) {
    sparse_bitmap_builder_skip_bitmap(cdna_text->sparse_bitmap_builder);
  } else {
    sparse_bitmap_builder_add_bitmap(cdna_text->sparse_bitmap_builder,cdna_text->layer_2);
  }
  // Clear Bitmaps
  cdna_text->layer_0 = 0;
  cdna_text->layer_1 = 0;
  cdna_text->layer_2 = 0;
  cdna_text->position_mod64 = 0;
}
void cdna_bitwise_text_builder_add_char(
    cdna_bitwise_text_builder_t* const cdna_text,
    const uint8_t enc_char) {
  // Write new character
  switch (enc_char) {
  case ENC_DNA_CHAR_A: /* 000 */
    // cdna_text->layer_2 |= CDNA_BITWISE_WRITE_LAST_MASK;
    // cdna_text->layer_1 |= CDNA_BITWISE_WRITE_LAST_MASK;
    // cdna_text->layer_0 |= CDNA_BITWISE_WRITE_LAST_MASK;
    break;
  case ENC_DNA_CHAR_C: /* 001 */
    // cdna_text->layer_2 |= CDNA_BITWISE_WRITE_LAST_MASK;
    // cdna_text->layer_1 |= CDNA_BITWISE_WRITE_LAST_MASK;
    cdna_text->layer_0 |= CDNA_BITWISE_WRITE_LAST_MASK;
    break;
  case ENC_DNA_CHAR_G: /* 010 */
    // cdna_text->layer_2 |= CDNA_BITWISE_WRITE_LAST_MASK;
    cdna_text->layer_1 |= CDNA_BITWISE_WRITE_LAST_MASK;
    // cdna_text->layer_0 |= CDNA_BITWISE_WRITE_LAST_MASK;
    break;
  case ENC_DNA_CHAR_T: /* 011 */
    // cdna_text->layer_2 |= CDNA_BITWISE_WRITE_LAST_MASK;
    cdna_text->layer_1 |= CDNA_BITWISE_WRITE_LAST_MASK;
    cdna_text->layer_0 |= CDNA_BITWISE_WRITE_LAST_MASK;
    break;
  case ENC_DNA_CHAR_N: /* 100 */
    cdna_text->layer_2 |= CDNA_BITWISE_WRITE_LAST_MASK;
    // cdna_text->layer_1 |= CDNA_BITWISE_WRITE_LAST_MASK;
    // cdna_text->layer_0 |= CDNA_BITWISE_WRITE_LAST_MASK;
    break;
  case ENC_DNA_CHAR_SEP: /* 101 */
    cdna_text->layer_2 |= CDNA_BITWISE_WRITE_LAST_MASK;
    // cdna_text->layer_1 |= CDNA_BITWISE_WRITE_LAST_MASK;
    cdna_text->layer_0 |= CDNA_BITWISE_WRITE_LAST_MASK;
    break;
  case ENC_DNA_CHAR_JUMP: /* 110 */
    // cdna_text->layer_2 |= CDNA_BITWISE_WRITE_LAST_MASK;
    cdna_text->layer_1 |= CDNA_BITWISE_WRITE_LAST_MASK;
    cdna_text->layer_0 |= CDNA_BITWISE_WRITE_LAST_MASK;
    break;
  default:
    gem_fatal_error(CDNA_BITWISE_NOT_VALID_CHARACTER,enc_char);
    break;
  }
  // Inc position
  ++(cdna_text->position_mod64);
  // Dump (if needed)
  if (cdna_text->position_mod64==CDNA_BITWISE_BLOCK_CHARS) {
    cdna_bitwise_text_builder_flush(cdna_text);
  } else {
    // Shift to allocate next character
    cdna_text->layer_0 >>= 1;
    cdna_text->layer_1 >>= 1;
    cdna_text->layer_2 >>= 1;
  }
}
void cdna_bitwise_text_builder_close(cdna_bitwise_text_builder_t* const cdna_text) {
  if (cdna_text->position_mod64 > 0) {
    // Padding
    cdna_text->layer_2 |= CDNA_BITWISE_WRITE_LAST_MASK;
    cdna_text->layer_1 |= CDNA_BITWISE_WRITE_LAST_MASK;
    cdna_text->layer_0 |= CDNA_BITWISE_WRITE_LAST_MASK;
    ++(cdna_text->position_mod64);
    while (cdna_text->position_mod64!=CDNA_BITWISE_BLOCK_CHARS) {
      cdna_text->layer_0 >>= 1;
      cdna_text->layer_1 >>= 1;
      cdna_text->layer_2 >>= 1;
      cdna_text->layer_2 |= CDNA_BITWISE_WRITE_LAST_MASK;
      cdna_text->layer_1 |= CDNA_BITWISE_WRITE_LAST_MASK;
      cdna_text->layer_0 |= CDNA_BITWISE_WRITE_LAST_MASK;
      ++(cdna_text->position_mod64);
    }
    // Flush
    cdna_bitwise_text_builder_flush(cdna_text);
  }
  // Write Sparse 3rd-layer
  sparse_bitmap_builder_write(cdna_text->file_manager,cdna_text->sparse_bitmap_builder);
}
void cdna_bitwise_text_builder_delete(cdna_bitwise_text_builder_t* const cdna_text) {
  // Sparse-Bitmap
  sparse_bitmap_builder_delete(cdna_text->sparse_bitmap_builder);
  // Handler
  mm_free(cdna_text);
}
/*
 * Accessors
 */
uint64_t cdna_bitwise_text_get_size(cdna_bitwise_text_t* const cdna_text) {
  const uint64_t sparse_bitmap_size = sparse_bitmap_get_size(cdna_text->sparse_bitmap);
  return cdna_text->text_size + sparse_bitmap_size;
}
/*
 * Display
 */
void cdna_bitwise_text_print(FILE* const stream,cdna_bitwise_text_t* const cdna_text) {
  tab_fprintf(stream,"[GEM]>Compacted DNA-text\n");
  tab_fprintf(stream,"  => Architecture CDNA.3b.2bm64.xl\n");
  tab_fprintf(stream,"  => Text.length %"PRIu64"\n",cdna_text->text_length);
  const uint64_t sparse_bitmap_size = sparse_bitmap_get_size(cdna_text->sparse_bitmap);
  const uint64_t total_size = cdna_text->text_size + sparse_bitmap_size;
  tab_fprintf(stream,"  => Text.size %"PRIu64" MB (100 %)\n",total_size);
  tab_fprintf(stream,"    => Text.size %"PRIu64" MB (%2.3f%%)\n",cdna_text->text_size,PERCENTAGE(cdna_text->text_size,total_size));
  tab_fprintf(stream,"    => SpaceBitmap.size %"PRIu64" MB (%2.3f%%)\n",sparse_bitmap_size,PERCENTAGE(sparse_bitmap_size,total_size));
  sparse_bitmap_print(stream,cdna_text->sparse_bitmap,false);
  // Flush
  fflush(stream);
}


