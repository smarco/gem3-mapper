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

#ifndef CDNA_BITWISE_TEXT_H_
#define CDNA_BITWISE_TEXT_H_

#include "utils/essentials.h"
#include "utils/sparse_bitmap.h"

/*
 * Compact DNA Bitwise Text
 */
typedef struct {
  /* Text */
  uint64_t* text;       // Interleaved less significant bit-layers
  uint64_t text_length; // Total text length (number of characters)
  uint64_t text_size;   // Text size (Bytes)
  /* Third bit-layer */
  sparse_bitmap_t* sparse_bitmap; // Sparse storage for the most significant bit-layer
  /* MM */
  mm_t* mm_text;
} cdna_bitwise_text_t;
// Compact DNA Bitwise Text Iterator
typedef struct {
  /* Compact DNA Text */
  cdna_bitwise_text_t* cdna_text;
  /* Current CDNA-Block */
  uint64_t layer_0;
  uint64_t layer_1;
  uint64_t layer_2;
  uint64_t position_mod64;
} cdna_bitwise_text_iterator_t;
// Compact DNA Bitwise Text Builder
typedef struct {
  /* File Manager */
  fm_t* file_manager;
  /* Third bit-layer */
  sparse_bitmap_builder_t* sparse_bitmap_builder;
  /* Current CDNA-Block */
  uint64_t layer_0;
  uint64_t layer_1;
  uint64_t layer_2;
  uint64_t position_mod64;
} cdna_bitwise_text_builder_t;

/*
 * CDNA Bitwise Text (Loader/Setup)
 */
cdna_bitwise_text_t* cdna_bitwise_text_read(fm_t* const file_manager);
cdna_bitwise_text_t* cdna_bitwise_text_read_mem(mm_t* const memory_manager);
void cdna_bitwise_text_delete(cdna_bitwise_text_t* const cdna_text);

// Iterator
// TODO void cdna_bitwise_text_iterator_new(
//    cdna_bitwise_text_iterator_t* const cdna_bitwise_text_iterator,
//    cdna_bitwise_text_t* const cdna_text,
//    const uint64_t position,
//    const traversal_direction_t text_traversal);
// TODO char cdna_bitwise_text_iterator_get_char(
//    cdna_bitwise_text_iterator_t* const cdna_bitwise_text_iterator);
// TODO uint8_t cdna_bitwise_text_iterator_get_enc(
//    cdna_bitwise_text_iterator_t* const cdna_bitwise_text_iterator);

/*
 * CDNA Bitwise Text (Builder)
 */
cdna_bitwise_text_builder_t* cdna_bitwise_text_builder_new(
    fm_t* const file_manager,
    const uint64_t text_length,
    mm_slab_t* const mm_slab);
void cdna_bitwise_text_builder_delete(
    cdna_bitwise_text_builder_t* const cdna_text);
void cdna_bitwise_text_builder_add_char(
    cdna_bitwise_text_builder_t* const cdna_text,
    const uint8_t enc_char);
void cdna_bitwise_text_builder_close(
    cdna_bitwise_text_builder_t* const cdna_text);

/*
 * Accessors
 */
uint64_t cdna_bitwise_text_get_size(cdna_bitwise_text_t* const cdna_text);

/*
 * Display
 */
void cdna_bitwise_text_print(FILE* const stream,cdna_bitwise_text_t* const cdna_text);

/*
 * Errors
 */
#define GEM_ERROR_CDNA_BITWISE_NOT_VALID_CHARACTER "Compact Bitwise DNA-Text Builder. Invalid character provided ASCII='%d'"

#endif /* CDNA_BITWISE_TEXT_H_ */
