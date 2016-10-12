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

#ifndef CDNA_TEXT_H_
#define CDNA_TEXT_H_

#include "utils/essentials.h"
#include "utils/segmented_vector.h"
#include "text/cdna_bitwise_text.h"

/*
 * Compact DNA Text
 */
typedef struct {
  /* Text */
  svector_t* text;         // Compacted 3bits-character text (21 chars in each 8Bytes)
  uint64_t text_length;    // Total text length (number of characters)
  /* Iterator */
  svector_iterator_t text_iterator;
  /* Internals (Builder) */
  uint64_t position_mod21;
  uint64_t current_word;
  /* MM */
  mm_slab_t* mm_slab;
} cdna_text_t;
typedef struct {
  /* Text */
  cdna_text_t* cdna_text;
  svector_iterator_t text_iterator;
  bool eoi;
  /* Internals (Iterator) */
  uint64_t position;
  uint64_t position_mod21;
  uint64_t current_word;
} cdna_text_iterator_t;
typedef struct {
  /* Text(s) */
  uint64_t num_texts;
  /* Iterators */
  cdna_text_iterator_t* iterators;
  uint64_t* limits;
  /* Current State */
  uint64_t text_position;
  uint64_t iterator_num;
  bool eoi;
} cdna_text_iterator_hub_t;

extern const uint64_t cdna_text_block_mask_left[256];

/*
 * CDNA Text (Builder incorporated)
 */
cdna_text_t* cdna_text_new(mm_slab_t* const mm_slab);
void cdna_text_delete(cdna_text_t* const cdna_text);
void cdna_text_add_char(cdna_text_t* const cdna_text,const uint8_t enc_char);
void cdna_text_close(cdna_text_t* const cdna_text);

void cdna_text_write_as_bitwise_text(fm_t* const file_manager,cdna_text_t* const cdna_text);

/*
 * CDNA Accessors
 */
uint64_t cdna_text_get_length(cdna_text_t* const cdna_text);

/*
 * CDNA Text Iterator
 */
void cdna_text_iterator_init(
    cdna_text_iterator_t* const iterator,
    cdna_text_t* const cdna_text,
    const uint64_t position);
bool cdna_text_iterator_eoi(cdna_text_iterator_t* const iterator);
uint8_t cdna_text_iterator_get_char_encoded(cdna_text_iterator_t* const iterator);
void cdna_text_iterator_next_char(cdna_text_iterator_t* const iterator);

/*
 * CDNA Text Reverse Iterator
 *   Note that:
 *    - Only works with an open cdna_text (cdna_text_close(.) not called)
 *    - Traverses the text from last added character until the first
 */
void cdna_text_reverse_iterator_init(
    cdna_text_iterator_t* const iterator,
    cdna_text_t* const cdna_text,
    const uint64_t position);
bool cdna_text_reverse_iterator_eoi(cdna_text_iterator_t* const iterator);
uint8_t cdna_text_reverse_iterator_get_char_encoded(cdna_text_iterator_t* const iterator);
void cdna_text_reverse_iterator_next_char(cdna_text_iterator_t* const iterator);

/*
 * CDNA Text Compare
 */
void cdna_text_block_iterator_init(
    cdna_text_iterator_t* const iterator,
    cdna_text_t* const cdna_text,
    const uint64_t position);
void cdna_text_block_iterator_next_block(cdna_text_iterator_t* const iterator);
uint64_t cdna_text_block_iterator_get_block(cdna_text_iterator_t* const iterator);
int cdna_text_block_cmp(cdna_text_iterator_t* const iterator_a,cdna_text_iterator_t* const iterator_b);

/*
 * CDNA Text (Iterator-HUB)
 */
cdna_text_iterator_hub_t* cdna_text_iterator_hub_new(
    cdna_text_t** cdna_texts,
    const uint64_t num_texts,
    const uint64_t starting_position);
void cdna_text_iterator_hub_delete(cdna_text_iterator_hub_t* const iterator_hub);
bool cdna_text_iterator_hub_eoi(cdna_text_iterator_hub_t* const iterator_hub);
uint8_t cdna_text_iterator_hub_get_char_encoded(cdna_text_iterator_hub_t* const iterator_hub);
void cdna_text_iterator_hub_next_char(cdna_text_iterator_hub_t* const iterator_hub);

/*
 * Display
 */
void cdna_text_print(FILE* const stream,cdna_text_t* const cdna_text);

#endif /* CDNA_TEXT_H_ */
