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
#include "text/dna_text.h"

/*
 * Tables/Conversions Implementation
 */
const bool dna_canonical_table[256] =
{
    [0 ... 255] = false,
    ['A'] = true, ['C'] = true, ['G'] = true, ['T'] = true,
    ['a'] = true, ['c'] = true, ['g'] = true, ['t'] = true,
};
const bool dna_canonical_encoded_table[DNA_EXT_RANGE] =
{
  [ENC_DNA_CHAR_A] = true,
  [ENC_DNA_CHAR_C] = true,
  [ENC_DNA_CHAR_G] = true,
  [ENC_DNA_CHAR_T] = true,
  [ENC_DNA_CHAR_N] = false,
  [ENC_DNA_CHAR_SEP] = false,
  [ENC_DNA_CHAR_JUMP] = false
};
const bool dna_table[256] =
{
    [0 ... 255] = false,
    ['A'] = true, ['C'] = true, ['G'] = true, ['T'] = true, ['N'] = true,
    ['a'] = true, ['c'] = true, ['g'] = true, ['t'] = true, ['n'] = true
};
const bool dna_encoded_table[DNA_EXT_RANGE] =
{
  [ENC_DNA_CHAR_A] = true,
  [ENC_DNA_CHAR_C] = true,
  [ENC_DNA_CHAR_G] = true,
  [ENC_DNA_CHAR_T] = true,
  [ENC_DNA_CHAR_N] = true,
  [ENC_DNA_CHAR_SEP] = false,
  [ENC_DNA_CHAR_JUMP] = false
};
const bool extended_dna_table[256] =
{
    [0 ... 255] = false,
    ['A'] = true, ['C'] = true, ['G'] = true, ['T'] = true, ['N'] = true,
    ['a'] = true, ['c'] = true, ['g'] = true, ['t'] = true, ['n'] = true,
    [DNA_CHAR_SEP] = true,
    [DNA_CHAR_JUMP] = true
};
const bool extended_dna_encoded_table[DNA_EXT_RANGE] =
{
  [ENC_DNA_CHAR_A] = true,
  [ENC_DNA_CHAR_C] = true,
  [ENC_DNA_CHAR_G] = true,
  [ENC_DNA_CHAR_T] = true,
  [ENC_DNA_CHAR_N] = true,
  [ENC_DNA_CHAR_SEP] = true,
  [ENC_DNA_CHAR_JUMP] = true
};
const bool unmasked_dna_table[256] =
{
    [0 ... 255] = false,
    ['A'] = true, ['C'] = true, ['G'] = true, ['T'] = true
};
const bool iupac_code_table[256] =
{
  [0 ... 255] = false,
  /* Upper case */
  ['A'] = true, ['C'] = true, ['G'] = true, ['T'] = true,
  ['N'] = true, ['R'] = true, ['N'] = true, ['D'] = true,
  ['E'] = true, ['H'] = true, ['I'] = true, ['L'] = true,
  ['K'] = true, ['M'] = true, ['F'] = true, ['P'] = true,
  ['S'] = true, ['W'] = true, ['Y'] = true, ['V'] = true, ['Q'] = true,
  ['B'] = true, ['Z'] = true, ['X'] = true, ['U'] = true, ['R'] = true,
  /* Lower case */
  ['a'] = true, ['c'] = true, ['g'] = true, ['t'] = true,
  ['n'] = true, ['r'] = true, ['n'] = true, ['d'] = true,
  ['e'] = true, ['h'] = true, ['i'] = true, ['l'] = true,
  ['k'] = true, ['m'] = true, ['f'] = true, ['p'] = true,
  ['s'] = true, ['w'] = true, ['y'] = true, ['v'] = true, ['q'] = true,
  ['b'] = true, ['z'] = true, ['x'] = true, ['u'] = true, ['r'] = true
};
const char dna_normalized_table[256] =
{
    [0 ... 255] = 'N',
    ['A'] = 'A', ['C'] = 'C', ['G'] = 'G', ['T'] = 'T',
    ['a'] = 'a', ['c'] = 'c', ['g'] = 'g', ['t'] = 't',
};
const char dna_strictly_normalized_table[256] =
{
    [0 ... 255] = 'N',
    ['A'] = 'A', ['C'] = 'C', ['G'] = 'G', ['T'] = 'T',
};
const char dna_complement_table[256] =
{
  [0 ... 255] = '~',
  ['A'] = 'T', ['C'] = 'G', ['G'] = 'C',  ['T'] = 'A', ['N'] = 'N',
  ['a'] = 'T', ['c'] = 'G', ['g'] = 'C',  ['t'] = 'A', ['n'] = 'N',
  [DNA_CHAR_SEP] = DNA_CHAR_SEP,
  [DNA_CHAR_JUMP] = DNA_CHAR_JUMP
};
const uint8_t dna_encoded_complement_table[DNA_EXT_RANGE] =
{
  [ENC_DNA_CHAR_A] = ENC_DNA_CHAR_T,
  [ENC_DNA_CHAR_C] = ENC_DNA_CHAR_G,
  [ENC_DNA_CHAR_G] = ENC_DNA_CHAR_C,
  [ENC_DNA_CHAR_T] = ENC_DNA_CHAR_A,
  [ENC_DNA_CHAR_N] = ENC_DNA_CHAR_N,
  [ENC_DNA_CHAR_SEP] = ENC_DNA_CHAR_SEP,
  [ENC_DNA_CHAR_JUMP] = ENC_DNA_CHAR_JUMP
};
const uint8_t dna_encode_table[256] =
{
  [0 ... 255] = 4,
  ['A'] = 0, ['C'] = 1, ['G'] = 2,  ['T'] = 3, ['N'] = 4,
  ['a'] = 0, ['c'] = 1, ['g'] = 2,  ['t'] = 3, ['n'] = 4,
  [DNA_CHAR_SEP] = ENC_DNA_CHAR_SEP,
  [DNA_CHAR_JUMP] = ENC_DNA_CHAR_JUMP
};
const char dna_decode_table[DNA_EXT_RANGE] =
{
  [ENC_DNA_CHAR_A] = DNA_CHAR_A,
  [ENC_DNA_CHAR_C] = DNA_CHAR_C,
  [ENC_DNA_CHAR_G] = DNA_CHAR_G,
  [ENC_DNA_CHAR_T] = DNA_CHAR_T,
  [ENC_DNA_CHAR_N] = DNA_CHAR_N,
  [ENC_DNA_CHAR_SEP] = DNA_CHAR_SEP,
  [ENC_DNA_CHAR_JUMP] = DNA_CHAR_JUMP
};
const char dna_bisulfite_C2T_table[256] =
{
  [0 ... 255] = '~',
  ['A'] = 'A', ['C'] = 'T', ['G'] = 'G',  ['T'] = 'T', ['N'] = 'N',
  ['a'] = 'A', ['c'] = 'T', ['g'] = 'G',  ['t'] = 'T', ['n'] = 'N',
  [DNA_CHAR_SEP] = DNA_CHAR_SEP,
  [DNA_CHAR_JUMP] = DNA_CHAR_JUMP
};
const char dna_bisulfite_G2A_table[256] =
{
  [0 ... 255] = '~',
  ['A'] = 'A', ['C'] = 'C', ['G'] = 'A',  ['T'] = 'T', ['N'] = 'N',
  ['a'] = 'A', ['c'] = 'C', ['g'] = 'A',  ['t'] = 'T', ['n'] = 'N',
  [DNA_CHAR_SEP] = DNA_CHAR_SEP,
  [DNA_CHAR_JUMP] = DNA_CHAR_JUMP
};
const uint8_t dna_encoded_bisulfite_C2T_table[DNA_EXT_RANGE] =
{
  [ENC_DNA_CHAR_A] = ENC_DNA_CHAR_A,
  [ENC_DNA_CHAR_C] = ENC_DNA_CHAR_T,
  [ENC_DNA_CHAR_G] = ENC_DNA_CHAR_G,
  [ENC_DNA_CHAR_T] = ENC_DNA_CHAR_T,
  [ENC_DNA_CHAR_N] = ENC_DNA_CHAR_N,
  [ENC_DNA_CHAR_SEP] = ENC_DNA_CHAR_SEP,
  [ENC_DNA_CHAR_JUMP] = ENC_DNA_CHAR_JUMP
};
const uint8_t dna_encoded_bisulfite_G2A_table[DNA_EXT_RANGE] =
{
  [ENC_DNA_CHAR_A] = ENC_DNA_CHAR_A,
  [ENC_DNA_CHAR_C] = ENC_DNA_CHAR_C,
  [ENC_DNA_CHAR_G] = ENC_DNA_CHAR_A,
  [ENC_DNA_CHAR_T] = ENC_DNA_CHAR_T,
  [ENC_DNA_CHAR_N] = ENC_DNA_CHAR_N,
  [ENC_DNA_CHAR_SEP] = ENC_DNA_CHAR_SEP,
  [ENC_DNA_CHAR_JUMP] = ENC_DNA_CHAR_JUMP
};
const uint8_t dna_encoded_colorspace_table[DNA_EXT_RANGE][DNA_EXT_RANGE] = {
   /* A */  {ENC_DNA_CHAR_A,   ENC_DNA_CHAR_C,   ENC_DNA_CHAR_G,   ENC_DNA_CHAR_T,   ENC_DNA_CHAR_SEP, ENC_DNA_CHAR_N},
   /* C */  {ENC_DNA_CHAR_C,   ENC_DNA_CHAR_A,   ENC_DNA_CHAR_T,   ENC_DNA_CHAR_G,   ENC_DNA_CHAR_SEP, ENC_DNA_CHAR_N},
   /* G */  {ENC_DNA_CHAR_G,   ENC_DNA_CHAR_T,   ENC_DNA_CHAR_A,   ENC_DNA_CHAR_C,   ENC_DNA_CHAR_SEP, ENC_DNA_CHAR_N},
   /* T */  {ENC_DNA_CHAR_T,   ENC_DNA_CHAR_G,   ENC_DNA_CHAR_C,   ENC_DNA_CHAR_A,   ENC_DNA_CHAR_SEP, ENC_DNA_CHAR_N},
  /* SEP */ {ENC_DNA_CHAR_SEP, ENC_DNA_CHAR_SEP, ENC_DNA_CHAR_SEP, ENC_DNA_CHAR_SEP, ENC_DNA_CHAR_SEP, ENC_DNA_CHAR_SEP},
   /* N */  {ENC_DNA_CHAR_N,   ENC_DNA_CHAR_N,   ENC_DNA_CHAR_N,   ENC_DNA_CHAR_N,   ENC_DNA_CHAR_SEP, ENC_DNA_CHAR_N}
};
/*
 * DNA-Text Model & Version
 */
#define DNA_TEXT_MODEL_NO  5002ull
/*
 * DNA-Text
 */
typedef enum { dna_text_raw = 0, dna_text_compact=UINT64_MAX } dna_text_type;
struct _dna_text_t {
  /* Metadata */
  dna_text_type type;
  /* Text */
  uint8_t* buffer;
  uint8_t* text;
  uint64_t length;
  uint64_t allocated;
  /* MM */
  mm_t* mm_text;
  bool mm_extern;
};
/*
 * Setup/Loader
 */
dna_text_t* dna_text_read_mem(mm_t* const memory_manager) {
  // Alloc
  dna_text_t* const dna_text = mm_alloc(dna_text_t);
  // Read header
  const uint64_t dna_text_model_no = mm_read_uint64(memory_manager);
  gem_cond_error(dna_text_model_no!=DNA_TEXT_MODEL_NO,
      DNA_TEXT_WRONG_MODEL_NO,dna_text_model_no,(uint64_t)DNA_TEXT_MODEL_NO);
  dna_text->type = mm_read_uint64(memory_manager);
  dna_text->length = mm_read_uint64(memory_manager);
  dna_text->allocated = dna_text->length;
  // Read Text
  dna_text->mm_extern = true;
  dna_text->mm_text = NULL;
  dna_text->buffer = mm_read_mem(memory_manager,dna_text->length*UINT8_SIZE);
  dna_text->text = dna_text->buffer;
  // Return
  return dna_text;
}
void dna_text_delete(dna_text_t* const dna_text) {
  if (dna_text->mm_text!=NULL) {
    mm_bulk_free(dna_text->mm_text);
  } else if (!dna_text->mm_extern) {
    mm_free(dna_text->buffer);
  }
  mm_free(dna_text);
}
/*
 * Builder
 */
dna_text_t* dna_text_new(const uint64_t dna_text_length) {
  dna_text_t* const dna_text = mm_alloc(dna_text_t);
  dna_text->type = dna_text_raw;
  dna_text->allocated = dna_text_length;
  dna_text->buffer = mm_calloc(dna_text->allocated,uint8_t,false);
  dna_text->text = dna_text->buffer;
  dna_text->length = 0;
  dna_text->mm_text = NULL;
  dna_text->mm_extern = false;
  return dna_text;
}
dna_text_t* dna_text_padded_new(
    const uint64_t dna_text_length,
    const uint64_t init_padding,
    const uint64_t end_padding) {
  dna_text_t* const dna_text = mm_alloc(dna_text_t);
  dna_text->type = dna_text_raw;
  dna_text->allocated = dna_text_length+init_padding+end_padding;
  dna_text->buffer = mm_calloc(dna_text->allocated,uint8_t,false);
  dna_text->text = dna_text->buffer + init_padding;
  dna_text->length = 0;
  dna_text->mm_text = NULL;
  dna_text->mm_extern = false;
  return dna_text;
}
// Builder writer
void dna_text_write_chunk(
    fm_t* const output_file_manager,
    dna_text_t* const dna_text,
    const uint64_t chunk_length) {
  fm_write_uint64(output_file_manager,DNA_TEXT_MODEL_NO);
  fm_write_uint64(output_file_manager,dna_text->type);
  fm_write_uint64(output_file_manager,chunk_length);
  fm_write_mem(output_file_manager,dna_text->text,chunk_length*UINT8_SIZE);
}
void dna_text_write(
    fm_t* const output_file_manager,
    dna_text_t* const dna_text) {
  dna_text_write_chunk(output_file_manager,dna_text,dna_text->length);
}
/*
 * Accessors
 */
uint64_t dna_text_get_length(const dna_text_t* const dna_text) {
  return dna_text->length;
}
void dna_text_set_length(
    dna_text_t* const dna_text,
    const uint64_t length) {
  dna_text->length = length;
}
uint64_t dna_text_get_size(const dna_text_t* const dna_text) {
  return dna_text->length;
}
uint8_t dna_text_get_char(
    const dna_text_t* const dna_text,
    const uint64_t position) {
  gem_fatal_check(position >= dna_text->allocated,DNA_TEXT_OOR,position,dna_text->allocated);
  return dna_text->text[position];
}
void dna_text_set_char(
    const dna_text_t* const dna_text,
    const uint64_t position,
    const uint8_t enc_char) {
  gem_fatal_check(position >= dna_text->allocated,DNA_TEXT_OOR,position,dna_text->allocated);
  dna_text->text[position] = enc_char;
}
uint8_t* dna_text_get_text(const dna_text_t* const dna_text) {
  return dna_text->text;
}
uint8_t* dna_text_retrieve_sequence(
    const dna_text_t* const dna_text,
    const uint64_t position) {
  uint8_t* const sequence = dna_text->text+position;
  PREFETCH(sequence); // Prefetch text // TODO Hint later on (LLC)
  return dna_text->text+position;
}
/*
 * Utils
 */
strand_t dna_text_strand_get_complement(const strand_t strand) {
  return (strand==Forward ? Reverse : Forward);
}
/*
 * Display
 */
void dna_text_print(
    FILE* const stream,
    dna_text_t* const dna_text) {
  fprintf(stream,"[GEM]>DNA-text\n");
  switch (dna_text->type) {
    case dna_text_raw:
      fprintf(stream,"  => Architecture Plain.encoded\n");
      break;
    case dna_text_compact:
      fprintf(stream,"  => Architecture Compact.encoded\n");
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  fprintf(stream,"  => Text.Length %"PRIu64"\n",dna_text->length);
  fprintf(stream,"  => Text.Size %"PRIu64" MB\n",CONVERT_B_TO_MB(dna_text->length*UINT8_SIZE));
  fflush(stream); // Flush
}
void dna_text_print_content(
    FILE* const stream,
    dna_text_t* const dna_text) {
  const uint8_t* const enc_text = dna_text->text;
  const uint64_t text_length = dna_text->length;
  fwrite(enc_text,1,text_length,stream);
}
void dna_text_pretty_print_content(
    FILE* const stream,
    dna_text_t* const dna_text,
    const uint64_t width) {
  // Iterate over all indexed text
  const uint8_t* const enc_text = dna_text->text;
  const uint64_t text_length = dna_text->length;
  uint64_t i, imod=0;
  for (i=0;i<text_length;++i) {
    // Print each character
    fprintf(stream,"%c",dna_decode(enc_text[i]));
    // Print column-wise
    if (++imod==width) {
      fprintf(stream,"\n"); imod = 0;
    }
  }
  if (imod!=width) fprintf(stream,"\n");
}
// Buffer Display
void dna_buffer_print(
    FILE* const stream,
    const uint8_t* const dna_buffer,
    const uint64_t dna_buffer_length,
    const bool print_reverse) {
  int64_t i;
  if (print_reverse) {
    for (i=dna_buffer_length-1;i>=0;--i) {
      fprintf(stream,"%c",dna_decode(dna_buffer[i]));
    }
  } else {
    for (i=0;i<dna_buffer_length;++i) {
      fprintf(stream,"%c",dna_decode(dna_buffer[i]));
    }
  }
}

