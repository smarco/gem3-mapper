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
 * DESCRIPTION:
 *   Output buffer provides storage before final dump to output file
 */

#include "io/output_buffer.h"

/*
 * Setup
 */
output_buffer_t* output_buffer_new(const uint64_t output_buffer_size) {
  output_buffer_t* out_buffer = mm_alloc(output_buffer_t);
  out_buffer->buffer_mem = mm_calloc(output_buffer_size,char,false);
  output_buffer_clear(out_buffer);
  output_buffer_set_state(out_buffer,OUTPUT_BUFFER_FREE);
  return out_buffer;
}
void output_buffer_clear(output_buffer_t* const out_buffer) {
  // Clear State
  out_buffer->mayor_block_id=0;
  out_buffer->minor_block_id=0;
  out_buffer->is_final_block=true;
  // Clear Buffer
  out_buffer->buffer_used = 0;
  out_buffer->buffer_cursor = out_buffer->buffer_mem;
}
void output_buffer_delete(output_buffer_t* const out_buffer) {
  mm_free(out_buffer->buffer_mem);
  mm_free(out_buffer);
}
/*
 * Accessors
 */
void output_buffer_set_state(output_buffer_t* const out_buffer,const output_buffer_state_t buffer_state) {
  out_buffer->buffer_state=buffer_state;
}
output_buffer_state_t output_buffer_get_state(output_buffer_t* const out_buffer) {
  return out_buffer->buffer_state;
}
void output_buffer_set_incomplete(output_buffer_t* const out_buffer) {
  out_buffer->is_final_block=false;
}
uint64_t output_buffer_get_used(output_buffer_t* const out_buffer) {
  return out_buffer->buffer_used;
}
char* output_buffer_get_buffer(output_buffer_t* const out_buffer) {
  return out_buffer->buffer_mem;
}
/*
 * Fast-printer functions
 */
void bprintf_uint64(output_buffer_t* const out_buffer,const uint64_t number) {
  const int chars_printed = integer_to_ascii(out_buffer->buffer_cursor,number);
  out_buffer->buffer_cursor += chars_printed;
  out_buffer->buffer_used += chars_printed;
}
void bprintf_int64(output_buffer_t* const out_buffer,const int64_t number) {
  if (number >= 0) {
    const int chars_printed = integer_to_ascii(out_buffer->buffer_cursor,number);
    out_buffer->buffer_cursor += chars_printed;
    out_buffer->buffer_used += chars_printed;
  } else {
    bprintf_char(out_buffer,'-');
    const int chars_printed = integer_to_ascii(out_buffer->buffer_cursor,(-number));
    out_buffer->buffer_cursor += chars_printed;
    out_buffer->buffer_used += chars_printed;
  }
}
void bprintf_char(output_buffer_t* const out_buffer,const char character) {
  *(out_buffer->buffer_cursor) = character;
  ++(out_buffer->buffer_cursor);
  ++(out_buffer->buffer_used);
}
void bprintf_buffer(
    output_buffer_t* const out_buffer,const int string_length,const char* const string) {
  memcpy(out_buffer->buffer_cursor,string,string_length);
  out_buffer->buffer_cursor += string_length;
  out_buffer->buffer_used += string_length;
}

