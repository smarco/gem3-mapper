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
 *   Input Buffer data structure provides a simple input buffer to store
 *   a single input-block of the input file and enable multiple readers
 *   to read its content (keeping track of the readers). Keeps track of
 *   line offsets at to quickly compose/fill buffered-input-files
 */

#include "io/input_buffer.h"

/*
 * Constants
 */
#define INPUT_BUFFER_LINE_OFFSETS_INIT 10000

/*
 * Setup
 */
input_buffer_t* input_buffer_new(const uint64_t buffer_size) {
  // Alloc
  input_buffer_t* const input_buffer = mm_alloc(input_buffer_t);
  // Buffer info
  input_buffer->buffer_id = 0;
  input_buffer->buffer_state = input_buffer_empty;
  input_buffer->num_readers = 0;
  // Buffer
  input_buffer->buffer = mm_malloc(buffer_size);
  input_buffer->buffer_size = 0;
  input_buffer->buffer_allocated = buffer_size;
  // Line Index
  input_buffer->line_lengths = vector_new(INPUT_BUFFER_LINE_OFFSETS_INIT,uint64_t);
  // Return
  return input_buffer;
}
void input_buffer_delete(input_buffer_t* const input_buffer) {
  mm_free(input_buffer->buffer);
  vector_delete(input_buffer->line_lengths);
  mm_free(input_buffer);
}
/*
 * Annotate lines
 */
void input_buffer_annotate_lines(input_buffer_t* const input_buffer) {
  // Clear index
  vector_t* const line_lengths = input_buffer->line_lengths;
  vector_clear(line_lengths);
  // Traverse buffer & annotate line offsets
  const uint64_t buffer_size = input_buffer->buffer_size;
  const char* const buffer = input_buffer->buffer;
  // Implementation based on memchr (better for large lines)
  const char* sentinel = buffer;
  const char* last_sentinel = buffer;
  const char* const end = sentinel + buffer_size;
  while ((sentinel = memchr(sentinel,EOL,end-sentinel))) {
    ++sentinel;
    vector_insert(line_lengths,sentinel-last_sentinel,uint64_t);
    last_sentinel = sentinel;
  }
  // Insert the length of the remaining chars (Possible zero length line)
  vector_insert(line_lengths,end-last_sentinel,uint64_t); // But not a line (not '\n' found)
}
uint64_t input_buffer_get_num_lines(input_buffer_t* const input_buffer) {
  return vector_get_used(input_buffer->line_lengths)-1;
}
