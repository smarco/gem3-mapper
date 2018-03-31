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

#ifndef INPUT_BUFFER_H_
#define INPUT_BUFFER_H_

#include "utils/essentials.h"

typedef enum {
  input_buffer_empty,
  input_buffer_processing,
  input_buffer_ready,
} input_buffer_state;
typedef struct {
  /* Buffer info */
  uint64_t buffer_id;                // Buffer ID
  input_buffer_state buffer_state;   // Buffer State
  uint64_t num_readers;              // Current number of readers
  /* Buffer */
  char* buffer;                      // Pointer to the buffer
  uint64_t buffer_size;              // Total bytes in buffer
  uint64_t buffer_allocated;         // Total bytes allocated
  /* Line Index */
  vector_t* line_lengths;            // Length of every line in the buffer (uint64_t)
} input_buffer_t;

/*
 * Setup
 */
input_buffer_t* input_buffer_new(const uint64_t buffer_size);
void input_buffer_delete(input_buffer_t* const input_buffer);

/*
 * Annotate lines
 */
void input_buffer_annotate_lines(input_buffer_t* const input_buffer);
uint64_t input_buffer_get_num_lines(input_buffer_t* const input_buffer);

#endif /* INPUT_BUFFER_H_ */
