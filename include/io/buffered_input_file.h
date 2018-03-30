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
 *   Buffered input file provides buffered-reading from an input file.
 *     Data reading is line-oriented
 *     Compliant with regular input_file
 *     Compliant with input_file_sliced (for fast block parallel-reading)
 *     Keeps track of input-block(s) ordering
 */

#ifndef BUFFERED_INPUT_FILE_H_
#define BUFFERED_INPUT_FILE_H_

#include "system/commons.h"
#include "io/input_file_sliced.h"
#include "io/buffered_output_file.h"

/*
 * Codes status
 */
#define INPUT_STATUS_OK 1
#define INPUT_STATUS_EOF 0
#define INPUT_STATUS_FAIL -1

/*
 * Buffered input file
 */
typedef struct {
  /* Input file */
  input_file_sliced_t* input_file_sliced;   // Source input file
  uint64_t prefered_read_size;              // Size to reload the buffer
  /* Block buffer info */
  uint32_t block_id;                        // Block ID
  uint64_t num_lines;                       // Total number of lines in buffer
  uint64_t current_line_no;                 // Current line no
  /* Input-file Buffers */
  vector_t* input_buffers;                  // Associated input-buffers (input_buffer_t*)
  uint64_t input_buffer_next;               // Next input-buffer offset
  uint64_t input_first_buffer_offset;       // First-Buffer beginning offset
  uint64_t input_first_buffer_line_begin;   // First-Buffer beginning line
  uint64_t input_last_buffer_line_end;      // Last-Buffer end line
  /* Current input-file buffer */
  input_buffer_t* current_buffer;           // Current buffer
  char* current_buffer_sentinel;            // Current buffer sentinel
  uint64_t current_buffer_line_no;          // Current buffer line number
  uint64_t current_buffer_line_max;         // Current buffer max. line number
  /* Attached output buffer */
  buffered_output_file_t* attached_buffered_output_file;  // Attach output buffer file
} buffered_input_file_t;

/*
 * Setup
 */
buffered_input_file_t* buffered_input_file_new(
    input_file_sliced_t* const input_file_sliced,
    const uint64_t prefered_read_size);
void buffered_input_file_close(buffered_input_file_t* const buffered_input);

/*
 * Accessors
 */
char* buffered_input_file_get_file_name(buffered_input_file_t* const buffered_input);
uint32_t buffered_input_file_get_block_id(buffered_input_file_t* const buffered_input);
uint64_t buffered_input_file_get_num_lines(buffered_input_file_t* const buffered_input);
uint64_t buffered_input_file_get_current_line_num(buffered_input_file_t* const buffered_input);
bool buffered_input_file_eob(buffered_input_file_t* const buffered_input);

/*
 * Line Reader
 */
int buffered_input_file_get_line(
    buffered_input_file_t* const buffered_input,
    string_t* const input_line);

/*
 * Buffer Reload
 */
uint64_t buffered_input_file_reload(
    buffered_input_file_t* const buffered_input,
    const uint64_t forced_read_lines);
void buffered_input_file_attach_buffered_output(
    buffered_input_file_t* const buffered_input_file,
    buffered_output_file_t* const buffered_output_file);

#endif /* BUFFERED_INPUT_FILE_H_ */
