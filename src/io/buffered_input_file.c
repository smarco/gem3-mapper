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

#include "io/buffered_input_file.h"
#include "mapper/mapper_profile.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PLOW

/*
 * Constants
 */
#define BUFFERED_INPUT_FILE_NUM_BUFFERS_INIT 5

/*
 * Setup
 */
buffered_input_file_t* buffered_input_file_new(
    input_file_sliced_t* const input_file_sliced,
    const uint64_t prefered_read_size) {
  // Alloc
  buffered_input_file_t* const buffered_input = mm_alloc(buffered_input_file_t);
  // Input file
  buffered_input->input_file_sliced = input_file_sliced;
  buffered_input->prefered_read_size = prefered_read_size;
  // Block buffer info
  buffered_input->block_id = 0;
  buffered_input->num_lines = 0;
  buffered_input->current_line_no = 0;
  // Input buffers
  buffered_input->input_buffers = vector_new(BUFFERED_INPUT_FILE_NUM_BUFFERS_INIT,input_buffer_t*);
  buffered_input->input_buffer_next = 0;
  buffered_input->input_first_buffer_line_begin = 0;
  buffered_input->input_last_buffer_line_end = 0;
  // Current input-file buffer
  buffered_input->current_buffer = NULL;
  buffered_input->current_buffer_line_no = 0;
  buffered_input->current_buffer_line_max = 0;
  // Attached output buffer
  buffered_input->attached_buffered_output_file = NULL;
  // Return
  return buffered_input;
}
void buffered_input_file_close(buffered_input_file_t* const buffered_input) {
  vector_delete(buffered_input->input_buffers);
  mm_free(buffered_input);
}
/*
 * Accessors
 */
char* buffered_input_file_get_file_name(buffered_input_file_t* const buffered_input) {
  return input_file_sliced_get_file_name(buffered_input->input_file_sliced);
}
uint32_t buffered_input_file_get_block_id(buffered_input_file_t* const buffered_input) {
  return buffered_input->block_id;
}
uint64_t buffered_input_file_get_num_lines(buffered_input_file_t* const buffered_input) {
  return buffered_input->num_lines;
}
uint64_t buffered_input_file_get_current_line_num(buffered_input_file_t* const buffered_input) {
  return buffered_input->current_line_no;
}
bool buffered_input_file_eob(buffered_input_file_t* const buffered_input) {
  return buffered_input->current_buffer_line_no >= buffered_input->current_buffer_line_max &&
         buffered_input->input_buffer_next >= vector_get_used(buffered_input->input_buffers);
}
/*
 * Input-File Buffer Reader
 */
void buffered_input_file_read_buffer(
    input_file_sliced_t* const input_file_sliced,
    buffered_input_file_t* const buffered_input,
    const uint64_t forced_read_lines) {
  // Return exhausted input-buffers
  input_file_sliced_discard_exhausted_buffers(
      buffered_input->input_file_sliced,buffered_input->input_buffers);
  // Process input-buffers (read file & annotate line lengths)
  input_file_sliced_process(buffered_input->input_file_sliced);
  // Read a new buffer-block
  MUTEX_BEGIN_SECTION(input_file_sliced->input_read_lines_mutex);
  buffered_input->block_id = input_file_sliced_get_next_id(input_file_sliced);
  buffered_input->current_line_no = input_file_sliced->current_input_line;
  vector_clear(buffered_input->input_buffers);
  buffered_input->input_buffer_next = 0;
  buffered_input->input_first_buffer_offset = input_file_sliced->current_buffer_offset;
  buffered_input->input_first_buffer_line_begin = input_file_sliced->current_buffer_line;
  // Read Buffer-block
  uint64_t total_read_lines = 0, total_read_size = 0;
  bool buffer_filled = false;
  while (!buffer_filled) {
    // Get current input-buffer
    input_buffer_t* const input_buffer = input_file_sliced_input_buffer_get_current(input_file_sliced);
    if (input_buffer==NULL) break; // EOF
    vector_insert(buffered_input->input_buffers,input_buffer,input_buffer_t*);
    // Read lines
    buffer_filled = input_file_sliced_read_lines(
        input_file_sliced,input_buffer,buffered_input->prefered_read_size,
        forced_read_lines,&total_read_lines,&total_read_size,
        &buffered_input->input_last_buffer_line_end);
  }
  buffered_input->num_lines = total_read_lines; // Set buffer lines
  input_file_sliced->current_input_line += total_read_lines; // Update input line
  MUTEX_END_SECTION(input_file_sliced->input_read_lines_mutex);
}
/*
 * Line Reader
 */
int buffered_input_file_load_next_chunk(buffered_input_file_t* const buffered_input) {
  // Load next chunk
  vector_t* const input_buffers = buffered_input->input_buffers;
  const uint64_t num_input_buffers = vector_get_used(input_buffers);
  const uint64_t input_buffer_next = buffered_input->input_buffer_next;
  if (input_buffer_next < num_input_buffers) {
    // Buffer
    input_buffer_t* const input_buffer = *(vector_get_elm(input_buffers,input_buffer_next,input_buffer_t*));
    buffered_input->current_buffer = input_buffer;
    // Compute Offsets
    if (input_buffer_next == 0) {
      buffered_input->current_buffer_sentinel = input_buffer->buffer + buffered_input->input_first_buffer_offset;
      buffered_input->current_buffer_line_no = buffered_input->input_first_buffer_line_begin;
    } else {
      buffered_input->current_buffer_sentinel = input_buffer->buffer;
      buffered_input->current_buffer_line_no = 0;
    }
    if (input_buffer_next == num_input_buffers-1) {
      buffered_input->current_buffer_line_max = buffered_input->input_last_buffer_line_end;
    } else {
      buffered_input->current_buffer_line_max = input_buffer_get_num_lines(input_buffer);
    }
    // Next
    ++(buffered_input->input_buffer_next);
    return 1;
  } else {
    return 0; // EOB
  }
}
int buffered_input_file_get_line(
    buffered_input_file_t* const buffered_input,
    string_t* const input_line) {
  // Clear string
  string_clear(input_line);
  // Copy the remaining in the buffer (looking for the end-of-line)
  while (true) {
    // Delimit line
    char* const line = buffered_input->current_buffer_sentinel;
    const uint64_t* const lengths = vector_get_mem(buffered_input->current_buffer->line_lengths,uint64_t);
    const uint64_t line_length = lengths[buffered_input->current_buffer_line_no];
    // Append to line
    if (line_length > 0) string_right_append_buffer(input_line,line,line_length);
    // Update input buffer
    if (buffered_input->current_buffer_line_no >= buffered_input->current_buffer_line_max) {
      if (!buffered_input_file_load_next_chunk(buffered_input)) {
        string_append_char(input_line,EOL);
        string_append_eos(input_line);
        return string_get_length(input_line);
      }
    } else {
      buffered_input->current_buffer_sentinel += line_length;
      ++(buffered_input->current_buffer_line_no);
      ++(buffered_input->current_line_no);
      // Handle EOL
      const uint64_t input_line_length = string_get_length(input_line);
      if (input_line_length >= 2) {
        char* const input_line_buffer = string_get_buffer(input_line);
        if (input_line_buffer[input_line_length-2] == DOS_EOL) {
          input_line_buffer[input_line_length-2] = EOL;
          string_set_length(input_line,input_line_length-1);
        }
      }
      // Return
      return string_get_length(input_line);
    }
  }
}
/*
 * Utils
 */
uint64_t buffered_input_file_reload(
    buffered_input_file_t* const buffered_input,
    const uint64_t forced_read_lines) {
  PROFILE_START(GP_BUFFERED_INPUT_RELOAD__DUMP_ATTACHED,PROFILE_LEVEL);
  // Dump attached buffered-output
  if (buffered_input->attached_buffered_output_file != NULL) {
    buffered_output_file_dump_buffer(buffered_input->attached_buffered_output_file);
  }
  PROFILE_START(GP_BUFFERED_INPUT_RELOAD,PROFILE_LEVEL);
  // Read a new buffer-block
  buffered_input_file_read_buffer(buffered_input->input_file_sliced,buffered_input,forced_read_lines);
  if (buffered_input->num_lines==0) {
    PROFILE_STOP(GP_BUFFERED_INPUT_RELOAD,PROFILE_LEVEL);
    PROFILE_STOP(GP_BUFFERED_INPUT_RELOAD__DUMP_ATTACHED,PROFILE_LEVEL);
    return INPUT_STATUS_EOF; // EOF
  } else {
    buffered_input_file_load_next_chunk(buffered_input); // Load first chunk
  }
  PROFILE_STOP(GP_BUFFERED_INPUT_RELOAD,PROFILE_LEVEL);
  // Get output buffer (block ID)
  if (buffered_input->attached_buffered_output_file != NULL) {
    buffered_output_file_request_buffer(buffered_input->attached_buffered_output_file,buffered_input->block_id);
  }
  PROFILE_STOP(GP_BUFFERED_INPUT_RELOAD__DUMP_ATTACHED,PROFILE_LEVEL);
  return INPUT_STATUS_OK; // OK
}
void buffered_input_file_attach_buffered_output(
    buffered_input_file_t* const buffered_input_file,
    buffered_output_file_t* const buffered_output_file) {
  buffered_input_file->attached_buffered_output_file = buffered_output_file;
}
