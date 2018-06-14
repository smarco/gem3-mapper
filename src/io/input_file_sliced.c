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
 *   Input module allows reading of input files block-wise (Prepares input-buffers)
 */

#include "io/input_file_sliced.h"
#include "io/input_buffer.h"

/*
 * Setup
 */
void input_file_sliced_initialize(
    input_file_sliced_t* const input_file_sliced,
    const uint64_t input_num_blocks,
    const uint64_t block_size) {
  // Input-buffers
  input_file_sliced->total_input_buffers = input_num_blocks;
  input_file_sliced->num_buffers_empty = input_num_blocks;
  input_file_sliced->input_buffers = mm_calloc(input_num_blocks,input_buffer_t*,false);
  uint64_t i;
  for (i=0;i<input_num_blocks;++i) {
    input_file_sliced->input_buffers[i] = input_buffer_new(block_size);
  }
  // Input-buffers order
  input_file_sliced->processed_buffer_id = 0;
  input_file_sliced->current_buffer_id = 0;
  input_file_sliced->current_buffer = NULL;
  // Input-buffers access
  MUTEX_INIT(input_file_sliced->input_buffers_mutex);
  CV_INIT(input_file_sliced->input_buffers_ready_cond);
  CV_INIT(input_file_sliced->input_buffers_empty_cond);
  // Input-buffer request queue
  input_file_sliced->input_buffer_id = 0;
  input_file_sliced->current_buffer_offset = 0;
  input_file_sliced->current_buffer_line = 0;
  MUTEX_INIT(input_file_sliced->input_read_lines_mutex);
}
input_file_sliced_t* input_stream_sliced_open(
    FILE* stream,
    const uint64_t input_num_blocks,
    const uint64_t block_size) {
  // Allocate handler
  input_file_sliced_t* const input_file_sliced = mm_alloc(input_file_sliced_t);
  // Initialize
  input_file_sliced->file_manager = fm_open_FILE(stream,FM_READ);
  input_file_sliced_initialize(input_file_sliced,input_num_blocks,block_size);
  // Return
  return input_file_sliced;
}
input_file_sliced_t* input_gzip_stream_sliced_open(
    FILE* stream,
    const uint64_t input_num_blocks,
    const uint64_t block_size){
  // Allocate handler
  input_file_sliced_t* const input_file_sliced = mm_alloc(input_file_sliced_t);
  // Initialize
  input_file_sliced->file_manager = fm_open_gzFILE(stream,FM_READ);
  input_file_sliced_initialize(input_file_sliced,input_num_blocks,block_size);
  // Return
  return input_file_sliced;
}
input_file_sliced_t* input_bzip_stream_sliced_open(
    FILE* stream,
    const uint64_t input_num_blocks,
    const uint64_t block_size){
  // Allocate handler
  input_file_sliced_t* const input_file_sliced = mm_alloc(input_file_sliced_t);
  // Initialize
  input_file_sliced->file_manager = fm_open_bzFILE(stream,FM_READ);
  input_file_sliced_initialize(input_file_sliced,input_num_blocks,block_size);
  // Return
  return input_file_sliced;
}
input_file_sliced_t* input_file_sliced_open(
    char* const file_name,
    const uint64_t input_num_blocks,
    const uint64_t block_size){
  // Allocate handler
  input_file_sliced_t* const input_file_sliced = mm_alloc(input_file_sliced_t);
  // Initialize
  input_file_sliced->file_manager = fm_open_file(file_name,FM_READ);
  input_file_sliced_initialize(input_file_sliced,input_num_blocks,block_size);
  // Return
  return input_file_sliced;
}
input_file_sliced_t* input_file_sliced_popen(
    char* const file_name,
    const uint64_t input_num_blocks,
    const uint64_t block_size){
  // Allocate handler
  input_file_sliced_t* const input_file_sliced = mm_alloc(input_file_sliced_t);
  // Initialize
  input_file_sliced->file_manager = fm_open_popen(file_name,FM_READ);
  input_file_sliced_initialize(input_file_sliced,input_num_blocks,block_size);
  // Return
  return input_file_sliced;
}
void input_file_sliced_close(input_file_sliced_t* const input_file_sliced) {
  fm_close(input_file_sliced->file_manager);
  MUTEX_DESTROY(input_file_sliced->input_buffers_mutex);
  MUTEX_DESTROY(input_file_sliced->input_read_lines_mutex);
  // Buffers
  const uint64_t total_input_buffers = input_file_sliced->total_input_buffers;
  uint64_t i;
  for (i=0;i<total_input_buffers;++i) {
    input_buffer_delete(input_file_sliced->input_buffers[i]);
  }
  mm_free(input_file_sliced->input_buffers);
  // Handler
  mm_free(input_file_sliced);
}
/*
 * Accessors
 */
char* input_file_sliced_get_file_name(
    input_file_sliced_t* const input_file_sliced) {
  char* const file_name = fm_get_file_name(input_file_sliced->file_manager);
  return (file_name!=NULL) ? file_name : STREAM_FILE_NAME;
}
uint64_t input_file_sliced_get_next_id(
    input_file_sliced_t* const input_file_sliced) {
  const uint64_t assign_id = input_file_sliced->input_buffer_id;
  input_file_sliced->input_buffer_id = (input_file_sliced->input_buffer_id+1) % UINT32_MAX;
  return assign_id;
}
/*
 * Buffer Request
 */
input_buffer_t* input_file_sliced_request_buffer_empty(input_file_sliced_t* const input_file_sliced) {
  // Check total buffers empty
  while (input_file_sliced->num_buffers_empty == 0) {
    CV_WAIT(input_file_sliced->input_buffers_empty_cond,input_file_sliced->input_buffers_mutex);
  }
  // Get buffer empty
  const uint64_t total_input_buffers = input_file_sliced->total_input_buffers;
  input_buffer_t* input_buffer = NULL;
  uint64_t i;
  for (i=0;i<total_input_buffers;++i) {
    if (input_file_sliced->input_buffers[i]->buffer_state == input_buffer_empty) {
      input_buffer = input_file_sliced->input_buffers[i];
      break;
    }
  }
  --(input_file_sliced->num_buffers_empty);
  return input_buffer;
}
input_buffer_t* input_file_sliced_find_buffer(
    input_file_sliced_t* const input_file_sliced,
    const uint64_t buffer_id) {
  // Get next buffer
  const uint64_t total_input_buffers = input_file_sliced->total_input_buffers;
  input_buffer_t* input_buffer = NULL;
  uint64_t i;
  for (i=0;i<total_input_buffers;++i) {
    if (input_file_sliced->input_buffers[i]->buffer_id == buffer_id) {
      input_buffer = input_file_sliced->input_buffers[i];
      break;
    }
  }
  return input_buffer;
}
/*
 * Process Input-buffers
 */
void input_file_sliced_input_buffer_read(
    input_file_sliced_t* const input_file_sliced,
    input_buffer_t* const input_buffer) {
  input_buffer->buffer_size = fm_read_mem(input_file_sliced->file_manager,
      input_buffer->buffer,input_buffer->buffer_allocated);
  fm_prefetch_next(input_file_sliced->file_manager,input_buffer->buffer_allocated);
  PROF_ADD_COUNTER(GP_BUFFERED_INPUT_BUFFER_SIZE,input_buffer->buffer_size);
}
void input_file_sliced_input_buffer_process(
    input_file_sliced_t* const input_file_sliced,
    input_buffer_t* const input_buffer) {
  input_buffer->buffer_id = input_file_sliced->processed_buffer_id;
  input_buffer->buffer_state = input_buffer_processing;
  ++(input_file_sliced->processed_buffer_id); // Next Input-Block
  input_file_sliced_input_buffer_read(input_file_sliced,input_buffer);
  MUTEX_END_SECTION(input_file_sliced->input_buffers_mutex);
  // Annotate the lines in the buffer
  input_buffer_annotate_lines(input_buffer);
  // Store the annotated buffer
  MUTEX_BEGIN_SECTION(input_file_sliced->input_buffers_mutex);
  input_buffer->buffer_state = input_buffer_ready;
  if (input_buffer->buffer_id == input_file_sliced->current_buffer_id) {
    CV_BROADCAST(input_file_sliced->input_buffers_ready_cond);
  }
}
void input_file_sliced_process(
    input_file_sliced_t* const input_file_sliced) {
  // Check input buffers
  MUTEX_BEGIN_SECTION(input_file_sliced->input_buffers_mutex);
  while (true) {
    // Check the number of count-pending buffers
    if (input_file_sliced->num_buffers_empty == 0 || fm_eof(input_file_sliced->file_manager)) break;
    // Read into the next empty buffer
    input_buffer_t* const input_buffer = input_file_sliced_request_buffer_empty(input_file_sliced);
    input_file_sliced_input_buffer_process(input_file_sliced,input_buffer);
  }
  MUTEX_END_SECTION(input_file_sliced->input_buffers_mutex);
}
void input_file_sliced_discard_exhausted_buffers(
    input_file_sliced_t* const input_file_sliced,
    vector_t* const input_buffers) {
  const uint64_t num_buffers = vector_get_used(input_buffers);
  if (num_buffers==0) return; // Quick return
  // Set exhausted buffers as empty
  MUTEX_BEGIN_SECTION(input_file_sliced->input_buffers_mutex);
  const uint64_t current_buffer_id = input_file_sliced->current_buffer_id;
  uint64_t i;
  for (i=0;i<num_buffers;++i) {
    input_buffer_t* const input_buffer = *(vector_get_elm(input_buffers,i,input_buffer_t*));
    --(input_buffer->num_readers);
    if (input_buffer->num_readers==0 && input_buffer->buffer_id < current_buffer_id) {
      input_buffer->buffer_state = input_buffer_empty; // Set empty
      if (input_file_sliced->num_buffers_empty == 0) {
        // Send broadcast if any waiting for empty
        CV_BROADCAST(input_file_sliced->input_buffers_empty_cond);
      }
      ++(input_file_sliced->num_buffers_empty);
    }
  }
  MUTEX_END_SECTION(input_file_sliced->input_buffers_mutex);
}
/*
 * Read input-buffer lines
 */
bool input_file_sliced_read_lines(
    input_file_sliced_t* const input_file_sliced,
    input_buffer_t* const current_file_buffer,
    const uint64_t prefered_read_size,
    const uint64_t forced_read_lines,
    uint64_t* const total_read_lines,
    uint64_t* const total_read_size,
    uint64_t* const input_last_buffer_line_end) {
  // Parameters
  const uint64_t total_lines = input_buffer_get_num_lines(current_file_buffer);
  uint64_t* const line_lengths = vector_get_mem(current_file_buffer->line_lengths,uint64_t);
  // Get current buffer line
  uint64_t current_buffer_offset = input_file_sliced->current_buffer_offset;
  uint64_t current_buffer_line = input_file_sliced->current_buffer_line;
  uint64_t read_size = *total_read_size, read_lines = *total_read_lines;
  while (current_buffer_line < total_lines) {
    const uint64_t line_length = line_lengths[current_buffer_line];
    read_size += line_length;
    current_buffer_offset += line_length;
    ++read_lines;
    ++current_buffer_line;
    // Check amount of data read
    if ((forced_read_lines>0 && read_lines==forced_read_lines) ||
        (forced_read_lines==0 && read_lines%8==0 && read_size>=prefered_read_size)) {
      input_file_sliced->current_buffer_offset = current_buffer_offset;
      input_file_sliced->current_buffer_line = current_buffer_line;
      *input_last_buffer_line_end = input_file_sliced->current_buffer_line; // Set buffer limit
      *total_read_lines = read_lines;
      *total_read_size = read_size;
      return true;
    }
  }
  // Current line-buffer exhausted
  const uint64_t line_length = line_lengths[current_buffer_line]; // Last characters of the buffer
  read_size += line_length;
  current_buffer_offset += line_length;
  *input_last_buffer_line_end = current_buffer_line; // Set buffer limit
  *total_read_lines = read_lines;
  *total_read_size = read_size;
  // Next input-buffer
  input_file_sliced->current_buffer = NULL;
  ++(input_file_sliced->current_buffer_id);
  input_file_sliced->current_buffer_offset = 0;
  input_file_sliced->current_buffer_line = 0;
  // Check amount of data read
  return (forced_read_lines==0 && read_lines%8==0 && read_size>=prefered_read_size);
}
/*
 * Current Input-Buffer
 */
input_buffer_t* input_file_sliced_input_buffer_get_current(
    input_file_sliced_t* const input_file_sliced) {
  // Check current buffer
  input_buffer_t* input_buffer = NULL;
  if (input_file_sliced->current_buffer != NULL) {
    input_buffer = input_file_sliced->current_buffer;
    // Increase the number of readers
    ++(input_buffer->num_readers);
    // Return the buffer
    return input_buffer;
  }
  // Find next buffer ready
  MUTEX_BEGIN_SECTION(input_file_sliced->input_buffers_mutex);
  while (input_buffer == NULL) {
    // Find next buffer
    const uint64_t current_buffer_id = input_file_sliced->current_buffer_id;
    input_buffer = input_file_sliced_find_buffer(input_file_sliced,current_buffer_id);
    if (input_buffer != NULL) {
      // Buffer found
      switch (input_buffer->buffer_state) {
        case input_buffer_ready:
          input_file_sliced->current_buffer = input_buffer;
          break;
        case input_buffer_processing:
          // Another thread is annotating the buffer (just wait for it)
          CV_WAIT(input_file_sliced->input_buffers_ready_cond,input_file_sliced->input_buffers_mutex);
          input_file_sliced->current_buffer = input_buffer;
          break;
        default:
          GEM_INVALID_CASE();
          break;
      }
    } else {
      // EOF
      if (fm_eof(input_file_sliced->file_manager)) {
        MUTEX_END_SECTION(input_file_sliced->input_buffers_mutex);
        return NULL;
      }
      // Buffer not processed yet
      if (input_file_sliced->num_buffers_empty > 0) {
        // Read into the next empty buffer
        input_buffer = input_file_sliced_request_buffer_empty(input_file_sliced);
        input_file_sliced_input_buffer_process(input_file_sliced,input_buffer);
        input_file_sliced->current_buffer = input_buffer;
        break;
      } else {
        // No empty buffers (go sleep & wait)
        CV_WAIT(input_file_sliced->input_buffers_empty_cond,input_file_sliced->input_buffers_mutex);
      }
    }
  }
  MUTEX_END_SECTION(input_file_sliced->input_buffers_mutex);
  // Increase the number of readers
  ++(input_buffer->num_readers);
  // Return the buffer
  return input_buffer;
}

