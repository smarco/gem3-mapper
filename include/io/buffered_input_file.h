/*
 * PROJECT: GEMMapper
 * FILE: buffered_input_file.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef BUFFERED_INPUT_FILE_H_
#define BUFFERED_INPUT_FILE_H_

#include "system/commons.h"
#include "io/input_file_sliced.h"
#include "io/buffered_output_file.h"

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
