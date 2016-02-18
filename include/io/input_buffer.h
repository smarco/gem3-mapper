/*
 * PROJECT: GEMMapper
 * FILE: input_buffer.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef INPUT_BUFFER_H_
#define INPUT_BUFFER_H_

#include "utils/essentials.h"

typedef struct {
  uint32_t block_id;             // Block ID
  vector_t* block_buffer;        // Buffer
  char* cursor;                  // Pointer to the current char
  uint64_t lines_in_buffer;      // Actual number of lines in buffer
  uint64_t current_line_num;     // Current line no
} input_buffer_t;

/*
 * Setup
 */
input_buffer_t* input_buffer_new();
void input_buffer_clear(input_buffer_t* const input_buffer);
void input_buffer_delete(input_buffer_t* const input_buffer);

/*
 * Accessors
 */
char** input_buffer_get_cursor(input_buffer_t* const input_buffer);
uint64_t input_buffer_get_cursor_pos(input_buffer_t* const input_buffer);
bool input_buffer_eob(input_buffer_t* const input_buffer);

/*
 * Utils
 */
void input_buffer_skip_line(input_buffer_t* const input_buffer);

#endif /* INPUT_BUFFER_H_ */
