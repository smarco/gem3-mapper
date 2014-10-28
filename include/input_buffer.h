/*
 * PROJECT: GEMMapper
 * FILE: input_buffer.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef INPUT_BUFFER_H_
#define INPUT_BUFFER_H_

#include "essentials.h"

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
GEM_INLINE input_buffer_t* input_buffer_new();
GEM_INLINE void input_buffer_clear(input_buffer_t* const input_buffer);
GEM_INLINE void input_buffer_delete(input_buffer_t* const input_buffer);

/*
 * Accessors
 */
GEM_INLINE char** input_buffer_get_cursor(input_buffer_t* const input_buffer);
GEM_INLINE uint64_t input_buffer_get_cursor_pos(input_buffer_t* const input_buffer);
GEM_INLINE bool input_buffer_eob(input_buffer_t* const input_buffer);

/*
 * Utils
 */
GEM_INLINE void input_buffer_skip_line(input_buffer_t* const input_buffer);

#endif /* INPUT_BUFFER_H_ */
