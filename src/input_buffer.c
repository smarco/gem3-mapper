/*
 * PROJECT: GEMMapper
 * FILE: input_buffer.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "input_buffer.h"
#include "vector.h"

/*
 * Setup
 */
GEM_INLINE input_buffer_t* input_buffer_new() {
  // Allocate
  input_buffer_t* const input_buffer = mm_alloc(input_buffer_t);
  // Init
  input_buffer->block_buffer = vector_new(1,sizeof(uint8_t));
  input_buffer_clear(input_buffer);
  // Return
  return input_buffer;
}
GEM_INLINE void input_buffer_clear(input_buffer_t* const input_buffer) {
  input_buffer->block_id = UINT32_MAX;
  input_buffer->cursor = vector_get_mem(input_buffer->block_buffer,char);
  input_buffer->current_line_num = UINT64_MAX;
}
GEM_INLINE void input_buffer_delete(input_buffer_t* const input_buffer) {
  vector_delete(input_buffer->block_buffer);
  mm_free(input_buffer);
}
/*
 * Accessors
 */
GEM_INLINE char** input_buffer_get_cursor(input_buffer_t* const input_buffer) {
  return &(input_buffer->cursor);
}
GEM_INLINE uint64_t input_buffer_get_cursor_pos(input_buffer_t* const input_buffer) {
  GEM_CHECK_NULL(input_buffer->cursor);
  return input_buffer->cursor-vector_get_mem(input_buffer->block_buffer,char);
}
GEM_INLINE bool input_buffer_eob(input_buffer_t* const input_buffer) {
  return input_buffer_get_cursor_pos(input_buffer) >= vector_get_used(input_buffer->block_buffer);
}
/*
 * Utils
 */
GEM_INLINE void input_buffer_skip_line(input_buffer_t* const input_buffer) {
  if (!input_buffer_eob(input_buffer)) {
    while (input_buffer->cursor[0]!=EOL) {
      ++input_buffer->cursor;
    }
    input_buffer->cursor[0]=EOS;
    ++input_buffer->cursor;
    ++input_buffer->current_line_num;
  }
}

