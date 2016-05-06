/*
 * PROJECT: GEMMapper
 * FILE: input_buffer.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
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
  input_buffer->line_lengths = vector_new(INPUT_BUFFER_LINE_OFFSETS_INIT,uint32_t);
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
  uint32_t i, offset = 0;
  for (i=0;i<buffer_size;++i) {
    const char current_char = buffer[i];
    if (gem_expect_false(current_char==EOL)) {
      vector_insert(line_lengths,offset+1,uint32_t);
      offset = 0;
    } else {
      ++offset;
    }
  }
  // Insert the length of the remaining chars
  vector_insert(line_lengths,offset,uint32_t);
}
uint64_t input_buffer_get_num_lines(input_buffer_t* const input_buffer) {
  return vector_get_used(input_buffer->line_lengths)-1;
}
