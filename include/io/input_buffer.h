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

typedef enum {
  input_buffer_ready,
  input_buffer_empty,
  input_buffer_processing,
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
  vector_t* line_lengths;            // Length of every line in the buffer (uint32_t)
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
