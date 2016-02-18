/*
 * PROJECT: GEMMapper
 * FILE: output_buffer.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef OUTPUT_BUFFER_H_
#define OUTPUT_BUFFER_H_

#include "utils/essentials.h"

typedef enum { OUTPUT_BUFFER_FREE, OUTPUT_BUFFER_BUSY, OUTPUT_BUFFER_WRITE_PENDING } output_buffer_state_t;

typedef struct {
  /* Buffer state */
  uint32_t mayor_block_id; // Block mayor ID (for synchronization purposes)
  uint32_t minor_block_id; // Block minor ID (for segmented dumps)
  bool is_final_block;
  output_buffer_state_t buffer_state;
  /* Buffer */
  char* buffer_mem;        // Buffer Memory
  uint64_t buffer_used;    // Number of characters written
  char* buffer_cursor;     // Pointer to the current char
} output_buffer_t;

/*
 * Setup
 */
output_buffer_t* output_buffer_new(const uint64_t output_buffer_size);
void output_buffer_clear(output_buffer_t* const out_buffer);
void output_buffer_delete(output_buffer_t* const out_buffer);

/*
 * Accessors
 */
void output_buffer_set_state(output_buffer_t* const output_buffer,const output_buffer_state_t buffer_state);
output_buffer_state_t output_buffer_get_state(output_buffer_t* const output_buffer);
void output_buffer_set_incomplete(output_buffer_t* const output_buffer);
uint64_t output_buffer_get_used(output_buffer_t* const output_buffer);
char* output_buffer_get_buffer(output_buffer_t* const out_buffer);

/*
 * Fast-printer functions
 */
void bprintf_uint64(output_buffer_t* const out_buffer,const uint64_t number);
void bprintf_int64(output_buffer_t* const out_buffer,const int64_t number);
void bprintf_char(output_buffer_t* const out_buffer,const char character);
void bprintf_buffer(output_buffer_t* const out_buffer,const int string_length,const char* const string);

#endif /* OUTPUT_BUFFER_H_ */
