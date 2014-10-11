/*
 * PROJECT: GEMMapper
 * FILE: output_buffer.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

/*
 * TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
 *  Generalize minor/mayor to 3-levels
 *  At unsorted buffers extend to segment buffers
 */

#ifndef OUTPUT_BUFFER_H_
#define OUTPUT_BUFFER_H_

#include "essentials.h"

typedef enum { OUTPUT_BUFFER_FREE, OUTPUT_BUFFER_BUSY, OUTPUT_BUFFER_WRITE_PENDING } output_buffer_state_t;

typedef struct {
  /* Buffer state */
  uint32_t mayor_block_id; // Block mayor ID (for synchronization purposes)
  uint32_t minor_block_id; // Block minor ID (for segmented dumps)
  bool is_final_block;
  output_buffer_state_t buffer_state;
  /* Buffer */
  vector_t* buffer;
} output_buffer_t;

/*
 * Checkers
 */
#define OUTPUT_BUFFER_CHECK(output_buffer) GEM_CHECK_NULL(output_buffer)

/*
 * Setup
 */
GEM_INLINE output_buffer_t* output_buffer_new(const uint64_t output_buffer_size);
GEM_INLINE void output_buffer_clear(output_buffer_t* const out_buffer);
GEM_INLINE void output_buffer_delete(output_buffer_t* const out_buffer);

/*
 * Accessors
 */
GEM_INLINE void output_buffer_set_state(output_buffer_t* const output_buffer,const output_buffer_state_t buffer_state);
GEM_INLINE output_buffer_state_t output_buffer_get_state(output_buffer_t* const output_buffer);
GEM_INLINE void output_buffer_set_incomplete(output_buffer_t* const output_buffer);
GEM_INLINE uint64_t output_buffer_get_used(output_buffer_t* const output_buffer);
GEM_INLINE void output_buffer_reserve(output_buffer_t* const out_buffer,const uint64_t num_bytes);

/*
 * Adaptors
 */
GEM_INLINE char* output_buffer_to_char(output_buffer_t* const output_buffer);
GEM_INLINE vector_t* output_buffer_to_vchar(output_buffer_t* const output_buffer);

/*
 * Fast-printer functions
 */
GEM_INLINE void bprintf_uint64(output_buffer_t* const out_buffer,const uint64_t number);
GEM_INLINE void bprintf_int64(output_buffer_t* const out_buffer,const uint64_t number);
GEM_INLINE void bprintf_char(output_buffer_t* const out_buffer,const char character);
GEM_INLINE void bprintf_buffer(output_buffer_t* const out_buffer,const int string_length,const char* const string);

#endif /* OUTPUT_BUFFER_H_ */
