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
GEM_INLINE output_buffer_t* output_buffer_new(void);
GEM_INLINE void output_buffer_clear(output_buffer_t* const out_buffer);
GEM_INLINE void output_buffer_delete(output_buffer_t* const out_buffer);

/*
 * Accessors
 */
GEM_INLINE void output_buffer_set_state(output_buffer_t* const output_buffer,const output_buffer_state_t buffer_state);
GEM_INLINE output_buffer_state_t output_buffer_get_state(output_buffer_t* const output_buffer);
GEM_INLINE void output_buffer_set_incomplete(output_buffer_t* const output_buffer);
GEM_INLINE uint64_t output_buffer_get_used(output_buffer_t* const output_buffer);

/*
 * Adaptors
 */
GEM_INLINE char* output_buffer_to_char(output_buffer_t* const output_buffer);
GEM_INLINE vector_t* output_buffer_to_vchar(output_buffer_t* const output_buffer);

/*
 * Buffer printer
 */
GEM_INLINE int vbprintf(output_buffer_t* const output_buffer,const char *template,va_list v_args);
GEM_INLINE int bprintf(output_buffer_t* const output_buffer,const char *template,...);
// If you know how much memory is going to be used
GEM_INLINE int vbprintf_fixed(
    output_buffer_t* const output_buffer,const uint64_t expected_mem_usage,
    const char *template,va_list v_args);
GEM_INLINE int bprintf_fixed(
    output_buffer_t* const output_buffer,const uint64_t expected_mem_usage,const char *template,...);

#endif /* OUTPUT_BUFFER_H_ */
