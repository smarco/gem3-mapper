/*
 * PROJECT: GEMMapper
 * FILE: output_buffer.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "output_buffer.h"
#include "vector.h"
#define OUTPUT_BUFFER_INITIAL_SIZE BUFFER_SIZE_16M

/*
 * Setup
 */
GEM_INLINE output_buffer_t* output_buffer_new(void) {
  output_buffer_t* out_buffer = mm_alloc(output_buffer_t);
  out_buffer->buffer = vector_new(OUTPUT_BUFFER_INITIAL_SIZE,char);
  output_buffer_initiallize(out_buffer,OUTPUT_BUFFER_FREE);
  return out_buffer;
}
GEM_INLINE void output_buffer_clear(output_buffer_t* const out_buffer) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  out_buffer->mayor_block_id=UINT32_MAX;
  out_buffer->minor_block_id=0;
  out_buffer->is_final_block=true;
  vector_clear(out_buffer->buffer);
}
GEM_INLINE void output_buffer_initiallize(output_buffer_t* const out_buffer,const output_buffer_state_t buffer_state) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  output_buffer_clear(out_buffer);
  output_buffer_set_state(out_buffer,buffer_state);
}
GEM_INLINE void output_buffer_delete(output_buffer_t* const out_buffer) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  vector_delete(out_buffer->buffer);
  mm_free(out_buffer);
}
/*
 * Accessors
 */
GEM_INLINE void output_buffer_set_state(output_buffer_t* const out_buffer,const output_buffer_state_t buffer_state) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  out_buffer->buffer_state=buffer_state;
}
GEM_INLINE output_buffer_state_t output_buffer_get_state(output_buffer_t* const out_buffer) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  return out_buffer->buffer_state;
}
GEM_INLINE void output_buffer_set_partial_block(output_buffer_t* const out_buffer) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  out_buffer->is_final_block=false;
}
GEM_INLINE void output_buffer_set_mayor_block_id(output_buffer_t* const out_buffer,const uint32_t mayor_block_id) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  out_buffer->mayor_block_id=mayor_block_id;
}
GEM_INLINE uint32_t output_buffer_get_mayor_block_id(output_buffer_t* const out_buffer) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  return out_buffer->mayor_block_id;
}
GEM_INLINE void output_buffer_set_minor_block_id(output_buffer_t* const out_buffer,const uint32_t minor_block_id) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  out_buffer->minor_block_id=minor_block_id;
}
GEM_INLINE uint32_t output_buffer_get_minor_block_id(output_buffer_t* const out_buffer) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  return out_buffer->minor_block_id;
}
GEM_INLINE void output_buffer_inc_minor_block_id(output_buffer_t* const out_buffer) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  ++out_buffer->minor_block_id;
}
GEM_INLINE uint64_t output_buffer_get_used(output_buffer_t* const out_buffer) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  return vector_get_used(out_buffer->buffer);
}
/*
 * Adaptors
 */
GEM_INLINE char* output_buffer_to_char(output_buffer_t* const out_buffer) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  vector_insert(out_buffer->buffer,EOS,char);
  return vector_get_mem(out_buffer->buffer,char);
}
GEM_INLINE vector_t* output_buffer_to_vchar(output_buffer_t* const out_buffer) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  return out_buffer->buffer;
}
/*
 * Buffer printer
 */
GEM_INLINE int vbprintf(output_buffer_t* const out_buffer,const char *template,va_list v_args) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  GEM_CHECK_NULL(template);
  GEM_CHECK_NULL(v_args);
  // Calculate buffer size required to dump template{v_args}
  int64_t mem_required = calculate_memory_required_v(template,v_args);
  return vbprintf_fixed(out_buffer,mem_required,template,v_args);
}
GEM_INLINE int bprintf(output_buffer_t* const out_buffer,const char *template,...) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  GEM_CHECK_NULL(template);
  va_list v_args;
  va_start(v_args,template);
  const int chars_printed = vbprintf(out_buffer,template,v_args);
  va_end(v_args);
  return chars_printed;
}
GEM_INLINE int vbprintf_fixed(
    output_buffer_t* const out_buffer,const uint64_t expected_mem_usage,
    const char *template,va_list v_args) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  GEM_CHECK_NULL(template);
  GEM_CHECK_NULL(v_args);
  // Calculate buffer size required to dump template{v_args}
  vector_reserve_additional(out_buffer->buffer,expected_mem_usage);
  const int chars_printed = vsprintf(vector_get_free_elm(out_buffer->buffer,char),template,v_args);
  if (gem_expect_true(chars_printed>=0)) {
    vector_add_used(out_buffer->buffer,chars_printed);
  }
  return chars_printed;
}
GEM_INLINE int bprintf_fixed(
    output_buffer_t* const out_buffer,const uint64_t expected_mem_usage,const char *template,...) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  GEM_CHECK_NULL(template);
  va_list v_args;
  va_start(v_args,template);
  const int chars_printed = vbprintf_fixed(out_buffer,expected_mem_usage,template,v_args);
  va_end(v_args);
  return chars_printed;
}
