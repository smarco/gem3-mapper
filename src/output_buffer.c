/*
 * PROJECT: GEMMapper
 * FILE: output_buffer.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "output_buffer.h"
#include "vector.h"

/*
 * Setup
 */
GEM_INLINE output_buffer_t* output_buffer_new(const uint64_t output_buffer_size) {
  output_buffer_t* out_buffer = mm_alloc(output_buffer_t);
  out_buffer->buffer = vector_new(output_buffer_size,char);
  output_buffer_clear(out_buffer);
  output_buffer_set_state(out_buffer,OUTPUT_BUFFER_FREE);
  return out_buffer;
}
GEM_INLINE void output_buffer_clear(output_buffer_t* const out_buffer) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  out_buffer->mayor_block_id=0;
  out_buffer->minor_block_id=0;
  out_buffer->is_final_block=true;
  vector_clear(out_buffer->buffer);
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
GEM_INLINE void output_buffer_set_incomplete(output_buffer_t* const out_buffer) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  out_buffer->is_final_block=false;
}
GEM_INLINE uint64_t output_buffer_get_used(output_buffer_t* const out_buffer) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  return vector_get_used(out_buffer->buffer);
}
GEM_INLINE void output_buffer_reserve(output_buffer_t* const out_buffer,const uint64_t num_bytes) {
  vector_reserve_additional(out_buffer->buffer,num_bytes);
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
 * Buffer printer (Disabled for performance issues)
 */
//GEM_INLINE int vbprintf(output_buffer_t* const out_buffer,const char *template,va_list v_args) {
//  OUTPUT_BUFFER_CHECK(out_buffer);
//  GEM_CHECK_NULL(template);
//  GEM_CHECK_NULL(v_args);
//  // Calculate buffer size required to dump template{v_args}
//  int64_t mem_required = calculate_memory_required_v(template,v_args);
//  return vbprintf_fixed(out_buffer,mem_required,template,v_args);
//}
//GEM_INLINE int bprintf(output_buffer_t* const out_buffer,const char *template,...) {
//  OUTPUT_BUFFER_CHECK(out_buffer);
//  GEM_CHECK_NULL(template);
//  va_list v_args;
//  va_start(v_args,template);
//  const int chars_printed = vbprintf(out_buffer,template,v_args);
//  va_end(v_args);
//  return chars_printed;
//}
//GEM_INLINE int vbprintf_fixed(
//    output_buffer_t* const out_buffer,const uint64_t expected_mem_usage,
//    const char *template,va_list v_args) {
//  OUTPUT_BUFFER_CHECK(out_buffer);
//  GEM_CHECK_NULL(template);
//  GEM_CHECK_NULL(v_args);
//  // Calculate buffer size required to dump template{v_args}
//  vector_reserve_additional(out_buffer->buffer,expected_mem_usage);
//  const int chars_printed = vsprintf(vector_get_free_elm(out_buffer->buffer,char),template,v_args);
//  if (gem_expect_true(chars_printed>=0)) {
//    vector_add_used(out_buffer->buffer,chars_printed);
//  }
//  return chars_printed;
//}
//GEM_INLINE int bprintf_fixed(
//    output_buffer_t* const out_buffer,const uint64_t expected_mem_usage,const char *template,...) {
//  OUTPUT_BUFFER_CHECK(out_buffer);
//  GEM_CHECK_NULL(template);
//  va_list v_args;
//  va_start(v_args,template);
//  const int chars_printed = vbprintf_fixed(out_buffer,expected_mem_usage,template,v_args);
//  va_end(v_args);
//  return chars_printed;
//}
/*
 * Fast-printer functions
 */
GEM_INLINE void bprintf_uint64(output_buffer_t* const out_buffer,const uint64_t number) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  const int chars_printed = integer_to_ascii(vector_get_free_elm(out_buffer->buffer,char),number);
  vector_add_used(out_buffer->buffer,chars_printed);
}
GEM_INLINE void bprintf_int64(output_buffer_t* const out_buffer,const uint64_t number) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  const int chars_printed = integer_to_ascii(vector_get_free_elm(out_buffer->buffer,char),number);
  vector_add_used(out_buffer->buffer,chars_printed);
}
GEM_INLINE void bprintf_char(output_buffer_t* const out_buffer,const char character) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  *vector_get_free_elm(out_buffer->buffer,char) = character;
  vector_add_used(out_buffer->buffer,1);
}
GEM_INLINE void bprintf_buffer(
    output_buffer_t* const out_buffer,const int string_length,const char* const string) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  memcpy(vector_get_free_elm(out_buffer->buffer,char),string,string_length);
  vector_add_used(out_buffer->buffer,string_length);
}



