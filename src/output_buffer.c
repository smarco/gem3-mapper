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
output_buffer_t* output_buffer_new(const uint64_t output_buffer_size) {
  output_buffer_t* out_buffer = mm_alloc(output_buffer_t);
  out_buffer->buffer_mem = mm_calloc(output_buffer_size,char,false);
  output_buffer_clear(out_buffer);
  output_buffer_set_state(out_buffer,OUTPUT_BUFFER_FREE);
  return out_buffer;
}
void output_buffer_clear(output_buffer_t* const out_buffer) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  // Clear State
  out_buffer->mayor_block_id=0;
  out_buffer->minor_block_id=0;
  out_buffer->is_final_block=true;
  // Clear Buffer
  out_buffer->buffer_used = 0;
  out_buffer->buffer_cursor = out_buffer->buffer_mem;
}
void output_buffer_delete(output_buffer_t* const out_buffer) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  mm_free(out_buffer->buffer_mem);
  mm_free(out_buffer);
}
/*
 * Accessors
 */
void output_buffer_set_state(output_buffer_t* const out_buffer,const output_buffer_state_t buffer_state) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  out_buffer->buffer_state=buffer_state;
}
output_buffer_state_t output_buffer_get_state(output_buffer_t* const out_buffer) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  return out_buffer->buffer_state;
}
void output_buffer_set_incomplete(output_buffer_t* const out_buffer) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  out_buffer->is_final_block=false;
}
uint64_t output_buffer_get_used(output_buffer_t* const out_buffer) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  return out_buffer->buffer_used;
}
char* output_buffer_get_buffer(output_buffer_t* const out_buffer) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  return out_buffer->buffer_mem;
}
/*
 * Buffer printer (Disabled for performance issues)
 */
//int vbprintf(output_buffer_t* const out_buffer,const char *template,va_list v_args) {
//  OUTPUT_BUFFER_CHECK(out_buffer);
//  GEM_CHECK_NULL(template);
//  GEM_CHECK_NULL(v_args);
//  // Calculate buffer size required to dump template{v_args}
//  int64_t mem_required = calculate_memory_required_v(template,v_args);
//  return vbprintf_fixed(out_buffer,mem_required,template,v_args);
//}
//int bprintf(output_buffer_t* const out_buffer,const char *template,...) {
//  OUTPUT_BUFFER_CHECK(out_buffer);
//  GEM_CHECK_NULL(template);
//  va_list v_args;
//  va_start(v_args,template);
//  const int chars_printed = vbprintf(out_buffer,template,v_args);
//  va_end(v_args);
//  return chars_printed;
//}
//int vbprintf_fixed(
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
//int bprintf_fixed(
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
void bprintf_uint64(output_buffer_t* const out_buffer,const uint64_t number) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  const int chars_printed = integer_to_ascii(out_buffer->buffer_cursor,number);
  out_buffer->buffer_cursor += chars_printed;
  out_buffer->buffer_used += chars_printed;
}
void bprintf_int64(output_buffer_t* const out_buffer,const int64_t number) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  if (number >= 0) {
    const int chars_printed = integer_to_ascii(out_buffer->buffer_cursor,number);
    out_buffer->buffer_cursor += chars_printed;
    out_buffer->buffer_used += chars_printed;
  } else {
    bprintf_char(out_buffer,'-');
    const int chars_printed = integer_to_ascii(out_buffer->buffer_cursor,(-number));
    out_buffer->buffer_cursor += chars_printed;
    out_buffer->buffer_used += chars_printed;
  }
}
void bprintf_char(output_buffer_t* const out_buffer,const char character) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  *(out_buffer->buffer_cursor) = character;
  ++(out_buffer->buffer_cursor);
  ++(out_buffer->buffer_used);
}
void bprintf_buffer(
    output_buffer_t* const out_buffer,const int string_length,const char* const string) {
  OUTPUT_BUFFER_CHECK(out_buffer);
  memcpy(out_buffer->buffer_cursor,string,string_length);
  out_buffer->buffer_cursor += string_length;
  out_buffer->buffer_used += string_length;
}

