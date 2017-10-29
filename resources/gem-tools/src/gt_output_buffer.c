/*
 * PROJECT: GEM-Tools library
 * FILE: gt_output_buffer.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_output_buffer.h"

#define GT_OUTPUT_BUFFER_INITIAL_SIZE GT_BUFFER_SIZE_16M

/*
 * Setup
 */
GT_INLINE gt_output_buffer* gt_output_buffer_new(void) {
  gt_output_buffer* output_buffer = gt_alloc(gt_output_buffer);
  output_buffer->buffer=gt_vector_new(GT_OUTPUT_BUFFER_INITIAL_SIZE,sizeof(char));
  gt_output_buffer_initiallize(output_buffer,GT_OUTPUT_BUFFER_FREE);
  return output_buffer;
}
GT_INLINE void gt_output_buffer_clear(gt_output_buffer* const output_buffer) {
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  output_buffer->mayor_block_id=UINT32_MAX;
  output_buffer->minor_block_id=0;
  output_buffer->is_final_block=true;
  gt_vector_clear(output_buffer->buffer);
}
GT_INLINE void gt_output_buffer_initiallize(gt_output_buffer* const output_buffer,const gt_output_buffer_state buffer_state) {
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  gt_output_buffer_clear(output_buffer);
  gt_output_buffer_set_state(output_buffer,buffer_state);
}
GT_INLINE void gt_output_buffer_delete(gt_output_buffer* const output_buffer) {
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  gt_vector_delete(output_buffer->buffer);
  gt_free(output_buffer);
}

/*
 * Accessors
 */
GT_INLINE void gt_output_buffer_set_state(gt_output_buffer* const output_buffer,const gt_output_buffer_state buffer_state) {
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  output_buffer->buffer_state=buffer_state;
}
GT_INLINE gt_output_buffer_state gt_output_buffer_get_state(gt_output_buffer* const output_buffer) {
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  return output_buffer->buffer_state;
}
GT_INLINE void gt_output_buffer_set_partial_block(gt_output_buffer* const output_buffer) {
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  output_buffer->is_final_block=false;
}
GT_INLINE void gt_output_buffer_set_mayor_block_id(gt_output_buffer* const output_buffer,const uint32_t mayor_block_id) {
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  output_buffer->mayor_block_id=mayor_block_id;
}
GT_INLINE uint32_t gt_output_buffer_get_mayor_block_id(gt_output_buffer* const output_buffer) {
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  return output_buffer->mayor_block_id;
}
GT_INLINE void gt_output_buffer_set_minor_block_id(gt_output_buffer* const output_buffer,const uint32_t minor_block_id) {
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  output_buffer->minor_block_id=minor_block_id;
}
GT_INLINE uint32_t gt_output_buffer_get_minor_block_id(gt_output_buffer* const output_buffer) {
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  return output_buffer->minor_block_id;
}
GT_INLINE void gt_output_buffer_inc_minor_block_id(gt_output_buffer* const output_buffer) {
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  ++output_buffer->minor_block_id;
}
GT_INLINE uint64_t gt_output_buffer_get_used(gt_output_buffer* const output_buffer) {
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  return gt_vector_get_used(output_buffer->buffer);
}

/*
 * Adaptors
 */
GT_INLINE char* gt_output_buffer_to_char(gt_output_buffer* const output_buffer) {
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  gt_vector_insert(output_buffer->buffer,EOS,char);
  return gt_vector_get_mem(output_buffer->buffer,char);
}
GT_INLINE gt_vector* gt_output_buffer_to_vchar(gt_output_buffer* const output_buffer) {
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  return output_buffer->buffer;
}

/*
 * Buffer printer
 */
GT_INLINE gt_status gt_vbprintf(gt_output_buffer* const output_buffer,const char *template,va_list v_args) {
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  GT_NULL_CHECK(template); GT_NULL_CHECK(v_args);
  // Calculate buffer size required to dump template{v_args}
  int64_t mem_required = gt_calculate_memory_required_v(template,v_args);
  return gt_vbprintf_(output_buffer,mem_required,template,v_args);
}
GT_INLINE gt_status gt_bprintf(gt_output_buffer* const output_buffer,const char *template,...) {
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  GT_NULL_CHECK(template);
  va_list v_args;
  va_start(v_args,template);
  const gt_status chars_printed = gt_vbprintf(output_buffer,template,v_args);
  va_end(v_args);
  return chars_printed;
}
GT_INLINE gt_status gt_vbprintf_(
    gt_output_buffer* const output_buffer,const uint64_t expected_mem_usage,
    const char *template,va_list v_args) {
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  GT_NULL_CHECK(template); GT_NULL_CHECK(v_args);
  // Calculate buffer size required to dump template{v_args}
  gt_vector_reserve_additional(output_buffer->buffer,expected_mem_usage);
  const gt_status chars_printed=vsprintf(gt_vector_get_free_elm(output_buffer->buffer,char),template,v_args);
  if (gt_expect_true(chars_printed>=0)) {
    gt_vector_add_used(output_buffer->buffer,chars_printed);
  }
  return chars_printed;
}
GT_INLINE gt_status gt_bprintf_(
    gt_output_buffer* const output_buffer,const uint64_t expected_mem_usage,const char *template,...) {
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  GT_NULL_CHECK(template);
  va_list v_args;
  va_start(v_args,template);
  const gt_status chars_printed = gt_vbprintf_(output_buffer,expected_mem_usage,template,v_args);
  va_end(v_args);
  return chars_printed;
}
