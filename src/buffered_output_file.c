/*
 * PROJECT: GEMMapper
 * FILE: buffered_output_file.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "buffered_output_file.h"

#define BUFFERED_OUTPUT_FILE_FORCE_DUMP_SIZE BUFFER_SIZE_64M

/*
 * Setup
 */
buffered_output_file_t* buffered_output_file_new(output_file_t* const output_file) {
  OUTPUT_FILE_CHECK(output_file);
  buffered_output_file_t* buffered_output = mm_alloc(buffered_output_file_t);
  // Initialize the bof
  buffered_output->output_file = output_file;
  buffered_output->buffer = output_file_request_buffer(buffered_output->output_file);
  return buffered_output;
}
void buffered_output_file_close(buffered_output_file_t* const buffered_output) {
  BUFFERED_OUTPUT_FILE_CHECK(buffered_output);
  if (buffered_output->buffer != NULL) {
    buffered_output_file_dump(buffered_output);
    output_file_return_buffer(buffered_output->output_file,buffered_output->buffer);
  }
  mm_free(buffered_output);
}
/*
 * Accessors
 */
GEM_INLINE void buffered_output_file_get_block_ids(
    buffered_output_file_t* const buffered_output,uint32_t* const mayor_id,uint32_t* const minor_id) {
  BUFFERED_OUTPUT_FILE_CHECK(buffered_output);
  *mayor_id = output_buffer_get_mayor_block_id(buffered_output->buffer);
  *minor_id = output_buffer_get_minor_block_id(buffered_output->buffer);
}
GEM_INLINE void buffered_output_file_set_block_ids(
    buffered_output_file_t* const buffered_output,const uint32_t mayor_id,const uint32_t minor_id) {
  BUFFERED_OUTPUT_FILE_CHECK(buffered_output);
  output_buffer_set_mayor_block_id(buffered_output->buffer,mayor_id);
  output_buffer_set_minor_block_id(buffered_output->buffer,minor_id);
}
GEM_INLINE output_buffer_t* buffered_output_file_get_buffer(buffered_output_file_t* const buffered_output) {
  BUFFERED_OUTPUT_FILE_CHECK(buffered_output);
  return buffered_output->buffer;
}
GEM_INLINE void buffered_output_file_set_buffer(
    buffered_output_file_t* const buffered_output,output_buffer_t* const out_buffer) {
  BUFFERED_OUTPUT_FILE_CHECK(buffered_output);
  buffered_output->buffer = out_buffer;
}
/*
 * Buffered Output Dump
 */
GEM_INLINE void buffered_output_file_dump(buffered_output_file_t* const buffered_output) {
  BUFFERED_OUTPUT_FILE_CHECK(buffered_output);
  // Skip empty-buffers with no ID
  if (buffered_output->buffer->mayor_block_id==UINT32_MAX &&
      output_buffer_get_used(buffered_output->buffer) == 0) return;
  // Dump
  buffered_output->buffer = output_file_dump_buffer(buffered_output->output_file,buffered_output->buffer,true);
  gem_cond_fatal_error(buffered_output->buffer==NULL,BUFFER_DUMP);
}
GEM_INLINE void buffered_output_file_safety_dump(buffered_output_file_t* const buffered_output) {
  BUFFERED_OUTPUT_FILE_CHECK(buffered_output);
  output_buffer_set_partial_block(buffered_output->buffer);
  const uint32_t mayor_id = output_buffer_get_mayor_block_id(buffered_output->buffer);
  const uint32_t minor_id = output_buffer_get_minor_block_id(buffered_output->buffer);
  buffered_output->buffer = output_file_dump_buffer(buffered_output->output_file,buffered_output->buffer,false);
  gem_cond_fatal_error(buffered_output->buffer==NULL,BUFFER_SAFETY_DUMP);
  output_buffer_set_mayor_block_id(buffered_output->buffer,mayor_id);
  output_buffer_set_minor_block_id(buffered_output->buffer,minor_id+1);
}
/*
 * Buffered Output File Printers
 */
GEM_INLINE int vbofprintf(buffered_output_file_t* const buffered_output,const char *template,va_list v_args) {
  BUFFERED_OUTPUT_FILE_CHECK(buffered_output);
  GEM_CHECK_NULL(template);
  if (gem_expect_false(
      output_buffer_get_used(buffered_output->buffer) >= BUFFERED_OUTPUT_FILE_FORCE_DUMP_SIZE)) {
    buffered_output_file_safety_dump(buffered_output);
  }
  return vbprintf(buffered_output->buffer,template,v_args);
}
GEM_INLINE int bofprintf(buffered_output_file_t* const buffered_output,const char *template,...) {
  BUFFERED_OUTPUT_FILE_CHECK(buffered_output);
  GEM_CHECK_NULL(template);
  va_list v_args;
  va_start(v_args,template);
  const int chars_printed = vbofprintf(buffered_output,template,v_args);
  va_end(v_args);
  return chars_printed;
}
GEM_INLINE int vbofprintf_fixed(
    buffered_output_file_t* const buffered_output,
    const uint64_t expected_mem_usage,const char *template,va_list v_args) {
  BUFFERED_OUTPUT_FILE_CHECK(buffered_output);
  GEM_CHECK_ZERO(expected_mem_usage);
  GEM_CHECK_NULL(template);
  if (gem_expect_false(
      output_buffer_get_used(buffered_output->buffer) >= BUFFERED_OUTPUT_FILE_FORCE_DUMP_SIZE)) {
    buffered_output_file_safety_dump(buffered_output);
  }
  return vbprintf_fixed(buffered_output->buffer,expected_mem_usage,template,v_args);
}
GEM_INLINE int bofprintf_fixed(
    buffered_output_file_t* const buffered_output,
    const uint64_t expected_mem_usage,const char *template,...) {
  BUFFERED_OUTPUT_FILE_CHECK(buffered_output);
  GEM_CHECK_ZERO(expected_mem_usage);
  GEM_CHECK_NULL(template);
  va_list v_args;
  va_start(v_args,template);
  const int chars_printed = vbofprintf_fixed(buffered_output,expected_mem_usage,template,v_args);
  va_end(v_args);
  return chars_printed;
}

