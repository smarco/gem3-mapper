/*
 * PROJECT: GEMMapper
 * FILE: buffered_output_file.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "buffered_output_file.h"

#define BUFFERED_OUTPUT_FILE_FORCE_DUMP_SIZE BUFFER_SIZE_4M

/*
 * Setup
 */
buffered_output_file_t* buffered_output_file_new(output_file_t* const output_file) {
  OUTPUT_FILE_CHECK(output_file);
  buffered_output_file_t* buffered_output = mm_alloc(buffered_output_file_t);
  // Initialize the bof
  buffered_output->output_file = output_file;
  buffered_output->buffer = NULL;
  return buffered_output;
}
void buffered_output_file_close(buffered_output_file_t* const buffered_output) {
  BUFFERED_OUTPUT_FILE_CHECK(buffered_output);
  buffered_output_file_dump_buffer(buffered_output);
  mm_free(buffered_output);
}
/*
 * Utils
 */
GEM_INLINE void buffered_output_file_request_buffer(
    buffered_output_file_t* const buffered_output,const uint32_t block_id) {
  buffered_output->buffer = output_file_request_buffer(buffered_output->output_file,block_id);
}
GEM_INLINE void buffered_output_file_dump_buffer(buffered_output_file_t* const buffered_output) {
  BUFFERED_OUTPUT_FILE_CHECK(buffered_output);
  // Dump
  if (buffered_output->buffer != NULL) {
    output_file_return_buffer(buffered_output->output_file,buffered_output->buffer);
  }
  buffered_output->buffer = NULL;
}
GEM_INLINE void buffered_output_file_safety_dump_buffer(buffered_output_file_t* const buffered_output) {
  BUFFERED_OUTPUT_FILE_CHECK(buffered_output);
  buffered_output->buffer = output_file_request_buffer_extension(buffered_output->output_file,buffered_output->buffer);
  gem_cond_fatal_error(buffered_output->buffer==NULL,BUFFER_SAFETY_DUMP);
}
GEM_INLINE void buffered_output_file_reserve(buffered_output_file_t* const buffered_output,const uint64_t num_chars) {
  BUFFERED_OUTPUT_FILE_CHECK(buffered_output);
  const uint64_t buffer_used = output_buffer_get_used(buffered_output->buffer);
  if (gem_expect_false(buffer_used+num_chars >= buffered_output->output_file->buffer_size)) {
    buffered_output_file_safety_dump_buffer(buffered_output);
  }
  output_buffer_reserve(buffered_output->buffer,num_chars);
}
/*
 * Printers (Disabled for performance issues)
 */
//GEM_INLINE int vbofprintf(buffered_output_file_t* const buffered_output,const char *template,va_list v_args) {
//  BUFFERED_OUTPUT_FILE_CHECK(buffered_output);
//  GEM_CHECK_NULL(template);
//  if (gem_expect_false(output_buffer_get_used(buffered_output->buffer)+BUFFER_SIZE_1K*10
//      >= buffered_output->output_file->buffer_size)) {
//    buffered_output_file_safety_dump_buffer(buffered_output);
//  }
//  return vbprintf(buffered_output->buffer,template,v_args);
//}
//GEM_INLINE int bofprintf(buffered_output_file_t* const buffered_output,const char *template,...) {
//  BUFFERED_OUTPUT_FILE_CHECK(buffered_output);
//  GEM_CHECK_NULL(template);
//  va_list v_args;
//  va_start(v_args,template);
//  const int chars_printed = vbofprintf(buffered_output,template,v_args);
//  va_end(v_args);
//  return chars_printed;
//}
//GEM_INLINE int vbofprintf_fixed(
//    buffered_output_file_t* const buffered_output,
//    const uint64_t expected_mem_usage,const char *template,va_list v_args) {
//  BUFFERED_OUTPUT_FILE_CHECK(buffered_output);
//  GEM_CHECK_ZERO(expected_mem_usage);
//  GEM_CHECK_NULL(template);
//  if (gem_expect_false(output_buffer_get_used(buffered_output->buffer)+BUFFER_SIZE_1K*10
//      >= buffered_output->output_file->buffer_size)) {
//    buffered_output_file_safety_dump_buffer(buffered_output);
//  }
//  return vbprintf_fixed(buffered_output->buffer,expected_mem_usage,template,v_args);
//}
//GEM_INLINE int bofprintf_fixed(
//    buffered_output_file_t* const buffered_output,
//    const uint64_t expected_mem_usage,const char *template,...) {
//  BUFFERED_OUTPUT_FILE_CHECK(buffered_output);
//  GEM_CHECK_ZERO(expected_mem_usage);
//  GEM_CHECK_NULL(template);
//  va_list v_args;
//  va_start(v_args,template);
//  const int chars_printed = vbofprintf_fixed(buffered_output,expected_mem_usage,template,v_args);
//  va_end(v_args);
//  return chars_printed;
//}
/*
 * Fast-printer
 */
GEM_INLINE void bofprintf_uint64(buffered_output_file_t* const buffered_output,const uint64_t number) {
  bprintf_uint64(buffered_output->buffer,number);
}
GEM_INLINE void bofprintf_int64(buffered_output_file_t* const buffered_output,const uint64_t number) {
  bprintf_int64(buffered_output->buffer,number);
}
GEM_INLINE void bofprintf_char(buffered_output_file_t* const buffered_output,const char character) {
  bprintf_char(buffered_output->buffer,character);
}
GEM_INLINE void bofprintf_string(
    buffered_output_file_t* const buffered_output,const int string_length,const char* const string) {
  bprintf_buffer(buffered_output->buffer,string_length,string);
}

