/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   Buffered output file provides buffered-writing from an output file.
 *     Preserves input-block - output-block ordering
 *     Enables deferred block dump to output file (non-blocking)
 *     Lazy output buffer allocation
 */

#include "io/buffered_output_file.h"

/*
 * Error Messages
 */
#define GEM_ERROR_BUFFER_SAFETY_DUMP "Output buffer. Could not perform safety dump (no buffer extension)"
#define GEM_ERROR_BUFFER_RESERVE "Output buffer. Could not reserve buffer-memory for dump"

/*
 * Profile
 */
#define PROFILE_LEVEL PLOW

/*
 * Constants
 */
#define BUFFERED_OUTPUT_FILE_FORCE_DUMP_SIZE BUFFER_SIZE_4M

/*
 * Setup
 */
buffered_output_file_t* buffered_output_file_new(output_file_t* const output_file) {
  buffered_output_file_t* buffered_output = mm_alloc(buffered_output_file_t);
  // Initialize the bof
  buffered_output->output_file = output_file;
  buffered_output->buffer = NULL;
  return buffered_output;
}
void buffered_output_file_close(buffered_output_file_t* const buffered_output) {
  buffered_output_file_dump_buffer(buffered_output);
  mm_free(buffered_output);
}
/*
 * Utils
 */
void buffered_output_file_request_buffer(
    buffered_output_file_t* const buffered_output,
    const uint32_t block_id) {
  buffered_output->buffer = output_file_request_buffer(buffered_output->output_file,block_id);
}
void buffered_output_file_dump_buffer(buffered_output_file_t* const buffered_output) {
  PROFILE_START(GP_BUFFERED_OUTPUT_DUMP,PROFILE_LEVEL);
  // Dump
  if (buffered_output->buffer != NULL) {
    output_file_return_buffer(buffered_output->output_file,buffered_output->buffer);
  }
  buffered_output->buffer = NULL;
  PROFILE_STOP(GP_BUFFERED_OUTPUT_DUMP,PROFILE_LEVEL);
}
void buffered_output_file_safety_dump_buffer(buffered_output_file_t* const buffered_output) {
  buffered_output->buffer = output_file_request_buffer_extension(buffered_output->output_file,buffered_output->buffer);
  gem_cond_fatal_error(buffered_output->buffer==NULL,BUFFER_SAFETY_DUMP);
}
void buffered_output_file_reserve(
    buffered_output_file_t* const buffered_output,
    const uint64_t num_chars) {
  const uint64_t buffer_size = buffered_output->output_file->buffer_size;
  gem_cond_fatal_error(num_chars>buffer_size,BUFFER_RESERVE);
  const uint64_t buffer_used = output_buffer_get_used(buffered_output->buffer);
  if (gem_expect_false(buffer_used+num_chars > buffer_size)) {
    buffered_output_file_safety_dump_buffer(buffered_output);
  }
}
/*
 * Printers (Disabled for performance issues)
 */
//int vbofprintf(buffered_output_file_t* const buffered_output,const char *template,va_list v_args) {
//  GEM_CHECK_NULL(template);
//  if (gem_expect_false(output_buffer_get_used(buffered_output->buffer)+BUFFER_SIZE_1K*10
//      >= buffered_output->output_file->buffer_size)) {
//    buffered_output_file_safety_dump_buffer(buffered_output);
//  }
//  return vbprintf(buffered_output->buffer,template,v_args);
//}
//int bofprintf(buffered_output_file_t* const buffered_output,const char *template,...) {
//  GEM_CHECK_NULL(template);
//  va_list v_args;
//  va_start(v_args,template);
//  const int chars_printed = vbofprintf(buffered_output,template,v_args);
//  va_end(v_args);
//  return chars_printed;
//}
//int vbofprintf_fixed(
//    buffered_output_file_t* const buffered_output,
//    const uint64_t expected_mem_usage,const char *template,va_list v_args) {
//  GEM_CHECK_ZERO(expected_mem_usage);
//  GEM_CHECK_NULL(template);
//  if (gem_expect_false(output_buffer_get_used(buffered_output->buffer)+BUFFER_SIZE_1K*10
//      >= buffered_output->output_file->buffer_size)) {
//    buffered_output_file_safety_dump_buffer(buffered_output);
//  }
//  return vbprintf_fixed(buffered_output->buffer,expected_mem_usage,template,v_args);
//}
//int bofprintf_fixed(
//    buffered_output_file_t* const buffered_output,
//    const uint64_t expected_mem_usage,const char *template,...) {
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
void bofprintf_uint64(buffered_output_file_t* const buffered_output,const uint64_t number) {
  bprintf_uint64(buffered_output->buffer,number);
}
void bofprintf_int64(buffered_output_file_t* const buffered_output,const int64_t number) {
  bprintf_int64(buffered_output->buffer,number);
}
void bofprintf_char(buffered_output_file_t* const buffered_output,const char character) {
  bprintf_char(buffered_output->buffer,character);
}
void bofprintf_string(
    buffered_output_file_t* const buffered_output,
    const int string_length,
    const char* const string) {
  bprintf_buffer(buffered_output->buffer,string_length,string);
}

