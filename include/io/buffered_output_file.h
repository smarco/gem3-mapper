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

#ifndef BUFFERED_OUTPUT_FILE_H_
#define BUFFERED_OUTPUT_FILE_H_

#include "utils/essentials.h"
#include "io/output_buffer.h"
#include "io/output_file.h"

typedef struct {
  /* Output file */
  output_file_t* output_file;
  /* Output Buffer */
  output_buffer_t* buffer;
} buffered_output_file_t;

// Codes status
#define BUFFERED_OUTPUT_FILE_OK 0
#define BUFFERED_OUTPUT_FILE_FAIL -1

/*
 * Setup
 */
buffered_output_file_t* buffered_output_file_new(output_file_t* const output_file);
void buffered_output_file_close(buffered_output_file_t* const buffered_output);

/*
 * Utils
 */
void buffered_output_file_request_buffer(
    buffered_output_file_t* const buffered_output,
    const uint32_t block_id);
void buffered_output_file_dump_buffer(buffered_output_file_t* const buffered_output);
void buffered_output_file_safety_dump_buffer(buffered_output_file_t* const buffered_output);
void buffered_output_file_reserve(
    buffered_output_file_t* const buffered_output,
    const uint64_t num_chars);

/*
 * Fast-printer
 */
void bofprintf_uint64(buffered_output_file_t* const buffered_output,const uint64_t number);
void bofprintf_int64(buffered_output_file_t* const buffered_output,const int64_t number);
void bofprintf_char(buffered_output_file_t* const buffered_output,const char character);
void bofprintf_string(
    buffered_output_file_t* const buffered_output,
    const int string_length,
    const char* const string);
#define bofprintf_string_literal(buffered_output,string) bofprintf_string(buffered_output,sizeof(string)-1,string)

#endif /* BUFFERED_OUTPUT_FILE_H_ */
