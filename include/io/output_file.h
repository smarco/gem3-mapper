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
 *   Output file data structure enables writing to a file
 *     Supports internal output-buffer queue for ordered output
 *     Supports deferred output-buffer dump
 *     Supports compressed (BZIP/GZIP) output
 */

#ifndef OUTPUT_FILE_H_
#define OUTPUT_FILE_H_

#include "utils/essentials.h"
#include "io/output_buffer.h"

/*
 * Output format
 */
typedef enum { FASTA, MAP, SAM, FILE_FORMAT_UNKNOWN } file_format_t;

/*
 * Output file
 */
typedef struct {
  /* Output file */
  fm_t* file_manager;
  /* Output Buffers */
  uint64_t num_buffers;
  uint64_t buffer_size;
  output_buffer_t** buffer;
  uint64_t buffer_free;
  uint64_t buffer_write_pending;
  /* Block ID prioritization (for synchronization purposes) */
  pqueue_t* buffer_requests;
  uint32_t next_output_mayor_id;
  uint32_t next_output_minor_id;
  /* Mutexes */
  pthread_cond_t  request_buffer_cond;
  pthread_mutex_t output_file_mutex;
} output_file_t;

// Codes status
#define OUTPUT_FILE_OK 0
#define OUTPUT_FILE_FAIL -1

/*
 * Setup
 */
output_file_t* output_file_new(
    char* const file_name,
    const uint64_t max_output_buffers,
    const uint64_t buffer_size);
output_file_t* output_stream_new(
    FILE* const stream,
    const uint64_t max_output_buffers,
    const uint64_t buffer_size);
output_file_t* output_gzip_stream_new(
    FILE* const stream,
    const uint64_t max_output_buffers,
    const uint64_t buffer_size);
output_file_t* output_bzip_stream_new(
    FILE* const stream,
    const uint64_t max_output_buffers,
    const uint64_t buffer_size);
void output_file_close(output_file_t* const out_file);

/*
 * Utils
 */
output_buffer_t* output_file_request_buffer(
    output_file_t* const output_file,
    const uint64_t block_id);
output_buffer_t* output_file_request_buffer_extension(
    output_file_t* const output_file,
    output_buffer_t* const output_buffer);
void output_file_return_buffer(
    output_file_t* const output_file,
    output_buffer_t* const output_buffer);

/*
 * Output File Printers
 */
int vofprintf(output_file_t* const out_file,const char *template,va_list v_args);
int ofprintf(output_file_t* const out_file,const char *template,...);

/*
 * Error Messages
 */
#define GEM_ERROR_OUTPUT_FILE_INCONSISTENCY "Output file state inconsistent"

#endif /* OUTPUT_FILE_H_ */
