/*
 * PROJECT: GEMMapper
 * FILE: output_file.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef OUTPUT_FILE_H_
#define OUTPUT_FILE_H_

#include "utils/essentials.h"
#include "io/output_buffer.h"

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
