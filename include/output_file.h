/*
 * PROJECT: GEMMapper
 * FILE: output_file.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef OUTPUT_FILE_H_
#define OUTPUT_FILE_H_

#include "essentials.h"
#include "output_buffer.h"

#define MAX_OUTPUT_BUFFERS 25

typedef enum { SORTED_FILE, UNSORTED_FILE } output_file_type;
typedef struct {
  /* Output file */
  fm_t* file_manager;
  output_file_type file_type;
  /* Output Buffers */
  output_buffer_t* buffer[MAX_OUTPUT_BUFFERS];
  uint64_t buffer_busy;
  uint64_t buffer_write_pending;
  /* Block ID (for synchronization purposes) */
  uint32_t mayor_block_id;
  uint32_t minor_block_id;
  /* Mutexes */
  pthread_cond_t  out_buffer_cond;
  pthread_cond_t  out_write_cond;
  pthread_mutex_t out_file_mutex;
} output_file_t;

// Codes status
#define OUTPUT_FILE_OK 0
#define OUTPUT_FILE_FAIL -1

/*
 * Checkers
 */
#define OUTPUT_FILE_CHECK(output_file) \
  GEM_CHECK_NULL(output_file); \
  FM_CHECK((output_file)->file_manager); \
  GEM_CHECK_NULL((output_file)->buffer)

#define OUTPUT_FILE_CONSISTENCY_CHECK(output_file) \
  OUTPUT_FILE_CHECK(output_file); \
  gem_fatal_check( \
    output_file->buffer_busy>MAX_OUTPUT_BUFFERS || \
    output_file->buffer_write_pending>MAX_OUTPUT_BUFFERS,OUTPUT_FILE_INCONSISTENCY)

/*
 * Output File Setup
 */
output_file_t* output_stream_new(FILE* const stream,const output_file_type output_file_type);
output_file_t* output_gzip_stream_new(FILE* const stream,const output_file_type output_file_type);
output_file_t* output_bzip_stream_new(FILE* const stream,const output_file_type output_file_type);
output_file_t* output_file_new(char* const file_name,const output_file_type output_file_type);
void output_file_close(output_file_t* const out_file);

/*
 * Internal Buffers Accessors
 */
GEM_INLINE output_buffer_t* output_file_request_buffer(output_file_t* const output_file);
GEM_INLINE void output_file_return_buffer(
    output_file_t* const output_file,output_buffer_t* const output_buffer);
GEM_INLINE output_buffer_t* output_file_dump_buffer(
    output_file_t* const output_file,output_buffer_t* const output_buffer,const bool asynchronous);

/*
 * Error Messages
 */
#define GEM_ERROR_OUTPUT_FILE_INCONSISTENCY "Output file state inconsistent"

#endif /* OUTPUT_FILE_H_ */
