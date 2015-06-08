/*
 * PROJECT: GEM-Tools library
 * FILE: gt_output_file.h
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_OUTPUT_FILE_H_
#define GT_OUTPUT_FILE_H_

#include "gt_essentials.h"
#include "gt_output_buffer.h"

#define GT_MAX_OUTPUT_BUFFERS 25

typedef enum { SORTED_FILE, UNSORTED_FILE } gt_output_file_type;
typedef struct {
  /* Output file */
  gt_fm* file_manager;
  gt_output_file_type file_type;
  /* Output Buffers */
  gt_output_buffer* buffer[GT_MAX_OUTPUT_BUFFERS];
  uint64_t buffer_busy;
  uint64_t buffer_write_pending;
  /* Block ID (for synchronization purposes) */
  uint32_t mayor_block_id;
  uint32_t minor_block_id;
  /* Mutexes */
  pthread_cond_t  out_buffer_cond;
  pthread_cond_t  out_write_cond;
  pthread_mutex_t out_file_mutex;
} gt_output_file;

// Codes gt_status
#define GT_OUTPUT_FILE_OK 0
#define GT_OUTPUT_FILE_FAIL -1

/*
 * Checkers
 */
#define GT_OUTPUT_FILE_CHECK(output_file) \
  GT_NULL_CHECK(output_file); \
  GT_NULL_CHECK((output_file)->buffer); \
  GT_FM_CHECK((output_file)->file_manager)

#define GT_OUTPUT_FILE_CONSISTENCY_CHECK(output_file) \
  GT_OUTPUT_FILE_CHECK(output_file); \
  gt_fatal_check( \
    output_file->buffer_busy>GT_MAX_OUTPUT_BUFFERS || \
    output_file->buffer_write_pending>GT_MAX_OUTPUT_BUFFERS,OUTPUT_FILE_INCONSISTENCY)

/*
 * Output File Setup
 */
gt_output_file* gt_output_file_new(char* const file_name,const gt_output_file_type output_file_type);
gt_output_file* gt_output_stream_new(FILE* const stream,const gt_output_file_type output_file_type);
gt_output_file* gt_output_gzip_stream_new(FILE* const stream,const gt_output_file_type output_file_type);
gt_output_file* gt_output_bzip_stream_new(FILE* const stream,const gt_output_file_type output_file_type);
gt_output_file* gt_output_file_new(char* const file_name,const gt_output_file_type output_file_type);
void gt_output_file_close(gt_output_file* const output_file);

/*
 * Output File Printers
 */
GT_INLINE gt_status gt_vofprintf(gt_output_file* const output_file,const char *template,va_list v_args);
GT_INLINE gt_status gt_ofprintf(gt_output_file* const output_file,const char *template,...);

/*
 * Internal Buffers Accessors
 */
GT_INLINE gt_output_buffer* gt_output_file_request_buffer(gt_output_file* const output_file);
GT_INLINE void gt_output_file_release_buffer(
    gt_output_file* const output_file,gt_output_buffer* const output_buffer);
GT_INLINE gt_output_buffer* gt_output_file_dump_buffer(
    gt_output_file* const output_file,gt_output_buffer* const output_buffer,const bool asynchronous);

/*
 * Error Messages. Output File
 */
#define GT_ERROR_OUTPUT_FILE_INCONSISTENCY "Output file state inconsistent"
#define GT_ERROR_OUTPUT_FILE_FAIL_WRITE "Output file. Error writing to to file"
#define GT_ERROR_BUFFER_SAFETY_DUMP "Output buffer. Could not perform safety dump"

#endif /* GT_OUTPUT_FILE_H_ */
