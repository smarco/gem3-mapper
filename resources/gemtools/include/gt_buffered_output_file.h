/*
 * PROJECT: GEM-Tools library
 * FILE: gt_buffered_output_file.h
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_BUFFERED_OUTPUT_FILE_H_
#define GT_BUFFERED_OUTPUT_FILE_H_

#include "gt_essentials.h"
#include "gt_output_buffer.h"
#include "gt_output_file.h"

typedef struct {
  /* Output file */
  gt_output_file* output_file;
  /* Output Buffer */
  gt_output_buffer* buffer;
} gt_buffered_output_file;

// Codes gt_status
#define GT_BUFFERED_OUTPUT_FILE_OK 0
#define GT_BUFFERED_OUTPUT_FILE_FAIL -1

/*
 * Checkers
 */
#define GT_BUFFERED_OUTPUT_FILE_CHECK(buffered_output_file) \
  GT_NULL_CHECK(buffered_output_file); \
  GT_OUTPUT_FILE_CHECK(buffered_output_file->output_file); \
  GT_OUTPUT_BUFFER_CHECK(buffered_output_file->buffer)

/*
 * Setup
 */
gt_buffered_output_file* gt_buffered_output_file_new(gt_output_file* const output_file);
void gt_buffered_output_file_close(gt_buffered_output_file* const buffered_output_file);

/*
 * Accessors
 */
GT_INLINE void gt_buffered_output_file_get_block_ids(
    gt_buffered_output_file* const buffered_output_file,uint32_t* const mayor_id,uint32_t* const minor_id);
GT_INLINE void gt_buffered_output_file_set_block_ids(
    gt_buffered_output_file* const buffered_output_file,const uint32_t mayor_id,const uint32_t minor_id);
GT_INLINE gt_output_buffer* gt_buffered_output_file_get_buffer(gt_buffered_output_file* const buffered_output_file);
GT_INLINE void gt_buffered_output_file_set_buffer(
    gt_buffered_output_file* const buffered_output_file,gt_output_buffer* const output_buffer);

/*
 * Dump
 */
GT_INLINE void gt_buffered_output_file_dump(gt_buffered_output_file* const buffered_output_file);
GT_INLINE void gt_buffered_output_file_safety_dump(gt_buffered_output_file* const buffered_output_file);

/*
 * Buffered Output File Printers
 */
GT_INLINE gt_status gt_vbofprintf(gt_buffered_output_file* const buffered_output_file,const char *template,va_list v_args);
GT_INLINE gt_status gt_bofprintf(gt_buffered_output_file* const buffered_output_file,const char *template,...);

#endif /* GT_BUFFERED_OUTPUT_FILE_H_ */
