/*
 * PROJECT: GEM-Tools library
 * FILE: gt_buffered_input_file.h
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_BUFFERED_INPUT_FILE_H_
#define GT_BUFFERED_INPUT_FILE_H_

#include "gt_essentials.h"
#include "gt_input_file.h"
#include "gt_template.h"
#include "gt_buffered_output_file.h"

typedef struct {
  /* Input file */
  gt_input_file* input_file;
  /* Block buffer and cursors */
  uint32_t block_id;
  gt_vector* block_buffer;
  char* cursor;
  uint64_t lines_in_buffer;
  uint64_t current_line_num;
  /* Attached output buffer */
  gt_vector* attached_buffered_output_file; /* (gt_buffered_output_file*) */
} gt_buffered_input_file;

/*
 * Checkers
 */
#define GT_BUFFERED_INPUT_FILE_CHECK(buffered_input_file) \
  GT_NULL_CHECK(buffered_input_file); \
  GT_NULL_CHECK(buffered_input_file->input_file); \
  GT_NULL_CHECK(buffered_input_file->block_buffer); \
  GT_NULL_CHECK(buffered_input_file->cursor)

/*
 * Buffered Input File Handlers
 */
gt_buffered_input_file* gt_buffered_input_file_new(gt_input_file* const input_file);
gt_status gt_buffered_input_file_close(gt_buffered_input_file* const buffered_input_file);
GT_INLINE const char** const gt_buffered_input_file_get_text_line(gt_buffered_input_file* const buffered_input_file);
GT_INLINE uint64_t gt_buffered_input_file_get_cursor_pos(gt_buffered_input_file* const buffered_input_file);
GT_INLINE bool gt_buffered_input_file_eob(gt_buffered_input_file* const buffered_input_file);
GT_INLINE gt_status gt_buffered_input_file_get_block(
    gt_buffered_input_file* const buffered_input_file,const uint64_t num_lines);
GT_INLINE gt_status gt_buffered_input_file_add_lines_to_block(
    gt_buffered_input_file* const buffered_input_file,const uint64_t num_lines);

/*
 * BufferedInputFile Utils
 */
GT_INLINE void gt_buffered_input_file_skip_line(gt_buffered_input_file* const buffered_input);
GT_INLINE gt_status gt_buffered_input_file_reload(gt_buffered_input_file* const buffered_input,const uint64_t num_lines);

/*
 * Block Synchronization with Output
 */
GT_INLINE void gt_buffered_input_file_attach_buffered_output(
    gt_buffered_input_file* const buffered_input_file,gt_buffered_output_file* const buffered_output_file);
GT_INLINE void gt_buffered_input_file_dump_attached_buffers(gt_vector* const attached_buffered_output_file);
GT_INLINE void gt_buffered_input_file_set_id_attached_buffers(gt_vector* const attached_buffered_output_file,const uint64_t block_id);

#endif /* GT_BUFFERED_INPUT_FILE_H_ */
