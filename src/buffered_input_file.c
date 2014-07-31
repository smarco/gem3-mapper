/*
 * PROJECT: GEMMapper
 * FILE: buffered_input_file.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "buffered_input_file.h"

#define BMI_BUFFER_SIZE BUFFER_SIZE_4M
#define BMI_NUM_LINES   NUM_LINES_5K

/*
 * Buffered map file handlers
 */
buffered_input_file_t* buffered_input_file_new(input_file_t* const in_file) {
  INPUT_FILE_CHECK(in_file);
  buffered_input_file_t* buffered_input = mm_alloc(buffered_input_file_t);
  /* Input file */
  buffered_input->input_file = in_file;
  /* Block buffer and cursors */
  buffered_input->block_id = UINT32_MAX;
  buffered_input->block_buffer = vector_new(BMI_BUFFER_SIZE,sizeof(uint8_t));
  buffered_input->cursor = (char*)vector_get_mem(buffered_input->block_buffer,uint8_t);
  buffered_input->current_line_num = UINT64_MAX;
  /* Attached output buffer */
  buffered_input->attached_buffered_output_file = NULL;
  return buffered_input;
}
void buffered_input_file_close(buffered_input_file_t* const buffered_input) {
  BUFFERED_INPUT_FILE_CHECK(buffered_input);
  vector_delete(buffered_input->block_buffer);
  mm_free(buffered_input);
}
/*
 * Accessors
 */
GEM_INLINE char** const buffered_input_file_get_text_line(buffered_input_file_t* const buffered_input) {
  return ((char** const)(&buffered_input->cursor));
}
GEM_INLINE uint64_t buffered_input_file_get_cursor_pos(buffered_input_file_t* const buffered_input) {
  BUFFERED_INPUT_FILE_CHECK(buffered_input);
  GEM_CHECK_NULL(buffered_input->cursor);
  return buffered_input->cursor-vector_get_mem(buffered_input->block_buffer,char);
}
GEM_INLINE bool buffered_input_file_eob(buffered_input_file_t* const buffered_input_file) {
  BUFFERED_INPUT_FILE_CHECK(buffered_input_file);
  return buffered_input_file_get_cursor_pos(buffered_input_file) >= vector_get_used(buffered_input_file->block_buffer);
}
GEM_INLINE error_code_t buffered_input_file_get_lines_block(
    buffered_input_file_t* const buffered_input,const uint64_t num_lines) {
  BUFFERED_INPUT_FILE_CHECK(buffered_input);
  input_file_t* const in_file = buffered_input->input_file;
  // Read lines
  if (input_file_eof(in_file)) return INPUT_STATUS_EOF;
  input_file_lock(in_file);
  if (input_file_eof(in_file)) {
    input_file_unlock(in_file);
    return INPUT_STATUS_EOF;
  }
  buffered_input->block_id = input_file_get_next_id(in_file);
  buffered_input->current_line_num = in_file->processed_lines+1;
  buffered_input->lines_in_buffer =
      input_file_get_lines(in_file,buffered_input->block_buffer,
          gem_expect_true(num_lines)?num_lines:BMI_NUM_LINES);
  input_file_unlock(in_file);
  // Setup the block
  buffered_input->cursor = vector_get_mem(buffered_input->block_buffer,char);
  return buffered_input->lines_in_buffer;
}
GEM_INLINE error_code_t buffered_input_file_add_lines_to_block(
    buffered_input_file_t* const buffered_input,const uint64_t num_lines) {
  BUFFERED_INPUT_FILE_CHECK(buffered_input);
  input_file_t* const input_file = buffered_input->input_file;
  // Read lines
  if (input_file_eof(input_file)) return INPUT_STATUS_EOF;
  const uint64_t current_position =
      buffered_input->cursor - vector_get_mem(buffered_input->block_buffer,char);
  const uint64_t lines_added =
      input_file_get_lines(input_file,buffered_input->block_buffer,
          gem_expect_true(num_lines)?num_lines:BMI_NUM_LINES);
  buffered_input->lines_in_buffer += lines_added;
  buffered_input->cursor = vector_get_elm(buffered_input->block_buffer,current_position,char);
  return lines_added;
}
/*
 * BufferedInputFile Utils
 */
GEM_INLINE void buffered_input_file_skip_line(buffered_input_file_t* const buffered_input) {
  BUFFERED_INPUT_FILE_CHECK(buffered_input);
  if (!buffered_input_file_eob(buffered_input)) {
    while (buffered_input->cursor[0]!=EOL) {
      ++buffered_input->cursor;
    }
    buffered_input->cursor[0]=EOS;
    ++buffered_input->cursor;
    ++buffered_input->current_line_num;
  }
}
GEM_INLINE error_code_t buffered_input_file_reload(
    buffered_input_file_t* const buffered_input,const uint64_t num_lines) {
  BUFFERED_INPUT_FILE_CHECK(buffered_input);
  // Dump buffer if BOF it attached to Map-input, and get new out block (always FIRST)
  buffered_input_file_dump_attached_buffers(buffered_input);
  // Read new input block
  if (gem_expect_false(buffered_input_file_get_lines_block(buffered_input,num_lines)==0)) return INPUT_STATUS_EOF;
  // Assign block ID
  buffered_input_file_set_id_attached_buffers(buffered_input,buffered_input->block_id);
  // Return OK
  return INPUT_STATUS_OK;
}
/*
 * Block Synchronization with Output
 */
GEM_INLINE void buffered_input_file_attach_buffered_output(
    buffered_input_file_t* const buffered_input_file,buffered_output_file_t* const buffered_output_file) {
  BUFFERED_INPUT_FILE_CHECK(buffered_input_file);
  BUFFERED_OUTPUT_FILE_CHECK(buffered_output_file);
  buffered_input_file->attached_buffered_output_file = buffered_output_file;
}
GEM_INLINE void buffered_input_file_dump_attached_buffers(buffered_input_file_t* const buffered_input_file) {
  BUFFERED_INPUT_FILE_CHECK(buffered_input_file);
  if (buffered_input_file->attached_buffered_output_file!=NULL) {
    buffered_output_file_dump(buffered_input_file->attached_buffered_output_file);
  }
}
GEM_INLINE void buffered_input_file_set_id_attached_buffers(buffered_input_file_t* const buffered_input_file,const uint64_t block_id) {
  BUFFERED_INPUT_FILE_CHECK(buffered_input_file);
  if (buffered_input_file->attached_buffered_output_file!=NULL) {
    buffered_output_file_set_block_ids(buffered_input_file->attached_buffered_output_file,block_id,0);
  }
}

