/*
 * PROJECT: GEMMapper
 * FILE: buffered_input_file.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "buffered_input_file.h"

/*
 * Constants
 */
// Profile Level
#define PROFILE_LEVEL PLOW

/*
 * Buffered map file handlers
 */
buffered_input_file_t* buffered_input_file_new(input_file_t* const in_file,const uint64_t buffer_num_lines) {
  INPUT_FILE_CHECK(in_file);
  buffered_input_file_t* buffered_input = mm_alloc(buffered_input_file_t);
  /* Input file */
  buffered_input->input_file = in_file;
  /* Block buffer and cursors */
  buffered_input->buffer_num_lines = buffer_num_lines;
  buffered_input->input_buffer = input_buffer_new();
  /* Attached output buffer */
  buffered_input->attached_buffered_output_file = NULL;
  return buffered_input;
}
void buffered_input_file_close(buffered_input_file_t* const buffered_input) {
  BUFFERED_INPUT_FILE_CHECK(buffered_input);
  input_buffer_delete(buffered_input->input_buffer);
  mm_free(buffered_input);
}
/*
 * Accessors
 */
char** buffered_input_file_get_text_line(buffered_input_file_t* const buffered_input) {
  BUFFERED_INPUT_FILE_CHECK(buffered_input);
  return input_buffer_get_cursor(buffered_input->input_buffer);
}
uint64_t buffered_input_file_get_cursor_pos(buffered_input_file_t* const buffered_input) {
  BUFFERED_INPUT_FILE_CHECK(buffered_input);
  return input_buffer_get_cursor_pos(buffered_input->input_buffer);
}
uint64_t buffered_input_file_get_block_id(buffered_input_file_t* const buffered_input) {
  BUFFERED_INPUT_FILE_CHECK(buffered_input);
  return buffered_input->input_buffer->block_id;
}
bool buffered_input_file_eob(buffered_input_file_t* const buffered_input) {
  BUFFERED_INPUT_FILE_CHECK(buffered_input);
  return input_buffer_eob(buffered_input->input_buffer);
}
void buffered_input_file_attach_buffered_output(
    buffered_input_file_t* const buffered_input_file,buffered_output_file_t* const buffered_output_file) {
  BUFFERED_INPUT_FILE_CHECK(buffered_input_file);
  BUFFERED_OUTPUT_FILE_CHECK(buffered_output_file);
  buffered_input_file->attached_buffered_output_file = buffered_output_file;
}
/*
 * Utils
 */
uint64_t buffered_input_file_reload(buffered_input_file_t* const buffered_input) {
  BUFFERED_INPUT_FILE_CHECK(buffered_input);
  // Reload Buffer
  return input_file_reload_buffer(buffered_input->input_file,
      &(buffered_input->input_buffer),buffered_input->buffer_num_lines);
}
//#include "libittnotify.h"
//__itt_resume();
//__itt_pause();
uint64_t buffered_input_file_reload__dump_attached(buffered_input_file_t* const buffered_input) {
  BUFFERED_INPUT_FILE_CHECK(buffered_input);
  PROFILE_START(GP_BUFFERED_INPUT_RELOAD__DUMP_ATTACHED,PROFILE_LEVEL);
  // Dump buffer
  if (buffered_input->attached_buffered_output_file!=NULL) {
    buffered_output_file_dump_buffer(buffered_input->attached_buffered_output_file);
  }
  // Read new input block
  PROFILE_START(GP_BUFFERED_INPUT_RELOAD,PROFILE_LEVEL);
  if (gem_expect_false(buffered_input_file_reload(buffered_input)==0)) {
    PROFILE_STOP(GP_BUFFERED_INPUT_RELOAD,PROFILE_LEVEL);
    PROFILE_STOP(GP_BUFFERED_INPUT_RELOAD__DUMP_ATTACHED,PROFILE_LEVEL);
    return INPUT_STATUS_EOF;
  }
  PROFILE_STOP(GP_BUFFERED_INPUT_RELOAD,PROFILE_LEVEL);
  // Get output buffer (block ID)
  if (buffered_input->attached_buffered_output_file!=NULL) {
    buffered_output_file_request_buffer(
        buffered_input->attached_buffered_output_file,buffered_input->input_buffer->block_id);
  }
  PROFILE_STOP(GP_BUFFERED_INPUT_RELOAD__DUMP_ATTACHED,PROFILE_LEVEL);
  return INPUT_STATUS_OK; // Return OK
}

