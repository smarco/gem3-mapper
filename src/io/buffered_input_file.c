/*
 * PROJECT: GEMMapper
 * FILE: buffered_input_file.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "io/buffered_input_file.h"
#include "mapper/mapper_profile.h"

/*
 * Errors
 */
#define GEM_ERROR_BUFFERED_INPUT_BUFFER_LINES "Reading Input File ('%s'). Premature end of file (expected end of record)"

/*
 * Constants
 */
// Profile Level
#define PROFILE_LEVEL PLOW

/*
 * Buffered map file handlers
 */
buffered_input_file_t* buffered_input_file_new(input_file_t* const in_file,const uint64_t reload_buffer_size) {
  buffered_input_file_t* buffered_input = mm_alloc(buffered_input_file_t);
  /* Input file */
  buffered_input->input_file = in_file;
  /* Block buffer and cursors */
  buffered_input->reload_buffer_size = reload_buffer_size;
  buffered_input->input_buffer = input_buffer_new();
  /* Attached output buffer */
  buffered_input->attached_buffered_output_file = NULL;
  return buffered_input;
}
void buffered_input_file_close(buffered_input_file_t* const buffered_input) {
  input_buffer_delete(buffered_input->input_buffer);
  mm_free(buffered_input);
}
/*
 * Accessors
 */
char** buffered_input_file_get_text_line(buffered_input_file_t* const buffered_input) {
  return input_buffer_get_cursor(buffered_input->input_buffer);
}
uint64_t buffered_input_file_get_cursor_pos(buffered_input_file_t* const buffered_input) {
  return input_buffer_get_cursor_pos(buffered_input->input_buffer);
}
uint64_t buffered_input_file_get_block_id(buffered_input_file_t* const buffered_input) {
  return buffered_input->input_buffer->block_id;
}
bool buffered_input_file_eob(buffered_input_file_t* const buffered_input) {
  return input_buffer_eob(buffered_input->input_buffer);
}
void buffered_input_file_attach_buffered_output(
    buffered_input_file_t* const buffered_input_file,
    buffered_output_file_t* const buffered_output_file) {
  buffered_input_file->attached_buffered_output_file = buffered_output_file;
}
/*
 * Utils
 */
uint64_t buffered_input_file_reload(
    buffered_input_file_t* const buffered_input,
    const uint64_t min_lines) {
  // Reload Buffer
  const uint64_t lines_in_buffer = input_file_reload_buffer(
      buffered_input->input_file,&(buffered_input->input_buffer),
      min_lines,buffered_input->reload_buffer_size);
  gem_cond_fatal_error(lines_in_buffer%8!=0,
      BUFFERED_INPUT_BUFFER_LINES,input_file_get_file_name(buffered_input->input_file));
  return lines_in_buffer;
}
//#include "libittnotify.h"
//__itt_resume();
//__itt_pause();
uint64_t buffered_input_file_reload__dump_attached(
    buffered_input_file_t* const buffered_input,
    const uint64_t min_lines) {
  PROFILE_START(GP_BUFFERED_INPUT_RELOAD__DUMP_ATTACHED,PROFILE_LEVEL);
  // Dump buffer
  if (buffered_input->attached_buffered_output_file!=NULL) {
    buffered_output_file_dump_buffer(buffered_input->attached_buffered_output_file);
  }
  // Read new input block
  PROFILE_START(GP_BUFFERED_INPUT_RELOAD,PROFILE_LEVEL);
  if (gem_expect_false(buffered_input_file_reload(buffered_input,min_lines)==0)) {
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

