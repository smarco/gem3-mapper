/*
 * PROJECT: GEMMapper
 * FILE: buffered_input_file.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef BUFFERED_INPUT_FILE_H_
#define BUFFERED_INPUT_FILE_H_

#include "commons.h"
#include "input_file.h"
#include "buffered_output_file.h"

/*
 * Checker
 */
#define BUFFERED_INPUT_FILE_CHECK(buffered_input) GEM_CHECK_NULL((buffered_input)); \

typedef struct {
  /* Input file */
  input_file_t* input_file;        // Source input file
  /* Block buffer and cursors */
  input_buffer_t* input_buffer;    // Input Buffer
  uint64_t buffer_num_lines;       // Number of lines to reload the buffer
  /* Attached output buffer */
  buffered_output_file_t* attached_buffered_output_file;  // Attach output buffer file
} buffered_input_file_t;

/*
 * Buffered Input File Handlers
 */
buffered_input_file_t* buffered_input_file_new(input_file_t* const in_file,const uint64_t buffer_num_lines);
void buffered_input_file_close(buffered_input_file_t* const buffered_input);

/*
 * Accessors
 */
GEM_INLINE char** const buffered_input_file_get_text_line(buffered_input_file_t* const buffered_input);
GEM_INLINE uint64_t buffered_input_file_get_cursor_pos(buffered_input_file_t* const buffered_input);
GEM_INLINE uint64_t buffered_input_file_get_block_id(buffered_input_file_t* const buffered_input);
GEM_INLINE bool buffered_input_file_eob(buffered_input_file_t* const buffered_input);
GEM_INLINE void buffered_input_file_attach_buffered_output(
    buffered_input_file_t* const buffered_input_file,buffered_output_file_t* const buffered_output_file);

/*
 * Utils
 */
GEM_INLINE void buffered_input_file_skip_line(buffered_input_file_t* const buffered_input);
GEM_INLINE uint64_t buffered_input_file_reload(buffered_input_file_t* const buffered_input);
GEM_INLINE uint64_t buffered_input_file_reload__dump_attached(buffered_input_file_t* const buffered_input);

#endif /* BUFFERED_INPUT_FILE_H_ */
