/*
 * PROJECT: GEMMapper
 * FILE: input_file.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef INPUT_FILE_SLICED_H_
#define INPUT_FILE_SLICED_H_

#include "utils/essentials.h"
#include "io/input_buffer.h"

/*
 * Codes status
 */
#define INPUT_STATUS_OK 1
#define INPUT_STATUS_EOF 0
#define INPUT_STATUS_FAIL -1

/*
 * Input file
 */
typedef struct {
  /* File manager */
  fm_t* file_manager;
  /* Input-buffers */
  input_buffer_t** input_buffers;             // Buffers
  uint64_t total_input_buffers;               // Total number of buffers allocated
  uint64_t num_buffers_empty;                 // Total buffers empty
  /* Input-buffers order */
  uint64_t processed_buffer_id;               // Next buffer to be processed
  uint64_t current_buffer_id;                 // Current buffer being read
  input_buffer_t* current_buffer;             // Current input-buffer
  /* Input-buffers access */
  pthread_mutex_t input_buffers_mutex;        // Mutex
  pthread_cond_t input_buffers_ready_cond;    // CV
  pthread_cond_t input_buffers_empty_cond;    // CV
  /* Input-buffer request queue */
  uint64_t input_buffer_id;                   // ID generator
  uint64_t current_buffer_line;               // Current buffer line number
  uint64_t current_buffer_offset;             // Current buffer line number
  uint64_t current_input_line;                // Current input line number
  pthread_mutex_t input_read_lines_mutex;     // Mutex
} input_file_sliced_t;

/*
 * Setup
 */
input_file_sliced_t* input_stream_sliced_open(
    FILE* stream,
    const uint64_t input_num_blocks,
    const uint64_t block_size);
input_file_sliced_t* input_gzip_stream_sliced_open(
    FILE* stream,
    const uint64_t input_num_blocks,
    const uint64_t block_size);
input_file_sliced_t* input_bzip_stream_sliced_open(
    FILE* stream,
    const uint64_t input_num_blocks,
    const uint64_t block_size);
input_file_sliced_t* input_file_sliced_open(
    char* const file_name,
    const uint64_t input_num_blocks,
    const uint64_t block_size);
void input_file_sliced_close(
    input_file_sliced_t* const input_file_sliced);

/*
 * Accessors
 */
char* input_file_sliced_get_file_name(input_file_sliced_t* const input_file_sliced);
uint64_t input_file_sliced_get_next_id(input_file_sliced_t* const input_file_sliced);

/*
 * Process Input-buffers
 */
void input_file_sliced_return(
    input_file_sliced_t* const input_file_sliced,
    vector_t* const input_buffers);
void input_file_sliced_process(input_file_sliced_t* const input_file_sliced);
void input_file_sliced_discard_exhausted_buffers(
    input_file_sliced_t* const input_file_sliced,
    vector_t* const input_buffers);

/*
 * Read input-buffer lines
 */
bool input_file_sliced_eob(
    input_file_sliced_t* const input_file_sliced,
    input_buffer_t* const current_file_buffer);
bool input_file_sliced_read_lines(
    input_file_sliced_t* const input_file_sliced,
    input_buffer_t* const current_file_buffer,
    const uint64_t prefered_read_size,
    const uint64_t forced_read_lines,
    uint64_t* const total_read_lines,
    uint64_t* const total_read_size);

/*
 * Current Input-Buffer
 */
input_buffer_t* input_file_sliced_input_buffer_get_current(input_file_sliced_t* const input_file_sliced);
void input_file_sliced_input_buffer_next(input_file_sliced_t* const input_file_sliced);

#endif /* INPUT_FILE_SLICED_H_ */
