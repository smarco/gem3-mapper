/*
 * PROJECT: GEMMapper
 * FILE: input_file.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef INPUT_FILE_H_
#define INPUT_FILE_H_

#include "utils/essentials.h"
#include "io/input_buffer.h"

/*
 * Codes status
 */
#define INPUT_STATUS_OK 1
#define INPUT_STATUS_EOF 0
#define INPUT_STATUS_FAIL -1

/*
 * File format
 */
typedef enum { FASTA, MAP, SAM, FILE_FORMAT_UNKNOWN } file_format_t;

/*
 * Input file
 */
typedef struct {
  /* File manager */
  fm_t* file_manager;
  bool mmaped;
  mm_t* memory_manager;
  /* Internal Buffer */
  uint64_t buffer_allocated;
  uint8_t* file_buffer;
  uint64_t buffer_size;
  uint64_t buffer_begin;
  uint64_t buffer_pos;
  uint64_t global_pos;
  uint64_t processed_lines;
} input_file_t;

/*
 * Setup
 */
input_file_t* input_stream_open(FILE* stream,const uint64_t input_buffer_size);
input_file_t* input_gzip_stream_open(FILE* stream,const uint64_t input_buffer_size);
input_file_t* input_bzip_stream_open(FILE* stream,const uint64_t input_buffer_size);
input_file_t* input_file_open(
    char* const file_name,
    const uint64_t input_buffer_size,
    const bool mmap_file);
void input_file_rewind(input_file_t* const input_file);
void input_file_close(input_file_t* const input_file);

/*
 * Accessors
 */
uint8_t input_file_get_current_char(input_file_t* const input_file);
uint8_t input_file_get_char_at(input_file_t* const input_file,const uint64_t position_in_buffer);

char* input_file_get_file_name(input_file_t* const input_file);
char* input_file_get_nonull_file_name(input_file_t* const input_file);
uint64_t input_file_get_size(input_file_t* const input_file);
uint64_t input_file_get_current_line(input_file_t* const input_file);

bool input_file_eof(input_file_t* const input_file);
void input_file_lock(input_file_t* const input_file);
void input_file_unlock(input_file_t* const input_file);

/*
 * Basic Buffer Functions
 */
uint64_t input_file_fill_buffer(input_file_t* const input_file);
uint64_t input_file_dump_to_buffer(input_file_t* const input_file,vector_t* const buffer_dst);
bool input_file_check_buffer(input_file_t* const input_file);
bool input_file_check_buffer__dump(input_file_t* const input_file,vector_t* const buffer_dst);
bool input_file_next_char(input_file_t* const input_file);
bool input_file_next_char__dump(input_file_t* const input_file,vector_t* const buffer_dst);

/*
 * Basic Line Functions
 */
void input_file_skip_eol(input_file_t* const input_file);
void input_file_skip_eol__dump(input_file_t* const input_file,vector_t* const buffer_dst);
uint64_t input_file_next_line(input_file_t* const input_file,vector_t* const buffer_dst);

/*
 * Line Readers (thread-unsafe, must call mutex functions before)
 */
uint64_t input_file_get_lines(
    input_file_t* const input_file,
    vector_t* buffer_dst,
    const uint64_t num_lines);

/*
 * Display
 */
#define PRI_input_file "s:%"PRIu64
#define PRI_input_file_content(input_file) input_file_get_file_name(input_file),input_file_get_current_line(input_file)

#endif /* INPUT_FILE_H_ */
