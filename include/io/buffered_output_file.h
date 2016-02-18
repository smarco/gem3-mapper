/*
 * PROJECT: GEMMapper
 * FILE: buffered_output_file.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef BUFFERED_OUTPUT_FILE_H_
#define BUFFERED_OUTPUT_FILE_H_

#include "utils/essentials.h"
#include "io/output_buffer.h"
#include "io/output_file.h"

typedef struct {
  /* Output file */
  output_file_t* output_file;
  /* Output Buffer */
  output_buffer_t* buffer;
} buffered_output_file_t;

// Codes status
#define BUFFERED_OUTPUT_FILE_OK 0
#define BUFFERED_OUTPUT_FILE_FAIL -1

/*
 * Setup
 */
buffered_output_file_t* buffered_output_file_new(output_file_t* const output_file);
void buffered_output_file_close(buffered_output_file_t* const buffered_output);

/*
 * Utils
 */
void buffered_output_file_request_buffer(
    buffered_output_file_t* const buffered_output,
    const uint32_t block_id);
void buffered_output_file_dump_buffer(buffered_output_file_t* const buffered_output);
void buffered_output_file_safety_dump_buffer(buffered_output_file_t* const buffered_output);
void buffered_output_file_reserve(
    buffered_output_file_t* const buffered_output,
    const uint64_t num_chars);

/*
 * Fast-printer
 */
void bofprintf_uint64(buffered_output_file_t* const buffered_output,const uint64_t number);
void bofprintf_int64(buffered_output_file_t* const buffered_output,const int64_t number);
void bofprintf_char(buffered_output_file_t* const buffered_output,const char character);
void bofprintf_string(
    buffered_output_file_t* const buffered_output,
    const int string_length,
    const char* const string);
#define bofprintf_string_literal(buffered_output,string) bofprintf_string(buffered_output,sizeof(string)-1,string)

#endif /* BUFFERED_OUTPUT_FILE_H_ */
