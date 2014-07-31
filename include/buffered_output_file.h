/*
 * PROJECT: GEMMapper
 * FILE: buffered_output_file.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef BUFFERED_OUTPUT_FILE_H_
#define BUFFERED_OUTPUT_FILE_H_

#include "essentials.h"
#include "output_buffer.h"
#include "output_file.h"

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
 * Checkers
 */
#define BUFFERED_OUTPUT_FILE_CHECK(buffered_output_file) \
  GEM_CHECK_NULL(buffered_output_file); \
  OUTPUT_FILE_CHECK((buffered_output_file)->output_file); \
  OUTPUT_BUFFER_CHECK((buffered_output_file)->buffer)

/*
 * Setup
 */
buffered_output_file_t* buffered_output_file_new(output_file_t* const output_file);
void buffered_output_file_close(buffered_output_file_t* const buffered_output);

/*
 * Accessors
 */
GEM_INLINE void buffered_output_file_get_block_ids(
    buffered_output_file_t* const buffered_output,uint32_t* const mayor_id,uint32_t* const minor_id);
GEM_INLINE void buffered_output_file_set_block_ids(
    buffered_output_file_t* const buffered_output,const uint32_t mayor_id,const uint32_t minor_id);
GEM_INLINE output_buffer_t* buffered_output_file_get_buffer(buffered_output_file_t* const buffered_output_file);
GEM_INLINE void buffered_output_file_set_buffer(
    buffered_output_file_t* const buffered_output,output_buffer_t* const out_buffer);

/*
 * Dump
 */
GEM_INLINE void buffered_output_file_dump(buffered_output_file_t* const buffered_output);
GEM_INLINE void buffered_output_file_safety_dump(buffered_output_file_t* const buffered_output);

/*
 * Buffered Output File Printers
 */
GEM_INLINE int vbofprintf(buffered_output_file_t* const buffered_output,const char *template,va_list v_args);
GEM_INLINE int bofprintf(buffered_output_file_t* const buffered_output,const char *template,...);
GEM_INLINE int vbofprintf_fixed(
    buffered_output_file_t* const buffered_output,
    const uint64_t expected_mem_usage,const char *template,va_list v_args);
GEM_INLINE int bofprintf_fixed(
    buffered_output_file_t* const buffered_output,
    const uint64_t expected_mem_usage,const char *template,...);

/*
 * Error Messages
 */
#define GEM_ERROR_BUFFER_DUMP "Output buffer. Could not perform dump"
#define GEM_ERROR_BUFFER_SAFETY_DUMP "Output buffer. Could not perform safety dump"


#endif /* BUFFERED_OUTPUT_FILE_H_ */
