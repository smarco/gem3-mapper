/*
 * PROJECT: GEMMapper
 * FILE: input_file.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef INPUT_FILE_H_
#define INPUT_FILE_H_

#include "essentials.h"

/*
 * Codes status
 */
#define INPUT_STATUS_OK 1
#define INPUT_STATUS_EOF 0
#define INPUT_STATUS_FAIL -1

/*
 * File I/O constants
 */
#define INPUT_BUFFER_SIZE BUFFER_SIZE_64M

/*
 * Checkers
 */
#define INPUT_FILE_CHECK(input_file) \
  GEM_CHECK_NULL(input_file); \
  if (input_file->file_manager) {FM_CHECK(input_file->file_manager);} \
  if (input_file->memory_manager) {MM_CHECK(input_file->memory_manager);}

/*
 * FASTQ/FASTA/MULTIFASTA File Attribute
 */
typedef enum { F_FASTA, F_FASTQ } fasta_file_format_t;
typedef struct {
  fasta_file_format_t format;
} fasta_file_attributes_t;

/*
 * Input file
 */
typedef enum { FASTA, MAP, SAM, FILE_FORMAT_UNKNOWN } file_format_t;
typedef struct {
  /* File manager */
  fm_t* file_manager;
  bool mmaped;
  mm_t* memory_manager;
  pthread_mutex_t input_mutex;
  /* File format */
  file_format_t file_format;
  union {
    fasta_file_format_t fasta;
    // map_file_format map_type;
    // sam_headers sam_headers;
  };
  /* Buffer (for synch purposes) */
  uint8_t* file_buffer;
  uint64_t buffer_size;
  uint64_t buffer_begin;
  uint64_t buffer_pos;
  uint64_t global_pos;
  uint64_t processed_lines;
  /* ID generator */
  uint64_t processed_id;
} input_file_t;

/*
 * Basic I/O functions
 */
input_file_t* input_stream_open(FILE* stream);
input_file_t* input_gzip_stream_open(FILE* stream);
input_file_t* input_bzip_stream_open(FILE* stream);
input_file_t* input_file_open(char* const file_name,const bool mmap_file);
void input_file_rewind(input_file_t* const input_file);
void input_file_close(input_file_t* const input_file);

/* Format detection */
file_format_t input_file_detect_file_format(input_file_t* const input_file);

/*
 * Accessors
 */
GEM_INLINE uint8_t input_file_get_current_char(input_file_t* const input_file);
GEM_INLINE uint8_t input_file_get_char_at(input_file_t* const input_file,const uint64_t position_in_buffer);
GEM_INLINE uint64_t input_file_get_next_id(input_file_t* const input_file);

GEM_INLINE char* input_file_get_file_name(input_file_t* const input_file);
GEM_INLINE char* input_file_get_nonull_file_name(input_file_t* const input_file);
GEM_INLINE uint64_t input_file_get_size(input_file_t* const input_file);
GEM_INLINE uint64_t input_file_get_current_line(input_file_t* const input_file);

GEM_INLINE bool input_file_eof(input_file_t* const input_file);
GEM_INLINE void input_file_lock(input_file_t* const input_file);
GEM_INLINE void input_file_unlock(input_file_t* const input_file);

/*
 * Basic Buffer Functions
 */
GEM_INLINE uint64_t input_file_fill_buffer(input_file_t* const input_file);
GEM_INLINE uint64_t input_file_dump_to_buffer(input_file_t* const input_file,vector_t* const buffer_dst);
GEM_INLINE bool input_file_check_buffer(input_file_t* const input_file);
GEM_INLINE bool input_file_check_buffer__dump(input_file_t* const input_file,vector_t* const buffer_dst);
GEM_INLINE bool input_file_next_char(input_file_t* const input_file);
GEM_INLINE bool input_file_next_char__dump(input_file_t* const input_file,vector_t* const buffer_dst);

/*
 * Basic Line Functions
 */
GEM_INLINE void input_file_skip_eol(input_file_t* const input_file);
GEM_INLINE void input_file_skip_eol__dump(input_file_t* const input_file,vector_t* const buffer_dst);
GEM_INLINE uint64_t input_file_next_line(input_file_t* const input_file,vector_t* const buffer_dst);

/*
 * Line Readers (thread-unsafe, must call mutex functions before)
 */
GEM_INLINE uint64_t input_file_add_lines(
    input_file_t* const input_file,vector_t* buffer_dst,const uint64_t num_lines);
GEM_INLINE uint64_t input_file_get_line(input_file_t* const input_file,vector_t* buffer_dst);
GEM_INLINE uint64_t input_file_get_lines(
    input_file_t* const input_file,vector_t* buffer_dst,const uint64_t num_lines);

/*
 * Printers
 */
#define PRI_input_file "s:%"PRIu64
#define PRI_input_file_content(input_file) input_file_get_file_name(input_file),input_file_get_current_line(input_file)

/*
 * Error Messages
 */
#define GEM_ERROR_FILE_FORMAT "Could not determine file format of '%s'"

#endif /* INPUT_FILE_H_ */
