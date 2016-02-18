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
 * Input Buffer Scheme
 */
#define INPUT_FILE_READING_THREAD

/*
 * Codes status
 */
#define INPUT_STATUS_OK 1
#define INPUT_STATUS_EOF 0
#define INPUT_STATUS_FAIL -1

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
  /* File format */
  file_format_t file_format;
  union {
    fasta_file_format_t fasta;
    // map_file_format map_type;
    // sam_headers sam_headers;
  };
  /* Internal Buffer */
  uint64_t buffer_allocated;
  uint8_t* file_buffer;
  uint64_t buffer_size;
  uint64_t buffer_begin;
  uint64_t buffer_pos;
  uint64_t global_pos;
  uint64_t processed_lines;
  /* Input Buffer Queue */
  pthread_mutex_t input_mutex;           // Mutex
  pthread_cond_t requested_buffer_cond;  // CV (Reload input buffer)
  input_buffer_t* input_buffer;          // Input-Buffers
  /* ID generator */
  uint64_t processed_id;
} input_file_t;

/*
 * Basic I/O functions
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
uint64_t input_file_get_next_id(input_file_t* const input_file);

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
 * Line Readers (thread-unsafe)
 */
uint64_t input_file_add_lines(
    input_file_t* const input_file,
    vector_t* buffer_dst,
    const uint64_t num_lines);
uint64_t input_file_get_lines(
    input_file_t* const input_file,
    vector_t* buffer_dst,
    const uint64_t num_lines);

/*
 * Buffer reader (thread-safe)
 */
uint64_t input_file_reload_buffer(
    input_file_t* const input_file,
    input_buffer_t** const input_buffer,
    const uint64_t num_lines);

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
