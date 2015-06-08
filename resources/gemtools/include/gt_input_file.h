/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_file.h
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_INPUT_FILE_H_
#define GT_INPUT_FILE_H_

#include "gt_essentials.h"
#include "gt_attributes.h"
#include "gt_sam_attributes.h"

/*
 * Codes status
 */
#define GT_INPUT_STATUS_OK 1
#define GT_INPUT_STATUS_EOF 0
#define GT_INPUT_STATUS_FAIL -1

/*
 * Checkers
 */
#define GT_INPUT_FILE_CHECK(input_file) \
  GT_NULL_CHECK(input_file); \
  if (input_file->file_manager) {GT_FM_CHECK(input_file->file_manager);} \
  if (input_file->memory_manager) {GT_MM_CHECK(input_file->memory_manager);}

/*
 * GT Input file
 */
typedef enum { FASTA, MAP, SAM, FILE_FORMAT_UNKNOWN } gt_file_format;
typedef enum { STREAM, REGULAR_FILE, MAPPED_FILE, GZIPPED_FILE, BZIPPED_FILE } gt_file_type;
typedef struct {
  /* File manager */
  gt_fm* file_manager;
  bool mmaped;
  gt_mm* memory_manager;
  pthread_mutex_t input_mutex;
  /* File format */
  gt_file_format file_format;
  union {
    gt_map_file_format map_type;
    gt_fasta_file_format fasta_type;
    gt_sam_headers sam_headers;
  };
  /* Auxiliary Buffer (for synch purposes) */
  uint8_t* file_buffer;
  uint64_t buffer_size;
  uint64_t buffer_begin;
  uint64_t buffer_pos;
  uint64_t global_pos;
  uint64_t processed_lines;
  /* ID generator */
  uint64_t processed_id;
} gt_input_file;

/*
 * Basic I/O functions
 */
gt_input_file* gt_input_stream_open(FILE* stream);
gt_input_file* gt_input_gzip_stream_open(FILE* stream);
gt_input_file* gt_input_bzip_stream_open(FILE* stream);
gt_input_file* gt_input_file_open(char* const file_name,const bool mmap_file);
void gt_input_file_close(gt_input_file* const input_file);

/* Format detection */
gt_file_format gt_input_file_detect_file_format(gt_input_file* const input_file);

/*
 * Accessors
 */
GT_INLINE uint8_t gt_input_file_get_current_char(gt_input_file* const input_file);
GT_INLINE uint8_t gt_input_file_get_char_at(gt_input_file* const input_file,const uint64_t position_in_buffer);
GT_INLINE uint64_t gt_input_file_get_next_id(gt_input_file* const input_file);
GT_INLINE char* gt_input_file_get_file_name(gt_input_file* const input_file);

GT_INLINE bool gt_input_file_eof(gt_input_file* const input_file);
GT_INLINE void gt_input_file_lock(gt_input_file* const input_file);
GT_INLINE void gt_input_file_unlock(gt_input_file* const input_file);

/*
 * Basic Buffer Functions
 */
GT_INLINE size_t gt_input_file_fill_buffer(gt_input_file* const input_file);
GT_INLINE size_t gt_input_file_dump_to_buffer(gt_input_file* const input_file,gt_vector* const buffer_dst);
GT_INLINE void gt_input_file_check_buffer(gt_input_file* const input_file);
GT_INLINE void gt_input_file_check_buffer__dump(gt_input_file* const input_file,gt_vector* const buffer_dst);
GT_INLINE void gt_input_file_next_char(gt_input_file* const input_file);
GT_INLINE void gt_input_file_next_char__dump(gt_input_file* const input_file,gt_vector* const buffer_dst);

/*
 * Basic Line Functions
 */
GT_INLINE void gt_input_file_skip_eol(gt_input_file* const input_file,gt_vector* const buffer_dst);
GT_INLINE uint64_t gt_input_file_next_line(gt_input_file* const input_file,gt_vector* const buffer_dst);
GT_INLINE size_t gt_input_file_next_record(
    gt_input_file* const input_file,gt_vector* const buffer_dst,gt_string* const first_field,
    uint64_t* const num_spaces,uint64_t* const num_tabs);
GT_INLINE bool gt_input_file_next_record_cmp_first_field(gt_input_file* const input_file,gt_string* const first_field);

/*
 * Line Readers (thread-unsafe, must call mutex functions before)
 */
GT_INLINE uint64_t gt_input_file_add_lines(
    gt_input_file* const input_file,gt_vector* buffer_dst,const uint64_t num_lines);
GT_INLINE uint64_t gt_input_file_get_lines(
    gt_input_file* const input_file,gt_vector* buffer_dst,const uint64_t num_lines);
GT_INLINE size_t gt_input_file_next_record(
    gt_input_file* const input_file,gt_vector* const buffer_dst,gt_string* const first_field,
    uint64_t* const num_blocks,uint64_t* const num_tabs);

#endif /* GT_INPUT_FILE_H_ */
