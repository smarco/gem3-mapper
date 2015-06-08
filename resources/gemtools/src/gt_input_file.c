/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_file.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_input_file.h"

// File I/O constants
#define GT_INPUT_BUFFER_SIZE GT_BUFFER_SIZE_64M

/*
 * Basic I/O functions
 */
GT_INLINE void gt_input_file_initialize(gt_input_file* input_file) {
  GT_NULL_CHECK(input_file);
  // Input file
  input_file->mmaped = false;
  input_file->memory_manager = NULL;
  gt_cond_fatal_error__perror(pthread_mutex_init(&input_file->input_mutex,NULL),SYS_MUTEX_INIT);
  // Auxiliary Buffer (for synch purposes)
  input_file->file_buffer = gt_malloc(GT_INPUT_BUFFER_SIZE);
  input_file->buffer_size = 0;
  input_file->buffer_begin = 0;
  input_file->buffer_pos = 0;
  input_file->global_pos = 0;
  input_file->processed_lines = 0;
  // ID generator
  input_file->processed_id = 0;
  // Detect file format
  input_file->file_format = FILE_FORMAT_UNKNOWN;
  gt_input_file_detect_file_format(input_file);
}
gt_input_file* gt_input_stream_open(FILE* stream) {
  GT_NULL_CHECK(stream);
  // Allocate handler
  gt_input_file* input_file = gt_alloc(gt_input_file);
  // Create file manager
  input_file->file_manager = gt_fm_open_FILE(stream,GT_FM_READ);
  // Initialize input file
  gt_input_file_initialize(input_file);
  return input_file;
}
gt_input_file* gt_input_gzip_stream_open(FILE* stream) {
  GT_NULL_CHECK(stream);
  // Allocate handler
  gt_input_file* input_file = gt_alloc(gt_input_file);
  // Create file manager
  input_file->file_manager = gt_fm_open_gzFILE(stream,GT_FM_READ);
  // Initialize input file
  gt_input_file_initialize(input_file);
  return input_file;
}
gt_input_file* gt_input_bzip_stream_open(FILE* stream) {
  GT_NULL_CHECK(stream);
  // Allocate handler
  gt_input_file* input_file = gt_alloc(gt_input_file);
  // Create file manager
  input_file->file_manager = gt_fm_open_bzFILE(stream,GT_FM_READ);
  // Initialize input file
  gt_input_file_initialize(input_file);
  return input_file;
}
gt_input_file* gt_input_file_open(char* const file_name,const bool mmap_file) {
  GT_NULL_CHECK(file_name);
  // Allocate handler
  gt_input_file* input_file = gt_alloc(gt_input_file);
  // Prepare File
  if (!mmap_file) {
    input_file->file_manager = gt_fm_open_file(file_name,GT_FM_READ);
    gt_input_file_initialize(input_file);
  } else {
    // Check it can be mapped
    struct stat stat_info;
    gt_stat(file_name,&stat_info);
    if (S_ISREG(stat_info.st_mode)) {
      // Regular file
      input_file->file_manager = NULL;
      input_file->mmaped = true;
      input_file->memory_manager = gt_mm_bulk_mmap_file(file_name,GT_MM_READ_ONLY,false);
      GT_MUTEX_INIT(input_file->input_mutex);
      // Auxiliary Buffer (for synch purposes)
      input_file->file_buffer = gt_mm_get_base_mem(input_file->memory_manager);
      input_file->buffer_size = gt_mm_get_allocated(input_file->memory_manager);
      input_file->buffer_begin = 0;
      input_file->buffer_pos = 0;
      input_file->global_pos = 0;
      input_file->processed_lines = 0;
      // ID generator
      input_file->processed_id = 0;
      // Detect file format
      input_file->file_format = FILE_FORMAT_UNKNOWN;
      gt_input_file_detect_file_format(input_file);
    } else {
      // Cannot be mapped (sorry)
      input_file->file_manager = gt_fm_open_file(file_name,GT_FM_READ);
      gt_input_file_initialize(input_file);
    }
  }
  return input_file;
}
void gt_input_file_close(gt_input_file* const input_file) {
  GT_INPUT_FILE_CHECK(input_file);
  if (!input_file->mmaped) { // Regular file
    gt_fm_close(input_file->file_manager);
    gt_free(input_file->file_buffer);
  } else { // MMaped file
    gt_mm_free(input_file->memory_manager);
  }
}
/*
 * Accessors
 */
GT_INLINE uint8_t gt_input_file_get_current_char(gt_input_file* const input_file) {
  GT_INPUT_FILE_CHECK(input_file);
  return input_file->file_buffer[input_file->buffer_pos];
}
GT_INLINE uint8_t gt_input_file_get_char_at(gt_input_file* const input_file,const uint64_t position_in_buffer) {
  GT_INPUT_FILE_CHECK(input_file);
  return input_file->file_buffer[input_file->buffer_pos];
}
GT_INLINE uint64_t gt_input_file_get_next_id(gt_input_file* const input_file) {
  GT_INPUT_FILE_CHECK(input_file);
  const uint64_t id = input_file->processed_id;
  input_file->processed_id = (input_file->processed_id+1) % UINT32_MAX;
  return id;
}
GT_INLINE char* gt_input_file_get_file_name(gt_input_file* const input_file) {
  return (input_file->file_manager!=NULL) ?
      gt_fm_get_file_name(input_file->file_manager) : gt_mm_get_mfile_name(input_file->memory_manager);
}
GT_INLINE bool gt_input_file_eob(gt_input_file* const input_file) {
  return input_file->buffer_pos >= input_file->buffer_size;
}
GT_INLINE bool gt_input_file_eof(gt_input_file* const input_file) {
  return gt_input_file_eob(input_file) && (input_file->mmaped || gt_fm_eof(input_file->file_manager));
}
GT_INLINE void gt_input_file_lock(gt_input_file* const input_file) {
  GT_INPUT_FILE_CHECK(input_file);
  gt_cond_fatal_error(pthread_mutex_lock(&input_file->input_mutex),SYS_MUTEX);
}
GT_INLINE void gt_input_file_unlock(gt_input_file* const input_file) {
  GT_INPUT_FILE_CHECK(input_file);
  gt_cond_fatal_error(pthread_mutex_unlock(&input_file->input_mutex),SYS_MUTEX);
}
/*
 * Basic Buffer Functions
 */
GT_INLINE uint64_t gt_input_file_fill_buffer(gt_input_file* const input_file) {
  GT_INPUT_FILE_CHECK(input_file);
  input_file->global_pos += input_file->buffer_size;
  input_file->buffer_pos = 0;
  input_file->buffer_begin = 0;
  // Check EOF
  if (!gt_fm_eof(input_file->file_manager)) {
    input_file->buffer_size = gt_fm_read_mem(input_file->file_manager,input_file->file_buffer,GT_INPUT_BUFFER_SIZE);
    return input_file->buffer_size;
  } else {
    input_file->buffer_size = 0;
    return 0;
  }
}
GT_INLINE size_t gt_input_file_dump_to_buffer(gt_input_file* const input_file,gt_vector* const buffer_dst) {
  // FIXME: If mmap file, internal buffer is just pointers to mem (shortcut 2nd layer)
  GT_INPUT_FILE_CHECK(input_file);
  // Copy internal file buffer to buffer_dst
  const uint64_t chunk_size = input_file->buffer_pos-input_file->buffer_begin;
  if (gt_expect_false(chunk_size==0)) return 0;
  gt_vector_reserve_additional(buffer_dst,chunk_size);
  memcpy(gt_vector_get_mem(buffer_dst,uint8_t)+gt_vector_get_used(buffer_dst),
      input_file->file_buffer+input_file->buffer_begin,chunk_size);
  gt_vector_add_used(buffer_dst,chunk_size);
  // Update position
  input_file->buffer_begin=input_file->buffer_pos;
  // Return number of written bytes
  return chunk_size;
}
GT_INLINE void gt_input_file_check_buffer(gt_input_file* const input_file) {
  GT_INPUT_FILE_CHECK(input_file);
  if (gt_expect_false(gt_input_file_eob(input_file))) {
    gt_input_file_fill_buffer(input_file);
  }
}
GT_INLINE void gt_input_file_check_buffer__dump(gt_input_file* const input_file,gt_vector* const buffer_dst) {
  GT_INPUT_FILE_CHECK(input_file);
  if (gt_expect_false(gt_input_file_eob(input_file))) {
    if (gt_expect_true(buffer_dst!=NULL)) gt_input_file_dump_to_buffer(input_file,buffer_dst);
    gt_input_file_fill_buffer(input_file);
  }
}
GT_INLINE void gt_input_file_next_char(gt_input_file* const input_file) {
  GT_INPUT_FILE_CHECK(input_file);
  ++input_file->buffer_pos;
  gt_input_file_check_buffer(input_file);
}
GT_INLINE void gt_input_file_next_char__dump(gt_input_file* const input_file,gt_vector* const buffer_dst) {
  GT_INPUT_FILE_CHECK(input_file);
  ++input_file->buffer_pos;
  gt_input_file_check_buffer__dump(input_file,buffer_dst);
}
/*
 * Basic Line Functions
 */
GT_INLINE void gt_input_file_skip_eol(gt_input_file* const input_file,gt_vector* const buffer_dst) {
  GT_INPUT_FILE_CHECK(input_file);
  if (!gt_input_file_eof(input_file)) {
    if (gt_input_file_get_current_char(input_file)==DOS_EOL) {
      gt_input_file_next_char__dump(input_file,buffer_dst); /* Skip DOS_EOL */
      if (gt_expect_false(!gt_input_file_eof(input_file))) { // input_file_get_current_char(input_file)==EOL
        ++input_file->buffer_pos;
        if (gt_expect_true(buffer_dst!=NULL)) {
          gt_input_file_dump_to_buffer(input_file,buffer_dst);
          gt_vector_dec_used(buffer_dst);
          *gt_vector_get_last_elm(buffer_dst,char)=EOL;
        }
        if (gt_expect_false(gt_input_file_eob(input_file))) {
          gt_input_file_fill_buffer(input_file);
        }
      }
    } else { // gt_input_file_get_current_char(input_file)==EOL
      gt_input_file_next_char__dump(input_file,buffer_dst); /* Skip EOF*/
    }
  }
}
GT_INLINE uint64_t gt_input_file_next_line(gt_input_file* const input_file,gt_vector* const buffer_dst) {
  GT_INPUT_FILE_CHECK(input_file);
  GT_VECTOR_CHECK(buffer_dst);
  gt_input_file_check_buffer__dump(input_file,buffer_dst);
  if (gt_input_file_eof(input_file)) return GT_INPUT_STATUS_EOF;
  // Read line
  while (gt_expect_true(!gt_input_file_eof(input_file) &&
      gt_input_file_get_current_char(input_file)!=EOL &&
      gt_input_file_get_current_char(input_file)!=DOS_EOL)) {
    gt_input_file_next_char__dump(input_file,buffer_dst);
  }
  // Handle EOL
  gt_input_file_skip_eol(input_file,buffer_dst);
  return GT_INPUT_STATUS_OK;
}
/*
 * Line Readers (thread-unsafe, must call mutex functions before)
 */
GT_INLINE uint64_t gt_input_file_add_lines(
    gt_input_file* const input_file,gt_vector* buffer_dst,const uint64_t num_lines) {
  GT_INPUT_FILE_CHECK(input_file);
  GT_VECTOR_CHECK(buffer_dst);
  // Read lines
  uint64_t lines_read = 0;
  while (lines_read < num_lines && gt_input_file_next_line(input_file,buffer_dst)) {
    ++lines_read;
  }
  // Dump remaining content into the buffer
  gt_input_file_dump_to_buffer(input_file,buffer_dst);
  if (lines_read > 0 && *gt_vector_get_last_elm(buffer_dst,char) != EOL) {
    gt_vector_insert(buffer_dst,EOL,char);
  }
  input_file->processed_lines+=lines_read;
  return lines_read;
}
GT_INLINE uint64_t gt_input_file_get_lines(
    gt_input_file* const input_file,gt_vector* buffer_dst,const uint64_t num_lines) {
  GT_INPUT_FILE_CHECK(input_file);
  GT_VECTOR_CHECK(buffer_dst);
  // Clear dst buffer
  gt_vector_clear(buffer_dst);
  // Read lines
  return gt_input_file_add_lines(input_file,buffer_dst,num_lines);
}
GT_INLINE size_t gt_input_file_next_record(
    gt_input_file* const input_file,gt_vector* const buffer_dst,gt_string* const first_field,
    uint64_t* const num_blocks,uint64_t* const num_tabs) {
  GT_INPUT_FILE_CHECK(input_file);
  GT_VECTOR_CHECK(buffer_dst);
  gt_input_file_check_buffer__dump(input_file,buffer_dst);
  if (gt_input_file_eof(input_file)) return GT_INPUT_STATUS_EOF;
  // Read line
  uint64_t const begin_line_pos_at_file = input_file->buffer_pos;
  uint64_t const begin_line_pos_at_buffer = gt_vector_get_used(buffer_dst);
  uint64_t current_pfield = 0, length_first_field = 0;
  while (gt_expect_true(!gt_input_file_eof(input_file) &&
      gt_input_file_get_current_char(input_file)!=EOL &&
      gt_input_file_get_current_char(input_file)!=DOS_EOL)) {
    if (current_pfield==0) {
      if (gt_expect_false(gt_input_file_get_current_char(input_file)==TAB)) {
        ++current_pfield; ++(*num_tabs);
      } else {
        ++length_first_field;
      }
    } else if (current_pfield==1) {
      if (gt_expect_false(gt_input_file_get_current_char(input_file)==SPACE)) {
        ++(*num_blocks);
      } else if (gt_expect_false(gt_input_file_get_current_char(input_file)==TAB)) {
        ++current_pfield; ++(*num_tabs); ++(*num_blocks);
      }
    } else {
      if (gt_expect_false(gt_input_file_get_current_char(input_file)==TAB)) {
        ++current_pfield; ++(*num_tabs);
      }
    }
    gt_input_file_next_char__dump(input_file,buffer_dst);
  }
  // Handle EOL
  gt_input_file_skip_eol(input_file,buffer_dst);
  // Set first field (from the input_file_buffer or the buffer_dst)
  if (first_field) {
    char* first_field_begin;
    if (input_file->buffer_pos <= begin_line_pos_at_file) {
      gt_input_file_dump_to_buffer(input_file,buffer_dst); // Forced to dump to buffer
      first_field_begin = gt_vector_get_elm(buffer_dst,begin_line_pos_at_buffer,char);
    } else {
      first_field_begin = (char*)input_file->file_buffer + begin_line_pos_at_file;
    }
    gt_string_set_nstring(first_field,first_field_begin,length_first_field);
  }
  return GT_INPUT_STATUS_OK;
}
#define GT_INPUT_FILE_TEST_NEXT_CHAR(input_file,buffer_centinel) \
  ++buffer_centinel; \
  if (gt_expect_false(buffer_centinel >= input_file->buffer_size)) return true;
#define GT_INPUT_FILE_TEST_CURRENT_CHAR(input_file,buffer_centinel) input_file->file_buffer[buffer_centinel]
GT_INLINE bool gt_input_file_next_record_cmp_first_field(gt_input_file* const input_file,gt_string* const first_field) {
  GT_INPUT_FILE_CHECK(input_file);
  GT_STRING_CHECK(first_field);
  if (gt_expect_false(gt_input_file_eob(input_file) || gt_input_file_eof(input_file))) return true;
  // Read line
  char* const tag_begin = (char*)(input_file->file_buffer+input_file->buffer_pos);
  uint64_t buffer_centinel = input_file->buffer_pos;
  while (gt_expect_true(!gt_input_file_eof(input_file) &&
      GT_INPUT_FILE_TEST_CURRENT_CHAR(input_file,buffer_centinel)!=EOL &&
      GT_INPUT_FILE_TEST_CURRENT_CHAR(input_file,buffer_centinel)!=DOS_EOL)) {
    if (gt_expect_false(
        (GT_INPUT_FILE_TEST_CURRENT_CHAR(input_file,buffer_centinel)==SPACE ||
         GT_INPUT_FILE_TEST_CURRENT_CHAR(input_file,buffer_centinel)==TAB) )) {
      char* const tag_end = (char*)(input_file->file_buffer+buffer_centinel);
      uint64_t tag_lenth = tag_end-tag_begin;
      if (tag_lenth>2 && tag_begin[tag_lenth-2]==SLASH) tag_lenth-=2;
      if (first_field->length != tag_lenth) return false;
      return gt_strneq(first_field->buffer,tag_begin,tag_lenth);
    }
    GT_INPUT_FILE_TEST_NEXT_CHAR(input_file,buffer_centinel);
  }
  return true;
}
/*
 * Format detection (cascade of checkers)
 */
/* Forward declarations (gt_file_format_test_<FORMAT> in each logic module) */
GT_INLINE bool gt_input_file_test_fasta(
    gt_input_file* const input_file,gt_fasta_file_format* const fasta_file_format,const bool show_errors);
GT_INLINE bool gt_input_file_test_map(
    gt_input_file* const input_file,gt_map_file_format* const map_file_format,const bool show_errors);
GT_INLINE bool gt_input_file_test_sam(
    gt_input_file* const input_file,gt_sam_headers* const sam_headers,const bool show_errors);
/* */
gt_file_format gt_input_file_detect_file_format(gt_input_file* const input_file) {
  GT_INPUT_FILE_CHECK(input_file);
  if (input_file->file_format != FILE_FORMAT_UNKNOWN) return input_file->file_format;
  // Try to determine the file format
  gt_input_file_fill_buffer(input_file);
  // MAP test
  if (gt_input_file_test_map(input_file,&(input_file->map_type),false)) {
    input_file->file_format = MAP;
    return MAP;
  }
  // FASTA test
  if (gt_input_file_test_fasta(input_file,&(input_file->fasta_type),false)) {
    input_file->file_format = FASTA;
    return FASTA;
  }
  // SAM test
  if (gt_input_file_test_sam(input_file,&(input_file->sam_headers),false)) {
    input_file->file_format = SAM;
    return SAM;
  }
  // gt_error(FILE_FORMAT);
  return FILE_FORMAT_UNKNOWN;
}

