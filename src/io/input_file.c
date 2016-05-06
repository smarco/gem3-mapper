/*
 * PROJECT: GEMMapper
 * FILE: input_file.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "io/input_file.h"

/*
 * Setup
 */
void input_file_initialize(input_file_t* const input_file,const uint64_t input_buffer_size) {
  // Input file
  input_file->mmaped = false;
  input_file->memory_manager = NULL;
  // Auxiliary Buffer (for synch purposes)
  input_file->buffer_allocated = input_buffer_size;
  input_file->file_buffer = mm_malloc(input_buffer_size);
  input_file->buffer_size = 0;
  input_file->buffer_begin = 0;
  input_file->buffer_pos = 0;
  input_file->global_pos = 0;
  input_file->processed_lines = 0;
}
input_file_t* input_stream_open(FILE* stream,const uint64_t input_buffer_size) {
  // Allocate handler
  input_file_t* const input_file = mm_alloc(input_file_t);
  // Create file manager
  input_file->file_manager = fm_open_FILE(stream,FM_READ);
  // Initialize input file
  input_file_initialize(input_file,input_buffer_size);
  return input_file;
}
input_file_t* input_gzip_stream_open(FILE* stream,const uint64_t input_buffer_size) {
  // Allocate handler
  input_file_t* const input_file = mm_alloc(input_file_t);
  // Create file manager
  input_file->file_manager = fm_open_gzFILE(stream,FM_READ);
  // Initialize input file
  input_file_initialize(input_file,input_buffer_size);
  return input_file;
}
input_file_t* input_bzip_stream_open(FILE* stream,const uint64_t input_buffer_size) {
  // Allocate handler
  input_file_t* const input_file = mm_alloc(input_file_t);
  // Create file manager
  input_file->file_manager = fm_open_bzFILE(stream,FM_READ);
  // Initialize input file
  input_file_initialize(input_file,input_buffer_size);
  return input_file;
}
input_file_t* input_file_open(
    char* const file_name,
    const uint64_t input_buffer_size,
    const bool mmap_file) {
  // Allocate handler
  input_file_t* const input_file = mm_alloc(input_file_t);
  // Prepare File
  if (!mmap_file) {
    input_file->file_manager = fm_open_file(file_name,FM_READ);
    input_file_initialize(input_file,input_buffer_size);
  } else {
    // Check it can be mapped
    struct stat stat_info;
    gem_stat(file_name,&stat_info);
    if (S_ISREG(stat_info.st_mode)) {
      // Regular file
      input_file->file_manager = NULL;
      input_file->mmaped = true;
      input_file->memory_manager = mm_bulk_mmap_file(file_name,MM_READ_ONLY,false);
      // Auxiliary Buffer (for synch purposes)
      input_file->file_buffer = mm_get_base_mem(input_file->memory_manager);
      input_file->buffer_size = mm_get_allocated(input_file->memory_manager);
      input_file->buffer_begin = 0;
      input_file->buffer_pos = 0;
      input_file->global_pos = 0;
      input_file->processed_lines = 0;
    } else {
      // Cannot be mapped (sorry)
      input_file->file_manager = fm_open_file(file_name,FM_READ);
      input_file_initialize(input_file,input_buffer_size);
    }
  }
  return input_file;
}
void input_file_rewind(input_file_t* const input_file) {
  if (input_file->mmaped) {
    input_file->buffer_size = mm_get_allocated(input_file->memory_manager);
  } else {
    fm_seek(input_file->file_manager,0);
    input_file->buffer_size = 0;
  }
  // Auxiliary Buffer
  input_file->buffer_begin = 0;
  input_file->buffer_pos = 0;
  input_file->global_pos = 0;
  input_file->processed_lines = 0;
}
void input_file_close(input_file_t* const input_file) {
  if (!input_file->mmaped) { // Regular file
    fm_close(input_file->file_manager);
    mm_free(input_file->file_buffer);
  } else { // MMaped file
    mm_bulk_free(input_file->memory_manager);
  }
  mm_free(input_file);
}
/*
 * Accessors
 */
uint8_t input_file_get_current_char(input_file_t* const input_file) {
  return input_file->file_buffer[input_file->buffer_pos];
}
uint8_t input_file_get_char_at(input_file_t* const input_file,const uint64_t position_in_buffer) {
  return input_file->file_buffer[input_file->buffer_pos];
}
char* input_file_get_nonull_file_name(input_file_t* const input_file) {
  return (input_file->file_manager!=NULL) ?
      fm_get_file_name(input_file->file_manager) : mm_get_mfile_name(input_file->memory_manager);
}
char* input_file_get_file_name(input_file_t* const input_file) {
  char* const file_name = input_file_get_nonull_file_name(input_file);
  return (file_name!=NULL) ? file_name : STREAM_FILE_NAME;
}
uint64_t input_file_get_size(input_file_t* const input_file) {
  return (input_file->mmaped) ?
    mm_get_allocated(input_file->memory_manager):
    fm_get_file_size(input_file->file_manager);
}
uint64_t input_file_get_current_line(input_file_t* const input_file) {
  return input_file->processed_lines;
}
bool input_file_eob(input_file_t* const input_file) {
  return input_file->buffer_pos >= input_file->buffer_size;
}
bool input_file_eof(input_file_t* const input_file) {
  return input_file_eob(input_file) && (input_file->mmaped || fm_eof(input_file->file_manager));
}
/*
 * Basic Buffer Functions
 */
uint64_t input_file_fill_buffer(input_file_t* const input_file) {
  input_file->global_pos += input_file->buffer_size;
  input_file->buffer_pos = 0;
  input_file->buffer_begin = 0;
  // Check EOF
  if (!fm_eof(input_file->file_manager)) {
    input_file->buffer_size = fm_read_mem(input_file->file_manager,input_file->file_buffer,input_file->buffer_allocated);
    fm_prefetch_next(input_file->file_manager,input_file->buffer_allocated);
    return input_file->buffer_size;
  } else {
    input_file->buffer_size = 0;
    return 0;
  }
}
uint64_t input_file_dump_to_buffer(input_file_t* const input_file,vector_t* const buffer_dst) {
  // Copy internal file buffer to buffer_dst
  const uint64_t chunk_size = input_file->buffer_pos-input_file->buffer_begin;
  if (gem_expect_false(chunk_size==0)) return 0;
  vector_reserve_additional(buffer_dst,chunk_size);
  memcpy(vector_get_mem(buffer_dst,uint8_t)+vector_get_used(buffer_dst),
      input_file->file_buffer+input_file->buffer_begin,chunk_size);
  vector_add_used(buffer_dst,chunk_size);
  // Update position
  input_file->buffer_begin=input_file->buffer_pos;
  // Return number of written bytes
  return chunk_size;
}
bool input_file_check_buffer(input_file_t* const input_file) {
  if (gem_expect_false(input_file_eob(input_file))) {
    // Returns if any character has been read
    return (input_file_fill_buffer(input_file) > 0);
  } else {
    return true;
  }
}
bool input_file_check_buffer__dump(input_file_t* const input_file,vector_t* const buffer_dst) {
  if (gem_expect_false(input_file_eob(input_file))) {
    if (gem_expect_true(buffer_dst!=NULL)) input_file_dump_to_buffer(input_file,buffer_dst);
    return (input_file_fill_buffer(input_file) > 0);
  } else {
    return true;
  }
}
bool input_file_next_char(input_file_t* const input_file) {
  ++input_file->buffer_pos;
  // Returns if there is any character in buffer
  return input_file_check_buffer(input_file);
}
bool input_file_next_char__dump(input_file_t* const input_file,vector_t* const buffer_dst) {
  ++input_file->buffer_pos;
  return input_file_check_buffer__dump(input_file,buffer_dst);
}
/*
 * Basic Line Functions
 */
void input_file_skip_eol(input_file_t* const input_file) {
  if (!input_file_eof(input_file)) {
    if (input_file_get_current_char(input_file)==DOS_EOL) {
      input_file_next_char(input_file); /* Skip DOS_EOL */
      if (gem_expect_false(!input_file_eof(input_file))) { // input_file_get_current_char(input_file)==EOL
        ++input_file->buffer_pos;
        if (gem_expect_false(input_file->buffer_pos >= input_file->buffer_size)) {
          input_file_fill_buffer(input_file);
        }
      }
    } else { // input_file_get_current_char(input_file)==EOL
      input_file_next_char(input_file); /* Skip EOF */
    }
  }
  ++(input_file->processed_lines);
}
void input_file_skip_eol__dump(input_file_t* const input_file,vector_t* const buffer_dst) {
  if (!input_file_eof(input_file)) {
    if (input_file_get_current_char(input_file)==DOS_EOL) {
      input_file_next_char__dump(input_file,buffer_dst); /* Skip DOS_EOL */
      if (gem_expect_false(!input_file_eof(input_file))) { // input_file_get_current_char(input_file)==EOL
        ++input_file->buffer_pos;
        if (gem_expect_true(buffer_dst!=NULL)) {
          input_file_dump_to_buffer(input_file,buffer_dst);
          vector_dec_used(buffer_dst);
          *vector_get_last_elm(buffer_dst,char)=EOL;
        }
        if (gem_expect_false(input_file->buffer_pos >= input_file->buffer_size)) {
          input_file_fill_buffer(input_file);
        }
      }
    } else { // input_file_get_current_char(input_file)==EOL
      input_file_next_char__dump(input_file,buffer_dst); /* Skip EOF*/
    }
  }
  ++(input_file->processed_lines);
}
uint64_t input_file_next_line(input_file_t* const input_file,vector_t* const buffer_dst) {
  input_file_check_buffer__dump(input_file,buffer_dst);
  if (input_file_eof(input_file)) return INPUT_STATUS_EOF;
  // Read line [OPT]
  if (gem_expect_false(input_file->mmaped)) { // MMAPED
    // Read line
    while (gem_expect_true(!input_file_eof(input_file) &&
        input_file_get_current_char(input_file)!=EOL &&
        input_file_get_current_char(input_file)!=DOS_EOL)) {
      input_file_next_char__dump(input_file,buffer_dst);
    }
  } else { // Not MMAPED
    uint64_t buffer_pos = input_file->buffer_pos;
    // Read line
    while (gem_expect_true((buffer_pos < input_file->buffer_size) || !fm_eof(input_file->file_manager) )) {
      const uint8_t current_char = input_file->file_buffer[buffer_pos];
      if (gem_expect_false(IS_ANY_EOL(current_char))) break;
      ++buffer_pos;
      if (buffer_pos >= input_file->buffer_size) {
        // Copy internal file buffer to buffer_dst
        const uint64_t chunk_size = buffer_pos-input_file->buffer_begin;
        if (gem_expect_false(chunk_size!=0)) {
          vector_reserve_additional(buffer_dst,chunk_size);
          memcpy(vector_get_mem(buffer_dst,uint8_t)+vector_get_used(buffer_dst),
              input_file->file_buffer+input_file->buffer_begin,chunk_size);
          vector_add_used(buffer_dst,chunk_size);
          // Update position
          input_file->buffer_begin=buffer_pos;
        }
        input_file->global_pos += input_file->buffer_size;
        buffer_pos = 0;
        input_file->buffer_begin = 0;
        // Check EOF
        if (!fm_eof(input_file->file_manager)) {
          input_file->buffer_size = fm_read_mem(input_file->file_manager,
              input_file->file_buffer,input_file->buffer_allocated);
          fm_prefetch_next(input_file->file_manager,input_file->buffer_allocated);
        } else {
          input_file->buffer_size = 0;
        }
      }
    }
    input_file->buffer_pos = buffer_pos;
  }
  // Handle EOL
  input_file_skip_eol__dump(input_file,buffer_dst);
  return INPUT_STATUS_OK;
}
/*
 * Line Readers (thread-unsafe, must call mutex functions before)
 */
uint64_t input_file_get_lines(
    input_file_t* const input_file,
    vector_t* buffer_dst,
    const uint64_t num_lines) {
  // Clear dst buffer
  vector_clear(buffer_dst);
  // Read lines
  uint64_t lines_read = 0;
  while (lines_read < num_lines && input_file_next_line(input_file,buffer_dst)) {
    ++lines_read;
  }
  // Dump remaining content into the buffer
  input_file_dump_to_buffer(input_file,buffer_dst);
  if (lines_read > 0 && *vector_get_last_elm(buffer_dst,char) != EOL) {
    vector_insert(buffer_dst,EOL,char);
  }
  return lines_read;
}

