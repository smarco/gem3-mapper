/*
 * PROJECT: GEMMapper
 * FILE: output_file.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "output_file.h"

/*
 * Setup
 */
GEM_INLINE void output_file_init_buffers(output_file_t* const out_file) {
  /* Output Buffers */
  uint64_t i;
  for (i=0;i<MAX_OUTPUT_BUFFERS;++i) {
    out_file->buffer[i]=NULL;
  }
  out_file->buffer_busy=0;
  out_file->buffer_write_pending=0;
  /* Block ID (for synchronization purposes) */
  out_file->mayor_block_id=0;
  out_file->minor_block_id=0;
  /* Mutexes */
  CV_INIT(out_file->out_buffer_cond);
  CV_INIT(out_file->out_write_cond);
  MUTEX_INIT(out_file->out_file_mutex);
}
output_file_t* output_file_new(char* const file_name,const output_file_type output_file_type) {
  GEM_CHECK_NULL(file_name);
  output_file_t* out_file = mm_alloc(output_file_t);
  // Output file
  out_file->file_manager = fm_open_file(file_name,FM_WRITE);
  out_file->file_type = output_file_type;
  // Setup buffers
  output_file_init_buffers(out_file);
  return out_file;
}
output_file_t* output_stream_new(FILE* const stream,const output_file_type output_file_type) {
  GEM_CHECK_NULL(stream);
  output_file_t* out_file = mm_alloc(output_file_t);
  // Output file
  out_file->file_manager = fm_open_FILE(stream,FM_WRITE);
  out_file->file_type = output_file_type;
  // Setup buffers
  output_file_init_buffers(out_file);
  return out_file;
}
output_file_t* output_gzip_stream_new(FILE* const stream,const output_file_type output_file_type) {
  GEM_CHECK_NULL(stream);
  output_file_t* out_file = mm_alloc(output_file_t);
  // Output file
  out_file->file_manager = fm_open_gzFILE(stream,FM_WRITE);
  out_file->file_type = output_file_type;
  // Setup buffers
  output_file_init_buffers(out_file);
  return out_file;
}
output_file_t* output_bzip_stream_new(FILE* const stream,const output_file_type output_file_type) {
  GEM_CHECK_NULL(stream);
  output_file_t* out_file = mm_alloc(output_file_t);
  // Output file
  out_file->file_manager = fm_open_bzFILE(stream,FM_WRITE);
  out_file->file_type = output_file_type;
  // Setup buffers
  output_file_init_buffers(out_file);
  return out_file;
}
void output_file_close(output_file_t* const out_file) {
  OUTPUT_FILE_CONSISTENCY_CHECK(out_file);
  // Output file
  fm_close(out_file->file_manager);
  // Delete allocated buffers
  uint64_t i;
  for (i=0;i<MAX_OUTPUT_BUFFERS && out_file->buffer[i]!=NULL;++i) {
    output_buffer_delete(out_file->buffer[i]);
  }
  // Free mutex/CV
  CV_DESTROY(out_file->out_buffer_cond);
  CV_DESTROY(out_file->out_write_cond);
  MUTEX_DESTROY(out_file->out_file_mutex);
  // Free handler
  mm_free(out_file);
}
/*
 * Output File Printers
 */
GEM_INLINE int vofprintf(output_file_t* const out_file,const char *template,va_list v_args) {
  OUTPUT_FILE_CHECK(out_file);
  GEM_CHECK_NULL(template);
  int num_bytes;
  MUTEX_BEGIN_SECTION(out_file->out_file_mutex)
  {
    num_bytes = vfmprintf(out_file->file_manager,template,v_args);
  }
  MUTEX_END_SECTION(out_file->out_file_mutex);
  return num_bytes;
}
GEM_INLINE int ofprintf(output_file_t* const out_file,const char *template,...) {
  OUTPUT_FILE_CHECK(out_file);
  GEM_CHECK_NULL(template);
  va_list v_args;
  va_start(v_args,template);
  const int num_bytes = vofprintf(out_file,template,v_args);
  va_end(v_args);
  return num_bytes;
}
/*
 * Internal Buffers Accessors
 */
GEM_INLINE output_buffer_t* output_file_get_buffer(output_file_t* const out_file) {
  OUTPUT_FILE_CONSISTENCY_CHECK(out_file);
  // Conditional guard. Wait till there is any free buffer left
  while (out_file->buffer_busy==MAX_OUTPUT_BUFFERS) {
    CV_WAIT(out_file->out_buffer_cond,out_file->out_file_mutex);
  }
  // There is at least one free buffer. Get it!
  uint64_t i;
  for (i=0;i<MAX_OUTPUT_BUFFERS&&out_file->buffer[i]!=NULL;++i) {
    if (output_buffer_get_state(out_file->buffer[i])==OUTPUT_BUFFER_FREE) break;
  }
  if (out_file->buffer[i]==NULL) {
    out_file->buffer[i] = output_buffer_new();
  }
  ++out_file->buffer_busy;
  output_buffer_initiallize(out_file->buffer[i],OUTPUT_BUFFER_BUSY);
  return out_file->buffer[i];
}
GEM_INLINE output_buffer_t* output_file_request_buffer(output_file_t* const out_file) {
  OUTPUT_FILE_CONSISTENCY_CHECK(out_file);
  output_buffer_t* fresh_buffer;
  MUTEX_BEGIN_SECTION(out_file->out_file_mutex)
  {
    fresh_buffer = output_file_get_buffer(out_file);
  }
  MUTEX_END_SECTION(out_file->out_file_mutex);
  return fresh_buffer;
}
GEM_INLINE void output_file_release_buffer(
    output_file_t* const out_file,output_buffer_t* const output_buffer) {
  OUTPUT_FILE_CONSISTENCY_CHECK(out_file);
  OUTPUT_BUFFER_CHECK(output_buffer);
  // Broadcast. Wake up sleepy.
  if (out_file->buffer_busy==MAX_OUTPUT_BUFFERS) {
    CV_BROADCAST(out_file->out_buffer_cond);
  }
  // Free buffer
  output_buffer_set_state(output_buffer,OUTPUT_BUFFER_FREE);
  --out_file->buffer_busy;
  output_buffer_clear(output_buffer);
}
GEM_INLINE void output_file_return_buffer(
    output_file_t* const out_file,output_buffer_t* const output_buffer) {
  OUTPUT_FILE_CONSISTENCY_CHECK(out_file);
  OUTPUT_BUFFER_CHECK(output_buffer);
  MUTEX_BEGIN_SECTION(out_file->out_file_mutex)
  {
    output_file_release_buffer(out_file,output_buffer);
  }
  MUTEX_END_SECTION(out_file->out_file_mutex);
}
GEM_INLINE output_buffer_t* output_file_write_buffer(
    output_file_t* const out_file,output_buffer_t* const output_buffer) {
  OUTPUT_FILE_CONSISTENCY_CHECK(out_file);
  if (output_buffer_get_used(output_buffer) > 0) {
    vector_t* const vbuffer = output_buffer_to_vchar(output_buffer);
    MUTEX_BEGIN_SECTION(out_file->out_file_mutex)
    {
      fm_write_mem(out_file->file_manager,vector_get_mem(vbuffer,char),vector_get_used(vbuffer));
    }
    MUTEX_END_SECTION(out_file->out_file_mutex);
  }
  output_buffer_initiallize(output_buffer,OUTPUT_BUFFER_BUSY);
  return output_buffer;
}
GEM_INLINE output_buffer_t* output_file_sorted_write_buffer_asynchronous(
    output_file_t* const out_file,output_buffer_t* output_buffer,const bool asynchronous) {
  OUTPUT_FILE_CONSISTENCY_CHECK(out_file);
  // Set the block buffer as write pending and set the victim
  bool victim;
  uint32_t mayor_block_id, minor_block_id;
  MUTEX_BEGIN_SECTION(out_file->out_file_mutex)
  {
    victim = (out_file->mayor_block_id==output_buffer_get_mayor_block_id(output_buffer) &&
              out_file->minor_block_id==output_buffer_get_minor_block_id(output_buffer));
    while (!asynchronous && !victim) {
      CV_WAIT(out_file->out_write_cond,out_file->out_file_mutex);
      victim = (out_file->mayor_block_id==output_buffer_get_mayor_block_id(output_buffer) &&
                out_file->minor_block_id==output_buffer_get_minor_block_id(output_buffer));
    }
    // Set the buffer as write pending
    ++out_file->buffer_write_pending;
    output_buffer_set_state(output_buffer,OUTPUT_BUFFER_WRITE_PENDING);
    if (!victim) { // Enqueue the buffer and continue (someone else will do the writing)
      output_buffer = output_file_get_buffer(out_file);
      MUTEX_END_SECTION(out_file->out_file_mutex);
      return output_buffer;
    } else {
      mayor_block_id = out_file->mayor_block_id;
      minor_block_id = out_file->minor_block_id;
    }
  }
  MUTEX_END_SECTION(out_file->out_file_mutex);
  // I'm the victim, I will output as much as I can
  do {
    // Write the current buffer
    if (output_buffer_get_used(output_buffer) > 0) {
      vector_t* const vbuffer = output_buffer_to_vchar(output_buffer);
      fm_write_mem(out_file->file_manager,vector_get_mem(vbuffer,char),vector_get_used(vbuffer));
    }
    // Update buffers' state
    MUTEX_BEGIN_SECTION(out_file->out_file_mutex) {
      // Decrement write-pending blocks and update next block ID (mayorID,minorID)
      --out_file->buffer_write_pending;
      if (output_buffer->is_final_block) {
        ++mayor_block_id;
        minor_block_id = 0;
      } else {
        ++minor_block_id;
      }
      // Search for the next block buffer in order
      bool next_found = false;
      if (out_file->buffer_write_pending>0) {
        uint64_t i;
        for (i=0;i<MAX_OUTPUT_BUFFERS&&out_file->buffer[i]!=NULL;++i) {
          if (mayor_block_id==output_buffer_get_mayor_block_id(out_file->buffer[i]) &&
              minor_block_id==output_buffer_get_minor_block_id(out_file->buffer[i])) {
            if (output_buffer_get_state(out_file->buffer[i])!=OUTPUT_BUFFER_WRITE_PENDING) {
              break; // Cannot dump a busy buffer
            } else {
              // I'm still the victim, free the current buffer and output the new one
              output_file_release_buffer(out_file,output_buffer);
              next_found = true;
              output_buffer = out_file->buffer[i];
              break;
            }
          }
        }
      }
      if (!next_found) {
        // Fine, I'm done, let's get out of here ASAP
        out_file->mayor_block_id = mayor_block_id;
        out_file->minor_block_id = minor_block_id;
        output_buffer_initiallize(output_buffer,OUTPUT_BUFFER_BUSY);
        CV_BROADCAST(out_file->out_write_cond);
        MUTEX_END_SECTION(out_file->out_file_mutex);
        return output_buffer;
      }
    } MUTEX_END_SECTION(out_file->out_file_mutex);
  } while (true);
}
GEM_INLINE output_buffer_t* output_file_dump_buffer(
    output_file_t* const out_file,output_buffer_t* const output_buffer,const bool asynchronous) {
  OUTPUT_FILE_CONSISTENCY_CHECK(out_file);
  switch (out_file->file_type) {
    case SORTED_FILE:
      return output_file_sorted_write_buffer_asynchronous(out_file,output_buffer,asynchronous);
      break;
    case UNSORTED_FILE:
      return output_file_write_buffer(out_file,output_buffer);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
