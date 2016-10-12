/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   Output file data structure enables writing to a file
 *     Supports internal output-buffer queue for ordered output
 *     Supports deferred output-buffer dump
 *     Supports compressed (BZIP/GZIP) output
 */

#include "io/output_file.h"

/*
 * Setup
 */
void output_file_init_buffers(
    output_file_t* const output_file,
    const uint64_t max_output_buffers,
    const uint64_t buffer_size) {
  /* Output Buffers */
  uint64_t i;
  output_file->buffer_size = buffer_size;
  output_file->num_buffers = max_output_buffers;
  output_file->buffer = mm_calloc(max_output_buffers,output_buffer_t*,true);
  for (i=0;i<output_file->num_buffers;++i) {
    output_file->buffer[i]=NULL;
  }
  output_file->buffer_free=max_output_buffers;
  output_file->buffer_write_pending=0;
  /* Block ID prioritization (for synchronization purposes) */
  output_file->buffer_requests = pqueue_new(max_output_buffers);
  output_file->next_output_mayor_id=0;
  output_file->next_output_minor_id=0;
  /* Mutexes */
  CV_INIT(output_file->request_buffer_cond);
  MUTEX_INIT(output_file->output_file_mutex);
}
output_file_t* output_file_new(
    char* const file_name,
    const uint64_t max_output_buffers,
    const uint64_t buffer_size) {
  output_file_t* output_file = mm_alloc(output_file_t);
  // Output file
  output_file->file_manager = fm_open_file(file_name,FM_WRITE);
  // Setup buffers
  output_file_init_buffers(output_file,max_output_buffers,buffer_size);
  return output_file;
}
output_file_t* output_stream_new(
    FILE* const stream,
    const uint64_t max_output_buffers,
    const uint64_t buffer_size) {
  output_file_t* output_file = mm_alloc(output_file_t);
  // Output file
  output_file->file_manager = fm_open_FILE(stream,FM_WRITE);
  // Setup buffers
  output_file_init_buffers(output_file,max_output_buffers,buffer_size);
  return output_file;
}
output_file_t* output_gzip_stream_new(
    FILE* const stream,
    const uint64_t max_output_buffers,
    const uint64_t buffer_size) {
  output_file_t* output_file = mm_alloc(output_file_t);
  // Output file
  output_file->file_manager = fm_open_gzFILE(stream,FM_WRITE);
  // Setup buffers
  output_file_init_buffers(output_file,max_output_buffers,buffer_size);
  return output_file;
}
output_file_t* output_bzip_stream_new(
    FILE* const stream,
    const uint64_t max_output_buffers,
    const uint64_t buffer_size) {
  output_file_t* output_file = mm_alloc(output_file_t);
  // Output file
  output_file->file_manager = fm_open_bzFILE(stream,FM_WRITE);
  // Setup buffers
  output_file_init_buffers(output_file,max_output_buffers,buffer_size);
  return output_file;
}
void output_file_close(output_file_t* const output_file) {
  // Output file
  fm_close(output_file->file_manager);
  // Delete allocated buffers
  uint64_t i;
  for (i=0;i<output_file->num_buffers && output_file->buffer[i]!=NULL;++i) {
    output_buffer_delete(output_file->buffer[i]);
  }
  mm_free(output_file->buffer);
  // Prioritization
  pqueue_delete(output_file->buffer_requests);
  // Free mutex/CV
  CV_DESTROY(output_file->request_buffer_cond);
  MUTEX_DESTROY(output_file->output_file_mutex);
  // Free handler
  mm_free(output_file);
}
/*
 * Conditions
 */
bool output_file_serve_buffer_cond(
    output_file_t* const output_file,
    const uint64_t block_id) {
  // 0. We need a buffer free to be served [Security check to make it more robust]
  if (output_file->buffer_free==0) return false;
  // 1. Serve if request is the next in-order (next_request_id)
  const bool next_in_order = output_file->next_output_mayor_id==block_id;
  if (next_in_order) return true;
  // 2.1 Serve if we have at least one more buffer left to serve the next in-order (next_request_id)
  // 2.2 and the request is the next in the priority queue
  if (output_file->buffer_free>1 && pqueue_top_priority(output_file->buffer_requests)==block_id) {
    return true;
  } else {
    return false;
  }
}
bool output_file_eligible_request_cond(output_file_t* const output_file) {
  return !pqueue_is_empty(output_file->buffer_requests) && output_file->buffer_free>0;
}
bool output_file_write_victim_cond(
    output_file_t* const output_file,
    output_buffer_t* const output_buffer) {
  return (output_file->next_output_mayor_id==output_buffer->mayor_block_id &&
          output_file->next_output_minor_id==output_buffer->minor_block_id);
}
/*
 * Accessors
 */
output_buffer_t* output_file_get_free_buffer(output_file_t* const output_file) {
  // There is at least one free buffer. Get it!
  GEM_INTERNAL_CHECK(output_file->buffer_free > 0,"Output file. Cannot serve; no free buffers");
  uint64_t i;
  for (i=0;i<output_file->num_buffers && output_file->buffer[i]!=NULL;++i) {
    if (output_buffer_get_state(output_file->buffer[i])==OUTPUT_BUFFER_FREE) {
      return output_file->buffer[i];
    }
  }
  // Lazy allocation
  if (i>=output_file->num_buffers) {
    GEM_INTERNAL_CHECK(i<output_file->num_buffers,"Output file. Could not find free buffer");
  }
  GEM_INTERNAL_CHECK(i<output_file->num_buffers,"Output file. Could not find free buffer");
  if (output_file->buffer[i]==NULL) output_file->buffer[i] = output_buffer_new(output_file->buffer_size);
  // Return
  return output_file->buffer[i];
}
output_buffer_t* output_file_get_next_buffer_to_write(output_file_t* const output_file) {
  if (output_file->buffer_write_pending > 0) {
    const uint32_t mayor_block_id = output_file->next_output_mayor_id;
    const uint32_t minor_block_id = output_file->next_output_minor_id;
    uint64_t i;
    for (i=0;i<output_file->num_buffers && output_file->buffer[i]!=NULL;++i) {
      if (output_buffer_get_state(output_file->buffer[i]) == OUTPUT_BUFFER_WRITE_PENDING &&
          mayor_block_id==output_file->buffer[i]->mayor_block_id &&
          minor_block_id==output_file->buffer[i]->minor_block_id) {
        return output_file->buffer[i];
      }
    }
  }
  return NULL;
}
/*
 * Utils
 */
void output_file_print_buffers(FILE* stream,output_file_t* const output_file) {
  uint64_t i;
  for (i=0;i<output_file->num_buffers;++i) {
    if (output_file->buffer[i] == NULL) break;
    fprintf(stream,"[(%u,%u),%s,%s]",
        output_file->buffer[i]->mayor_block_id,
        output_file->buffer[i]->minor_block_id,
        output_file->buffer[i]->is_final_block ? "F":"I",
        output_file->buffer[i]->buffer_state==OUTPUT_BUFFER_FREE ? "Free" :
            (output_file->buffer[i]->buffer_state==OUTPUT_BUFFER_BUSY ? "Busy" : "PendW"));
  }
  fprintf(stream,"\n");
}
void output_file_release_buffer(
    output_file_t* const output_file,
    output_buffer_t* const output_buffer) {
  // Free buffer
  output_buffer_set_state(output_buffer,OUTPUT_BUFFER_FREE);
  output_buffer_clear(output_buffer);
  ++output_file->buffer_free; // Inc number of free buffers
}
output_buffer_t* output_file_request_buffer(
    output_file_t* const output_file,
    const uint64_t block_id) {
  PROF_INC_COUNTER(GP_OUTPUT_BUFFER_REQUESTS);
  output_buffer_t* output_buffer = NULL;
  MUTEX_BEGIN_SECTION(output_file->output_file_mutex)
  {
    // Add request to queue
    pqueue_push(output_file->buffer_requests,NULL,block_id);
    while (!output_file_serve_buffer_cond(output_file,block_id)) {
      PROF_INC_COUNTER(GP_OUTPUT_BUFFER_REQUESTS_STALLS);
      PROF_BLOCK() {
        if (output_file->buffer_free==0) {
          PROF_INC_COUNTER(GP_OUTPUT_BUFFER_REQUESTS_STALLS_BUSY);
        } else {
          PROF_INC_COUNTER(GP_OUTPUT_BUFFER_REQUESTS_STALLS_NOT_PRIORITY);
        }
      }
      CV_WAIT(output_file->request_buffer_cond,output_file->output_file_mutex);
    }
    // Get a free buffer & set state busy
    output_buffer = output_file_get_free_buffer(output_file);
    output_buffer_clear(output_buffer);
    output_buffer_set_state(output_buffer,OUTPUT_BUFFER_BUSY);
    output_buffer->mayor_block_id = block_id;
    output_buffer->minor_block_id = 0;
    --output_file->buffer_free; // Dec number of free buffers
    // Update next in-order (next_request_id) & priority-queue
    pqueue_pop_(output_file->buffer_requests);
    // Broadcast if any there are requests possible eligible
    if (output_file_eligible_request_cond(output_file)) {
      CV_BROADCAST(output_file->request_buffer_cond);
    }
  }
  MUTEX_END_SECTION(output_file->output_file_mutex);
  // Return buffer
  return output_buffer;
}
output_buffer_t* output_file_request_buffer_extension(
    output_file_t* const output_file,
    output_buffer_t* const output_buffer) {
  PROF_INC_COUNTER(GP_OUTPUT_BUFFER_EXTENSIONS);
  // Set current output-buffer as incomplete
  const uint64_t mayor_block_id = output_buffer->mayor_block_id;
  const uint64_t minor_block_id = output_buffer->minor_block_id;
  output_buffer_set_incomplete(output_buffer);
  // Output current output-buffer
  output_file_return_buffer(output_file,output_buffer);
  // Request a new one
  output_buffer_t* const output_buffer_extension =
      output_file_request_buffer(output_file,mayor_block_id);
  output_buffer_extension->minor_block_id = minor_block_id+1;
  return output_buffer_extension;
}
void output_file_next_block_id(
    output_file_t* const output_file,
    output_buffer_t* const output_buffer) {
  if (output_buffer->is_final_block) {
    ++output_file->next_output_mayor_id;
    output_file->next_output_minor_id = 0;
  } else {
    ++output_file->next_output_minor_id;
  }
}
void output_file_return_buffer(
    output_file_t* const output_file,
    output_buffer_t* output_buffer) {
  // Set the block buffer as write pending and set the victim
  MUTEX_BEGIN_SECTION(output_file->output_file_mutex)
  {
    // Set the buffer as write pending
    ++output_file->buffer_write_pending;
    output_buffer_set_state(output_buffer,OUTPUT_BUFFER_WRITE_PENDING);
    if (!output_file_write_victim_cond(output_file,output_buffer)) {
      // Enqueue the buffer and continue (someone else will do the writing)
      MUTEX_END_SECTION(output_file->output_file_mutex);
      return;
    }
  }
  MUTEX_END_SECTION(output_file->output_file_mutex);
  // I'm the victim (I will output as much as I can)
  bool keep_on_writing = true;
  do {
    // Write the current buffer
    if (output_buffer_get_used(output_buffer) > 0) {
      PROF_START_TIMER(GP_OUTPUT_WRITE_BUFFER);
      const uint64_t buffer_used = output_buffer_get_used(output_buffer);
      const char* const buffer = output_buffer_get_buffer(output_buffer);
      fm_write_mem(output_file->file_manager,buffer,buffer_used);
      PROF_ADD_COUNTER(GP_OUTPUT_BYTES_WRITTEN,buffer_used);
      PROF_STOP_TIMER(GP_OUTPUT_WRITE_BUFFER);
    }
    // Update buffers' state
    MUTEX_BEGIN_SECTION(output_file->output_file_mutex) {
      // Decrement write-pending blocks
      --output_file->buffer_write_pending;
      // Update next block ID (mayorID,minorID)
      output_file_next_block_id(output_file,output_buffer);
      // Release the current buffer
      output_file_release_buffer(output_file,output_buffer);
      // Search for the next output-buffer (in order)
      output_buffer_t* const next_output_buffer = output_file_get_next_buffer_to_write(output_file);
      if (next_output_buffer==NULL) {
        // Fine, I'm done, let's get out of here ASAP
        keep_on_writing = false;
        // Broadcast. Wake up sleepy.
        if (output_file_eligible_request_cond(output_file)) {
          CV_BROADCAST(output_file->request_buffer_cond);
        }
      } else {
        output_buffer = next_output_buffer;
      }
    } MUTEX_END_SECTION(output_file->output_file_mutex);
  } while (keep_on_writing);
}
/*
 * Output File Printers
 */
int vofprintf(output_file_t* const out_file,const char *template,va_list v_args) {
  int num_bytes;
  MUTEX_BEGIN_SECTION(out_file->output_file_mutex)
  {
    num_bytes = vfmprintf(out_file->file_manager,template,v_args);
  }
  MUTEX_END_SECTION(out_file->output_file_mutex);
  return num_bytes;
}
int ofprintf(output_file_t* const out_file,const char *template,...) {
  va_list v_args;
  va_start(v_args,template);
  const int num_bytes = vofprintf(out_file,template,v_args);
  va_end(v_args);
  return num_bytes;
}

