/*
 * PROJECT: GEMMapper
 * FILE: mm_stack.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *     - StackMemory
 *         Memory allocator that serves memory in an stack fashion.
 *         All memory requested in successive calls to mm_stack_alloc() it's hold in memory.
 *         Once mm_stack_free() is called, all memory requests are freed (all at once).
 *         No individual memory free is possible with this type of memory.
 */
// TODO
// Implement DEBUG/SECURITY mechanism with which we can say if memory has been
// touched beyond limits (memsetting mem to sth and the checking after pop/push)

#include "system/mm_stack.h"

/*
 * Debug
 */
#define MM_STACK_LOG false
//#define MM_STACK_DEBUG

/*
 * Errors
 */
#define GEM_ERROR_MM_STACK_LARGE_MEM "Stack Allocator. Allocating large chunk of memory " \
                                     "(%"PRIu64" MBytes requested)[mm_slab_unit = %"PRIu64"MB]"

/*
 * Constants
 */
#define MM_STACK_INITIAL_SEGMENTS           10
#define MM_STACK_INITIAL_SEGMENTS_ALLOCATED  1
#define MM_STACK_INITIAL_STATES             10
#define MM_STACK_INITIAL_MALLOC_REQUESTS    10

/*
 * Segment handling
 */
void mm_stack_segment_allocate(mm_stack_t* const mm_stack,mm_stack_segment_t* const stack_segment) {
  stack_segment->slab_unit = mm_slab_request(mm_stack->mm_slab);
  stack_segment->memory = stack_segment->slab_unit->memory;
  stack_segment->memory_available = mm_stack->segment_size;
}
void mm_stack_segment_reset(mm_stack_t* const mm_stack,mm_stack_segment_t* const stack_segment) {
  stack_segment->memory = stack_segment->slab_unit->memory;
  stack_segment->memory_available = mm_stack->segment_size;
}
void mm_stack_segment_free(mm_stack_t* const mm_stack,mm_stack_segment_t* const stack_segment) {
  mm_slab_put(mm_stack->mm_slab,stack_segment->slab_unit);
}
/*
 * Setup
 */
mm_stack_t* mm_stack_new(mm_slab_t* const mm_slab) {
  static int no = 0;
  // Allocate handler
  mm_stack_t* const mm_stack = mm_alloc(mm_stack_t);
  mm_stack->id = no++;
  gem_cond_log(MM_STACK_LOG,"[GEM]> mm_stack(%"PRIu64").new()",mm_stack->id);
  // Initialize slab
  mm_stack->mm_slab = mm_slab;
  mm_stack->state = vector_new(MM_STACK_INITIAL_STATES,mm_stack_state_t); // Initialize stack state & dimensions
  mm_stack->segment_size = mm_slab_get_slab_size(mm_slab);
  mm_stack->current_segment = 0;
  // Init segments
  mm_stack->segments = vector_new(MM_STACK_INITIAL_SEGMENTS,mm_stack_segment_t);
  vector_set_used(mm_stack->segments,MM_STACK_INITIAL_SEGMENTS_ALLOCATED);
  VECTOR_ITERATE(mm_stack->segments,stack_segment,position,mm_stack_segment_t) {
    mm_stack_segment_allocate(mm_stack,stack_segment);
  }
  // Init malloc-requests
  mm_stack->malloc_requests = vector_new(MM_STACK_INITIAL_MALLOC_REQUESTS,void*);
  // Return
  return mm_stack;
}
void mm_stack_reap_segments(mm_stack_t* const mm_stack,const uint64_t resident_segments) {
  mm_slab_lock(mm_stack->mm_slab);
  VECTOR_ITERATE_OFFSET(mm_stack->segments,stack_segment,position,resident_segments,mm_stack_segment_t) {
    mm_stack_segment_free(mm_stack,stack_segment);
  }
  mm_slab_unlock(mm_stack->mm_slab);
  vector_set_used(mm_stack->segments,resident_segments); // Set segment available
}
void mm_stack_clear(mm_stack_t* const mm_stack) {
  // Clear all states
  vector_clear(mm_stack->state);
  // Clear first Segment
  mm_stack_segment_reset(mm_stack,vector_get_elm(mm_stack->segments,0,mm_stack_segment_t));
  mm_stack->current_segment = 0; // Set current segment
}
void mm_stack_delete(mm_stack_t* const mm_stack) {
  // Delete Slab
  mm_slab_lock(mm_stack->mm_slab);
  VECTOR_ITERATE(mm_stack->segments,stack_segment,p1,mm_stack_segment_t) {
    mm_stack_segment_free(mm_stack,stack_segment); // Return all slabs
  }
  mm_slab_unlock(mm_stack->mm_slab);
  vector_delete(mm_stack->state); // Free stack state
  // Delete Segments vector
  vector_delete(mm_stack->segments);
  // Delete malloc-requests
  VECTOR_ITERATE(mm_stack->malloc_requests,malloc_request,p2,void*) {
    mm_free(*malloc_request); // Free remaining requests
  }
  vector_delete(mm_stack->malloc_requests);
  // Free handler
  mm_free(mm_stack);
}
/*
 * State
 */
void mm_stack_push_state(mm_stack_t* const mm_stack) {
  // Allocate new state
  mm_stack_state_t* state;
  vector_alloc_new(mm_stack->state,mm_stack_state_t,state);
  // Setup state
  mm_stack_segment_t* const current_segment = vector_get_elm(
      mm_stack->segments,mm_stack->current_segment,mm_stack_segment_t);
  state->memory = current_segment->memory;
  state->memory_available = current_segment->memory_available;
  state->current_segment = mm_stack->current_segment;
  state->malloc_requests = vector_get_used(mm_stack->malloc_requests);
}
void mm_stack_pop_state(mm_stack_t* const mm_stack) {
  // Pop state
  mm_stack_state_t* const state = vector_get_last_elm(mm_stack->state,mm_stack_state_t);
  vector_dec_used(mm_stack->state);
  // Restore state
  mm_stack_segment_t* const current_segment = vector_get_elm(
      mm_stack->segments,state->current_segment,mm_stack_segment_t);
  current_segment->memory = state->memory; // Memory
  current_segment->memory_available = state->memory_available; // Memory available
  mm_stack->current_segment = state->current_segment; // Last segment
  // Restore malloc requests
  VECTOR_ITERATE_OFFSET(mm_stack->malloc_requests,malloc_request,n,state->malloc_requests,void*) {
    mm_free(*malloc_request); // Free requested
  }
  vector_set_used(mm_stack->malloc_requests,state->malloc_requests); // Set segment available
}
/*
 * Align stack memory
 */
void mm_stack_skip_align(mm_stack_t* const mm_stack,const uint64_t num_bytes) {
  GEM_CHECK_ZERO(num_bytes);
  if (gem_expect_true(num_bytes > 1)) {
    // Get last stack segment
    mm_stack_segment_t* const current_segment =
        vector_get_elm(mm_stack->segments,mm_stack->current_segment,mm_stack_segment_t);
    // Get current memory position
    void* const current_memory = current_segment->memory + (num_bytes-1);
    // Calculate the number of bytes to skip
    const uint64_t num_bytes_to_skip = (num_bytes-1) - (MM_CAST_ADDR(current_memory)%num_bytes);
    mm_stack_memory_allocate(mm_stack,num_bytes_to_skip,false);
  }
}
/*
 * Add new segment
 */
mm_stack_segment_t* mm_stack_add_segment(mm_stack_t* const mm_stack) {
  mm_stack_segment_t* stack_segment;
  // Increment last segment
  const uint64_t total_segments = vector_get_used(mm_stack->segments);
  ++(mm_stack->current_segment);
  if (mm_stack->current_segment < total_segments) {
    stack_segment = vector_get_elm(mm_stack->segments,mm_stack->current_segment,mm_stack_segment_t);
  } else {
    // Add new segment
    vector_reserve_additional(mm_stack->segments,1);
    vector_inc_used(mm_stack->segments);
    stack_segment = vector_get_last_elm(mm_stack->segments,mm_stack_segment_t);
    // Init segment
    mm_stack_segment_allocate(mm_stack,stack_segment);
    gem_cond_log(MM_STACK_LOG,"[GEM]> mm_stack(%"PRIu64").addSegment(%"PRIu64" x %"PRIu64" MB)",
        mm_stack->id,vector_get_used(mm_stack->segments),CONVERT_B_TO_MB(mm_stack->segment_size));
  }
  // Clear segment
  mm_stack_segment_reset(mm_stack,stack_segment);
  return stack_segment;
}
/*
 * Allocate/Deallocate
 */
#ifdef MM_STACK_DEBUG
void* mm_stack_memory_allocate(
    mm_stack_t* const mm_stack,
    const uint64_t num_bytes,
    const bool zero_mem) {
  // Issue malloc request
  if (num_bytes > BUFFER_SIZE_1G) {
    gem_warn(MM_STACK_LARGE_MEM,CONVERT_B_TO_MB(num_bytes),CONVERT_B_TO_MB(mm_stack->segment_size));
  }
  void** memory;
  vector_alloc_new(mm_stack->malloc_requests,void*,memory);
  *memory = mm_malloc_(1,num_bytes,zero_mem,0);
  // Return memory
  return *memory;
}
#else
void* mm_stack_memory_allocate(
    mm_stack_t* const mm_stack,
    const uint64_t num_bytes,
    const bool zero_mem) {
  // Get last stack segment
  mm_stack_segment_t* current_segment = vector_get_elm(mm_stack->segments,mm_stack->current_segment,mm_stack_segment_t);
  // Check if there is enough free memory in the segment
  if (gem_expect_false(num_bytes > current_segment->memory_available)) {
    // Check we can fit the request into a slab unit
    if (num_bytes <= mm_stack->segment_size) {
      // Use new segment
      current_segment = mm_stack_add_segment(mm_stack);
    } else {
      // Issue malloc request
      if (num_bytes > BUFFER_SIZE_1G) {
        gem_warn(MM_STACK_LARGE_MEM,CONVERT_B_TO_MB(num_bytes),CONVERT_B_TO_MB(mm_stack->segment_size));
      }
      void** memory;
      vector_alloc_new(mm_stack->malloc_requests,void*,memory);
      *memory = mm_malloc_(1,num_bytes,zero_mem,0);
      // Return memory
      return *memory;
    }
  }
  // Serve Memory Request
  void* const memory = current_segment->memory;
  current_segment->memory += num_bytes;
  current_segment->memory_available -= num_bytes;
  if (gem_expect_false(zero_mem)) memset(memory,0,num_bytes); // Set zero
  // Return memory
  return memory;
}
#endif

