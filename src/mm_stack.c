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

#include "mm_stack.h"

#define MM_STACK_INITIAL_SEGMENTS 10
#define MM_STACK_INITIAL_SEGMENTS_ALLOCATED 1
#define MM_STACK_INITIAL_STATES 10

/*
 * Stack Segment
 */
GEM_INLINE void mm_stack_segment_allocate(mm_stack_t* const mm_stack,mm_stack_segment_t* const stack_segment) {
#ifdef MM_STACK_DEBUG
    stack_segment->slab_unit = (mm_slab_unit_t*) malloc(mm_stack->segment_size);
    stack_segment->memory = (void*) stack_segment->slab_unit;
#else
    stack_segment->slab_unit = mm_slab_request(mm_stack->mm_slab);
    stack_segment->memory = stack_segment->slab_unit->memory;
#endif
    stack_segment->memory_available = mm_stack->segment_size;
}
GEM_INLINE void mm_stack_segment_reset(mm_stack_t* const mm_stack,mm_stack_segment_t* const stack_segment) {
#ifdef MM_STACK_DEBUG
  stack_segment->memory = (void*) stack_segment->slab_unit;
  stack_segment->memory_available = mm_stack->segment_size;
#else
  stack_segment->memory = stack_segment->slab_unit->memory;
  stack_segment->memory_available = mm_stack->segment_size;
#endif
}
GEM_INLINE void mm_stack_segment_free(mm_stack_t* const mm_stack,mm_stack_segment_t* const stack_segment) {
#ifdef MM_STACK_DEBUG
    free((void*)stack_segment->slab_unit);
#else
    mm_slab_put(mm_stack->mm_slab,stack_segment->slab_unit);
#endif
}
/*
 * Setup
 */
GEM_INLINE mm_stack_t* mm_stack_new(mm_slab_t* const mm_slab) {
  MM_SLAB_CHECK(mm_slab);
  // Allocate handler
  mm_stack_t* const mm_stack = mm_alloc(mm_stack_t);
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
  // Return
  return mm_stack;
}
GEM_INLINE void mm_stack_delete(mm_stack_t* const mm_stack) {
  MM_STACK_CHECK(mm_stack);
  // Return all slabs
  mm_slab_lock(mm_stack->mm_slab);
  VECTOR_ITERATE(mm_stack->segments,stack_segment,position,mm_stack_segment_t) {
    mm_stack_segment_free(mm_stack,stack_segment);
  }
  mm_slab_unlock(mm_stack->mm_slab);
  // Free segments' vector
  vector_delete(mm_stack->segments);
  // Free stack state
  vector_delete(mm_stack->state);
  // Free handler
  mm_free(mm_stack);
}
/*
 * State
 */
GEM_INLINE void mm_stack_push_state(mm_stack_t* const mm_stack) {
  // Allocate new state
  mm_stack_state_t* state;
  vector_alloc_new(mm_stack->state,mm_stack_state_t,state);
  // Setup state
  mm_stack_segment_t* const current_segment = vector_get_elm(mm_stack->segments,mm_stack->current_segment,mm_stack_segment_t);
  state->current_segment = mm_stack->current_segment;
  state->memory = current_segment->memory;
  state->memory_available = current_segment->memory_available;
}
GEM_INLINE void mm_stack_pop_state(mm_stack_t* const mm_stack,const bool reap_segments) {
  // Pop state
  mm_stack_state_t* const state = vector_get_last_elm(mm_stack->state,mm_stack_state_t);
  vector_dec_used(mm_stack->state);
  // Restore state
  mm_stack->current_segment = state->current_segment; // Last segment
  mm_stack_segment_t* const current_segment = vector_get_elm(mm_stack->segments,state->current_segment,mm_stack_segment_t);
  current_segment->memory = state->memory; // Memory
  current_segment->memory_available = state->memory_available; // Memory available
  // Reap non-resident segments
  if (gem_expect_false(reap_segments)) {
    mm_slab_lock(mm_stack->mm_slab);
    VECTOR_ITERATE_OFFSET(mm_stack->segments,stack_segment,position,state->current_segment+1,mm_stack_segment_t) {
      mm_stack_segment_free(mm_stack,stack_segment);
    }
    mm_slab_unlock(mm_stack->mm_slab);
    vector_set_used(mm_stack->segments,state->current_segment+1); // Set segment available
  }
}
/*
 * Allocators
 */
GEM_INLINE mm_stack_segment_t* mm_stack_add_segment(mm_stack_t* const mm_stack) {
  mm_stack_segment_t* stack_segment;
  // Increment last segment
  ++(mm_stack->current_segment);
  if (mm_stack->current_segment < vector_get_used(mm_stack->segments)) {
    stack_segment = vector_get_elm(mm_stack->segments,mm_stack->current_segment,mm_stack_segment_t);
  } else {
    // Add new segment
    vector_reserve_additional(mm_stack->segments,1);
    vector_inc_used(mm_stack->segments);
    stack_segment = vector_get_last_elm(mm_stack->segments,mm_stack_segment_t);
    // Init segment
    mm_stack_segment_allocate(mm_stack,stack_segment);
  }
  // Clear segment
  mm_stack_segment_reset(mm_stack,stack_segment);
  return stack_segment;
}
GEM_INLINE void* mm_stack_memory_allocate(mm_stack_t* const mm_stack,const uint64_t num_bytes,const bool zero_mem) {
  MM_STACK_CHECK(mm_stack);
  // Get last stack segment
  mm_stack_segment_t* current_segment = vector_get_elm(mm_stack->segments,mm_stack->current_segment,mm_stack_segment_t);
  // Check if there is enough free memory in the segment
  if (gem_expect_false(num_bytes > current_segment->memory_available)) {
    // Check we can fit the request into a slab unit
    gem_cond_fatal_error(num_bytes > mm_stack->segment_size,
        MM_STACK_CANNOT_ALLOCATE,num_bytes,mm_stack->segment_size);
    current_segment = mm_stack_add_segment(mm_stack); // Use new segment
  }
  // Serve Memory Request
  current_segment->memory_available -= num_bytes;
  void* const memory = current_segment->memory;
  current_segment->memory += num_bytes;
  // Check zero mem
  if (gem_expect_false(zero_mem)) memset(memory,0,num_bytes);
  // Return memory
  return memory;
}
GEM_INLINE void mm_stack_free(mm_stack_t* const mm_stack) {
  MM_STACK_CHECK(mm_stack);
  // Clear all states
  vector_clear(mm_stack->state);
  // Clear first Segment
  mm_stack_segment_reset(mm_stack,vector_get_elm(mm_stack->segments,0,mm_stack_segment_t));
  // Reap non-resident segments
  if (vector_get_used(mm_stack->segments) > 1) {
    mm_slab_lock(mm_stack->mm_slab);
    VECTOR_ITERATE_OFFSET(mm_stack->segments,stack_segment,position,1,mm_stack_segment_t) {
      mm_stack_segment_free(mm_stack,stack_segment);
    }
    mm_slab_unlock(mm_stack->mm_slab);
  }
  vector_set_used(mm_stack->segments,1); // Set used to 1
  mm_stack->current_segment = 0; // Set current segment
}

