/*
 * PROJECT: GEMMapper
 * FILE: mm_stack.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *     - StackMemory
 *         Memory allocator that serves memory in an stack fashion.
 *         All memory requested in successive calls to mm_stack_alloc() it's hold in memory.
 *         Once mm_stack_free() is called, all memory requests are freed (all at once).
 *         No individual memory free is possible with this type of memory.
 */

#ifndef MM_STACK_H_
#define MM_STACK_H_

#include "commons.h"
#include "errors.h"
#include "vector.h"
#include "fm.h"
#include "mm.h"
#include "mm_slab.h"

/*
 * Checkers
 */
#define MM_STACK_CHECK(mm_stack) GEM_CHECK_NULL(mm_stack)

/*
 * MemoryManager Stack-like
 */
typedef struct {
  mm_slab_unit_t* slab_unit;   // Slab-Unit
  void* memory;                // Pointer to free memory
  uint64_t memory_available;   // Total memory available
} mm_stack_segment_t;
typedef struct {
  uint64_t segment_pos;            // Last segments being used
  uint64_t memory_available;   // Total memory available (last segment)
} mm_stack_state_t;
typedef struct {
  /* Stack state(s) */
  vector_t* state;             // Vector of states (mm_stack_state_t)
  uint64_t segment_size;       // Total size of each memory segment
  uint64_t current_segment;    // Last segment being used
  /* Vector segments */
  vector_t* segments;          // Memory segments (mm_stack_segment_t)
  /* Slab allocator */
  mm_slab_t* mm_slab;          // Memory allocator
} mm_stack_t;

/*
 * Setup
 */
GEM_INLINE mm_stack_t* mm_stack_new(mm_slab_t* const mm_slab);
GEM_INLINE void mm_stack_delete(mm_stack_t* const mm_stack);

/*
 * State
 */
GEM_INLINE void mm_stack_push_state(mm_stack_t* const mm_stack);
GEM_INLINE void mm_stack_pop_state(mm_stack_t* const mm_stack,const bool reap_segments);

/*
 * Allocators
 */
GEM_INLINE void* mm_stack_memory_allocate(mm_stack_t* const mm_stack,const uint64_t num_bytes,const bool zero_mem);
GEM_INLINE void mm_stack_free(mm_stack_t* const mm_stack);

#define mm_stack_alloc(mm_stack,type)                         ((type*)mm_stack_memory_allocate(mm_stack,sizeof(type),false))
#define mm_stack_malloc(mm_stack,num_bytes)                   (       mm_stack_memory_allocate(mm_stack,num_bytes,false))
#define mm_stack_calloc(mm_stack,num_elements,type,clear_mem) ((type*)mm_stack_memory_allocate(mm_stack,num_elements*sizeof(type),clear_mem))

#define mm_stack_malloc_uint64(mm_stack) mm_stack_malloc(mm_stack,sizeof(uint64_t))
#define mm_stack_malloc_uint32(mm_stack) mm_stack_malloc(mm_stack,sizeof(uint32_t))
#define mm_stack_malloc_uint16(mm_stack) mm_stack_malloc(mm_stack,sizeof(uint16_t))
#define mm_stack_malloc_uint8(mm_stack)  mm_stack_malloc(mm_stack,sizeof(uint8_t))

/*
 * Errors
 */
#define GEM_ERROR_MM_STACK_CANNOT_ALLOCATE "Stack Allocator. Could not allocate memory (%"PRIu64" Bytes requested)[mm_slab_unit = %"PRIu64"B]"

#endif /* MM_STACK_H_ */
