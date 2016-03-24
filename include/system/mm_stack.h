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

#include "system/commons.h"
#include "system/errors.h"
#include "system/fm.h"
#include "system/mm.h"
#include "system/mm_slab.h"
#include "utils/vector.h"


/*
 *  MM-Stack State
 */
typedef struct {
  /* Stacked Segments */
  uint64_t current_segment;    // Last segment being used
  uint64_t memory_available;   // Total memory available (last segment)
  void* memory;                // Pointer memory
  /* Malloc Requests */
  uint64_t malloc_requests;    // Total malloc requests performed
} mm_stack_state_t;

/*
 * MM-Stack Segments
 */
typedef struct {
  /* Memory */
  void* memory;                // Pointer to free memory
  uint64_t memory_available;   // Total memory available
  /* Slab-Unit */
  mm_slab_unit_t* slab_unit;   // Slab-Unit
} mm_stack_segment_t;

/*
 * MM-Stack (MemoryManager Stack-like)
 */
typedef struct {
  /* Stack state(s) */
  uint64_t id;
  vector_t* state;             // Vector of states (mm_stack_state_t)
  uint64_t segment_size;       // Total size of each memory segment
  /* Stack Segments */
  uint64_t current_segment;    // Last segment being used
  vector_t* segments;          // Memory segments (mm_stack_segment_t)
  /* Malloc requests */
  vector_t* malloc_requests;   // malloc() requests (Non stacked)
  /* Slab allocator */
  mm_slab_t* mm_slab;          // Memory allocator
} mm_stack_t;

/*
 * Setup
 */
mm_stack_t* mm_stack_new(mm_slab_t* const mm_slab);
void mm_stack_reap_segments(mm_stack_t* const mm_stack,const uint64_t resident_segments);
void mm_stack_clear(mm_stack_t* const mm_stack);
void mm_stack_delete(mm_stack_t* const mm_stack);

/*
 * State
 */
void mm_stack_push_state(mm_stack_t* const mm_stack);
void mm_stack_pop_state(mm_stack_t* const mm_stack);

/*
 * Align stack memory
 */
void mm_stack_skip_align(mm_stack_t* const mm_stack,const uint64_t num_bytes);

/*
 * Allocators
 */
void* mm_stack_memory_allocate(
    mm_stack_t* const mm_stack,
    const uint64_t num_bytes,
    const bool zero_mem);
/*
#define mm_stack_alloc(mm_stack,type) \
  ((type*)mm_stack_memory_allocate(mm_stack,sizeof(type),false))
#define mm_stack_malloc(mm_stack,num_bytes) \
  (fprintf(stderr,">MM-Alloc %lu from %s\n",(uint64_t)(num_bytes),__func__), \
  (       mm_stack_memory_allocate(mm_stack,(num_bytes),false)))
#define mm_stack_calloc(mm_stack,num_elements,type,clear_mem) \
  (fprintf(stderr,">MM-Alloc %lu from %s\n",(uint64_t)((num_elements)*sizeof(type)),__func__), \
  ((type*)mm_stack_memory_allocate(mm_stack,(num_elements)*sizeof(type),clear_mem)))
*/

#define mm_stack_alloc(mm_stack,type)                         ((type*)mm_stack_memory_allocate(mm_stack,sizeof(type),false))
#define mm_stack_malloc(mm_stack,num_bytes)                   (       mm_stack_memory_allocate(mm_stack,(num_bytes),false))
#define mm_stack_calloc(mm_stack,num_elements,type,clear_mem) ((type*)mm_stack_memory_allocate(mm_stack,(num_elements)*sizeof(type),clear_mem))


#define mm_stack_malloc_uint64(mm_stack) mm_stack_malloc(mm_stack,sizeof(uint64_t))
#define mm_stack_malloc_uint32(mm_stack) mm_stack_malloc(mm_stack,sizeof(uint32_t))
#define mm_stack_malloc_uint16(mm_stack) mm_stack_malloc(mm_stack,sizeof(uint16_t))
#define mm_stack_malloc_uint8(mm_stack)  mm_stack_malloc(mm_stack,sizeof(uint8_t))

#endif /* MM_STACK_H_ */
