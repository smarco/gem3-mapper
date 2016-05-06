/*
 * PROJECT: GEMMapper
 * FILE: segmented_vector.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Memory segmented vector for generic type elements.
 *              Based on slabs for an efficient memory management
 */

#include "utils/segmented_vector.h"

#define SVECTOR_INITIAL_SEGMENTS 10
#define SVECTOR_INITIAL_SEGMENTS_ALLOCATED 1

/*
 * Vector Setup (Initialization & Allocation)
 */
svector_t* svector_new_(mm_slab_t* const mm_slab,const uint64_t element_size) {
  svector_t* const svector = mm_alloc(svector_t);
  // Initialize slab & dimensions
  svector->mm_slab = mm_slab;
  svector->segment_size = mm_slab_get_slab_size(mm_slab);
  svector->element_size = element_size;
  svector->elements_per_segment = svector->segment_size/svector->element_size;
  svector->min_resident_segments = 1;
  // Init segments
  svector->segments = vector_new(SVECTOR_INITIAL_SEGMENTS,vector_segment_t);
  vector_set_used(svector->segments,SVECTOR_INITIAL_SEGMENTS_ALLOCATED);
  VECTOR_ITERATE(svector->segments,segment,position,vector_segment_t) {
    segment->slab_unit = mm_slab_request(mm_slab);
    segment->memory = segment->slab_unit->memory;
  }
  // Reset & return
  svector->elements_used = 0;
  return svector;
}
void svector_delete(svector_t* const svector) {
  // Return all slabs
  mm_slab_lock(svector->mm_slab);
  VECTOR_ITERATE(svector->segments,segment,position,vector_segment_t) {
    mm_slab_put(svector->mm_slab,segment->slab_unit);
  }
  mm_slab_unlock(svector->mm_slab);
  // Free segments' vector
  vector_delete(svector->segments);
  // Free svector
  mm_free(svector);
}
void svector_clear(svector_t* const svector) {
  svector->elements_used = 0;
  svector_reap(svector); // Forced reap
}
void svector_reap(svector_t* const svector) {
  if (svector->segments->used > svector->min_resident_segments) {
    // Reap non-resident segments // TODO Register Number of reaps
    mm_slab_lock(svector->mm_slab);
    VECTOR_ITERATE_OFFSET(svector->segments,segment,position,svector->min_resident_segments,vector_segment_t) {
      mm_slab_put(svector->mm_slab,segment->slab_unit);
    }
    mm_slab_unlock(svector->mm_slab);
    vector_set_used(svector->segments,svector->min_resident_segments);
  }
}
/*
 * Accessors
 */
#define svector_get_elm_memory(svector,segment,segment_pos) \
    vector_get_elm(svector->segments,segment,vector_segment_t)->memory + segment_pos*svector->element_size
#define svector_get_location(svector,global_position,segment,segment_pos) \
  const uint64_t segment = global_position/svector->elements_per_segment; \
  const uint64_t segment_pos = global_position%svector->elements_per_segment
uint64_t svector_get_used(svector_t* const svector) {
  return svector->elements_used;
}
void* svector_get_elm(
    svector_t* const svector,
    const uint64_t position) {
  svector_get_location(svector,position,segment,segment_pos);
  return svector_get_elm_memory(svector,segment,segment_pos);
}
vector_segment_t* svector_add_segment(svector_t* const svector) {
  // Add segment
  vector_reserve_additional(svector->segments,1);
  vector_inc_used(svector->segments);
  vector_segment_t* const segment = vector_get_last_elm(svector->segments,vector_segment_t);
  // Init segment
  segment->slab_unit = mm_slab_request(svector->mm_slab);
  segment->memory = segment->slab_unit->memory;
  return segment;
}
void* svector_get_free_elm(svector_t* const svector) {
  svector_get_location(svector,svector->elements_used,num_segment,segment_pos);
  vector_segment_t* const segment = gem_expect_false(num_segment==vector_get_used(svector->segments)) ?
    svector_add_segment(svector) : vector_get_elm(svector->segments,num_segment,vector_segment_t);
  ++(svector->elements_used);
  return segment + segment_pos*svector->element_size;
}
char* svector_request_char_buffer(
    svector_t* const svector,
    uint64_t* const buffer_offset,
    const uint64_t length) {
  svector_get_location(svector,svector->elements_used,num_segment,segment_pos);
  // Check that there is enough space in the segment
  const uint64_t remaining = svector->segment_size-segment_pos;
  const uint64_t total_length = length+1;
  uint64_t segment_offset = 0;
  vector_segment_t* segment;
  if (gem_expect_false(num_segment >= vector_get_used(svector->segments))) {
    gem_cond_fatal_error(total_length>svector->segment_size,SVECTOR_INSERT_CHAR_BUFFER_TOO_LONG,total_length);
    segment = svector_add_segment(svector); // Add new segment
  } else if (gem_expect_false(remaining < total_length)) {
    gem_cond_fatal_error(total_length>svector->segment_size,SVECTOR_INSERT_CHAR_BUFFER_TOO_LONG,total_length);
    segment = svector_add_segment(svector); // Add new segment
    svector->elements_used += remaining;
  } else {
    segment = vector_get_elm(svector->segments,num_segment,vector_segment_t); // Add to last segment
    segment_offset = segment_pos; // Return
  }
  // Return buffer + offset
  if (buffer_offset) *buffer_offset = svector->elements_used; // Set offset
  svector->elements_used += total_length; // Update used
  return segment->memory + segment_offset; // Return
}
char* svector_insert_char_buffer(
    svector_t* const svector,
    uint64_t* const buffer_offset,
    const char* const buffer,
    const uint64_t length) {
  // Request char buffer
  char* const dst_memory = svector_request_char_buffer(svector,buffer_offset,length);
  // Copy buffer
  memcpy(dst_memory,buffer,length);
  dst_memory[length] = '\0';
  return dst_memory;
}
/*
 * Writer
 */
void svector_write(
    fm_t* const file_manager,
    svector_t* const svector) {
  // Write all the segments
  uint64_t pending_elements = svector->elements_used;
  VECTOR_ITERATE(svector->segments,segment,position,vector_segment_t) {
    if (pending_elements > svector->elements_per_segment) {
      // Write the whole segment
      fm_write_mem(file_manager,segment->memory,svector->elements_per_segment*svector->element_size);
      pending_elements -= svector->elements_per_segment;
    } else {
      // Write the remaining part of the vector
      fm_write_mem(file_manager,segment->memory,pending_elements*svector->element_size);
      break;
    }
  }
}
/*
 * Display/Profile
 */
void svector_print(
    FILE* const stream,
    svector_t* const svector) {
  // TODO // TODO // TODO // TODO
  // TODO // TODO // TODO // TODO
  // TODO // TODO // TODO // TODO
  // TODO // TODO // TODO // TODO
  // TODO // TODO // TODO // TODO
  // TODO // TODO // TODO // TODO
}
void svector_record_stats(svector_t* const svector) {
  // TODO // TODO // TODO // TODO
  // TODO // TODO // TODO // TODO
  // TODO // TODO // TODO // TODO
  // TODO // TODO // TODO // TODO
  // TODO // TODO // TODO // TODO
  // TODO // TODO // TODO // TODO
}
void svector_display_stats(
    FILE* const stream,
    svector_t* const svector) {
  // TODO // TODO // TODO // TODO
  // TODO // TODO // TODO // TODO
  // TODO // TODO // TODO // TODO
  // TODO // TODO // TODO // TODO
  // TODO // TODO // TODO // TODO
  // TODO // TODO // TODO // TODO
}
/*
 * Iterators
 */
#define svector_set_iterator_location(iterator,svector,position,segment,segment_pos) \
  svector_get_location(svector,position,segment,segment_pos); \
  iterator->current_segment = segment; \
  iterator->global_position = position; \
  iterator->local_position = segment_pos
void svector_iterator_new(
    svector_iterator_t* const iterator,
    svector_t* const svector,
    const svector_iterator_type iterator_type,
    const uint64_t init_position) {
  // Init
  iterator->svector = svector;
  iterator->element_size = svector->element_size;
  iterator->elements_per_segment = svector->elements_per_segment;
  iterator->elements_used = svector->elements_used;
  iterator->type = iterator_type;
  // Seek
  switch (iterator->type) {
    case SVECTOR_READ_ITERATOR:
      svector_read_iterator_seek(iterator,init_position);
      break;
    case SVECTOR_WRITE_ITERATOR:
    {
      /*
       * Writes always from the beginning (simplifies the code and for the moment seeking in writing mode is no needed)
       */
      iterator->current_segment = 0;
      iterator->global_position = 0;
      iterator->local_position = 0;
      iterator->memory = (gem_expect_false(vector_get_used(svector->segments)==0)) ?
        svector_add_segment(svector)->memory :
        vector_get_elm(svector->segments,0,vector_segment_t)->memory;
      break;
    }
    default:
      GEM_INVALID_CASE();
      break;
  }
}
void svector_read_iterator_seek(svector_iterator_t* const iterator,const uint64_t init_position) {
  svector_t* const svector = iterator->svector;
  // Locate
  svector_set_iterator_location(iterator,svector,init_position,segment,segment_pos);
  iterator->eoi = (init_position>=svector->elements_used);
  if (!iterator->eoi) iterator->memory = svector_get_elm_memory(svector,segment,segment_pos);
}
void* svector_iterator_get_elm(svector_iterator_t* const iterator) {
  return iterator->memory;
}
// Reading
bool svector_read_iterator_eoi(svector_iterator_t* const iterator) {
  return iterator->eoi;
}
void svector_read_iterator_next(svector_iterator_t* const iterator) {
  // Inc position
  ++(iterator->local_position);
  ++(iterator->global_position);
  // Check boundary condition
  if (gem_expect_true(iterator->global_position < iterator->elements_used)) {
    if (gem_expect_true(iterator->local_position < iterator->elements_per_segment)) {
      // Regular next
      iterator->memory += iterator->element_size;
    } else {
      // Next segment
      ++(iterator->current_segment);
      iterator->local_position = 0;
      iterator->memory = vector_get_elm(iterator->svector->segments,iterator->current_segment,vector_segment_t)->memory;
    }
  } else {
    iterator->eoi = true;
  }
}
// Writing
void svector_write_iterator_next(svector_iterator_t* const iterator) {
  // Inc position
  ++(iterator->local_position);
  // Check boundary condition
  if (gem_expect_true(iterator->local_position < iterator->elements_per_segment)) {
    // Regular next
    iterator->memory += iterator->element_size;
  } else {
    // Next segment
    ++(iterator->current_segment);
    iterator->local_position = 0;
    iterator->memory = (gem_expect_false(iterator->current_segment < vector_get_used(iterator->svector->segments))) ?
        vector_get_elm(iterator->svector->segments,iterator->current_segment,vector_segment_t)->memory :
        svector_add_segment(iterator->svector)->memory;
  }
  ++(iterator->svector->elements_used);
}
