/*
 * PROJECT: GEMMapper
 * FILE: text_collection.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "text_collection.h"

#define TEXT_COLLECTION_INIT_TEXT_BLOCKS  1000
#define TEXT_COLLECTION_INIT_TRACES       1000
#define TEXT_COLLECTION_INIT_TRACE_BLOCKS 1000

/*
 * Setup
 */
GEM_INLINE text_collection_t* text_collection_new() {
  // Allocate
  text_collection_t* const text_collection = mm_alloc(text_collection_t);
  // Request slabs
  mm_slab_t* const mm_slab_text = mm_pool_get_slab(mm_pool_8MB);
  // Text-Blocks
  text_collection->text_block_vector = vector_new(TEXT_COLLECTION_INIT_TEXT_BLOCKS,text_block_t);
  text_collection->text_buffer = svector_new(mm_slab_text,uint8_t);
  // Traces (Sequence of Text-Blocks)
  text_collection->trace_vector = vector_new(TEXT_COLLECTION_INIT_TRACES,text_trace_t);
  text_collection->trace_blocks_vector = vector_new(TEXT_COLLECTION_INIT_TRACE_BLOCKS,text_block_t*);
  // Return
  return text_collection;
}
GEM_INLINE void text_collection_clear(text_collection_t* const text_collection) {
  TEXT_COLLECTION_CHECK(text_collection);
  // Clear Text-Blocks & Traces // TODO: Reap vectors. vector->num_init_elements based.
  vector_clear(text_collection->text_block_vector);
  svector_clear(text_collection->text_buffer);
  vector_clear(text_collection->trace_vector);
  vector_clear(text_collection->trace_blocks_vector);
}
GEM_INLINE void text_collection_delete(text_collection_t* const text_collection) {
  TEXT_COLLECTION_CHECK(text_collection);
  // Delete Text-Blocks & Traces
  vector_delete(text_collection->text_block_vector);
  svector_delete(text_collection->text_buffer);
  vector_delete(text_collection->trace_vector);
  vector_delete(text_collection->trace_blocks_vector);
  // Delete handler
  mm_free(text_collection);
}
/*
 * Accessors
 */
// [Text-Block]
GEM_INLINE text_block_t* text_collection_new_block(text_collection_t* const text_collection) {
  TEXT_COLLECTION_CHECK(text_collection);
  // Allocate new text-block
  text_block_t* text_block;
  vector_alloc_new(text_collection->text_block_vector,text_block_t,text_block);
  return text_block;
}
GEM_INLINE text_block_t* text_collection_get_block(
    text_collection_t* const text_collection,const uint64_t text_block_offset) {
  TEXT_COLLECTION_CHECK(text_collection);
  return vector_get_elm(text_collection->text_block_vector,text_block_offset,text_block_t);
}
// [Text-Trace]
GEM_INLINE text_trace_t* text_collection_new_trace(text_collection_t* const text_collection,const uint64_t num_trace_blocks) {
  TEXT_COLLECTION_CHECK(text_collection);
  // Allocate new trace-block
  text_trace_t* text_trace;
  vector_alloc_new(text_collection->trace_vector,text_trace_t,text_trace);
  text_trace->trace_blocks_offset = vector_get_used(text_collection->trace_blocks_vector);
  text_trace->num_trace_blocks = num_trace_blocks;
  // Allocate @num_text_blocks
  vector_reserve_additional(text_collection->trace_blocks_vector,num_trace_blocks);
  // Return trace
  return text_trace;
}
GEM_INLINE uint64_t text_collection_get_num_traces(text_collection_t* const text_collection) {
  TEXT_COLLECTION_CHECK(text_collection);
  return vector_get_used(text_collection->trace_vector);
}
GEM_INLINE text_trace_t* text_collection_get_trace(text_collection_t* const text_collection,const uint64_t trace_offset) {
  TEXT_COLLECTION_CHECK(text_collection);
  return vector_get_elm(text_collection->trace_vector,trace_offset,text_trace_t);
}
GEM_INLINE uint64_t text_collection_get_trace_text_length(text_trace_t* const text_trace) {
  // TODO
  GEM_NOT_IMPLEMENTED();
  return 0;
}
GEM_INLINE uint64_t text_collection_get_trace_num_blocks(text_trace_t* const text_trace) {
  TEXT_TRACE_CHECK(text_trace);
  return text_trace->num_trace_blocks;
}
GEM_INLINE text_block_t** text_collection_get_trace_block(
    text_collection_t* const text_collection,text_trace_t* const text_trace,const uint64_t trace_block_offset) {
  TEXT_COLLECTION_CHECK(text_collection);
  gem_fatal_check(trace_block_offset>=text_trace->num_trace_blocks,TEXT_COLLECTION_BLOCK_NUM_OUT_OF_TRACE,
      trace_block_offset,text_trace->trace_blocks_offset,text_trace->trace_blocks_offset+text_trace->num_trace_blocks);
  return vector_get_elm(text_collection->trace_blocks_vector,text_trace->trace_blocks_offset+trace_block_offset,text_block_t*);
}
/*
 * Text-Trace Iterator
 */
GEM_INLINE void text_trace_iterator_new(
    text_trace_iterator_t* const text_iterator,text_collection_t* const text_collection,
    const uint64_t trace_offset,const traversal_direction_t traversal_direction) {
  GEM_CHECK_NULL(text_iterator);
  text_iterator->text_collection = text_collection;
  text_iterator->traversal_direction = traversal_direction;
  // Get Trace
  text_trace_t* const text_trace = text_collection_get_trace(text_collection,trace_offset);
  text_iterator->trace_begin_block_offset = text_trace->trace_blocks_offset;
  text_iterator->trace_end_block_offset = text_trace->trace_blocks_offset+text_trace->num_trace_blocks-1;
  // Locate
  if (traversal_direction==traversal_forward) {
    text_iterator->eoi = false;
    text_iterator->trace_current_block_offset = text_iterator->trace_begin_block_offset;
    text_iterator->current_position = 0;
    // Get current Text-Block
    text_block_t* const text_block = text_collection_get_block(text_collection,text_iterator->trace_current_block_offset);
    text_iterator->text = text_block->text;
    text_iterator->text_length = text_block->length;
  } else {
    text_iterator->eoi = false;
    text_iterator->trace_current_block_offset = text_iterator->trace_end_block_offset;
    text_iterator->current_position = 0;
    // Get current Text-Block
    text_block_t* const text_block = text_collection_get_block(text_collection,text_iterator->trace_current_block_offset);
    text_iterator->text = text_block->text;
    text_iterator->text_length = text_block->length;
    text_iterator->current_position = text_block->length-1;
  }
}
GEM_INLINE bool text_trace_iterator_eoi(text_trace_iterator_t* const text_iterator) {
  TEXT_TRACE_ITERATOR_CHECK(text_iterator);
  return text_iterator->eoi;
}
GEM_INLINE uint8_t text_trace_iterator_get_char(text_trace_iterator_t* const text_iterator) {
  TEXT_TRACE_ITERATOR_CHECK(text_iterator);
  return text_iterator->text[text_iterator->current_position];
}
GEM_INLINE void text_trace_iterator_next_char(text_trace_iterator_t* const text_iterator) {
  TEXT_TRACE_ITERATOR_CHECK(text_iterator);
  gem_fatal_check(text_iterator->traversal_direction!=traversal_forward,TEXT_TRACE_ITERATOR_WRONG_DIRECTION);
  ++(text_iterator->current_position);
  // Check end of text-block
  if (gem_expect_false(text_iterator->current_position >= text_iterator->text_length)) {
    // Check end of trace
    if (gem_expect_false(text_iterator->trace_current_block_offset >= text_iterator->trace_end_block_offset)) {
      text_iterator->eoi = true;
    } else {
      ++(text_iterator->trace_current_block_offset);
      // Get current Text-Block
      text_block_t* const text_block = text_collection_get_block(text_iterator->text_collection,text_iterator->trace_current_block_offset);
      text_iterator->text = text_block->text;
      text_iterator->text_length = text_block->length;
      text_iterator->current_position = 0;
    }
  }
}
GEM_INLINE void text_trace_iterator_previous_char(text_trace_iterator_t* const text_iterator) {
  TEXT_TRACE_ITERATOR_CHECK(text_iterator);
  gem_fatal_check(text_iterator->traversal_direction!=traversal_backward,TEXT_TRACE_ITERATOR_WRONG_DIRECTION);
  // Check end of text-block
  if (gem_expect_false(text_iterator->current_position==0)) {
    // Check end of trace
    if (gem_expect_false(text_iterator->trace_current_block_offset == text_iterator->trace_begin_block_offset)) {
      text_iterator->eoi = true;
    } else {
      --(text_iterator->trace_current_block_offset);
      // Get current Text-Block
      text_block_t* const text_block = text_collection_get_block(text_iterator->text_collection,text_iterator->trace_current_block_offset);
      text_iterator->text = text_block->text;
      text_iterator->text_length = text_block->length;
      text_iterator->current_position = text_block->length-1;
    }
  } else {
    --(text_iterator->current_position);
  }
}
