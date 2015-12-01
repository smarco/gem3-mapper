/*
 * PROJECT: GEMMapper
 * FILE: text_collection.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "text_collection.h"

#define TEXT_COLLECTION_INIT_TEXT_TRACES  1000
#define TEXT_COLLECTION_INIT_TRACE_BLOCKS 1000

/*
 * Setup
 */
void text_collection_init(text_collection_t* const text_collection) {
  // Text-Blocks
  text_collection->text_traces = vector_new(TEXT_COLLECTION_INIT_TEXT_TRACES,text_trace_t);
  // mm_slab_t* const mm_slab_text = mm_pool_get_slab(mm_pool_8MB); // Request slabs // TODO Graph
  // text_collection->text_buffer = svector_new(mm_slab_text,uint8_t); // TODO Graph
  // Traces (Sequence of Text-Blocks)
  text_collection->trace_block_vector = vector_new(TEXT_COLLECTION_INIT_TRACE_BLOCKS,trace_block_t);
}
void text_collection_clear(text_collection_t* const text_collection) {
  TEXT_COLLECTION_CHECK(text_collection);
  // Clear Text-Blocks & Traces
  vector_clear(text_collection->text_traces);
  // svector_clear(text_collection->text_buffer); // TODO Graph
  vector_clear(text_collection->trace_block_vector);
}
void text_collection_destroy(text_collection_t* const text_collection) {
  TEXT_COLLECTION_CHECK(text_collection);
  // Delete Text-Blocks & Traces
  vector_delete(text_collection->text_traces);
  // svector_delete(text_collection->text_buffer); // TODO Graph
  vector_delete(text_collection->trace_block_vector);
}
/*
 * Accessors
 */
// [Text-Trace]
uint64_t text_collection_new_trace(const text_collection_t* const text_collection) {
  TEXT_COLLECTION_CHECK(text_collection);
  // Allocate new text-trace
  const uint64_t text_trace_offset = vector_get_used(text_collection->text_traces);
  vector_reserve_additional(text_collection->text_traces,1);
  vector_add_used(text_collection->text_traces,1);
  // Return text_trace_offset
  return text_trace_offset;
}
text_trace_t* text_collection_get_trace(
    const text_collection_t* const text_collection,const uint64_t text_trace_offset) {
  TEXT_COLLECTION_CHECK(text_collection);
  return vector_get_elm(text_collection->text_traces,text_trace_offset,text_trace_t);
}
uint64_t text_collection_get_num_traces(const text_collection_t* const text_collection) {
  TEXT_COLLECTION_CHECK(text_collection);
  return vector_get_used(text_collection->text_traces);
}
// [Trace-Blocks]
uint64_t text_collection_allocate_trace_blocks(
    const text_collection_t* const text_collection,const uint64_t num_trace_blocks) {
  TEXT_COLLECTION_CHECK(text_collection);
  // Allocate new trace-blocks
  const uint64_t trace_blocks_offset = vector_get_used(text_collection->trace_block_vector);
  vector_reserve_additional(text_collection->trace_block_vector,num_trace_blocks);
  vector_add_used(text_collection->trace_block_vector,num_trace_blocks);
  // Return trace_blocks_offset
  return trace_blocks_offset;
}
trace_block_t* text_collection_get_trace_block(
    const text_collection_t* const text_collection,const uint64_t trace_block_offset) {
  TEXT_COLLECTION_CHECK(text_collection);
  return vector_get_elm(text_collection->trace_block_vector,trace_block_offset,trace_block_t);
}

