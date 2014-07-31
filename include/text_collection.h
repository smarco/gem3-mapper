/*
 * PROJECT: GEMMapper
 * FILE: text_collection.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef TEXT_COLLECTION_H_
#define TEXT_COLLECTION_H_

#include "essentials.h"

/*
 * Checkers
 */
#define TEXT_BLOCK_CHECK(text_block) GEM_CHECK_NULL(text_block)
#define TEXT_TRACE_CHECK(text_trace) GEM_CHECK_NULL(text_trace)
#define TEXT_COLLECTION_CHECK(text_collection) GEM_CHECK_NULL(text_collection)
#define TEXT_TRACE_ITERATOR_CHECK(iterator) GEM_CHECK_NULL(iterator)

/*
 * Text Collection
 */
typedef struct {
  /* Location */
  uint64_t position;               // Position at index // FIXME: Dubious to need it
  /* Encoded Text */
  uint8_t* text;                   // Encoded text
  uint64_t length;                 // Text-Block length
} text_block_t;
typedef struct {
  uint64_t trace_blocks_offset;
  uint64_t num_trace_blocks;
} text_trace_t;
typedef struct {
  /* Text-Blocks */
  vector_t* text_block_vector;              // Text-Blocks (text_block_t)
  svector_t* text_buffer;                   // Encoded Text-Blocks Buffer (uint8_t)
  /* Traces (Sequence of Text-Blocks) */
  vector_t* trace_vector;                   // Traces (text_trace_t)
  vector_t* trace_blocks_vector;            // Traces-Buffer (Text-Blocks Pointers) (text_block_t*)
} text_collection_t;

typedef struct {
  /*
   * Raw storage of all reads for Myers-CUDA API // TODO
   *   -- Based on trace_offsets for compatibility with
   */
} text_collection_raw_t;

/*
 * Text Iterator
 */
typedef struct {
  /* Iterator info */
  text_collection_t* text_collection;
  traversal_direction_t traversal_direction;
  /* Current State */
  bool eoi;
  /* Current Trace */
  uint64_t trace_current_block_offset;
  uint64_t trace_begin_block_offset;
  uint64_t trace_end_block_offset;
  /* Current Block */
  uint64_t current_position;
  uint8_t* text;
  uint64_t text_length;
} text_trace_iterator_t;

/*
 * Setup
 */
GEM_INLINE text_collection_t* text_collection_new();
GEM_INLINE void text_collection_clear(text_collection_t* const text_collection);
GEM_INLINE void text_collection_delete(text_collection_t* const text_collection);

/*
 * Accessors
 */
// [Text-Block]
GEM_INLINE text_block_t* text_collection_new_block(text_collection_t* const text_collection);
GEM_INLINE text_block_t* text_collection_get_block(
    text_collection_t* const text_collection,const uint64_t text_block_offset);
// [Text-Trace]
GEM_INLINE text_trace_t* text_collection_new_trace(text_collection_t* const text_collection,const uint64_t num_trace_blocks);
GEM_INLINE uint64_t text_collection_get_num_traces(text_collection_t* const text_collection);
GEM_INLINE text_trace_t* text_collection_get_trace(text_collection_t* const text_collection,const uint64_t trace_offset);
GEM_INLINE uint64_t text_collection_get_trace_text_length(text_trace_t* const text_trace);
GEM_INLINE uint64_t text_collection_get_trace_num_blocks(text_trace_t* const text_trace);
GEM_INLINE text_block_t** text_collection_get_trace_block(text_collection_t* const text_collection,text_trace_t* const text_trace,const uint64_t trace_block_offset);
/*
 * Text-Trace Iterator
 */
GEM_INLINE void text_trace_iterator_new(
    text_trace_iterator_t* const text_iterator,text_collection_t* const text_collection,
    const uint64_t trace_offset,const traversal_direction_t traversal_direction);

GEM_INLINE bool text_trace_iterator_eoi(text_trace_iterator_t* const text_iterator);
GEM_INLINE uint8_t text_trace_iterator_get_char(text_trace_iterator_t* const text_iterator);

GEM_INLINE void text_trace_iterator_next_char(text_trace_iterator_t* const text_iterator);
GEM_INLINE void text_trace_iterator_previous_char(text_trace_iterator_t* const text_iterator);

/*
 * Error Messages
 */
#define GEM_ERROR_TEXT_COLLECTION_BLOCK_NUM_OUT_OF_TRACE "Text-Collection. Requested Text-Block (%lu) out of trace [%lu,%lu)"
#define GEM_ERROR_TEXT_TRACE_ITERATOR_WRONG_DIRECTION "Text-Trace iterator. Wrong traversal direction"

#endif /* TEXT_COLLECTION_H_ */
