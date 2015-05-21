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
#define TEXT_TRACE_CHECK(text_trace) GEM_CHECK_NULL(text_trace)
#define TEXT_COLLECTION_CHECK(text_collection) GEM_CHECK_NULL(text_collection)

/*
 * Text Collection
 */
typedef struct {
  /* Encoded Text */
  uint8_t* text;                // Encoded text
  uint64_t length;              // Text length
  /* Trace */
  uint64_t trace_blocks_offset;
  uint64_t trace_length;
} text_trace_t;
typedef struct {
  uint64_t position;
  uint64_t length;
} trace_block_t;
typedef struct {
  /* Text-Blocks */
  vector_t* text_traces;        // Text-Traces (text_trace_t)
  // svector_t* text_buffer;       // Encoded Text Buffer (uint8_t) // TODO Graph
  vector_t* trace_block_vector; // Traces (trace_block_t)
} text_collection_t;

/*
 * Setup
 */
GEM_INLINE void text_collection_init(text_collection_t* const text_collection);
GEM_INLINE void text_collection_clear(text_collection_t* const text_collection);
GEM_INLINE void text_collection_destroy(text_collection_t* const text_collection);

/*
 * Accessors
 */
// [Text-Block]
GEM_INLINE uint64_t text_collection_new_trace(const text_collection_t* const text_collection);
GEM_INLINE text_trace_t* text_collection_get_trace(
    const text_collection_t* const text_collection,const uint64_t text_trace_offset);
GEM_INLINE uint64_t text_collection_get_num_traces(const text_collection_t* const text_collection);
// [Text-Trace]
GEM_INLINE uint64_t text_collection_allocate_trace_blocks(
    const text_collection_t* const text_collection,const uint64_t num_trace_blocks);
GEM_INLINE trace_block_t* text_collection_get_trace_block(
    const text_collection_t* const text_collection,const uint64_t trace_block_offset);

/*
 * Error Messages
 */
//#define GEM_ERROR_TEXT_COLLECTION_ ""

#endif /* TEXT_COLLECTION_H_ */
