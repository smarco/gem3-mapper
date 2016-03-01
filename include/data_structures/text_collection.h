/*
 * PROJECT: GEMMapper
 * FILE: text_collection.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef TEXT_COLLECTION_H_
#define TEXT_COLLECTION_H_

#include "utils/essentials.h"

/*
 * Text Collection
 */
typedef struct {
  /* Encoded Text    */
  uint8_t* text;                // Encoded text
  uint64_t text_length;         // Text length
  /* RL-Encoded Text */
  uint8_t* rl_text;             // RL-Encoded Text
  uint64_t rl_text_length;      // RL-Encoded Text length
  uint32_t* rl_runs_acc;        // Length of each run (accumulated)
} text_trace_t;
typedef struct {
  /* Text-Traces */
  vector_t* text_traces;        // Text-Traces (text_trace_t)
} text_collection_t;

/*
 * Setup
 */
void text_collection_init(text_collection_t* const text_collection);
void text_collection_clear(text_collection_t* const text_collection);
void text_collection_destroy(text_collection_t* const text_collection);

/*
 * Accessors
 */
// [Text-Block]
uint64_t text_collection_new_trace(const text_collection_t* const text_collection);
text_trace_t* text_collection_get_trace(
    const text_collection_t* const text_collection,
    const uint64_t text_trace_offset);
uint64_t text_collection_get_num_traces(const text_collection_t* const text_collection);

#endif /* TEXT_COLLECTION_H_ */
