/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
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
  /* Text Traces */
  vector_t* text_traces;        // Text-Traces (text_trace_t)
  /* MM */
  mm_stack_t* mm_text;          // MM-Stack
} text_collection_t;

/*
 * Setup
 */
void text_collection_init(text_collection_t* const text_collection);
void text_collection_clear(text_collection_t* const text_collection);
void text_collection_destroy(text_collection_t* const text_collection);

void text_collection_inject_mm(
    text_collection_t* const text_collection,
    mm_stack_t* const mm_text);

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
