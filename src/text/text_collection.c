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

#include "text/text_collection.h"

/*
 * Constants
 */
#define TEXT_COLLECTION_INIT_TEXT_TRACES  1000

/*
 * Setup
 */
void text_collection_init(text_collection_t* const text_collection) {
  text_collection->text_traces = vector_new(TEXT_COLLECTION_INIT_TEXT_TRACES,text_trace_t);
}
void text_collection_clear(text_collection_t* const text_collection) {
  vector_clear(text_collection->text_traces);
}
void text_collection_destroy(text_collection_t* const text_collection) {
  vector_delete(text_collection->text_traces);
}
void text_collection_inject_mm(
    text_collection_t* const text_collection,
    mm_stack_t* const mm_text) {
  text_collection->mm_text = mm_text;
}
/*
 * Accessors
 */
uint64_t text_collection_new_trace(const text_collection_t* const text_collection) {
  // Allocate new text-trace
  const uint64_t text_trace_offset = vector_get_used(text_collection->text_traces);
  vector_reserve_additional(text_collection->text_traces,1);
  vector_add_used(text_collection->text_traces,1);
  // Return text_trace_offset
  return text_trace_offset;
}
text_trace_t* text_collection_get_trace(
    const text_collection_t* const text_collection,
    const uint64_t text_trace_offset) {
  return vector_get_elm(text_collection->text_traces,text_trace_offset,text_trace_t);
}
uint64_t text_collection_get_num_traces(const text_collection_t* const text_collection) {
  return vector_get_used(text_collection->text_traces);
}
