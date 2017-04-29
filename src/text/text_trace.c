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

#include "text/text_trace.h"
#include "text/dna_text.h"

/*
 * Setup
 */
void text_trace_destroy(
    text_trace_t* const text_trace,
    mm_allocator_t* const mm_allocator) {
  if (text_trace->text != NULL) {
    if (text_trace->text_allocated) mm_allocator_free(mm_allocator,text_trace->text);
    if (text_trace->text_padded_allocated) mm_allocator_free(mm_allocator,text_trace->text_padded);
    if (text_trace->rl_text != NULL) mm_allocator_free(mm_allocator,text_trace->rl_text);
    if (text_trace->rl_runs_acc != NULL) mm_allocator_free(mm_allocator,text_trace->rl_runs_acc);
  }
}

/*
 * Padding
 */
void text_trace_compose_padded_text(
    text_trace_t* const text_trace,
    const uint64_t key_trim_left,
    const uint64_t key_trim_right,
    mm_allocator_t* const mm_allocator) {
  // Check trims
  if (key_trim_left > 0 || key_trim_right > 0) {
    // Allocate
    const uint64_t text_padded_length = text_trace->text_length + key_trim_left + key_trim_right;
    text_trace->text_padded = mm_allocator_calloc(mm_allocator,text_padded_length,uint8_t,false);
    text_trace->text_padded_allocated = true;
    text_trace->text_padded_length = text_padded_length;
    // Add left-trim
    uint64_t i;
    for (i=0;i<key_trim_left;++i) text_trace->text_padded[i] = ENC_DNA_CHAR_N;
    // Copy original text
    memcpy(text_trace->text_padded+i,text_trace->text,text_trace->text_length);
    i += text_trace->text_length;
    // Add rigth-trim
    for (;i<text_padded_length;++i) text_trace->text_padded[i] = ENC_DNA_CHAR_N;
  }
  text_trace->text_padded_left = key_trim_left;
  text_trace->text_padded_right = key_trim_right;
}
