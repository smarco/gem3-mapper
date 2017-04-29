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

#ifndef TEXT_TRACE_H_
#define TEXT_TRACE_H_

#include "utils/essentials.h"
#include "text/dna_text.h"

/*
 * Text Collection
 */
typedef struct {
  /* Location */
  uint64_t position;            // Text position
  strand_t strand;              // Strand
  /* Encoded Text */
  uint8_t* text;                // Encoded text
  bool text_allocated;          // MM-Allocated text
  uint64_t text_length;         // Text length
  /* Padded Text */
  uint8_t* text_padded;         // Encoded text
  bool text_padded_allocated;   // MM-Allocated text-padded
  uint64_t text_padded_left;    // Text left padding
  uint64_t text_padded_right;   // Text right padding
  uint64_t text_padded_length;  // Text length
  /* RL-Encoded Text */
  uint8_t* rl_text;             // RL-Encoded Text
  uint64_t rl_text_length;      // RL-Encoded Text length
  uint32_t* rl_runs_acc;        // Length of each run (accumulated)
} text_trace_t;

/*
 * Setup
 */
void text_trace_destroy(
    text_trace_t* const text_trace,
    mm_allocator_t* const mm_allocator);

/*
 * Padding
 */
void text_trace_compose_padded_text(
    text_trace_t* const text_trace,
    const uint64_t key_trim_left,
    const uint64_t key_trim_right,
    mm_allocator_t* const mm_allocator);

#endif /* TEXT_TRACE_H_ */
