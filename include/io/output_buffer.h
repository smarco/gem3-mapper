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
 * DESCRIPTION:
 *   Output buffer provides storage before final dump to output file
 */

#ifndef OUTPUT_BUFFER_H_
#define OUTPUT_BUFFER_H_

#include "utils/essentials.h"

typedef enum { OUTPUT_BUFFER_FREE, OUTPUT_BUFFER_BUSY, OUTPUT_BUFFER_WRITE_PENDING } output_buffer_state_t;

typedef struct {
  /* Buffer state */
  uint32_t mayor_block_id; // Block mayor ID (for synchronization purposes)
  uint32_t minor_block_id; // Block minor ID (for segmented dumps)
  bool is_final_block;
  output_buffer_state_t buffer_state;
  /* Buffer */
  char* buffer_mem;        // Buffer Memory
  uint64_t buffer_used;    // Number of characters written
  char* buffer_cursor;     // Pointer to the current char
} output_buffer_t;

/*
 * Setup
 */
output_buffer_t* output_buffer_new(const uint64_t output_buffer_size);
void output_buffer_clear(output_buffer_t* const out_buffer);
void output_buffer_delete(output_buffer_t* const out_buffer);

/*
 * Accessors
 */
void output_buffer_set_state(output_buffer_t* const output_buffer,const output_buffer_state_t buffer_state);
output_buffer_state_t output_buffer_get_state(output_buffer_t* const output_buffer);
void output_buffer_set_incomplete(output_buffer_t* const output_buffer);
uint64_t output_buffer_get_used(output_buffer_t* const output_buffer);
char* output_buffer_get_buffer(output_buffer_t* const out_buffer);

/*
 * Fast-printer functions
 */
void bprintf_uint64(output_buffer_t* const out_buffer,const uint64_t number);
void bprintf_int64(output_buffer_t* const out_buffer,const int64_t number);
void bprintf_char(output_buffer_t* const out_buffer,const char character);
void bprintf_buffer(output_buffer_t* const out_buffer,const int string_length,const char* const string);

#endif /* OUTPUT_BUFFER_H_ */
