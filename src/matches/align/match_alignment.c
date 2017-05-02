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

#include "matches/align/match_alignment.h"
#include "matches/matches_cigar.h"
#include "align/alignment.h"
#include "text/dna_text.h"

/*
 * Check
 */
bool match_alignment_check(
    FILE* const stream,
    match_alignment_t* const match_alignment,
    uint8_t* const key,
    const uint64_t key_length,
    uint8_t* const text,
    const uint64_t text_length,
    vector_t* const cigar_vector,
    const bool verbose,
    mm_allocator_t* const mm_allocator) {
  return alignment_check(
      stream,key,key_length,
      text+match_alignment->match_text_offset,
      match_alignment->effective_length,
      cigar_vector,match_alignment->cigar_offset,
      match_alignment->cigar_length,false,verbose);
}
/*
 * Display
 */
void match_alignment_print_pretty(
    FILE* const stream,
    match_alignment_t* const match_alignment,
    uint8_t* const key,
    const uint64_t key_length,
    uint8_t* const text,
    const uint64_t text_length,
    vector_t* const cigar_vector,
    mm_allocator_t* const mm_allocator) {
  return alignment_print_pretty(
      stream,key,key_length,
      text+match_alignment->match_text_offset,
      match_alignment->effective_length,
      cigar_vector,match_alignment->cigar_offset,
      match_alignment->cigar_length,mm_allocator);
}

