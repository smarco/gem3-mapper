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

#include "matches/align/match_alignment_region_rl.h"
#include "matches/matches_cigar.h"
#include "archive/archive_text_rl.h"

/*
 * RL-Translate CIGAR-Match
 */
void match_alignment_region_rl_translate_cigar(
    cigar_element_t** const cigar_buffer,
    const uint64_t rl_match_length,
    const uint32_t* const rl_key_runs,
    const uint32_t* const rl_text_runs,
    uint64_t rl_key_pos,
    uint64_t rl_text_pos,
    const bool left_gap_alignment) {
  // Traverse all matching characters
  uint64_t i;
  for (i=0;i<rl_match_length;++i) {
    // Check Run Length
    const uint64_t rl_key = archive_text_rl_get_run_length(rl_key_runs,rl_key_pos);
    const uint64_t rl_text = archive_text_rl_get_run_length(rl_text_runs,rl_text_pos);
    // Compute the CIGAR-translated
    cigar_t cigar_op;
    uint64_t match_length, indel_length;
    if (rl_key > rl_text) {
      cigar_op = cigar_del;
      indel_length = rl_key - rl_text;
      match_length = rl_text;
    } else {
      cigar_op = cigar_ins;
      indel_length = rl_text - rl_key;
      match_length = rl_key;
    }
    // Add CIGAR-translated
    if (left_gap_alignment) {
      if (indel_length > 0) {
        matches_cigar_buffer_add_cigar_element__attr(cigar_buffer,cigar_op,indel_length,cigar_attr_homopolymer);
      }
      matches_cigar_buffer_add_cigar_element(cigar_buffer,cigar_match,match_length);
    } else {
      matches_cigar_buffer_add_cigar_element(cigar_buffer,cigar_match,match_length);
      if (indel_length > 0) {
        matches_cigar_buffer_add_cigar_element__attr(cigar_buffer,cigar_op,indel_length,cigar_attr_homopolymer);
      }
    }
    // Next
    ++rl_key_pos;
    ++rl_text_pos;
  }
}
/*
 * RL-Translate CIGAR
 */
void match_alignment_region_rl_translate(
    match_alignment_region_t* const match_alignment_region,
    uint32_t* const rl_key_runs_acc,
    uint32_t* const rl_text_runs_acc,
    const bool left_gap_alignment,
    vector_t* const cigar_vector) {
  // Parameters
  const uint64_t region_key_begin = match_alignment_region_get_key_begin(match_alignment_region);
  const uint64_t region_key_end = match_alignment_region_get_key_end(match_alignment_region);
  const uint64_t region_text_begin = match_alignment_region_get_text_begin(match_alignment_region);
  // Allocate Translated CIGAR
  const uint64_t rl_match_length = region_key_end-region_key_begin;
  const uint64_t cigar_offset = vector_get_used(cigar_vector);
  vector_reserve_additional(cigar_vector,2*rl_match_length);
  cigar_element_t* cigar_buffer = vector_get_free_elm(cigar_vector,cigar_element_t); // Sentinel
  cigar_element_t* const cigar_buffer_base = cigar_buffer;
  cigar_buffer->type = cigar_null; // Trick
  // Translate all match-CIGAR element
  match_alignment_region_rl_translate_cigar(
      &cigar_buffer,rl_match_length,
      rl_key_runs_acc,rl_text_runs_acc,
      region_key_begin,region_text_begin,left_gap_alignment);
  // Set CIGAR buffer used
  if (cigar_buffer->type!=cigar_null) ++(cigar_buffer);
  const uint64_t num_cigar_elements = cigar_buffer - cigar_buffer_base;
  vector_add_used(cigar_vector,num_cigar_elements);
  // Setup translated CIGAR
  match_alignment_region_set_cigar_length(match_alignment_region,num_cigar_elements);
  match_alignment_region_set_cigar_buffer_offset(match_alignment_region,cigar_offset);
}
