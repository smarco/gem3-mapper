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
#include "matches/align/match_align_normalize.h"
#include "archive/archive_text_rl.h"

/*
 * Normalize Alignment
 */
bool match_align_normalize_cigar_trim(
    const cigar_element_t* const cigar_element,
    const uint64_t min_context_length,
    uint64_t* const trim_length,
    uint64_t* const match_position) {
  switch (cigar_element->type) {
    case cigar_match:
      // Small Match
      if (cigar_element->length < min_context_length) {
        *trim_length += cigar_element->length;
        if (match_position!=NULL) *match_position += cigar_element->length;
        return true;
      }
    break;
    case cigar_mismatch:
      ++(*trim_length);
      if (match_position!=NULL) ++(*match_position);
      return true;
      break;
    case cigar_del:
      *trim_length += cigar_element->length;
      return true;
      break;
    case cigar_ins:
      if (match_position!=NULL) *match_position += cigar_element->length;
      return true;
      break;
    default:
      return false;
      break;
  }
  return false;
}
void match_align_normalize_cigar(
    match_trace_t* const match_trace,
    search_parameters_t* const search_parameters,
    const bool left_gap_alignment,
    vector_t* const cigar_vector) {
  // Parameters
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  const uint64_t min_context_length = search_parameters->cigar_curation_min_end_context_nominal;
  const uint64_t cigar_length = match_alignment->cigar_length;
  cigar_element_t* const cigar_buffer =
      vector_get_elm(cigar_vector,match_alignment->cigar_offset,cigar_element_t);
  uint64_t normalized_cigar_length = 0, indel_length = 0, i = 0, j = 0;
  // Trim the beginning of the read
  while (j < cigar_length) {
    cigar_element_t* const cigar_element = cigar_buffer + j;
    if (!match_align_normalize_cigar_trim(
        cigar_element,min_context_length,
        &indel_length,&match_alignment->match_position)) {
      break;
    }
    // Trim
    cigar_element->type = cigar_del;
    cigar_element->attributes = cigar_attr_trim;
    cigar_element->length = indel_length;
    i = j++;
  }
  // Traverse all CIGAR elements
  while (i < cigar_length) {
    cigar_element_t* cigar_element = cigar_buffer + i;
    switch (cigar_element->type) {
      case cigar_mismatch:
        cigar_buffer[normalized_cigar_length++] = *cigar_element;
        ++i;
        break;
      case cigar_match:
        if (i+1<cigar_length && cigar_buffer[i+1].type == cigar_match) {
          cigar_buffer[i+1].length += cigar_element->length;
        } else {
          cigar_buffer[normalized_cigar_length++] = *cigar_element;
        }
        ++i;
        break;
      case cigar_del:
      case cigar_ins: {
        cigar_element_t accum_del = { .type = cigar_del, .attributes = cigar_attr_none, .length = 0 };
        cigar_element_t accum_ins = { .type = cigar_ins, .attributes = cigar_attr_none, .length = 0 };
        if (cigar_element->type==cigar_del && cigar_element->attributes==cigar_attr_trim) {
          accum_del.attributes = cigar_attr_trim;
        }
        // Compact all deletions/insertions in a row
        while (i<cigar_length) {
          cigar_element = cigar_buffer + i;
          if (cigar_element->type == cigar_del) {
            accum_del.length += cigar_element->length;
            ++i;
          } else if (cigar_element->type == cigar_ins) {
            accum_ins.length += cigar_element->length;
            ++i;
          } else {
            break;
          }
        }
        // Copy normalized Ins/Del
//        if (accum_del.length==1 && accum_ins.length==1) {
//          cigar_buffer[normalized_cigar_length].type = cigar_mismatch;
//          cigar_buffer[normalized_cigar_length].mismatch = ;
//        } else {
          if (left_gap_alignment) {
            if (accum_del.length > 0) cigar_buffer[normalized_cigar_length++] = accum_del;
            if (accum_ins.length > 0) cigar_buffer[normalized_cigar_length++] = accum_ins;
          } else {
            if (accum_ins.length > 0) cigar_buffer[normalized_cigar_length++] = accum_ins;
            if (accum_del.length > 0) cigar_buffer[normalized_cigar_length++] = accum_del;
          }
//        }
        break;
      }
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
  // Handle last CIGAR element(s)
  if (normalized_cigar_length > 0) {
    // Check match or deletion at the end
    cigar_element_t* last_cigar_element = cigar_buffer + (normalized_cigar_length-1);
    indel_length = 0;
    if (match_align_normalize_cigar_trim(last_cigar_element,min_context_length,&indel_length,NULL)) {
      // Chain all the mismatches & deletions at the end
      while (last_cigar_element > cigar_buffer) {
        --last_cigar_element;
        if (!match_align_normalize_cigar_trim(last_cigar_element,min_context_length,&indel_length,NULL)) break;
        --normalized_cigar_length;
      }
      // Merge all of them
      last_cigar_element = cigar_buffer + (normalized_cigar_length-1);
      last_cigar_element->type = cigar_del;
      last_cigar_element->attributes = cigar_attr_trim;
      last_cigar_element->length = indel_length;
      if (indel_length==0) --normalized_cigar_length;
    }
  }
  // Set normalized-CIGAR length
  match_alignment->cigar_length = normalized_cigar_length;
}
/*
 * SWG Normalize CIGAR & Adjust Position (Translate RL-CIGAR if required)
 */
void match_align_normalize(
    matches_t* const matches,
    match_trace_t* const match_trace,
    search_parameters_t* const search_parameters) {
  // Parameters
  text_trace_t* const text_trace = match_trace->text_trace;
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  vector_t* const cigar_vector = matches->cigar_vector;
  //  // DEBUG
  //  match_align_normalize_cigar(match_trace,cigar_vector,align_parameters);
  //  match_alignment->match_text_offset = match_alignment->match_position - align_input->text_position;
  //  //  fprintf(stderr,"[GEM]>RL-Read\n");
  //  //  dna_buffer_print(stderr,align_input->key,align_input->key_length,false);
  //  //  fprintf(stderr,"[GEM]>RL-Text\n");
  //  //  dna_buffer_print(stderr,align_input->text,align_input->text_length,false);
  //  fprintf(stderr,"[GEM]>RL-Alignment\n");
  //  match_alignment_print_pretty(stderr,&match_trace->match_alignment,
  //      matches->cigar_vector,align_input->key,align_input->key_length,
  //      align_input->text+match_trace->match_alignment.match_text_offset,
  //      align_input->text_length,mm_allocator);
  // Normalize alignment
  const bool left_gap_alignment = (text_trace->strand == Forward);
  match_align_normalize_cigar(match_trace,search_parameters,left_gap_alignment,cigar_vector);
  // Adjust
  match_alignment->match_text_offset = match_alignment->match_position - text_trace->position;
  match_trace->text = text_trace->text + match_alignment->match_text_offset;
}
