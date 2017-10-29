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

#include "matches/match_trace.h"

/*
 * Accessors
 */
cigar_element_t* match_trace_get_cigar_buffer(
    const match_trace_t* const match_trace,
    vector_t* const cigar_vector) {
  return vector_get_elm(cigar_vector,match_trace->match_alignment.cigar_offset,cigar_element_t);
}
uint64_t match_trace_get_cigar_length(const match_trace_t* const match_trace) {
  return match_trace->match_alignment.cigar_length;
}
uint64_t match_trace_get_event_distance(const match_trace_t* const match_trace) {
  return match_trace->event_distance;
}
/*
 * Compare
 */
int match_trace_cmp_swg_score(const match_trace_t** const _a,const match_trace_t** const _b) {
  const match_trace_t* const a = *_a;
  const match_trace_t* const b = *_b;
  const int distance_swg = (int)b->swg_score - (int)a->swg_score;
  if (distance_swg) return distance_swg;
  const int distance_event = (int)a->event_distance - (int)b->event_distance;
  if (distance_event) return distance_event;
  const int distance_edit = (int)a->edit_distance - (int)b->edit_distance;
  if (distance_edit) return distance_edit;
  if (a->error_quality > b->error_quality) return -1;
  if (a->error_quality < b->error_quality) return  1;
  // Untie using position (helps to stabilize & cmp results)
  if (a->sequence_name < b->sequence_name) return -1;
  if (a->sequence_name > b->sequence_name) return  1;
  if (a->text_position < b->text_position) return -1;
  if (a->text_position > b->text_position) return  1;
  return 0;
}
int match_trace_cmp_genomic_position(const match_trace_t** const _a,const match_trace_t** const _b) {
  const match_trace_t* const a = *_a;
  const match_trace_t* const b = *_b;
  const int cmp_name = gem_strcmp(a->sequence_name,b->sequence_name);
  if (cmp_name!=0) return cmp_name;
  if (a->text_position < b->text_position) return -1;
  if (a->text_position > b->text_position) return  1;
  return 0;
}
int matche_trace_cigar_cmp(
    vector_t* const cigar_vector_match0,
    match_trace_t* const match0,
    vector_t* const cigar_vector_match1,
    match_trace_t* const match1) {
  return matches_cigar_cmp(
      cigar_vector_match0,match0->match_alignment.cigar_offset,match0->match_alignment.cigar_length,
      cigar_vector_match1,match1->match_alignment.cigar_offset,match1->match_alignment.cigar_length);
}
