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

#ifndef MATCH_TRACE_H_
#define MATCH_TRACE_H_

#include "utils/essentials.h"
#include "align/pattern/pattern.h"
#include "text/text_trace.h"
#include "matches/matches_cigar.h"

/*
 * Matches Type
 */
typedef enum {
  match_type_regular,           // Regular Match
  match_type_local,             // Local Match (product of local alignment)
  match_type_extended           // Extended Match (product of extending PE)
} match_type;

/*
 * Match (Trace-Match)
 */
typedef struct {
  /* Type */
  match_type type;                   // Match type
  /* Location */
  char* sequence_name;               // Sequence name (After decoding.Eg Chr1)
  strand_t strand;                   // Mapping Strand
  bs_strand_t bs_strand;             // Bisulfite Strand
  uint64_t text_position;            // Position of the match in the text. Local text (Eg wrt Chr1)
  /* Reference-Text */
  text_trace_t* text_trace;          // Text-trace
  uint8_t* text;                     // Pointer to the matching-text
  uint64_t text_length;              // Length of the matching-text
  /* Distance/Score */
  uint64_t event_distance;           // Distance
  uint64_t edit_distance;            // Edit-Distance
  int32_t swg_score;                 // SWG Distance/Score
  float error_quality;               // Average errors quality
  uint8_t mapq_score;                // MAPQ Score
  /* Alignment */
  match_alignment_t match_alignment; // Match Alignment (CIGAR + ...)
  void* match_scaffold;              // Supporting Scaffolding
} match_trace_t;

/*
 * Accessors
 */
cigar_element_t* match_trace_get_cigar_buffer(
    const match_trace_t* const match_trace,
    vector_t* const cigar_vector);
uint64_t match_trace_get_cigar_length(
    const match_trace_t* const match_trace);
uint64_t match_trace_get_event_distance(
    const match_trace_t* const match_trace);

/*
 * Compare
 */
int match_trace_cmp_swg_score(
    const match_trace_t** const _a,
    const match_trace_t** const _b);
int match_trace_cmp_genomic_position(
    const match_trace_t** const _a,
    const match_trace_t** const _b);
int matche_trace_cigar_cmp(
    vector_t* const cigar_vector_match0,
    match_trace_t* const match0,
    vector_t* const cigar_vector_match1,
    match_trace_t* const match1);

#endif /* MATCH_TRACE_H_ */
