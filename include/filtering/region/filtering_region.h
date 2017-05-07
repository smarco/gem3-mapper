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
 *   Filtering region data structure holds all the required information
 *   of a candidate region of the index that can potentially align against
 *   the search pattern
 */

#ifndef FILTERING_REGION_H_
#define FILTERING_REGION_H_

#include "utils/essentials.h"
#include "text/sequence.h"
#include "align/pattern/pattern.h"
#include "archive/archive_text.h"
#include "matches/scaffold/match_scaffold.h"

/*
 * Filtering Region Status
 */
typedef enum {
  filtering_region_pending=0,               // On-hold state (unverified)
  filtering_region_unverified=1,            // Region still unverified
  filtering_region_verified_discarded=2,    // Region verified but too distant
  filtering_region_accepted=3,              // Region accepted by pre-filters
  filtering_region_accepted_subdominant=4   // Region accepted by pre-filters but sub-dominant (by ranking)
} filtering_region_status_t;
extern const char* filtering_region_status_label[5];
/*
 * Filtering Region
 */
typedef struct {
  /* State */
  filtering_region_status_t status;         // Filtering Region Status
  /* Source Region Offset */
  uint64_t text_source_region_offset;       // Text-Offset to the begin of the source-region
  uint64_t key_source_region_offset;        // Key-Offset to the begin of the source-region
  /* Text */
  text_trace_t text_trace;                  // Text-trace
  uint64_t text_begin_position;             // Region effective begin position (adjusted to error boundaries)
  uint64_t text_end_position;               // Region effective end position (adjusted to error boundaries)
  /* Key */
  uint64_t key_trimmed_length;              // Key trimmed length
  uint64_t key_trim_left;                   // Key left trim  (due to position correction wrt the reference limits)
  uint64_t key_trim_right;                  // Key right trim (due to position correction wrt the reference limits)
  bool key_trimmed;                         // Key has been trimmed (due to dangling ends; projected to the text-candidate)
  /* Alignment */
  alignment_t alignment;                    // Filtering Region Alignment
  match_scaffold_t match_scaffold;          // Alignment-regions supporting the filtering region (Scaffolding)
} filtering_region_t;

/*
 * Filtering Region Auxiliary Structures
 */
typedef struct {
  /* Location */
  uint64_t begin_position;
  uint64_t end_position;
} verified_region_t;

/*
 * Retrieve filtering region text-candidate
 */
void filtering_region_retrieve_text(
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    archive_text_t* const archive_text,
    mm_allocator_t* const mm_allocator);

/*
 * Filtering Region Key Trims
 */
void filtering_region_compute_key_trims(
    filtering_region_t* const filtering_region,
    pattern_t* const pattern);

/*
 * Display
 */
void filtering_region_print(
    FILE* const stream,
    filtering_region_t* const region,
    const bool print_region_text,
    const bool print_alignment_regions,
    const bool print_alignment);

#endif /* FILTERING_REGION_H_ */
