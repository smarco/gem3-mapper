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
 *   Output module provides SAM output generation
 */

#ifndef OUTPUT_SAM_H_
#define OUTPUT_SAM_H_

#include "utils/essentials.h"
#include "io/buffered_output_file.h"
#include "archive/search/archive_search.h"
#include "text/sequence.h"
#include "matches/matches.h"
#include "matches/paired_matches.h"

/*
 * SAM Parameters
 */
typedef struct {
  /* Header & RG */
  char *read_group_header;
  string_t *read_group_id;
  /* Read & Qualities */
  bool omit_secondary_read__qualities;
  /* CIGAR */
  bool print_mismatches;
  /* XA */
  bool compact_xa;
  /* Bisulfite */
  bool bisulfite_output;
  /* GEM compatibility */
  bool print_gem_fields;
} output_sam_parameters_t;

/*
 * Setup
 */
void output_sam_parameters_set_defaults(
    output_sam_parameters_t* const sam_parameters);
void output_sam_parse_read_group_header(
    char* const read_group_buffer,
    output_sam_parameters_t* const sam_parameters);

/*
 * SAM Headers
 */
void output_sam_print_header(
    output_file_t* const output_file,
    archive_t* const archive,
    output_sam_parameters_t* const sam_parameters,
    int argc,
    char** argv,
    char* const gem_version);
			

/*
 * SAM output SE
 */
void output_sam_single_end_matches(
    buffered_output_file_t* const buffered_output_file,
    archive_search_t* const archive_search,
    matches_t* const matches,
    const output_sam_parameters_t* const output_sam_parameters);

/*
 * SAM output PE
 */
void output_sam_paired_end_matches(
    buffered_output_file_t* const buffered_output_file,
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches,
    const output_sam_parameters_t* const output_sam_parameters);

#endif /* OUTPUT_SAM_H_ */
