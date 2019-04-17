/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *  Copyright (c) 2011-2017 by Simon Heath  <simon.heath@gmail.com>
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
 *            Simon Heath <simon.heath@gmail.com>
 * DESCRIPTION:
 *   Mapper module encapsulates and provides accessors to all
 *   the parameters used by the mapper
 */

#ifndef RESTRICTION_TEXT_H_
#define RESTRICTION_TEXT_H_

#include "utils/essentials.h"
#include "matches/match_trace.h"
#include "archive/archive.h"
#include "align/pattern/pattern.h"
#include "text/restriction_locate.h"

typedef struct {
  string_t restriction_site;
  uint64_t search_pattern[3][DNA_EXT_RANGE];  // For searching in read
  uint64_t reference_pattern[DNA_EXT_RANGE];  // For searching in reference
  uint64_t search_pattern_mask;
  uint64_t reference_pattern_mask;
  uint64_t search_pattern_len;
  uint64_t reference_pattern_len;
  char* restriction_enzyme;
  uint64_t cut_site_index;
  uint64_t search_pattern_split_offset;   // Offset from end of search pattern to first character after split point
  bool interaction_site;
} restriction_site_t;

typedef struct {
	restriction_site_locator_t * restriction_site_locator;
	archive_t* archive;
	vector_t* restriction_sites;
	pthread_mutex_t mutex;
	pthread_t restriction_site_writer;
	uint64_t num_tags;
	uint64_t next_tag_idx;
	pthread_t *threads;
	FILE *output_handle;
} restriction_site_locator_builder_t;

restriction_site_t *restriction_new(char * const);
void restriction_site_delete(restriction_site_t * const);
void restriction_create_interaction_junctions(vector_t * const restriction_sites);
void find_restriction_site_matches(
		pattern_t const * const pattern,
		vector_t const * const restriction_sites,
		vector_t * const restriction_hits,
		bisulfite_conversion_t const bisulfite_conversion);
restriction_site_locator_builder_t * const restriction_text_init_locator(
		archive_t* const archive,
		vector_t* const restriction_sites,
		restriction_site_locator_t * const restriction_site_locator,
		const char * const output_sites_name,
		const bool verbose_user,
		const uint64_t num_threads);

void restriction_match_trace_locate(
		match_trace_t * const match_trace,
		restriction_site_locator_t const * const restriction_site_locator);



#endif
