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
 *   Archive-Search main module & support data structures
 */

#ifndef ARCHIVE_SEARCH_H_
#define ARCHIVE_SEARCH_H_

#include "utils/essentials.h"
#include "archive/archive.h"
#include "approximate_search/approximate_search.h"
#include "text/sequence.h"
#include "matches/classify/matches_classify.h"

/*
 * Archive Search State
 */
typedef enum {
  archive_search_pe_state_begin = 0,                     // Beginning of the search
  archive_search_pe_state_search_end1 = 1,               // Generate candidates for end1
  archive_search_pe_state_search_end2 = 2,               // Generate candidates for end2
  archive_search_pe_state_recovery = 3,                  // Recover by extension when needed
  archive_search_pe_state_find_pairs = 4,                // Cross-link matches from both ends
  archive_search_pe_state_end = 5                        // End of the current workflow
} archive_search_pe_state_t;
extern const char* archive_search_pe_state_label[7];

/*
 * Archive Search
 */
typedef struct {
  /* Archive */
  archive_t* archive;                        // Archive
  /* Archive Paired-End-Search (Only end/1 used in PE search) */
  archive_search_pe_state_t pe_search_state; // Search State
  bool pair_searched;                        // Paired search performed
  bool pair_extended;                        // Paired extension performed
  bool pair_extended_shortcut;               // Paired extension performed (to shortcut)
  /* Parameters */
  search_parameters_t search_parameters;     // Search parameters
  sequence_t sequence;                       // Input Sequence
  bool buffered_search;                      // Buffered Search
  /* Parameters-BS */
  string_t* bs_original_sequence;            // Bisulfite original sequence before conversion
  sequence_end_t bs_sequence_end;            // Bisulfite sequence end
  /* Approximate Search */
  approximate_search_t approximate_search;   // Approximate Search State
  /* Stats */
  mapper_stats_t* mapper_stats;              // Mapping statistics
  /* MM */
  mm_stack_t* mm_stack;                      // MM-Stack
} archive_search_t;

/*
 * Setup
 */
void archive_search_init(
    archive_search_t* const archive_search,
    archive_t* const archive,
    search_parameters_t* const search_parameters,
    const bool buffered_search,
    mm_stack_t* const mm_stack);
void archive_search_destroy(archive_search_t* const archive_search);

void archive_search_se_new(
    archive_t* const archive,
    search_parameters_t* const search_parameters,
    const bool buffered_search,
    mm_stack_t* const mm_stack,
    archive_search_t** const archive_search);
void archive_search_pe_new(
    archive_t* const archive,
    search_parameters_t* const search_parameters,
    const bool buffered_search,
    mm_stack_t* const mm_stack,
    archive_search_t** const archive_search_end1,
    archive_search_t** const archive_search_end2);
void archive_search_reset(archive_search_t* const archive_search);
void archive_search_delete(archive_search_t* const archive_search);

/*
 * Accessors
 */
sequence_t* archive_search_get_sequence(const archive_search_t* const archive_search);
bool archive_search_finished(const archive_search_t* const archive_search);

uint64_t archive_search_get_num_regions_profile(const archive_search_t* const archive_search);
uint64_t archive_search_get_num_decode_candidates(const archive_search_t* const archive_search);
uint64_t archive_search_get_num_verify_candidates(const archive_search_t* const archive_search);

/*
 * Errors
 */
#define GEM_ERROR_ARCHIVE_SEARCH_INDEX_COMPLEMENT_REQUIRED "Archive Search. Explicit indexed complement required"

#endif /* ARCHIVE_SEARCH_H_ */
