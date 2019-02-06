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
 *   Archive-Search handlers module encapsulates all required
 *   objects & memory-dispatchers an archive-search needs to operate.
 *   Once a new archive-search is created, all dependencies
 *   and handlers are injected into it
 */

#ifndef ARCHIVE_SEARCH_HANDLERS_H_
#define ARCHIVE_SEARCH_HANDLERS_H_

#include "utils/essentials.h"
#include "archive/search/archive_search.h"
#include "filtering/candidates/filtering_candidates.h"
#include "neighborhood_search/nsearch_schedule.h"
#include "mapper/mapper_stats.h"

/*
 * Archive Search Handlers
 */
typedef struct {
  /* Archive */
  archive_t* archive;
  /* Filtering candidates */
  filtering_candidates_t filtering_candidates_end1;
  filtering_candidates_t filtering_candidates_end2;
  /* Neighborhood Search */
  nsearch_schedule_t nsearch_schedule;
  /* Stats */
  mapper_stats_t* mapper_stats;         // Mapping Statistics
  /* MM */
  mm_slab_t* mm_slab;                   // MM-Slab
  mm_allocator_t* mm_allocator;         // MM-Allocator
} archive_search_handlers_t;

/*
 * Archive Search Handlers
 */
archive_search_handlers_t* archive_search_handlers_new(archive_t* const archive);
void archive_search_handlers_clear(archive_search_handlers_t* const handlers);
void archive_search_handlers_delete(archive_search_handlers_t* const handlers);

/*
 * Archive Search Handlers Injection (Support Data Structures)
 */
void archive_search_handlers_prepare_se(
    archive_search_t* const archive_search,
    sequence_t* const sequence,
    bisulfite_conversion_t const bisulfite_conversion,
    archive_search_handlers_t* const handlers);
void archive_search_handlers_prepare_pe(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    sequence_t* const sequence_end1,
    sequence_t* const sequence_end2,
    bisulfite_conversion_t const bisulfite_conversion_end1,
    bisulfite_conversion_t const bisulfite_conversion_end2,
    archive_search_handlers_t* const handlers);

#endif /* ARCHIVE_SEARCH_HANDLERS_H_ */
