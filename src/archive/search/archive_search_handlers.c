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

#include "archive/search/archive_search_handlers.h"

/*
 * Archive Search Handlers
 */
archive_search_handlers_t* archive_search_handlers_new(archive_t* const archive) {
  // Allocate
  archive_search_handlers_t* const handlers = mm_alloc(archive_search_handlers_t);
  // Archive
  handlers->archive = archive;
  // Filtering candidates
  filtering_candidates_init(&handlers->filtering_candidates_end1);
  filtering_candidates_init(&handlers->filtering_candidates_end2);
  // Stats
  handlers->mapper_stats = mapper_stats_new();
  // MM
  handlers->mm_slab = mm_slab_new_(BUFFER_SIZE_8M,BUFFER_SIZE_8M,MM_UNLIMITED_MEM);
  handlers->mm_allocator = mm_allocator_new(handlers->mm_slab);
  // Return
  return handlers;
}
void archive_search_handlers_clear(archive_search_handlers_t* const handlers) {
  filtering_candidates_clear(&handlers->filtering_candidates_end1,false);
  filtering_candidates_clear(&handlers->filtering_candidates_end2,false);
  mm_allocator_clear(handlers->mm_allocator);
}
void archive_search_handlers_delete(archive_search_handlers_t* const handlers) {
  filtering_candidates_destroy(&handlers->filtering_candidates_end1,false);
  filtering_candidates_destroy(&handlers->filtering_candidates_end2,false);
  mapper_stats_delete(handlers->mapper_stats);
  mm_allocator_delete(handlers->mm_allocator);
  //mm_slab_delete(handlers->mm_slab);
  mm_free(handlers);
}
/*
 * Archive Search Handlers Injection (Support Data Structures)
 */
void archive_search_handlers_prepare_se(
    archive_search_t* const archive_search,
    sequence_t* const sequence,
    archive_search_handlers_t* const handlers) {
  // Inject Handlers
  archive_search_inject_handlers(
      archive_search,handlers->archive,
      &handlers->filtering_candidates_end1,
      &handlers->nsearch_schedule,
      handlers->mapper_stats,handlers->mm_allocator);
  // Prepare sequence
  archive_search_prepare_sequence(archive_search,sequence);
}
void archive_search_handlers_prepare_pe(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    sequence_t* const sequence_end1,
    sequence_t* const sequence_end2,
    archive_search_handlers_t* const handlers) {
  // Inject Handlers
  archive_search_inject_handlers(
      archive_search_end1,handlers->archive,
      &handlers->filtering_candidates_end1,
      &handlers->nsearch_schedule,
      handlers->mapper_stats,handlers->mm_allocator);
  archive_search_inject_handlers(
      archive_search_end2,handlers->archive,
      &handlers->filtering_candidates_end2,
      &handlers->nsearch_schedule,
      NULL,handlers->mm_allocator);
  // Prepare sequences
  archive_search_prepare_sequence(archive_search_end1,sequence_end1);
  archive_search_prepare_sequence(archive_search_end2,sequence_end2);
}
