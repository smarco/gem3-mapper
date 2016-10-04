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
archive_search_handlers_t* archive_search_handlers_new(void) {
  // Allocate
  archive_search_handlers_t* const handlers = mm_alloc(archive_search_handlers_t);
  // Filtering candidates
  filtering_candidates_mm_init(&handlers->filtering_candidates_mm);
  filtering_candidates_init(&handlers->filtering_candidates_end1);
  filtering_candidates_init(&handlers->filtering_candidates_end2);
  // Stats
  handlers->mapper_stats = mapper_stats_new();
  // MM
  handlers->mm_slab = mm_slab_new_(BUFFER_SIZE_8M,BUFFER_SIZE_8M,MM_UNLIMITED_MEM,"archive_search_handlers.8MB");
  handlers->mm_stack = mm_stack_new(handlers->mm_slab);
  // Return
  return handlers;
}
void archive_search_handlers_clear(archive_search_handlers_t* const handlers) {
  filtering_candidates_clear(&handlers->filtering_candidates_end1);
  filtering_candidates_clear(&handlers->filtering_candidates_end2);
  filtering_candidates_mm_clear(&handlers->filtering_candidates_mm);
  mm_stack_clear(handlers->mm_stack);
}
void archive_search_handlers_delete(archive_search_handlers_t* const handlers) {
  filtering_candidates_destroy(&handlers->filtering_candidates_end1);
  filtering_candidates_destroy(&handlers->filtering_candidates_end2);
  filtering_candidates_mm_destroy(&handlers->filtering_candidates_mm);
  mapper_stats_delete(handlers->mapper_stats);
  mm_stack_delete(handlers->mm_stack);
  mm_slab_delete(handlers->mm_slab);
  mm_free(handlers);
}
/*
 * Handlers Injection (Support Data Structures)
 */
void archive_search_handlers_inject_se(
    archive_search_t* const archive_search,
    archive_search_handlers_t* const handlers) {
  // Archive Search
  archive_search->mapper_stats = handlers->mapper_stats;
  archive_search->mm_stack = handlers->mm_stack;
  // Approximate Search
  approximate_search_inject_handlers(
      &archive_search->approximate_search,archive_search->archive,
      &archive_search->search_parameters,&handlers->filtering_candidates_end1,
      &handlers->filtering_candidates_mm,NULL,&handlers->nsearch_schedule,
      handlers->mm_stack,handlers->mm_stack);
}
void archive_search_handlers_inject_pe(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    archive_search_handlers_t* const handlers) {
  // Archive Search
  archive_search_end1->mapper_stats = handlers->mapper_stats;
  archive_search_end2->mapper_stats = NULL;
  archive_search_end1->mm_stack = handlers->mm_stack;
  archive_search_end2->mm_stack = handlers->mm_stack;
  // Approximate Search
  approximate_search_inject_handlers(
      &archive_search_end1->approximate_search,archive_search_end1->archive,
      &archive_search_end1->search_parameters,&handlers->filtering_candidates_end1,
      &handlers->filtering_candidates_mm,NULL,&handlers->nsearch_schedule,
      handlers->mm_stack,handlers->mm_stack);
  approximate_search_inject_handlers(
      &archive_search_end2->approximate_search,archive_search_end2->archive,
      &archive_search_end2->search_parameters,&handlers->filtering_candidates_end2,
      &handlers->filtering_candidates_mm,NULL,&handlers->nsearch_schedule,
      handlers->mm_stack,handlers->mm_stack);
}
