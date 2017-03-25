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
 *   Archive-search cache module provides bulk allocation of multiple
 *   archive-search objects, bulk initialization & re-use of them
 *   (avoiding allocations/initialization overheads)
 */

#include "archive/search/archive_search_cache.h"
#include "archive/search/archive_search_se.h"
#include "archive/search/archive_search_pe.h"
#include "archive/search/archive_select.h"

/*
 * Constants
 */
#define ARCHIVE_SEARCH_CACHE_INIT_SIZE 1000

/*
 * Setup
 */
archive_search_cache_t* archive_search_cache_new(
    search_parameters_t* const search_parameters) {
  // Alloc
  archive_search_cache_t* const archive_search_cache = mm_alloc(archive_search_cache_t);
  // Initialize cache
  archive_search_cache->archive_search_cache = vector_new(ARCHIVE_SEARCH_CACHE_INIT_SIZE,archive_search_t*);
  // Search Parameters
  archive_search_cache->search_parameters = search_parameters;
  // Return
  return archive_search_cache;
}
void archive_search_cache_delete(archive_search_cache_t* const archive_search_cache) {
  // Delete all archive_search_t objects in cache
  VECTOR_ITERATE(archive_search_cache->archive_search_cache,archive_search,n,archive_search_t*) {
    archive_search_delete(*archive_search);
  }
  // Free handlers
  vector_delete(archive_search_cache->archive_search_cache);
  mm_free(archive_search_cache);
}
/*
 * Allocate/Free
 */
void archive_search_cache_se_alloc(
    archive_search_cache_t* const archive_search_cache,
    archive_search_t** const archive_search) {
  if (vector_get_used(archive_search_cache->archive_search_cache) > 0) {
    // Get from the cache an already prepared archive-search
    *archive_search = *vector_get_last_elm(archive_search_cache->archive_search_cache,archive_search_t*);
    vector_dec_used(archive_search_cache->archive_search_cache);
  } else {
    // Allocate a new archive-search
    archive_search_se_new(archive_search_cache->search_parameters,true,archive_search);
  }
}
void archive_search_cache_pe_alloc(
    archive_search_cache_t* const archive_search_cache,
    archive_search_t** const archive_search_end1,
    archive_search_t** const archive_search_end2) {
  // Allocate
  if (vector_get_used(archive_search_cache->archive_search_cache) > 1) {
    *archive_search_end1 = *vector_get_last_elm(archive_search_cache->archive_search_cache,archive_search_t*);
    vector_dec_used(archive_search_cache->archive_search_cache);
    *archive_search_end2 = *vector_get_last_elm(archive_search_cache->archive_search_cache,archive_search_t*);
    vector_dec_used(archive_search_cache->archive_search_cache);
  } else {
    archive_search_pe_new(archive_search_cache->search_parameters,true,archive_search_end1,archive_search_end2);
  }
}
void archive_search_cache_free(
    archive_search_cache_t* const archive_search_cache,
    archive_search_t* const archive_search) {
  // Return it to the cache
  vector_insert(archive_search_cache->archive_search_cache,archive_search,archive_search_t*);
}
