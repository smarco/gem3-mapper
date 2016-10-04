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
    archive_t* const archive,
    search_parameters_t* const search_parameters) {
  // Alloc
  archive_search_cache_t* const archive_search_cache = mm_alloc(archive_search_cache_t);
  // Initialize cache
  archive_search_cache->archive_search_cache = vector_new(ARCHIVE_SEARCH_CACHE_INIT_SIZE,archive_search_t*);
  archive_search_cache->archive = archive;
  archive_search_cache->search_parameters = search_parameters;
  // MM
  archive_search_cache->mm_slab = mm_slab_new_(BUFFER_SIZE_8M,BUFFER_SIZE_8M,MM_UNLIMITED_MEM,"archive_search_cache.8MB");
  archive_search_cache->mm_stack = mm_stack_new(archive_search_cache->mm_slab);
  // Return
  return archive_search_cache;
}
void archive_search_cache_delete(archive_search_cache_t* const archive_search_cache) {
  // Delete all archive_search_t objects in cache
  VECTOR_ITERATE(archive_search_cache->archive_search_cache,archive_search_ptr,n,archive_search_t*) {
    archive_search_destroy(*archive_search_ptr);
  }
  mm_stack_delete(archive_search_cache->mm_stack);
  mm_slab_delete(archive_search_cache->mm_slab);
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
    *archive_search = mm_stack_alloc(archive_search_cache->mm_stack,archive_search_t); // Allocate handler
    archive_search_init(
        *archive_search,archive_search_cache->archive,
        archive_search_cache->search_parameters,true,
        archive_search_cache->mm_stack); // Prepare
    archive_select_configure_se(*archive_search); // Select align
  }
}
void archive_search_cache_pe_alloc(
    archive_search_cache_t* const archive_search_cache,
    archive_search_t** const archive_search_end1,
    archive_search_t** const archive_search_end2) {
  // Allocate End/1
  if (vector_get_used(archive_search_cache->archive_search_cache) > 0) {
    *archive_search_end1 = *vector_get_last_elm(archive_search_cache->archive_search_cache,archive_search_t*);
    vector_dec_used(archive_search_cache->archive_search_cache);
  } else {
    *archive_search_end1 = mm_stack_alloc(archive_search_cache->mm_stack,archive_search_t); // Allocate handler
    archive_search_init(
        *archive_search_end1,archive_search_cache->archive,
        archive_search_cache->search_parameters,true,
        archive_search_cache->mm_stack); // Prepare
    archive_select_configure_pe(*archive_search_end1); // Select align
  }
  // Allocate End/2
  if (vector_get_used(archive_search_cache->archive_search_cache) > 0) {
    *archive_search_end2 = *vector_get_last_elm(archive_search_cache->archive_search_cache,archive_search_t*);
    vector_dec_used(archive_search_cache->archive_search_cache);
  } else {
    *archive_search_end2 = mm_stack_alloc(archive_search_cache->mm_stack,archive_search_t); // Allocate handler
    archive_search_init(
        *archive_search_end2,archive_search_cache->archive,
        archive_search_cache->search_parameters,true,
        archive_search_cache->mm_stack); // Prepare
    archive_select_configure_pe(*archive_search_end2); // Select align
  }
}
void archive_search_cache_free(
    archive_search_cache_t* const archive_search_cache,
    archive_search_t* const archive_search) {
  // Add it to the cache
  vector_insert(archive_search_cache->archive_search_cache,archive_search,archive_search_t*);
}
