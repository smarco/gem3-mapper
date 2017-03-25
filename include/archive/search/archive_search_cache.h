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

#ifndef ARCHIVE_SEARCH_CACHE_H_
#define ARCHIVE_SEARCH_CACHE_H_

#include "utils/essentials.h"
#include "archive/search/archive_search.h"

/*
 * Archive Search Cache
 */
typedef struct {
  /* Archive-Search Cache */
  vector_t* archive_search_cache;         // Already allocated & configured Slab (archive_search_t*)
  /* Search Parameters */
  search_parameters_t* search_parameters;
} archive_search_cache_t;

/*
 * Setup
 */
archive_search_cache_t* archive_search_cache_new(search_parameters_t* const search_parameters);
void archive_search_cache_delete(archive_search_cache_t* const archive_search_cache);

/*
 * Allocate/Free
 */
void archive_search_cache_se_alloc(
    archive_search_cache_t* const archive_search_cache,
    archive_search_t** const archive_search);
void archive_search_cache_pe_alloc(
    archive_search_cache_t* const archive_search_cache,
    archive_search_t** const archive_search_end1,
    archive_search_t** const archive_search_end2);
void archive_search_cache_free(
    archive_search_cache_t* const archive_search_cache,
    archive_search_t* const archive_search);

#endif /* ARCHIVE_SEARCH_CACHE_H_ */
