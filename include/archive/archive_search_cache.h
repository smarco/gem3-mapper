/*
 * PROJECT: GEMMapper
 * FILE: archive_search_cache.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef ARCHIVE_SEARCH_CACHE_H_
#define ARCHIVE_SEARCH_CACHE_H_

#include "utils/essentials.h"
#include "archive/archive_search.h"
#include "mapper/mapper.h"

/*
 * Archive Search Cache
 */
typedef struct {
  mapper_parameters_t* mapper_parameters; // Mapper parameters
  vector_t* archive_search_cache;         // Already allocated & configured Slab (archive_search_t*)
} archive_search_cache_t;

/*
 * Setup
 */
archive_search_cache_t* archive_search_cache_new(mapper_parameters_t* const mapper_parameters);
void archive_search_cache_delete(archive_search_cache_t* const archive_search_cache);

/*
 * Allocate/Free
 */
archive_search_t* archive_search_cache_alloc(archive_search_cache_t* const archive_search_cache);
void archive_search_cache_free(
    archive_search_cache_t* const archive_search_cache,
    archive_search_t* const archive_search);

#endif /* ARCHIVE_SEARCH_CACHE_H_ */
