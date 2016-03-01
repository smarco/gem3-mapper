/*
 * PROJECT: GEMMapper
 * FILE: archive_search_cache.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive/archive_search_cache.h"
#include "archive/archive_search_se.h"
#include "archive/archive_search_pe.h"

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
  // Return
  return archive_search_cache;
}
void archive_search_cache_delete(archive_search_cache_t* const archive_search_cache) {
  // Delete all archive_search_t objects in cache
  VECTOR_ITERATE(archive_search_cache->archive_search_cache,archive_search_ptr,n,archive_search_t*) {
    archive_search_delete(*archive_search_ptr);
  }
  // Free handlers
  vector_delete(archive_search_cache->archive_search_cache);
  mm_free(archive_search_cache);
}
/*
 * Allocate/Free
 */
archive_search_t* archive_search_cache_alloc(archive_search_cache_t* const archive_search_cache) {
  archive_search_t* archive_search = NULL;
  if (vector_get_used(archive_search_cache->archive_search_cache)>0) {
    // Get from cache already prepared archive_search_t
    archive_search = *vector_get_last_elm(archive_search_cache->archive_search_cache,archive_search_t*);
    vector_dec_used(archive_search_cache->archive_search_cache);
  } else {
    // Allocate new one
    archive_search_se_new(archive_search_cache->archive,archive_search_cache->search_parameters,&archive_search);
  }
  // Return
  return archive_search;
}
void archive_search_cache_free(
    archive_search_cache_t* const archive_search_cache,
    archive_search_t* const archive_search) {
  // Add it to the cache
  vector_insert(archive_search_cache->archive_search_cache,archive_search,archive_search_t*);
}
