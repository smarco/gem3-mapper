/*
 * PROJECT: GEMMapper
 * FILE: archive_search_cache.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive/archive_search_cache.h"
#include "archive/archive_search_se.h"
#include "archive/archive_search_pe.h"
#include "archive/archive_select.h"

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
  archive_search_cache->mm_stack = mm_stack_new(mm_pool_get_slab(mm_pool_2MB));
  // Return
  return archive_search_cache;
}
void archive_search_cache_delete(archive_search_cache_t* const archive_search_cache) {
  // Delete all archive_search_t objects in cache
  VECTOR_ITERATE(archive_search_cache->archive_search_cache,archive_search_ptr,n,archive_search_t*) {
    archive_search_destroy(*archive_search_ptr);
  }
  mm_stack_delete(archive_search_cache->mm_stack);
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
