/*
 * PROJECT: GEMMapper
 * FILE: archive_search_group.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive_search_group.h"
#include "bpm_align_gpu.h"

/*
 * Constants
 */
#define MAPPER_CUDA_ARCHIVE_SEARCH_CACHE_INIT_SIZE 1000

/*
 * Archive Search Cache
 */
typedef struct {
  mapper_parameters_t* mapper_parameters; // Mapper parameters
  vector_t* archive_search_cache;         // Already allocated & configured Slab (archive_search_t*)
} archive_search_cache_t;
/*
 * Archive-Search groups
 */
typedef struct {
  /* BPM-GPU candidates buffer */
  bpm_gpu_buffer_t* bpm_gpu_buffer;      // BPM-Buffer
  /* Archive searches */
  vector_t* archive_searches;            // Vector of archive-searches (archive_search_t*)
} archive_search_group_buffer_t;
struct _archive_search_group_t {
  /* Configuration */
  uint64_t hint_patterns_per_search;                   // Number of patterns per search-group
  /* Archive-search group buffers */
  archive_search_group_buffer_t* search_group_buffers; // Search-group buffers
  uint64_t total_search_groups;                        // Number of search-groups
  /* Archive-search group iterator */
  uint64_t current_search_group_idx;                   // Current search-group index
  archive_search_group_buffer_t* current_search_group; // Current search-group
  uint64_t current_archive_search_idx;                 // Current archive-search index
  archive_search_t** current_archive_search;           // Current archive-search
  uint64_t total_archive_searches;                     // Total archive-searches
  /* Archive-search cache */
  archive_search_cache_t* archive_search_cache;        // Archive search cache
  /* MM */
  mm_search_t* mm_search;
};
/*
 * Archive Search Cache
 */
GEM_INLINE archive_search_cache_t* archive_search_cache_new(mapper_parameters_t* const mapper_parameters) {
  // Alloc
  archive_search_cache_t* const archive_search_cache = mm_alloc(archive_search_cache_t);
  // Initialize cache
  archive_search_cache->mapper_parameters = mapper_parameters;
  archive_search_cache->archive_search_cache = vector_new(MAPPER_CUDA_ARCHIVE_SEARCH_CACHE_INIT_SIZE,archive_search_t*);
  // Return
  return archive_search_cache;
}
GEM_INLINE void archive_search_cache_delete(archive_search_cache_t* const archive_search_cache) {
  // Delete all archive_search_t objects in cache
  VECTOR_ITERATE(archive_search_cache->archive_search_cache,archive_search_ptr,n,archive_search_t*) {
    archive_search_delete(*archive_search_ptr);
  }
  // Free handlers
  vector_delete(archive_search_cache->archive_search_cache);
  mm_free(archive_search_cache);
}
GEM_INLINE archive_search_t* archive_search_cache_alloc(
    archive_search_cache_t* const archive_search_cache,mm_search_t* const mm_search) {
  archive_search_t* archive_search = NULL;
  if (vector_get_used(archive_search_cache->archive_search_cache)>0) {
    // Get from cache already prepared archive_search_t
    archive_search = *vector_get_last_elm(archive_search_cache->archive_search_cache,archive_search_t*);
    vector_dec_used(archive_search_cache->archive_search_cache);
  } else {
    // Allocate new one
    archive_search = archive_search_new(
        archive_search_cache->mapper_parameters->archive,
        &archive_search_cache->mapper_parameters->search_parameters,
        &archive_search_cache->mapper_parameters->select_parameters);
  }
  // Init archive search
  archive_search_configure(archive_search,mm_search);
  // Return
  return archive_search;
}
GEM_INLINE void archive_search_cache_free(
    archive_search_cache_t* const archive_search_cache,archive_search_t* const archive_search) {
  // Add it to the cache
  vector_insert(archive_search_cache->archive_search_cache,archive_search,archive_search_t*);
}
/*
 * Archive-search group
 */
GEM_INLINE archive_search_group_t* archive_search_group_new(
    mapper_parameters_t* const mapper_parameters,bpm_gpu_buffer_t* const bpm_gpu_buffers,const uint64_t total_search_groups) {
  // Alloc
  archive_search_group_t* archive_search_group = mm_alloc(archive_search_group_t);
  // Initialize
  archive_search_group->total_search_groups = total_search_groups;
  const uint64_t hint_patterns_per_search = (archive_is_indexed_complement(mapper_parameters->archive)) ? 1 : 2;
  archive_search_group->hint_patterns_per_search = hint_patterns_per_search;
  const uint64_t num_initial_searches = DIV_CEIL(bpm_gpu_buffer_get_max_queries(bpm_gpu_buffers),hint_patterns_per_search);
  archive_search_group->search_group_buffers = mm_calloc(total_search_groups,archive_search_group_buffer_t,true);
  uint64_t i;
  for (i=0;i<total_search_groups;++i) {
    archive_search_group->search_group_buffers[i].archive_searches = vector_new(num_initial_searches,archive_search_t*);
    archive_search_group->search_group_buffers[i].bpm_gpu_buffer = bpm_gpu_buffers + i;
  }
  archive_search_group->current_search_group_idx = 0;
  archive_search_group->current_search_group = archive_search_group->search_group_buffers;
  archive_search_group->mm_search = mm_search_new(mm_pool_get_slab(mm_pool_2MB));
  archive_search_group->archive_search_cache = archive_search_cache_new(mapper_parameters);
  // Return
  return archive_search_group;
}
GEM_INLINE void archive_search_group_clear(archive_search_group_t* const archive_search_group) {
  // Free all archive searches
  uint64_t i;
  for (i=0;i<archive_search_group->total_search_groups;++i) {
    vector_t* const archive_searches = archive_search_group->search_group_buffers[i].archive_searches;
    VECTOR_ITERATE(archive_searches,archive_search,j,archive_search_t*) {
      archive_search_cache_free(archive_search_group->archive_search_cache,*archive_search);
    }
    vector_clear(archive_searches);
    bpm_gpu_buffer_clear(archive_search_group->search_group_buffers[i].bpm_gpu_buffer);
  }
  archive_search_group->current_search_group_idx = 0;
  archive_search_group->current_search_group = archive_search_group->search_group_buffers;
  mm_search_clear(archive_search_group->mm_search);
}
GEM_INLINE void archive_search_group_delete(archive_search_group_t* const archive_search_group) {
  uint64_t i;
  for (i=0;i<archive_search_group->total_search_groups;++i) {
    vector_t* const archive_searches = archive_search_group->search_group_buffers[i].archive_searches;
    VECTOR_ITERATE(archive_searches,archive_search,j,archive_search_t*) {
      archive_search_delete(*archive_search);
    }
    vector_delete(archive_search_group->search_group_buffers[i].archive_searches);
  }
  mm_free(archive_search_group->search_group_buffers);
  archive_search_cache_delete(archive_search_group->archive_search_cache);
  mm_search_delete(archive_search_group->mm_search);
  mm_free(archive_search_group);
}
GEM_INLINE archive_search_t* archive_search_group_alloc(archive_search_group_t* const archive_search_group) {
  return archive_search_cache_alloc(archive_search_group->archive_search_cache,archive_search_group->mm_search);
}
GEM_INLINE bool archive_search_group_is_almost_full(archive_search_group_t* const archive_search_group) {
  return (archive_search_group->current_search_group_idx == archive_search_group->total_search_groups) &&
      bpm_gpu_buffer_almost_full(archive_search_group->current_search_group->bpm_gpu_buffer);
  // TODO Remove ME
}
GEM_INLINE bool archive_search_group_is_empty(archive_search_group_t* const archive_search_group) {
  return vector_is_empty(archive_search_group->current_search_group->archive_searches);
}
GEM_INLINE bool archive_search_group_fits_in_buffer(
    archive_search_group_t* const archive_search_group,archive_search_t* const archive_search) {
  // Calculate query dimensions
  const bool index_complement = archive_is_indexed_complement(archive_search->archive);
  const uint64_t num_patterns = (index_complement) ? 1 : 2;
  const uint64_t total_pattern_length = (index_complement) ?
      sequence_get_length(&archive_search->sequence) : 2*sequence_get_length(&archive_search->sequence);
  const uint64_t total_candidates = archive_search_get_search_canditates(archive_search);
  // Check if current search fits in buffer
  return bpm_gpu_buffer_fits_in_buffer(
      archive_search_group->current_search_group->bpm_gpu_buffer,
      num_patterns,total_pattern_length,total_candidates);
}
GEM_INLINE bool archive_search_group_add_search(
    archive_search_group_t* const archive_search_group,archive_search_t* const archive_search) {
  // Check if it fits in current buffer
  while (!archive_search_group_fits_in_buffer(archive_search_group,archive_search)) {
    gem_cond_fatal_error(vector_is_empty(archive_search_group->current_search_group->archive_searches),ARCHIVE_SEARCH_GROUP_QUERY_TOO_BIG);
    // Change group
    if (archive_search_group->current_search_group_idx < archive_search_group->total_search_groups) {
      // Send the current group to verification
      bpm_gpu_buffer_send(archive_search_group->current_search_group->bpm_gpu_buffer);
      // Next group
      ++(archive_search_group->current_search_group_idx);
      ++(archive_search_group->current_search_group);
    } else {
      return false;
    }
  }
  // Copy the candidates to the buffer
  archive_search_copy_candidates(archive_search,archive_search_group->current_search_group->bpm_gpu_buffer);
  // Add the search to the current group
  vector_insert(archive_search_group->current_search_group->archive_searches,archive_search,archive_search_t*);
  // Return ok
  return true;
}
GEM_INLINE void archive_search_group_retrieve_begin(archive_search_group_t* const archive_search_group) {
  // Send the current group to verification
  PROF_ADD_COUNTER(GP_ARCHIVE_SEARCH_GROUP_BUFFERS_USED,archive_search_group->current_search_group_idx);
  bpm_gpu_buffer_send(archive_search_group->current_search_group->bpm_gpu_buffer);
  // Initialize the iterator
  archive_search_group->current_search_group_idx = 0;
  archive_search_group->current_search_group = archive_search_group->search_group_buffers;
  archive_search_group->current_archive_search_idx = 0;
  archive_search_group->current_archive_search = vector_get_mem(archive_search_group->current_search_group->archive_searches,archive_search_t*);
  archive_search_group->total_archive_searches = vector_get_used(archive_search_group->current_search_group->archive_searches);
  // Fetch first group
  if (archive_search_group->total_archive_searches > 0) {
    PROF_START(GP_ARCHIVE_SEARCH_RETRIEVE_CANDIDATES_DELAY);
    bpm_gpu_buffer_receive(archive_search_group->current_search_group->bpm_gpu_buffer);
    PROF_STOP(GP_ARCHIVE_SEARCH_RETRIEVE_CANDIDATES_DELAY);
  }
}
GEM_INLINE bool archive_search_group_get_search(
    archive_search_group_t* const archive_search_group,
    archive_search_t** const archive_search,bpm_gpu_buffer_t** const bpm_gpu_buffer) {
  // Check end-of-iteration
  if (archive_search_group->current_archive_search_idx==archive_search_group->total_archive_searches) {
    // Next group
    ++(archive_search_group->current_search_group_idx);
    ++(archive_search_group->current_search_group);
    if (archive_search_group->current_search_group_idx==archive_search_group->total_search_groups) return false;
    // Clear archive-search iterator
    archive_search_group->current_archive_search_idx = 0;
    archive_search_group->current_archive_search = vector_get_mem(archive_search_group->current_search_group->archive_searches,archive_search_t*);
    archive_search_group->total_archive_searches = vector_get_used(archive_search_group->current_search_group->archive_searches);
    if (archive_search_group->total_archive_searches==0) return false;
    // Fetch group
    PROF_START(GP_ARCHIVE_SEARCH_RETRIEVE_CANDIDATES_DELAY);
    bpm_gpu_buffer_receive(archive_search_group->current_search_group->bpm_gpu_buffer);
    PROF_STOP(GP_ARCHIVE_SEARCH_RETRIEVE_CANDIDATES_DELAY);
  }
  // Return next archive-search
  *archive_search = *archive_search_group->current_archive_search;
  *bpm_gpu_buffer = archive_search_group->current_search_group->bpm_gpu_buffer;
  // Next archive-search
  ++(archive_search_group->current_archive_search_idx);
  ++(archive_search_group->current_archive_search);
  // Return ok
  return true;
}

