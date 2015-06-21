/*
 * PROJECT: GEMMapper
 * FILE: archive_search_group.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive_search_group.h"
#include "archive_search_se.h"
#include "archive_search_pe.h"
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
GEM_INLINE archive_search_t* archive_search_cache_alloc(archive_search_cache_t* const archive_search_cache) {
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
  archive_search_group->hint_patterns_per_search = 1;
  const uint64_t num_initial_searches = bpm_gpu_buffer_get_max_queries(bpm_gpu_buffers);
  archive_search_group->search_group_buffers = mm_calloc(total_search_groups,archive_search_group_buffer_t,true);
  uint64_t i;
  for (i=0;i<total_search_groups;++i) {
    archive_search_group->search_group_buffers[i].archive_searches = vector_new(num_initial_searches,archive_search_t*);
    archive_search_group->search_group_buffers[i].bpm_gpu_buffer = bpm_gpu_buffers + i;
  }
  archive_search_group->current_search_group_idx = 0;
  archive_search_group->current_search_group = archive_search_group->search_group_buffers;
  archive_search_group->mm_search = mm_search_new(mm_pool_get_slab(mm_pool_32MB));
  archive_search_group->archive_search_cache = archive_search_cache_new(mapper_parameters);
  // Return
  return archive_search_group;
}
GEM_INLINE void archive_search_group_init_bpm_buffers(archive_search_group_t* const archive_search_group) {
  // Initialize
  uint64_t i;
  for (i=0;i<archive_search_group->total_search_groups;++i) {
	  bpm_gpu_init_buffer(archive_search_group->search_group_buffers[i].bpm_gpu_buffer);
  }
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
GEM_INLINE archive_search_t* archive_search_group_allocate(archive_search_group_t* const archive_search_group) {
  // Alloc
  archive_search_t* const archive_search = archive_search_cache_alloc(archive_search_group->archive_search_cache);
  // Init archive search
  archive_search_single_end_configure(archive_search,archive_search_group->mm_search);
  text_collection_clear(&archive_search_group->mm_search->text_collection); // Clear text-collection
  // Return
  return archive_search;
}
GEM_INLINE void archive_search_group_allocate_pe(
    archive_search_group_t* const archive_search_group,
    archive_search_t** const archive_search_end1,archive_search_t** const archive_search_end2) {
  // Alloca
  *archive_search_end1 = archive_search_cache_alloc(archive_search_group->archive_search_cache);
  *archive_search_end2 = archive_search_cache_alloc(archive_search_group->archive_search_cache);
  // Init archive search
  archive_search_paired_end_configure(*archive_search_end1,*archive_search_end2,archive_search_group->mm_search);
  text_collection_clear(&archive_search_group->mm_search->text_collection); // Clear text-collection
}
GEM_INLINE bool archive_search_group_is_empty(archive_search_group_t* const archive_search_group) {
  return vector_is_empty(archive_search_group->current_search_group->archive_searches);
}
GEM_INLINE bool archive_search_group_fits_in_buffer(
    archive_search_group_t* const archive_search_group,
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2) {
  bpm_gpu_buffer_t* const bpm_gpu_buffer = archive_search_group->current_search_group->bpm_gpu_buffer;
  // Compute dimensions
  uint64_t total_entries = 0,total_query_chunks = 0,total_candidate_chunks = 0;
  bpm_gpu_buffer_compute_dimensions(bpm_gpu_buffer,
      &archive_search_end1->forward_search_state.pattern,
      archive_search_get_search_canditates(archive_search_end1),
      &total_entries,&total_query_chunks,&total_candidate_chunks);
  if (archive_search_end2!=NULL) {
    bpm_gpu_buffer_compute_dimensions(bpm_gpu_buffer,
        &archive_search_end2->forward_search_state.pattern,
        archive_search_get_search_canditates(archive_search_end2),
        &total_entries,&total_query_chunks,&total_candidate_chunks);
  }
  // Check if current search fits in buffer
  return bpm_gpu_buffer_fits_in_buffer(bpm_gpu_buffer,total_entries,total_query_chunks,total_candidate_chunks);
}
GEM_INLINE bool archive_search_group_add_search(
    archive_search_group_t* const archive_search_group,archive_search_t* const archive_search) {
  // Check if it fits in current buffer
  while (!archive_search_group_fits_in_buffer(archive_search_group,archive_search,NULL)) {
    const bool empty_group = vector_is_empty(archive_search_group->current_search_group->archive_searches);
    gem_cond_fatal_error(empty_group,ARCHIVE_SEARCH_GROUP_QUERY_TOO_BIG);
    // Change group
    if (archive_search_group->current_search_group_idx < archive_search_group->total_search_groups - 1) {
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
GEM_INLINE bool archive_search_group_add_paired_search(
    archive_search_group_t* const archive_search_group,
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2) {
  // Check if it fits in current buffer
  while (!archive_search_group_fits_in_buffer(archive_search_group,archive_search_end1,archive_search_end2)) {
    const bool empty_group = vector_is_empty(archive_search_group->current_search_group->archive_searches);
    gem_cond_fatal_error(empty_group,ARCHIVE_SEARCH_GROUP_QUERY_TOO_BIG);
    // Change group
    if (archive_search_group->current_search_group_idx < archive_search_group->total_search_groups - 1) {
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
  archive_search_copy_candidates(archive_search_end1,archive_search_group->current_search_group->bpm_gpu_buffer);
  archive_search_copy_candidates(archive_search_end2,archive_search_group->current_search_group->bpm_gpu_buffer);
  // Add the search to the current group
  vector_insert(archive_search_group->current_search_group->archive_searches,archive_search_end1,archive_search_t*);
  vector_insert(archive_search_group->current_search_group->archive_searches,archive_search_end2,archive_search_t*);
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
    PROF_START(GP_ARCHIVE_SEARCH_GROUP_RETRIEVE_CANDIDATES_DELAY);
    bpm_gpu_buffer_receive(archive_search_group->current_search_group->bpm_gpu_buffer);
    PROF_STOP(GP_ARCHIVE_SEARCH_GROUP_RETRIEVE_CANDIDATES_DELAY);
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
    PROF_START(GP_ARCHIVE_SEARCH_GROUP_RETRIEVE_CANDIDATES_DELAY);
    bpm_gpu_buffer_receive(archive_search_group->current_search_group->bpm_gpu_buffer);
    PROF_STOP(GP_ARCHIVE_SEARCH_GROUP_RETRIEVE_CANDIDATES_DELAY);
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
GEM_INLINE bool archive_search_group_get_paired_search(
    archive_search_group_t* const archive_search_group,
    archive_search_t** const archive_search_end1,bpm_gpu_buffer_t** const bpm_gpu_buffer_end1,
    archive_search_t** const archive_search_end2,bpm_gpu_buffer_t** const bpm_gpu_buffer_end2) {
  bool success;
  // Get End/1
  success = archive_search_group_get_search(archive_search_group,archive_search_end1,bpm_gpu_buffer_end1);
  if (!success) return false;
  // Get End/2
  success = archive_search_group_get_search(archive_search_group,archive_search_end2,bpm_gpu_buffer_end2);
  gem_cond_fatal_error(!success,ARCHIVE_SEARCH_GROUP_UNPAIRED_QUERY);
  // Return ok
  return true;
}
GEM_INLINE text_collection_t* archive_search_group_get_text_collection(archive_search_group_t* const archive_search_group) {
  return &archive_search_group->mm_search->text_collection;
}

