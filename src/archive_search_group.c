/*
 * PROJECT: GEMMapper
 * FILE: archive_search_group.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive_search_group.h"
#include "bpm_align_gpu.h"

/*
 * Archive-search group
 */
GEM_INLINE void archive_search_group_add(
    archive_search_group_t* const archive_search_group,
    archive_search_t* const archive_search,const uint64_t results_buffer_offset) {
  archive_search_member_t* archive_search_member;
  vector_alloc_new(archive_search_group->archive_searches,archive_search_member_t,archive_search_member);
  archive_search_member->archive_search = archive_search;
  archive_search_member->results_buffer_offset = results_buffer_offset;
}
/*
 * Dispatcher
 */
GEM_INLINE void archive_search_group_init(
    archive_search_group_t* const search_group,
    bpm_gpu_buffer_t* const bpm_gpu_buffers,const uint64_t num_initial_searches) {
  // State
//  search_group->mayor_group_id = ;
//  search_group->minor_group_id = ;
  search_group->state = archive_search_group_free;
  // BPM-GPU candidates buffer
  search_group->bpm_gpu_buffer = bpm_gpu_buffers;
  // Archive searches
  search_group->mm_search = mm_search_new(mm_pool_get_slab(mm_pool_8MB));
  search_group->archive_searches = vector_new(num_initial_searches,archive_search_member_t);
//  search_group->archive_searches_attr = vector_new();
}
GEM_INLINE void archive_search_group_clear(archive_search_group_t* const search_group) {
  mm_search_clear(search_group->mm_search);
  vector_clear(search_group->archive_searches);
}
GEM_INLINE void archive_search_group_destroy(archive_search_group_t* const search_group) {
  mm_search_delete(search_group->mm_search);
  vector_delete(search_group->archive_searches);
// vector_delete(search_group->archive_searches_attr);
}
GEM_INLINE archive_search_group_dispatcher_t* archive_search_group_dispatcher_new(
    archive_t* const archive,const uint64_t num_search_groups,
    const uint64_t average_query_size,const uint64_t candidates_per_query) {
  // Allocate
  archive_search_group_dispatcher_t* const dispatcher = mm_alloc(archive_search_group_dispatcher_t);
  // Dispatcher State
  dispatcher->num_groups = num_search_groups;
  dispatcher->num_groups_free = num_search_groups;
  dispatcher->num_groups_generating = 0;
  dispatcher->num_groups_verifying = 0;
  dispatcher->num_groups_selecting = 0;
  /* Dispatcher Record */
  dispatcher->num_threads_generating = 0;
  // BPM-GPU Buffer
  dispatcher->bpm_gpu_buffer_collection =
      bpm_gpu_init(archive->enc_text,num_search_groups,average_query_size,candidates_per_query);
  // Search Groups
  bpm_gpu_buffer_t* const bpm_gpu_buffers = dispatcher->bpm_gpu_buffer_collection->bpm_gpu_buffers;
  const uint64_t query_per_search = (archive_is_indexed_complement(archive)) ? 1 : 2;
  uint64_t i;
  dispatcher->search_group = mm_calloc(num_search_groups,archive_search_group_t,true);
  for (i=0;i<num_search_groups;++i) {
    const uint64_t num_initial_searches = DIV_CEIL(bpm_gpu_buffer_get_max_queries(bpm_gpu_buffers+i),query_per_search);
    archive_search_group_init(dispatcher->search_group+i,bpm_gpu_buffers+i,num_initial_searches);
  }
  // Mutex/CV
  MUTEX_INIT(dispatcher->dispatcher_mutex);
  CV_INIT(dispatcher->groups_free_cond);
  CV_INIT(dispatcher->groups_verifying_cond);
  // Return
  return dispatcher;
}
GEM_INLINE void archive_search_group_dispatcher_delete(archive_search_group_dispatcher_t* const dispatcher) {
  // BPM-GPU Buffer
  bpm_gpu_destroy(dispatcher->bpm_gpu_buffer_collection);
  // Search Groups
  uint64_t i;
  archive_search_group_t* search_group = dispatcher->search_group;
  for (i=0;i<dispatcher->num_groups;++i,++search_group) {
    archive_search_group_destroy(search_group);
  }
  mm_free(dispatcher->search_group);
  // Mutex/CV
  MUTEX_DESTROY(dispatcher->dispatcher_mutex);
  CV_DESTROY(dispatcher->groups_free_cond);
  CV_DESTROY(dispatcher->groups_verifying_cond);
  // Free handler
  mm_free(dispatcher);
}
GEM_INLINE void archive_search_group_dispatcher_register_generating(
    archive_search_group_dispatcher_t* const dispatcher,const uint64_t num_threads) {
  MUTEX_BEGIN_SECTION(dispatcher->dispatcher_mutex) {
    dispatcher->num_threads_generating += num_threads;
  } MUTEX_END_SECTION(dispatcher->dispatcher_mutex);
}
GEM_INLINE void archive_search_group_dispatcher_deregister_generating(
    archive_search_group_dispatcher_t* const dispatcher,const uint64_t num_threads) {
  MUTEX_BEGIN_SECTION(dispatcher->dispatcher_mutex) {
    dispatcher->num_threads_generating -= num_threads;
    if (dispatcher->num_threads_generating == 0) {
      CV_BROADCAST(dispatcher->groups_verifying_cond);
    }
  } MUTEX_END_SECTION(dispatcher->dispatcher_mutex);
}
GEM_INLINE archive_search_group_t* archive_search_group_dispatcher_request_generating(
    archive_search_group_dispatcher_t* const dispatcher) {
  archive_search_group_t* search_group = NULL;
  MUTEX_BEGIN_SECTION(dispatcher->dispatcher_mutex) {
    // Wait for one free group
    while (dispatcher->num_groups_free==0) {
      CV_WAIT(dispatcher->groups_free_cond,dispatcher->dispatcher_mutex);
    }
    // Find a free group
    const uint64_t num_groups = dispatcher->num_groups;
    uint64_t i;
    for(i=0;i<num_groups;++i) {
      if (dispatcher->search_group[i].state==archive_search_group_free) {
        search_group = dispatcher->search_group + i;
        break;
      }
    }
    // Change state
    search_group->state = archive_search_group_generating_candidates;
    --(dispatcher->num_groups_free);
    ++(dispatcher->num_groups_generating);
  } MUTEX_END_SECTION(dispatcher->dispatcher_mutex);
  // Return
  return search_group;
}
GEM_INLINE void archive_search_group_dispatcher_return_generating(
    archive_search_group_dispatcher_t* const dispatcher,archive_search_group_t* const search_group) {
  MUTEX_BEGIN_SECTION(dispatcher->dispatcher_mutex) {
    // Send buffer
    bpm_gpu_buffer_send(search_group->bpm_gpu_buffer);
    // Check if broadcast needed
    if (dispatcher->num_groups_verifying==0) {
      CV_BROADCAST(dispatcher->groups_verifying_cond);
    }
    // Change state
    search_group->state = archive_search_group_verifying_candidates;
    --(dispatcher->num_groups_generating);
    ++(dispatcher->num_groups_verifying);
  } MUTEX_END_SECTION(dispatcher->dispatcher_mutex);
}
GEM_INLINE archive_search_group_t* archive_search_group_dispatcher_request_selecting(
    archive_search_group_dispatcher_t* const dispatcher) {
  archive_search_group_t* search_group = NULL;
  MUTEX_BEGIN_SECTION(dispatcher->dispatcher_mutex) {
    // Wait for one group being verified
    while (dispatcher->num_groups_verifying==0 && dispatcher->num_threads_generating>0) {
      CV_WAIT(dispatcher->groups_verifying_cond,dispatcher->dispatcher_mutex);
    }
    // Check exit condition (num_groups_verifying == 0 && num_threads_generating == 0)
    if (dispatcher->num_groups_verifying > 0) {
      // Find a group being verified
      const uint64_t num_groups = dispatcher->num_groups;
      uint64_t i;
      for(i=0;i<num_groups;++i) {
        if (dispatcher->search_group[i].state==archive_search_group_verifying_candidates) {
          search_group = dispatcher->search_group + i;
          break;
        }
      }
      // Change state
      search_group->state = archive_search_group_selecting_candidates;
      --(dispatcher->num_groups_verifying);
      ++(dispatcher->num_groups_selecting);
    }
  } MUTEX_END_SECTION(dispatcher->dispatcher_mutex);
  // Receive buffer
  if (search_group!=NULL) bpm_gpu_buffer_receive(search_group->bpm_gpu_buffer);
  // Return
  return search_group;
}
GEM_INLINE void archive_search_group_dispatcher_return_selecting(
    archive_search_group_dispatcher_t* const dispatcher,archive_search_group_t* const search_group) {
  MUTEX_BEGIN_SECTION(dispatcher->dispatcher_mutex) {
    // Check if broadcast needed
    if (dispatcher->num_groups_free==0) {
      CV_BROADCAST(dispatcher->groups_free_cond);
    }
    // Clear search_group
    archive_search_group_clear(search_group);
    // Change state
    search_group->state = archive_search_group_free;
    --(dispatcher->num_groups_selecting);
    ++(dispatcher->num_groups_free);
  } MUTEX_END_SECTION(dispatcher->dispatcher_mutex);
}
/*
 * Step-wise SE-Search
 */
GEM_INLINE void archive_search_generate_candidates(
    archive_search_t* const archive_search,mm_search_t* const mm_search) {
  ARCHIVE_SEARCH_CHECK(archive_search);
  // Prepare
  archive_search_prepare_sequence(archive_search,mm_search->mm_stack); // Prepare pattern(s)
  archive_search_clear(archive_search); // Clean Matches
  // Check mapping mode
  gem_cond_fatal_error(
      archive_search->search_actual_parameters.search_parameters->mapping_mode!=mapping_adaptive_filtering ||
      archive_search->search_actual_parameters.search_parameters->mapping_mode!=mapping_fixed_filtering,
      ARCHIVE_SEARCH_GROUP_MAPPING_MODE_NOT_SUPPORTED);
  // Run the search (FORWARD)
  approximate_search_t* const forward_asearch = archive_search->forward_search_state;
  forward_asearch->stop_search_stage = asearch_filtering; // Stop before filtering
  forward_asearch->search_strand = Forward; // Configure forward search
  approximate_search(forward_asearch,NULL,mm_search);
  if (archive_search->search_reverse) {
    // Run the search (REVERSE)
    approximate_search_t* const reverse_asearch = archive_search->reverse_search_state;
    reverse_asearch->stop_search_stage = asearch_filtering; // Stop before filtering
    reverse_asearch->search_strand = Reverse; // Configure reverse search
    approximate_search(reverse_asearch,NULL,mm_search);
  }
}
GEM_INLINE void archive_search_copy_candidates(
    archive_search_t* const archive_search,bpm_gpu_buffer_t* const bpm_gpu_buffer) {
  const archive_t* const archive = archive_search->archive;
  // Add candidates (FORWARD)
  approximate_search_t* const forward_asearch = archive_search->forward_search_state;
  forward_asearch->num_potential_candidates = filtering_candidates_add_to_bpm_buffer(
      &forward_asearch->filtering_candidates,archive->locator,archive->fm_index,archive->enc_text,
      &forward_asearch->pattern,forward_asearch->search_strand,forward_asearch->search_actual_parameters,bpm_gpu_buffer);
  if (archive_search->search_reverse) {
    // Add candidates (REVERSE)
    approximate_search_t* const reverse_asearch = archive_search->reverse_search_state;
    reverse_asearch->num_potential_candidates = filtering_candidates_add_to_bpm_buffer(
        &reverse_asearch->filtering_candidates,archive->locator,archive->fm_index,archive->enc_text,
        &reverse_asearch->pattern,reverse_asearch->search_strand,reverse_asearch->search_actual_parameters,bpm_gpu_buffer);
  }
}
GEM_INLINE void archive_search_select_candidates(
    archive_search_t* const archive_search,
    bpm_gpu_buffer_t* const bpm_gpu_buffer,const uint64_t results_buffer_offset) {
  GEM_NOT_IMPLEMENTED(); // TODO
}
