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
  /* Mapper parameters */
  mapper_parameters_t* mapper_parameters;
  /* Slab of archive_search_t */
  vector_t* archive_search_cache;  // Already allocated & configured (archive_search_t*)
  pthread_mutex_t mutex;           // Mutex to access the cache
} archive_search_cache_t;
typedef struct {
  /* Search */
  archive_search_t* archive_search;     // Archive search
  /* Candidates verified */
  uint64_t results_buffer_offset;       // Offset in the results vector
} search_group_member_t;
typedef enum {
  archive_search_group_free,
  archive_search_group_generating_candidates,
  archive_search_group_verifying_candidates,
  archive_search_group_selecting_candidates
} archive_search_group_state_t;
struct _search_group_t {
  /* Dispatcher (owner) */
  search_group_dispatcher_t* dispatcher; // Dispatcher
  /* State */
  uint32_t mayor_group_id;               // Group mayor ID
  uint32_t minor_group_id;               // Group minor ID
  archive_search_group_state_t state;    // State of the search group
  bool group_incomplete;                 // Is incomplete if multisearch-group and not the last group
  /* BPM-GPU candidates buffer */
  bpm_gpu_buffer_t* bpm_gpu_buffer;      // BPM-Buffer
  /* Archive searches */
  vector_t* search_group_members;      // Vector of search members (search_group_member_t)
  /* MM */
  mm_search_t* mm_search;                // Memory Managed Search
};
struct _search_group_dispatcher_t {
  /* Dispatcher State */
  uint64_t num_groups;                  // Total number of search-groups allocated
  uint64_t num_groups_free;             // Free groups (ready to be filled & used)
  uint64_t num_groups_generating;       // Groups dispatched for candidate generation
  uint64_t num_groups_verifying;        // Groups being verified (BPM-CUDA)
  uint64_t num_groups_selecting;        // Groups dispatched for candidate selection
  uint64_t num_threads_generating;      // Dispatcher Record of Generating-Threads
  uint32_t next_request_id;             // Block ID  (for synchronization purposes)
  pqueue_t* requests;                   // Priority queue for the request
  /* Search Groups */
  uint64_t hint_queries_per_search;              // Hint: Expected number of queries per search-group
  search_group_t* search_group;                  // Search Groups
  uint64_t search_group_used;                    // Number of Search Groups Activated
  archive_search_cache_t* archive_search_cache;  // Archive-search cache
  /* BPM-GPU Buffer*/
  bpm_gpu_buffer_collection_t* bpm_gpu_buffer_collection; // BPM Buffers
  /* Output File */
  output_file_t* output_file;           // Output File
  /* Mutex/CV */
  pthread_mutex_t dispatcher_mutex;
  pthread_cond_t groups_free_cond;
  pthread_cond_t groups_verifying_cond;
};

/*
 * Archive Search Cache
 */
GEM_INLINE archive_search_cache_t* archive_search_cache_new(mapper_parameters_t* const mapper_parameters) {
  // Alloc
  archive_search_cache_t* const archive_search_cache = mm_alloc(archive_search_cache_t);
  // Initialize cache
  archive_search_cache->mapper_parameters = mapper_parameters;
  archive_search_cache->archive_search_cache =
      vector_new(MAPPER_CUDA_ARCHIVE_SEARCH_CACHE_INIT_SIZE,archive_search_t*);
  MUTEX_INIT(archive_search_cache->mutex);
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
  MUTEX_DESTROY(archive_search_cache->mutex);
  mm_free(archive_search_cache);
}
GEM_INLINE archive_search_t* archive_search_cache_alloc(
    archive_search_cache_t* const archive_search_cache,mm_search_t* const mm_search) {
  archive_search_t* archive_search = NULL;
  MUTEX_BEGIN_SECTION(archive_search_cache->mutex) {
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
  } MUTEX_END_SECTION(archive_search_cache->mutex);
  // Init archive search
  archive_search_configure(archive_search,mm_search);
  // Return
  return archive_search;
}
GEM_INLINE void archive_search_cache_free(
    archive_search_cache_t* const archive_search_cache,archive_search_t* const archive_search) {
  MUTEX_BEGIN_SECTION(archive_search_cache->mutex) {
    // Add it to the cache
    vector_insert(archive_search_cache->archive_search_cache,archive_search,archive_search_t*);
  } MUTEX_END_SECTION(archive_search_cache->mutex);
}
/*
 * Archive-search group
 */
GEM_INLINE void search_group_clear(search_group_t* const search_group) {
  search_group->mayor_group_id = 0;
  search_group->minor_group_id = 0;
  search_group->group_incomplete = false;
  bpm_buffer_clear(search_group->bpm_gpu_buffer);
  vector_clear(search_group->search_group_members);
  mm_search_clear(search_group->mm_search);
}
GEM_INLINE void search_group_destroy(search_group_t* const search_group) {
  mm_search_delete(search_group->mm_search);
  vector_delete(search_group->search_group_members);
}
GEM_INLINE void search_group_set_incomplete(search_group_t* const search_group) {
  search_group->group_incomplete = true;
}
GEM_INLINE bool search_group_is_incomplete(search_group_t* const search_group) {
  return search_group->group_incomplete;
}
GEM_INLINE uint64_t search_group_get_group_id(search_group_t* const search_group) {
  return search_group->mayor_group_id;
}
GEM_INLINE bpm_gpu_buffer_t* search_group_get_bpm_buffer(search_group_t* const search_group) {
  return search_group->bpm_gpu_buffer;
}
GEM_INLINE uint64_t search_group_get_num_searches(search_group_t* const search_group) {
  return vector_get_used(search_group->search_group_members);
}
GEM_INLINE void search_group_get_search(
    search_group_t* const search_group,const uint64_t position,
    archive_search_t** const archive_search,uint64_t* const results_buffer_offset) {
  search_group_member_t* const search_group_member =
      vector_get_elm(search_group->search_group_members,position,search_group_member_t);
  *archive_search = search_group_member->archive_search;
  *results_buffer_offset = search_group_member->results_buffer_offset;
}
GEM_INLINE void search_group_add_search(
    search_group_t* const search_group,
    archive_search_t* const archive_search,const uint64_t results_buffer_offset) {
  search_group_member_t* search_group_member;
  vector_alloc_new(search_group->search_group_members,search_group_member_t,search_group_member);
  search_group_member->archive_search = archive_search;
  search_group_member->results_buffer_offset = results_buffer_offset;
}
// Archive Search Group Allocator (Cache)
GEM_INLINE archive_search_t* search_group_alloc(search_group_t* const archive_search_group) {
  return archive_search_cache_alloc(archive_search_group->dispatcher->archive_search_cache,archive_search_group->mm_search);
}
GEM_INLINE void search_group_release(
    search_group_t* const archive_search_group,archive_search_t* const archive_search) {
  archive_search_cache_free(archive_search_group->dispatcher->archive_search_cache,archive_search);
}
/*
 * Dispatcher
 */
GEM_INLINE void search_group_dispatcher_init_search_group(
    search_group_dispatcher_t* const dispatcher,const uint64_t search_group_position) {
  // Search group
  search_group_t* const search_group = dispatcher->search_group + search_group_position;
  // Dispatcher
  search_group->dispatcher = dispatcher;
  // State
  search_group->mayor_group_id = 0;
  search_group->minor_group_id = 0;
  search_group->state = archive_search_group_free;
  search_group->group_incomplete = false;
  // BPM-GPU candidates buffer
  bpm_gpu_buffer_t* const bpm_gpu_buffers = dispatcher->bpm_gpu_buffer_collection->bpm_gpu_buffers;
  search_group->bpm_gpu_buffer = bpm_gpu_buffers+search_group_position;
  // Archive searches
  const uint64_t num_initial_searches =
      DIV_CEIL(bpm_gpu_buffer_get_max_queries(search_group->bpm_gpu_buffer),dispatcher->hint_queries_per_search);
  search_group->search_group_members = vector_new(num_initial_searches,search_group_member_t);
  search_group->mm_search = mm_search_new(mm_pool_get_slab(mm_pool_2MB)); // FIXME 8MB
}
GEM_INLINE search_group_dispatcher_t* search_group_dispatcher_new(
    mapper_parameters_t* const mapper_parameters,
    archive_t* const archive,const uint64_t num_search_groups,
    const uint64_t average_query_size,const uint64_t candidates_per_query) {
  // Allocate
  search_group_dispatcher_t* const dispatcher = mm_alloc(search_group_dispatcher_t);
  // Dispatcher State
  dispatcher->num_groups = num_search_groups;
  dispatcher->num_groups_free = num_search_groups;
  dispatcher->num_groups_generating = 0;
  dispatcher->num_groups_verifying = 0;
  dispatcher->num_groups_selecting = 0;
  dispatcher->num_threads_generating = 0; // Dispatcher Record
  dispatcher->next_request_id = 0; // Next Block ID
  dispatcher->requests = pqueue_new(2*num_search_groups);
  // BPM-GPU Buffer
  dispatcher->bpm_gpu_buffer_collection =
      bpm_gpu_init(archive->enc_text,num_search_groups,average_query_size,candidates_per_query);
  // Archive Search Groups
  dispatcher->hint_queries_per_search = (archive_is_indexed_complement(archive)) ? 1 : 2;
  dispatcher->search_group = mm_calloc(num_search_groups,search_group_t,true);
  dispatcher->search_group_used = 0; // No search-group initialized
  // Archive Search Cache
  dispatcher->archive_search_cache = archive_search_cache_new(mapper_parameters);
  // Mutex/CV
  MUTEX_INIT(dispatcher->dispatcher_mutex);
  CV_INIT(dispatcher->groups_free_cond);
  CV_INIT(dispatcher->groups_verifying_cond);
  // Return
  return dispatcher;
}
GEM_INLINE void search_group_dispatcher_delete(search_group_dispatcher_t* const dispatcher) {
  // Dispatcher State
  pqueue_delete(dispatcher->requests);
  // BPM-GPU Buffer
  bpm_gpu_destroy(dispatcher->bpm_gpu_buffer_collection);
  // Archive Search Groups
  uint64_t i;
  search_group_t* search_group = dispatcher->search_group;
  for (i=0;i<dispatcher->search_group_used;++i,++search_group) {
    search_group_destroy(search_group);
  }
  mm_free(dispatcher->search_group);
  // Archive Search Cache
  archive_search_cache_delete(dispatcher->archive_search_cache);
  // Mutex/CV
  MUTEX_DESTROY(dispatcher->dispatcher_mutex);
  CV_DESTROY(dispatcher->groups_free_cond);
  CV_DESTROY(dispatcher->groups_verifying_cond);
  // Free handler
  mm_free(dispatcher);
}
/*
 * Register/Deregister Threads
 */
GEM_INLINE void search_group_dispatcher_register_generating(
    search_group_dispatcher_t* const dispatcher,const uint64_t num_threads) {
  MUTEX_BEGIN_SECTION(dispatcher->dispatcher_mutex) {
    dispatcher->num_threads_generating += num_threads;
  } MUTEX_END_SECTION(dispatcher->dispatcher_mutex);
}
GEM_INLINE void search_group_dispatcher_deregister_generating(
    search_group_dispatcher_t* const dispatcher,const uint64_t num_threads) {
  MUTEX_BEGIN_SECTION(dispatcher->dispatcher_mutex) {
    dispatcher->num_threads_generating -= num_threads;
    if (dispatcher->num_threads_generating == 0) {
      CV_BROADCAST(dispatcher->groups_verifying_cond);
    }
  } MUTEX_END_SECTION(dispatcher->dispatcher_mutex);
}
/*
 * Generating group
 */
GEM_INLINE bool search_group_dispatcher_generating_serve_group_cond(
    search_group_dispatcher_t* const dispatcher,const uint32_t request_id) {
  // 0. We need a free group to serve
  if (dispatcher->num_groups_free==0) return false;
  // 1. Serve if request is the next in-order (next_request_id)
  const bool next_in_order = dispatcher->next_request_id==request_id;
  if (next_in_order) return true;
  // 2.1 Serve if we have at least one more search-group left to serve the next in-order (next_request_id)
  // 2.2 and the request is the next in the priority queue
  if (dispatcher->num_groups_free>1 && pqueue_top_priority(dispatcher->requests)==request_id) {
    return true;
  } else {
    return false;
  }
}
GEM_INLINE bool search_group_dispatcher_generating_eligible_request_cond(search_group_dispatcher_t* const dispatcher) {
  return !pqueue_is_empty(dispatcher->requests) && dispatcher->num_groups_free>0;
}
GEM_INLINE search_group_t* search_group_dispatcher_generating_get_free_buffer(search_group_dispatcher_t* const dispatcher) {
  // Find a free group
  const uint64_t num_groups = dispatcher->num_groups;
  const uint64_t search_group_used = dispatcher->search_group_used;
  uint64_t i;
  for (i=0;i<search_group_used;++i) {
    if (dispatcher->search_group[i].state==archive_search_group_free) {
      return dispatcher->search_group + i;
    }
  }
  if (search_group_used < num_groups) {
    search_group_dispatcher_init_search_group(dispatcher,search_group_used);
    ++(dispatcher->search_group_used);
    return dispatcher->search_group + search_group_used;
  }
  gem_fatal_check_msg(true,"Archive search-group dispatcher. No free group could be found");
  return NULL;
}
GEM_INLINE void search_group_dispatcher_print_groups(
    FILE* stream,search_group_dispatcher_t* const dispatcher) {
  const uint64_t search_group_used = dispatcher->search_group_used;
  uint64_t i;
  fprintf(stream,"#{%lu,%lu,%lu,%lu}  ",
      dispatcher->num_groups_free,dispatcher->num_groups_generating,
      dispatcher->num_groups_verifying,dispatcher->num_groups_selecting);
  for (i=0;i<search_group_used;++i) {
    fprintf(stream,"[%s,%u,%u,%s] ",
        dispatcher->search_group[i].state==archive_search_group_free ? "Fre" :
        (dispatcher->search_group[i].state==archive_search_group_generating_candidates ? "Gen" :
        (dispatcher->search_group[i].state==archive_search_group_verifying_candidates ? "Ver" : "Sel")),
        dispatcher->search_group[i].mayor_group_id,dispatcher->search_group[i].minor_group_id,
        dispatcher->search_group[i].group_incomplete ? "I" : "C");
  }
}
GEM_INLINE search_group_t* search_group_dispatcher_request_generating(
    search_group_dispatcher_t* const dispatcher,const uint32_t request_id) {
  PROF_INC_COUNTER(GP_SGDISPATCHER_REQUESTS_GENERATING);
  search_group_t* search_group = NULL;
  MUTEX_BEGIN_SECTION(dispatcher->dispatcher_mutex) {
    // Add request to queue
    pqueue_push(dispatcher->requests,NULL,request_id);
    while (!search_group_dispatcher_generating_serve_group_cond(dispatcher,request_id)) {
      PROF_INC_COUNTER(GP_SGDISPATCHER_REQUESTS_GENERATING_STALLS);
      PROF_BLOCK() {
        if (dispatcher->num_groups_free==0) {
          PROF_INC_COUNTER(GP_SGDISPATCHER_REQUESTS_GENERATING_STALLS_BUSY);
        } else {
          PROF_INC_COUNTER(GP_SGDISPATCHER_REQUESTS_GENERATING_STALLS_NOT_PRIORITY);
        }
      }
      CV_WAIT(dispatcher->groups_free_cond,dispatcher->dispatcher_mutex);
    }
    // Get a free buffer
    search_group = search_group_dispatcher_generating_get_free_buffer(dispatcher);
    // Change state [generating_candidates]
    search_group->mayor_group_id = request_id;
    search_group->minor_group_id = 0;       // First of possibly multisearch-group
    search_group->group_incomplete = false; // Not a multisearch-group (yet)
    --(dispatcher->num_groups_free);
    ++(dispatcher->num_groups_generating);
    search_group->state = archive_search_group_generating_candidates;
    // Update next in-order & priority-queue
    pqueue_pop(dispatcher->requests,void);
    // Broadcast if any there are requests eligible (possible satisfiable)
    if (search_group_dispatcher_generating_eligible_request_cond(dispatcher)) {
      CV_BROADCAST(dispatcher->groups_free_cond);
    }
  } MUTEX_END_SECTION(dispatcher->dispatcher_mutex);
  // Return
  return search_group;
}
GEM_INLINE search_group_t* search_group_dispatcher_request_generating_extension(
    search_group_dispatcher_t* const dispatcher,search_group_t* const search_group) {
  PROF_INC_COUNTER(GP_SGDISPATCHER_REQUESTS_GENERATING_EXTENSION);
  // Store IDs and set incomplete
  const uint64_t mayor_group_id = search_group->mayor_group_id;
  const uint64_t minor_group_id = search_group->minor_group_id;
  search_group_set_incomplete(search_group);
  // Return current group
  search_group_dispatcher_return_generating(dispatcher,search_group);
  // Request new one
  search_group_t* const search_group_extension =
      search_group_dispatcher_request_generating(dispatcher,mayor_group_id);
  search_group_extension->mayor_group_id = mayor_group_id;
  search_group_extension->minor_group_id = minor_group_id+1;
  return search_group_extension;
}
GEM_INLINE void search_group_dispatcher_return_generating(
    search_group_dispatcher_t* const dispatcher,search_group_t* const search_group) {
  // Send buffer
  bpm_gpu_buffer_send(search_group->bpm_gpu_buffer);
  // Return Buffer to Dispatcher
  MUTEX_BEGIN_SECTION(dispatcher->dispatcher_mutex) {
    // Update Next ID
    if (!search_group->group_incomplete) ++dispatcher->next_request_id;
    // Broadcast (new group generating available)
    CV_BROADCAST(dispatcher->groups_verifying_cond);
    // Change state
    --(dispatcher->num_groups_generating);
    ++(dispatcher->num_groups_verifying);
    search_group->state = archive_search_group_verifying_candidates;
  } MUTEX_END_SECTION(dispatcher->dispatcher_mutex);
}
/*
 * Selecting group
 */
GEM_INLINE search_group_t* search_group_dispatcher_generating_get_new_verifying(
    search_group_dispatcher_t* const dispatcher) {
  // Find the lowest ID
  search_group_t* candidate = NULL;
  const uint64_t search_group_used = dispatcher->search_group_used;
  uint64_t i;
  for (i=0;i<search_group_used;++i) {
    if (dispatcher->search_group[i].state==archive_search_group_verifying_candidates &&
        dispatcher->search_group[i].minor_group_id==0) {
      if (candidate==NULL ||
           candidate->mayor_group_id> dispatcher->search_group[i].mayor_group_id ||
          (candidate->mayor_group_id==dispatcher->search_group[i].mayor_group_id &&
           candidate->minor_group_id> dispatcher->search_group[i].minor_group_id)) {
        candidate = dispatcher->search_group + i;
      }
    }
  }
  return candidate;
}
GEM_INLINE search_group_t* search_group_dispatcher_generating_get_next_part(
    search_group_dispatcher_t* const dispatcher,
    const uint64_t mayor_group_id,const uint64_t minor_group_id) {
  const uint64_t search_group_used = dispatcher->search_group_used;
  uint64_t i;
  for (i=0;i<search_group_used;++i) {
    if (dispatcher->search_group[i].state==archive_search_group_verifying_candidates &&
        dispatcher->search_group[i].mayor_group_id==mayor_group_id &&
        dispatcher->search_group[i].minor_group_id==minor_group_id) {
      return dispatcher->search_group + i;
    }
  }
  return NULL;
}
GEM_INLINE search_group_t* search_group_dispatcher_request_selecting_(
    search_group_dispatcher_t* const dispatcher,const bool multisearch_group,
    const uint64_t mayor_group_id,const uint64_t minor_group_id) {
  PROF_INC_COUNTER(GP_SGDISPATCHER_REQUESTS_SELECTING);
  search_group_t* search_group = NULL;
  MUTEX_BEGIN_SECTION(dispatcher->dispatcher_mutex) {
    while (true) {
      // Wait for one group being verified
      while (dispatcher->num_groups_verifying==0 && dispatcher->num_threads_generating>0) {
        PROF_INC_COUNTER(GP_SGDISPATCHER_REQUESTS_SELECTING_STALLS);
        PROF_INC_COUNTER(GP_SGDISPATCHER_REQUESTS_SELECTING_STALLS_IDLE);
        CV_WAIT(dispatcher->groups_verifying_cond,dispatcher->dispatcher_mutex);
      }
      // Check exit condition
      if (dispatcher->num_groups_verifying==0 && dispatcher->num_threads_generating==0) break;
      // Find a group being verified
      if (multisearch_group) {
        // Find the extension of a multisearch-group
        search_group = search_group_dispatcher_generating_get_next_part(dispatcher,mayor_group_id,minor_group_id);
        if (search_group!=NULL) break; // Found!
        // Multisearch-groups extensions)
        PROF_INC_COUNTER(GP_SGDISPATCHER_REQUESTS_SELECTING_STALLS);
        PROF_INC_COUNTER(GP_SGDISPATCHER_REQUESTS_SELECTING_STALLS_EXTENSION_NOT_READY);
        CV_WAIT(dispatcher->groups_verifying_cond,dispatcher->dispatcher_mutex);
      } else {
        // Find a single group (Not the extension of a multisearch-group; minor_group_id==0)
        search_group = search_group_dispatcher_generating_get_new_verifying(dispatcher);
        if (search_group!=NULL) break; // Found!
        // No group eligible (only multisearch-groups extensions)
        PROF_INC_COUNTER(GP_SGDISPATCHER_REQUESTS_SELECTING_STALLS);
        PROF_INC_COUNTER(GP_SGDISPATCHER_REQUESTS_SELECTING_STALLS_NO_SINGLE_GROUPS);
        CV_WAIT(dispatcher->groups_verifying_cond,dispatcher->dispatcher_mutex);
      }
    }
    // Change state
    if (search_group!=NULL) {
      --(dispatcher->num_groups_verifying);
      ++(dispatcher->num_groups_selecting);
      search_group->state = archive_search_group_selecting_candidates;
    }
  } MUTEX_END_SECTION(dispatcher->dispatcher_mutex);
  // Receive buffer
  if (search_group!=NULL) bpm_gpu_buffer_receive(search_group->bpm_gpu_buffer);
  // Return
  return search_group;
}
GEM_INLINE search_group_t* search_group_dispatcher_request_selecting(
    search_group_dispatcher_t* const dispatcher) {
  return search_group_dispatcher_request_selecting_(dispatcher,false,0,0);
}
GEM_INLINE search_group_t* search_group_dispatcher_request_selecting_next(
    search_group_dispatcher_t* const dispatcher,search_group_t* search_group) {
  PROF_INC_COUNTER(GP_SGDISPATCHER_REQUESTS_SELECTING_EXTENSION);
  // Store IDs & return the current search-group
  const uint64_t mayor_group_id = search_group->mayor_group_id;
  const uint64_t minor_group_id = search_group->minor_group_id+1;
  search_group_dispatcher_return_selecting(dispatcher,search_group);
  // Look-up the next search-group (of the multisearch-groups)
  return search_group_dispatcher_request_selecting_(dispatcher,true,mayor_group_id,minor_group_id);
}
GEM_INLINE void search_group_dispatcher_return_selecting(
    search_group_dispatcher_t* const dispatcher,search_group_t* const search_group) {
  MUTEX_BEGIN_SECTION(dispatcher->dispatcher_mutex) {
    // Check if broadcast needed
    if (!pqueue_is_empty(dispatcher->requests)) {
      CV_BROADCAST(dispatcher->groups_free_cond);
    }
    // Clear search_group
    search_group_clear(search_group);
    // Change state
    --(dispatcher->num_groups_selecting);
    ++(dispatcher->num_groups_free);
    search_group->state = archive_search_group_free;
  } MUTEX_END_SECTION(dispatcher->dispatcher_mutex);
}
/*
 * Step-wise SE-Search
 */
GEM_INLINE void archive_search_generate_candidates(archive_search_t* const archive_search) {
  ARCHIVE_SEARCH_CHECK(archive_search);
  PROF_START(GP_ARCHIVE_SEARCH_GENERATE_CANDIDATES);
  // Reset initial values (Prepare pattern(s), instantiate parameters values, ...)
  archive_search_reset(archive_search,sequence_get_length(&archive_search->sequence));
  // Check mapping mode
  gem_cond_fatal_error(
      archive_search->search_actual_parameters.search_parameters->mapping_mode!=mapping_adaptive_filtering,
      ARCHIVE_SEARCH_GROUP_MAPPING_MODE_NOT_SUPPORTED);
  // Run the search (FORWARD)
  approximate_search_t* const forward_asearch = &archive_search->forward_search_state;
  forward_asearch->stop_search_stage = asearch_filtering; // Stop before filtering
  forward_asearch->search_strand = Forward; // Configure forward search
  approximate_search(forward_asearch,NULL);
  if (archive_search->search_reverse) {
    // Run the search (REVERSE)
    approximate_search_t* const reverse_asearch = &archive_search->reverse_search_state;
    reverse_asearch->stop_search_stage = asearch_filtering; // Stop before filtering
    reverse_asearch->search_strand = Reverse; // Configure reverse search
    approximate_search(reverse_asearch,NULL);
  }
  PROF_STOP(GP_ARCHIVE_SEARCH_GENERATE_CANDIDATES);
}
GEM_INLINE void archive_search_copy_candidates(
    archive_search_t* const archive_search,bpm_gpu_buffer_t* const bpm_gpu_buffer) {
  PROF_START(GP_ARCHIVE_SEARCH_COPY_CANDIDATES);
  const archive_t* const archive = archive_search->archive;
  // Add candidates (FORWARD)
  approximate_search_t* const forward_asearch = &archive_search->forward_search_state;
  forward_asearch->num_potential_candidates = filtering_candidates_add_to_bpm_buffer(
      forward_asearch->filtering_candidates,archive->locator,archive->fm_index,archive->enc_text,
      &forward_asearch->pattern,forward_asearch->search_strand,forward_asearch->search_actual_parameters,bpm_gpu_buffer);
  if (archive_search->search_reverse) {
    // Add candidates (REVERSE)
    approximate_search_t* const reverse_asearch = &archive_search->reverse_search_state;
    reverse_asearch->num_potential_candidates = filtering_candidates_add_to_bpm_buffer(
        reverse_asearch->filtering_candidates,archive->locator,archive->fm_index,archive->enc_text,
        &reverse_asearch->pattern,reverse_asearch->search_strand,reverse_asearch->search_actual_parameters,bpm_gpu_buffer);
  }
  PROF_STOP(GP_ARCHIVE_SEARCH_COPY_CANDIDATES);
}
GEM_INLINE void archive_search_retrieve_candidates(
    archive_search_t* const archive_search,bpm_gpu_buffer_t* const bpm_gpu_buffer,
    const uint64_t results_buffer_offset,matches_t* const matches) {
  PROF_START(GP_ARCHIVE_SEARCH_RETRIEVE_CANDIDATES);
  // Verified candidates (FORWARD)
  approximate_search_t* const forward_asearch = &archive_search->forward_search_state;
  const uint64_t results_buffer_forward_top = results_buffer_offset+forward_asearch->num_potential_candidates;
  filtering_candidates_verify_from_bpm_buffer(
      archive_search->text_collection,archive_search->archive->enc_text,&forward_asearch->pattern,Forward,
      bpm_gpu_buffer,results_buffer_offset,results_buffer_forward_top,matches,archive_search->mm_stack);
  if (archive_search->search_reverse) {
    // Verified candidates (REVERSE)
    approximate_search_t* const reverse_asearch = &archive_search->reverse_search_state;
    const uint64_t results_buffer_reverse_top = results_buffer_forward_top+reverse_asearch->num_potential_candidates;
    filtering_candidates_verify_from_bpm_buffer(
        archive_search->text_collection,archive_search->archive->enc_text,&reverse_asearch->pattern,Reverse,
        bpm_gpu_buffer,results_buffer_forward_top,results_buffer_reverse_top,matches,archive_search->mm_stack);
  }
  PROF_STOP(GP_ARCHIVE_SEARCH_RETRIEVE_CANDIDATES);
}
