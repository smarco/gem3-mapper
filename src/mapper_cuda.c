/*
 * PROJECT: GEMMapper
 * FILE: mapper_se.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "mapper.h"
#include "mapper_cuda.h"


/*
 * Constants
 */
#define MAPPER_CUDA_NUM_LINES (2*4*NUM_LINES_5K)
#define MAPPER_CUDA_ARCHIVE_SEARCH_CACHE_INIT_SIZE 1000

/*
 * Archive Search Cache
 */
typedef struct {
  /* Mapper parameters */
  mapper_parameters_t* mapper_parameters;
  search_parameters_t search_parameters;
  select_parameters_t select_parameters;
  /* Slab of archive_search_t */
  vector_t* archive_search_cache;  // Already allocated & configured (archive_search_t*)
  pthread_mutex_t mutex;           // Mutex to access the cache
} archive_search_cache_t;

/*
 * Mapper-CUDA Search
 */
typedef enum { mapper_cuda_thread_generate_candidates, mapper_cuda_thread_select_candidates } mapper_cuda_thread_type;
typedef struct {
  // Thread Info
  uint64_t thread_id;
  pthread_t* thread_data;
  mapper_cuda_thread_type thread_type;
  /* Mapper parameters */
  mapper_parameters_t* mapper_parameters;
  mapper_cuda_parameters_t* mapper_cuda_parameters;
  /* Archive search */
  archive_search_cache_t* archive_search_cache;
  archive_search_group_dispatcher_t* search_group_dispatcher;
} mapper_cuda_search_t;
mapper_cuda_search_t* g_mapper_cuda_search; // Global searches on going

/*
 * CUDA Mapper parameters
 */
GEM_INLINE void mapper_cuda_parameters_set_defaults(mapper_cuda_parameters_t* const mapper_cuda_parameters) {
  /* I/O */
  mapper_cuda_parameters->output_file_type = UNSORTED_FILE; // FIXME to remove
  /* Single-end Alignment */
  /* Paired-end Alignment */
  /* Reporting */
  /* System */
  mapper_cuda_parameters->num_generating_threads = 4; // Total number of threads generating candidates
  mapper_cuda_parameters->num_selecting_threads = 4;  // TOtal number of threads selecting candidates
  /* Miscellaneous */
  /* Extras */
}
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
  // Configure Parameters
  mapper_configure_archive_search(mapper_parameters,
      &archive_search_cache->search_parameters,&archive_search_cache->select_parameters);
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
GEM_INLINE archive_search_t* archive_search_cache_alloc(archive_search_cache_t* const archive_search_cache) {
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
          &archive_search_cache->search_parameters,&archive_search_cache->select_parameters);
    }
  } MUTEX_END_SECTION(archive_search_cache->mutex);
  // Return
  archive_search_clear(archive_search);
  return archive_search;
}
GEM_INLINE void archive_search_cache_free(
    archive_search_cache_t* const archive_search_cache,archive_search_t* const archive_search) {
  // Clear archive search to save some space
  archive_search_clear(archive_search);
  // Add it to the cache
  MUTEX_BEGIN_SECTION(archive_search_cache->mutex) {
    vector_insert(archive_search_cache->archive_search_cache,archive_search,archive_search_t*);
  } MUTEX_END_SECTION(archive_search_cache->mutex);
}
/*
 * Utils
 */
GEM_INLINE void mapper_SE_CUDA_calculate_query_dimensions(
    archive_search_t* const archive_search,
    uint64_t* const num_patterns,uint64_t* const total_pattern_length,uint64_t* const total_candidates) {
  const bool index_complement = archive_is_indexed_complement(archive_search->archive);
  *num_patterns += (index_complement) ? 1 : 2;
  *total_pattern_length += (index_complement) ?
      sequence_get_length(archive_search->sequence) : 2*sequence_get_length(archive_search->sequence);
  *total_candidates += archive_search_get_num_potential_canditates(archive_search);
}
GEM_INLINE void mapper_SE_CUDA_calculate_query_dimensions_2a(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    uint64_t* const num_patterns,uint64_t* const total_pattern_length,uint64_t* const total_candidates) {
  // Init
  *num_patterns = 0;
  *total_pattern_length = 0;
  *total_candidates = 0;
  // Add dimensions for archive_search_end1
  mapper_SE_CUDA_calculate_query_dimensions(
      archive_search_end1,num_patterns,total_pattern_length,total_candidates);
  // Add dimensions for archive_search_end2
  if (archive_search_end2!=NULL) {
    mapper_SE_CUDA_calculate_query_dimensions(
        archive_search_end2,num_patterns,total_pattern_length,total_candidates);
  }
}
/*
 * SE-CUDA Mapper (sorted output)
 */
void* mapper_SE_CUDA_run_generate_candidates_so(mapper_cuda_search_t* const mapper_cuda_search) {
  GEM_NOT_IMPLEMENTED(); // TODO
  pthread_exit(0);
}
void* mapper_SE_CUDA_run_select_candidates_so(mapper_cuda_search_t* const mapper_cuda_search) {
  GEM_NOT_IMPLEMENTED(); // TODO
  pthread_exit(0);
}
/*
 * SE-CUDA Mapper (Unsorted output)
 */
void* mapper_SE_CUDA_run_generate_candidates_uo(mapper_cuda_search_t* const mapper_cuda_search) {
  // Archive search
  archive_search_cache_t* const archive_search_cache = mapper_cuda_search->archive_search_cache;
  archive_search_group_dispatcher_t* const dispatcher = mapper_cuda_search->search_group_dispatcher;

  // GEM-thread error handler
  gem_thread_register_id(mapper_cuda_search->thread_id+1);

  // Create new buffered reader/writer
  const mapper_parameters_t* const parameters = mapper_cuda_search->mapper_parameters;
  buffered_input_file_t* const buffered_fasta_input = buffered_input_file_new(parameters->input_file);

  // Get archive-search group for candidate generation
  archive_search_group_t* archive_search_group = archive_search_group_dispatcher_request_generating(dispatcher);
  archive_search_t* archive_search_end1 = NULL, *archive_search_end2 = NULL;
  uint64_t results_offset_end1, results_offset_end2;

  // FASTA/FASTQ reading loop
  while (true) {
    // Check the end_of_block (Reload buffer if needed)
    if (buffered_input_file_eob(buffered_fasta_input)) {
      if (buffered_input_file_get_lines_block(buffered_fasta_input,MAPPER_CUDA_NUM_LINES)==0) break;
    }
    // Obtain a clean archive-search (end/1)
    archive_search_end1 = archive_search_cache_alloc(archive_search_cache);
    // Parse Sequence (end/1)
    error_code_t error_code;
    if ((error_code=input_fasta_parse_sequence(
        buffered_fasta_input,archive_search_get_sequence(archive_search_end1),
        parameters->fastq_strictly_normalized,parameters->fastq_try_recovery))==INPUT_STATUS_EOF) break;
    if (error_code==INPUT_STATUS_FAIL) continue; // Skip this read & Restart
    // Parse Sequence (end/2)
    if (gem_expect_false(buffered_input_file_eob(buffered_fasta_input))) {
      // Only one read left
      archive_search_cache_free(archive_search_cache,archive_search_end2);
      archive_search_end2 = NULL;
    } else {
      // Obtain a clean Archive-Search (end/2)
      archive_search_end2 = archive_search_cache_alloc(archive_search_cache);
      error_code_t error_code;
      if ((error_code=input_fasta_parse_sequence(
          buffered_fasta_input,archive_search_get_sequence(archive_search_end2),
          parameters->fastq_strictly_normalized,parameters->fastq_try_recovery))==INPUT_STATUS_EOF) break;
      if (error_code==INPUT_STATUS_FAIL) {
        archive_search_cache_free(archive_search_cache,archive_search_end2);
        archive_search_end2 = NULL;
      }
    }
    // Generate Candidates (Search into the archive)
    archive_search_generate_candidates(archive_search_end1,archive_search_group->mm_search);
    if (gem_expect_true(archive_search_end2!=NULL)) {
      archive_search_generate_candidates(archive_search_end2,archive_search_group->mm_search);
    }
    // Calculate query dimensions
    uint64_t num_patterns, total_pattern_length, total_candidates;
    mapper_SE_CUDA_calculate_query_dimensions_2a(
        archive_search_end1,archive_search_end2,&num_patterns,&total_pattern_length,&total_candidates);
    // Check if current search fits in buffer
    const bool fits_in_buffer =
        bpm_gpu_buffer_fits_in_buffer(
            archive_search_group->bpm_gpu_buffer,num_patterns,total_pattern_length,total_candidates);
    if (!fits_in_buffer) { // Return search-group (closed) & request new one
      archive_search_group_dispatcher_return_generating(dispatcher,archive_search_group);
      archive_search_group = archive_search_group_dispatcher_request_generating(dispatcher);
    }
    // Put candidates in buffer
    results_offset_end1 = bpm_gpu_buffer_get_num_candidates(archive_search_group->bpm_gpu_buffer);
    archive_search_copy_candidates(archive_search_end1,archive_search_group->bpm_gpu_buffer);
    if (gem_expect_true(archive_search_end2!=NULL)) {
      results_offset_end2 = bpm_gpu_buffer_get_num_candidates(archive_search_group->bpm_gpu_buffer);
      archive_search_copy_candidates(archive_search_end2,archive_search_group->bpm_gpu_buffer);
    }
    // Add archive-search to group
    archive_search_group_add(archive_search_group,archive_search_end1,results_offset_end1);
    if (gem_expect_true(archive_search_end2!=NULL)) {
      archive_search_group_add(archive_search_group,archive_search_end2,results_offset_end2);
    }
  }
  // Return last search-group
  archive_search_group_dispatcher_return_generating(dispatcher,archive_search_group);
  // Deregister thread
  archive_search_group_dispatcher_deregister_generating(dispatcher,1);
  // Clean up & exit
  buffered_input_file_close(buffered_fasta_input);
  pthread_exit(0);
}
void* mapper_SE_CUDA_run_select_candidates_uo(mapper_cuda_search_t* const mapper_cuda_search) {
  // Archive search
  archive_search_cache_t* const archive_search_cache = mapper_cuda_search->archive_search_cache;
  archive_search_group_dispatcher_t* const dispatcher = mapper_cuda_search->search_group_dispatcher;

  // GEM-thread error handler
  gem_thread_register_id(mapper_cuda_search->thread_id+1);

  // Create new buffered writer
  const mapper_parameters_t* const parameters = mapper_cuda_search->mapper_parameters;
  buffered_output_file_t* buffered_output_file = buffered_output_file_new(parameters->output_file);
  matches_t* const matches = matches_new();

  // Archive search-group loop
  archive_search_group_t* archive_search_group;
  while ((archive_search_group=archive_search_group_dispatcher_request_selecting(dispatcher))!=NULL) {
    // Process all archive_search from the group
    VECTOR_ITERATE_CONST(archive_search_group->archive_searches,archive_search_member,n,archive_search_member_t) {
      archive_search_t* const archive_search = archive_search_member->archive_search;
      // Select candidates
      archive_search_select_candidates(
          archive_search,archive_search_group->bpm_gpu_buffer,archive_search_member->results_buffer_offset);
      // Select matches
      archive_select_matches(archive_search,matches,archive_search_group->mm_search);

      // archive_search_score_matches(archive_search); // TODO
      // Output matches
      mapper_SE_output_matches(
          parameters,buffered_output_file,
          archive_search_get_sequence(archive_search),matches);
      // Free archive_search
      archive_search_cache_free(archive_search_cache,archive_search);
    }
    // Return search-group
    archive_search_group_dispatcher_return_selecting(dispatcher,archive_search_group);
    // Dump buffer
    buffered_output_file_dump(buffered_output_file);
  }

  // Clean up & exit
  buffered_output_file_close(buffered_output_file);
  pthread_exit(0);
}
/*
 * SE-CUDA Mapper (Main RUN)
 */
GEM_INLINE void mapper_SE_CUDA_run(
    mapper_parameters_t* const mapper_parameters,
    mapper_cuda_parameters_t* const mapper_cuda_parameters) {
  // Check CUDA-Support
  if (!bpm_gpu_support()) GEM_CUDA_NOT_SUPPORTED();
  // Setup Profiling (Add master thread)
  const uint64_t total_threads =
      mapper_cuda_parameters->num_generating_threads +
      mapper_cuda_parameters->num_selecting_threads;
  PROF_NEW(total_threads+1);
  // Prepare archive-search group dispatcher & cache
  archive_search_group_dispatcher_t* const search_group_dispatcher =
      archive_search_group_dispatcher_new(mapper_parameters->archive);
  archive_search_cache_t* const archive_search_cache = archive_search_cache_new(mapper_parameters);
  // Prepare Mapper-CUDA searches
  mapper_cuda_search_t* const mapper_cuda_search =
      mm_malloc(total_threads*sizeof(mapper_cuda_search_t)); // Allocate mapper-CUDA searches
  g_mapper_cuda_search = mapper_cuda_search; // Set global searches for error reporting
  // Launch 'Generating Candidates' threads
  void* (*generate_candidates_thread_function)(void*) =
      (mapper_cuda_parameters->output_file_type==SORTED_FILE) ?
          (void* (*)(void*) ) mapper_SE_CUDA_run_generate_candidates_so :
          (void* (*)(void*) ) mapper_SE_CUDA_run_generate_candidates_so;
  void* (*select_candidates_thread_function)(void*) =
      (mapper_cuda_parameters->output_file_type==SORTED_FILE) ?
          (void* (*)(void*) ) mapper_SE_CUDA_run_select_candidates_uo :
          (void* (*)(void*) ) mapper_SE_CUDA_run_select_candidates_uo;
  uint64_t i = 0;
  archive_search_group_dispatcher_register_generating(
      search_group_dispatcher,mapper_cuda_parameters->num_generating_threads);
  for (;i<mapper_cuda_parameters->num_generating_threads;++i) {
    // Setup Thread
    mapper_cuda_search[i].thread_id = i;
    mapper_cuda_search[i].thread_data = mm_alloc(pthread_t);
    mapper_cuda_search[i].thread_type = mapper_cuda_thread_generate_candidates;
    mapper_cuda_search[i].mapper_parameters = mapper_parameters;
    mapper_cuda_search[i].mapper_cuda_parameters = mapper_cuda_parameters;
    mapper_cuda_search[i].archive_search_cache = archive_search_cache;
    mapper_cuda_search[i].search_group_dispatcher = search_group_dispatcher;
    // Launch thread
    gem_cond_fatal_error__perror(
        pthread_create(mapper_cuda_search[i].thread_data,0,
            generate_candidates_thread_function,(void* )(mapper_cuda_search + i)),SYS_THREAD_CREATE);
  }
  // Launch 'Selecting Candidates' threads
  for (;i<mapper_cuda_parameters->num_selecting_threads;++i) {
    // Setup Thread
    mapper_cuda_search[i].thread_id = i;
    mapper_cuda_search[i].thread_data = mm_alloc(pthread_t);
    mapper_cuda_search[i].thread_type = mapper_cuda_thread_select_candidates;
    mapper_cuda_search[i].mapper_parameters = mapper_parameters;
    mapper_cuda_search[i].mapper_cuda_parameters = mapper_cuda_parameters;
    mapper_cuda_search[i].archive_search_cache = archive_search_cache;
    mapper_cuda_search[i].search_group_dispatcher = search_group_dispatcher;
    // Launch thread
    gem_cond_fatal_error__perror(
        pthread_create(mapper_cuda_search[i].thread_data,0,
            select_candidates_thread_function,(void* )(mapper_cuda_search + i)),SYS_THREAD_CREATE);
  }
  // Join all threads
  for (i=0;i<total_threads;++i) {
    gem_cond_fatal_error__perror(pthread_join(*(mapper_cuda_search[i].thread_data),0),SYS_THREAD_JOIN);
    mm_free(mapper_cuda_search[i].thread_data);
  }
  // Clean up
  archive_search_group_dispatcher_delete(search_group_dispatcher); // Delete dispatcher
  archive_search_cache_delete(archive_search_cache); // Delete cache
  mm_free(mapper_cuda_search); // Delete mapper-CUDA searches
  PROF_DELETE(); // Clean-up Profiler
}
