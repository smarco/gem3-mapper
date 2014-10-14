/*
 * PROJECT: GEMMapper
 * FILE: mapper_se.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "mapper_cuda.h"

/*
 * Mapper-CUDA Search
 */
typedef enum { mapper_cuda_thread_generate_candidates, mapper_cuda_thread_select_candidates } mapper_cuda_thread_type;
typedef struct {
  /* Thread Info */
  uint64_t thread_id;
  pthread_t* thread_data;
  mapper_cuda_thread_type thread_type;
  /* Mapper parameters */
  mapper_parameters_t* mapper_parameters;
  mapper_cuda_parameters_t* cuda_parameters;
  /* Archive search */
  search_group_dispatcher_t* search_group_dispatcher;
  /* Ticker */
  ticker_t* ticker;
} mapper_cuda_search_t;

/*
 * Error report
 */
mapper_cuda_search_t* g_mapper_cuda_search; // Global searches on going

void mapper_cuda_error_report(FILE* stream) {
//
//  // TODO Implement 1 only call of this function, otherwise kill
//  // TODO Patch this with crazy checkers that nothing is null
//
//  const uint64_t threads_id = gem_thread_get_thread_id();
//  if (threads_id==0) {
//    fprintf(stream,"GEM::RunnigThread (threadID = MASTER)\n");
//  } else {
//    uint64_t i;
//    // Display Generating threads
//    const uint64_t num_generating_threads = g_mapper_cuda_search->cuda_parameters->num_generating_threads;
//    for (i=0;i<num_generating_threads;++i) {
//      mapper_cuda_search_t* const mapper_search = g_mapper_cuda_search + i;
//      // Display Input State
//      mapper_display_input_state();
//    }
//    // Display Selecting threads
//    const uint64_t num_selecting_threads = g_mapper_cuda_search->cuda_parameters->num_selecting_threads;
//    for (i=0;i<num_generating_threads;++i) {
//      mapper_cuda_search_t* const mapper_search = g_mapper_cuda_search + (num_generating_threads+i);
//
//    }
//
//
//    fprintf(stream,"GEM::RunnigThread (threadID = %lu)\n",mapper_search->thread_id);
//    // Dump FASTA/FASTQ read
//    const sequence_t* const sequence = archive_search_get_sequence(mapper_search->archive_search);
//    if (!string_is_null(&sequence->tag) && !string_is_null(&sequence->read)) {
//      const bool has_qualities = sequence_has_qualities(sequence);
//      char* const end_tag =
//          (sequence->attributes.end_info == PAIRED_END1) ? "/1" :
//        ( (sequence->attributes.end_info == PAIRED_END2) ? "/2" : "" );
//      fprintf(stream,"GEM::Sequence (File '%s' Line '%lu')\n",
//          input_file_get_file_name(mapper_search->buffered_fasta_input->input_file),
//          mapper_search->buffered_fasta_input->current_line_num - (has_qualities ? 4 : 2));
//      if (has_qualities) {
//        fprintf(stream,"@%"PRIs"%s\n%"PRIs"\n+\n%"PRIs"\n",
//            PRIs_content(&sequence->tag),end_tag,
//            PRIs_content(&sequence->read),
//            PRIs_content(&sequence->qualities));
//      } else {
//        fprintf(stream,">%"PRIs"%s\n%"PRIs"\n",
//            PRIs_content(&sequence->tag),end_tag,
//            PRIs_content(&sequence->read));
//      }
//    } else {
//      fprintf(stream,"GEM::Sequence <<Empty>>\n");
//    }
//    // TODO ... More useful info
//  }
}

/*
 * CUDA Mapper parameters
 */
GEM_INLINE void mapper_cuda_parameters_set_defaults(mapper_cuda_parameters_t* const mapper_cuda_parameters) {
  const uint64_t num_processors = system_get_num_processors();
  /* I/O */
  mapper_cuda_parameters->input_block_size = BUFFER_SIZE_64M; // FIXME: Not used
  mapper_cuda_parameters->buffer_num_lines = (2*4*NUM_LINES_20K); // 2l-Paired x 4l-FASTQRecord x 5K-BufferSize // FIXME: Not used
  mapper_cuda_parameters->output_buffer_size = BUFFER_SIZE_8M; // FIXME: Not used
  mapper_cuda_parameters->max_output_buffers = 3*num_processors;  // Lazy allocation // FIXME: Not used
  /* Single-end Alignment */
  /* Paired-end Alignment */
  /* Reporting */
  /* BPM Buffers */
  mapper_cuda_parameters->num_search_groups=3*num_processors;
  mapper_cuda_parameters->average_query_size=200;
  mapper_cuda_parameters->candidates_per_query=20;
  /* System */
  mapper_cuda_parameters->num_generating_threads=num_processors;
  mapper_cuda_parameters->num_selecting_threads=num_processors;
  /* Miscellaneous */
  /* Extras */
}
/*
 * I/O
 */
GEM_INLINE error_code_t mapper_SE_CUDA_parse_sequence(
    const mapper_parameters_t* const parameters,buffered_input_file_t* const buffered_fasta_input,
    search_group_dispatcher_t* const dispatcher,search_group_t** const search_group,
    archive_search_t** const archive_search) {
  // Check the end_of_block (Reload input-buffer if needed)
  if (buffered_input_file_eob(buffered_fasta_input)) {
    // Reload input-buffer
    if (buffered_input_file_get_lines_block(buffered_fasta_input)==0) return INPUT_STATUS_EOF;
    // Return search-group
    if (*search_group) search_group_dispatcher_return_generating(dispatcher,*search_group);
    // Request a new search-group (ID=buffer_in->block_id)
    *search_group = search_group_dispatcher_request_generating(dispatcher,buffered_fasta_input->block_id);
  } else if (bpm_gpu_buffer_almost_full(search_group_get_bpm_buffer(*search_group))) {
    // Request an extension before running out of space
    PROF_INC_COUNTER(GP_SGDISPATCHER_REQUESTS_GENERATING_EXTENSION_VOLUNTEER);
    *search_group = search_group_dispatcher_request_generating_extension(dispatcher,*search_group);
  }
  // Request a clean archive-search
  *archive_search = search_group_alloc(*search_group);
  // Parse Sequence
  error_code_t error_code;
  error_code=input_fasta_parse_sequence_(
      buffered_fasta_input,archive_search_get_sequence(*archive_search),
      parameters->fastq_strictly_normalized,parameters->fastq_try_recovery);
  gem_cond_fatal_error(error_code==INPUT_STATUS_FAIL,MAPPER_CUDA_ERROR_PARSING);
  // OK
  return INPUT_STATUS_OK;
}
/*
 * Utils
 */
GEM_INLINE void mapper_SE_CUDA_calculate_query_dimensions(
    archive_search_t* const archive_search,
    uint64_t* const num_patterns,uint64_t* const total_pattern_length,uint64_t* const total_candidates) {
  // Add dimensions for archive_search_end1
  const bool index_complement = archive_is_indexed_complement(archive_search->archive);
  *num_patterns = (index_complement) ? 1 : 2;
  *total_pattern_length = (index_complement) ?
      sequence_get_length(&archive_search->sequence) : 2*sequence_get_length(&archive_search->sequence);
  *total_candidates = archive_search_get_num_potential_canditates(archive_search);
}
/*
 * SE-CUDA Mapper (sorted output)
 */
void* mapper_SE_CUDA_run_generate_candidates(mapper_cuda_search_t* const mapper_cuda_search) {
  // Archive search
  search_group_dispatcher_t* const dispatcher = mapper_cuda_search->search_group_dispatcher;

  // GEM-thread error handler
  gem_thread_register_id(mapper_cuda_search->thread_id+1);

  // Create new buffered reader/writer
  const mapper_parameters_t* const parameters = mapper_cuda_search->mapper_parameters;
  const mapper_cuda_parameters_t* const cuda_parameters = mapper_cuda_search->cuda_parameters;
  buffered_input_file_t* const buffered_fasta_input =
      buffered_input_file_new(parameters->input_file,cuda_parameters->buffer_num_lines);

  // Archive Search-Group
  search_group_t* search_group = NULL;
  archive_search_t* archive_search = NULL;

  // FASTA/FASTQ reading loop
  error_code_t error_code;
  while (true) {
    // Parse Sequence
    error_code = mapper_SE_CUDA_parse_sequence(parameters,buffered_fasta_input,dispatcher,&search_group,&archive_search);
    gem_cond_fatal_error(error_code==INPUT_STATUS_FAIL,MAPPER_CUDA_ERROR_PARSING);
    if (error_code==INPUT_STATUS_EOF) break;

    // Generate Candidates (Search into the archive)
    archive_search_generate_candidates(archive_search);

    // Calculate query dimensions
    PROF_START(GP_MAPPER_CUDA_PUT_BPM_BUFFER);
    uint64_t num_patterns, total_pattern_length, total_candidates;
    mapper_SE_CUDA_calculate_query_dimensions(archive_search,&num_patterns,&total_pattern_length,&total_candidates);
    // Check if current search fits in buffer
    const bool fits_in_buffer = bpm_gpu_buffer_fits_in_buffer(
        search_group_get_bpm_buffer(search_group),num_patterns,total_pattern_length,total_candidates);
    if (!fits_in_buffer) { // Multiple BPM-Buffer
      search_group = search_group_dispatcher_request_generating_extension(dispatcher,search_group);
      search_group_configure(search_group,archive_search);
      archive_search_generate_candidates(archive_search); // Restate the search in the new group
    }

    // Put candidates in buffer
    const uint64_t results_offset = bpm_gpu_buffer_get_num_candidates(search_group_get_bpm_buffer(search_group));
    archive_search_copy_candidates(archive_search,search_group_get_bpm_buffer(search_group));

    // Add archive-search to group
    search_group_add_search(search_group,archive_search,results_offset);
    PROF_STOP(GP_MAPPER_CUDA_PUT_BPM_BUFFER);
  }
  // Return last search-group
  if (search_group) search_group_dispatcher_return_generating(dispatcher,search_group);
  // Deregister thread
  search_group_dispatcher_deregister_generating(dispatcher,1);
  // Clean up & exit
  buffered_input_file_close(buffered_fasta_input);
  pthread_exit(0);
}
GEM_INLINE void mapper_SE_CUDA_select_candidates(
    const mapper_parameters_t* const parameters,archive_search_t* const archive_search,
    bpm_gpu_buffer_t* const bpm_gpu_buffer,const uint64_t results_buffer_offset,
    matches_t* const matches,buffered_output_file_t* const buffered_output_file,
    ticker_t* const ticker,uint64_t* const reads_processed) {
  // Init
  matches_configure(matches,archive_search->text_collection); // TODO Factorice
  matches_clear(matches);
  // Select candidates
  archive_search_retrieve_candidates(
      archive_search,bpm_gpu_buffer,results_buffer_offset,matches);
  // Select matches
  archive_select_matches(archive_search,matches);
  // archive_search_score_matches(archive_search); // TODO
  // Output matches
  mapper_SE_output_matches(
      parameters,buffered_output_file,
      archive_search_get_sequence(archive_search),matches);
  // Update processed
  if (++(*reads_processed) == MAPPER_TICKER_STEP) {
    ticker_update_mutex(ticker,*reads_processed);
    *reads_processed=0;
  }
}
void* mapper_SE_CUDA_run_select_candidates(mapper_cuda_search_t* const mapper_cuda_search) {
  // Archive search
  search_group_dispatcher_t* const dispatcher = mapper_cuda_search->search_group_dispatcher;

  // GEM-thread error handler
  gem_thread_register_id(mapper_cuda_search->thread_id+1);

  // Create new buffered writer
  const mapper_parameters_t* const parameters = mapper_cuda_search->mapper_parameters;
  buffered_output_file_t* const buffered_output_file = buffered_output_file_new(parameters->output_file);
  matches_t* const matches = matches_new();

  // Archive search-group loop
  uint64_t i, reads_processed = 0, results_buffer_offset=0;
  search_group_t* search_group = NULL;
  archive_search_t* archive_search = NULL;
  while ((search_group=search_group_dispatcher_request_selecting(dispatcher))!=NULL) {
    // Request an output buffer
    buffered_output_file_request_buffer(buffered_output_file,search_group_get_group_id(search_group));
    // Iterate over all search-groups (over a possible multisearch-group)
    bool search_group_incomplete;
    do {
      const uint64_t num_searches = search_group_get_num_searches(search_group);
      for (i=0;i<num_searches;++i) {
        // Get next archive_search
        search_group_get_search(search_group,i,&archive_search,&results_buffer_offset);
        // Select candidates
        mapper_SE_CUDA_select_candidates(
            parameters,archive_search,search_group_get_bpm_buffer(search_group),results_buffer_offset,
            matches,buffered_output_file,mapper_cuda_search->ticker,&reads_processed);
        // Free archive_search
        search_group_release(search_group,archive_search);
      }
      // Return search-group
      search_group_incomplete = search_group_is_incomplete(search_group);
      if (search_group_incomplete) {
        search_group = search_group_dispatcher_request_selecting_next(dispatcher,search_group);
      } else {
        search_group_dispatcher_return_selecting(dispatcher,search_group);
      }
    } while (search_group_incomplete);
    // Dump buffer
    buffered_output_file_dump_buffer(buffered_output_file);
  }

  // Update processed
  ticker_update_mutex(mapper_cuda_search->ticker,reads_processed);
  // Clean & exit
  buffered_output_file_close(buffered_output_file);
  matches_delete(matches);
  pthread_exit(0);
}
/*
 * SE-CUDA Mapper (Main RUN)
 */
GEM_INLINE void mapper_SE_CUDA_run(
    mapper_parameters_t* const mapper_parameters,mapper_cuda_parameters_t* const cuda_parameters) {
  // Check CUDA-Support & parameters compliance
  if (!bpm_gpu_support()) GEM_CUDA_NOT_SUPPORTED();
  const uint64_t total_threads =
      cuda_parameters->num_generating_threads + cuda_parameters->num_selecting_threads;
  if (cuda_parameters->num_search_groups < cuda_parameters->num_generating_threads+1) {
    cuda_parameters->num_search_groups = cuda_parameters->num_generating_threads+1; // HINT on performance
  }
  // Prepare archive-search group dispatcher & cache
  search_group_dispatcher_t* const search_group_dispatcher = search_group_dispatcher_new(
      mapper_parameters,mapper_parameters->archive,cuda_parameters->num_search_groups,
      cuda_parameters->average_query_size,cuda_parameters->candidates_per_query);
  // Prepare output file (SAM headers)
  if (mapper_parameters->output_format==SAM) {
    output_sam_print_header(mapper_parameters->output_file,
        mapper_parameters->archive,mapper_parameters->argc,mapper_parameters->argv);
  }
  // Ticker
  ticker_t ticker;
  ticker_count_reset(&ticker,mapper_parameters->verbose_user,"Mapping Sequences",0,MAPPER_TICKER_STEP,true);
  ticker_add_process_label(&ticker,"#","sequences processed");
  ticker_add_finish_label(&ticker,"Total","sequences processed");
  ticker_mutex_enable(&ticker);
  // Prepare Mapper-CUDA searches
  mapper_cuda_search_t* const mapper_cuda_search =
      mm_malloc(total_threads*sizeof(mapper_cuda_search_t)); // Allocate mapper-CUDA searches
  // Set error-report function
  g_mapper_cuda_search = mapper_cuda_search;
  gem_error_set_report_function(mapper_cuda_error_report);
  /*
   * Launch 'Generating Candidates' threads
   */
  uint64_t i;
  search_group_dispatcher_register_generating(
      search_group_dispatcher,cuda_parameters->num_generating_threads);
  for (i=0;i<cuda_parameters->num_generating_threads;++i) {
    // Setup Thread
    mapper_cuda_search[i].thread_id = i;
    mapper_cuda_search[i].thread_data = mm_alloc(pthread_t);
    mapper_cuda_search[i].thread_type = mapper_cuda_thread_generate_candidates;
    mapper_cuda_search[i].mapper_parameters = mapper_parameters;
    mapper_cuda_search[i].cuda_parameters = cuda_parameters;
    mapper_cuda_search[i].search_group_dispatcher = search_group_dispatcher;
    // Launch thread
    gem_cond_fatal_error__perror(
        pthread_create(mapper_cuda_search[i].thread_data,0,
            (void* (*)(void*) )mapper_SE_CUDA_run_generate_candidates,
            (void* )(mapper_cuda_search + i)),SYS_THREAD_CREATE);
  }
  /*
   * Launch 'Selecting Candidates' threads
   */
  for (i=0;i<cuda_parameters->num_selecting_threads;++i) {
    // Setup Thread
    const uint64_t thread_idx = cuda_parameters->num_generating_threads + i;
    mapper_cuda_search[thread_idx].thread_id = thread_idx;
    mapper_cuda_search[thread_idx].thread_data = mm_alloc(pthread_t);
    mapper_cuda_search[thread_idx].thread_type = mapper_cuda_thread_select_candidates;
    mapper_cuda_search[thread_idx].mapper_parameters = mapper_parameters;
    mapper_cuda_search[thread_idx].cuda_parameters = cuda_parameters;
    mapper_cuda_search[thread_idx].search_group_dispatcher = search_group_dispatcher;
    mapper_cuda_search[thread_idx].ticker = &ticker;
    // Launch thread
    gem_cond_fatal_error__perror(
        pthread_create(mapper_cuda_search[thread_idx].thread_data,0,
            (void* (*)(void*) )mapper_SE_CUDA_run_select_candidates,
            (void* )(mapper_cuda_search + thread_idx)),SYS_THREAD_CREATE);
  }
  // Join all threads
  for (i=0;i<total_threads;++i) {
    gem_cond_fatal_error__perror(pthread_join(*(mapper_cuda_search[i].thread_data),0),SYS_THREAD_JOIN);
    mm_free(mapper_cuda_search[i].thread_data);
  }
  // Clean up
  ticker_finish(&ticker);
  ticker_mutex_cleanup(&ticker);
  search_group_dispatcher_delete(search_group_dispatcher); // Delete dispatcher
  mm_free(mapper_cuda_search); // Delete mapper-CUDA searches
}
