/*
 * PROJECT: GEMMapper
 * FILE: mapper_cuda.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "mapper_cuda.h"

/*
 * Mapper-CUDA Search
 */
typedef struct {
  /* Thread Info */
  uint64_t thread_id;
  pthread_t* thread_data;
  /* I/O */
  bool paired_end;
  buffered_input_file_t* buffered_fasta_input;
  buffered_input_file_t* buffered_fasta_input_end1;
  buffered_input_file_t* buffered_fasta_input_end2;
  /* Mapper parameters */
  mapper_parameters_t* mapper_parameters;
  /* Archive-Search group */
  archive_search_group_t* search_group;
  /* Ticker */
  ticker_t* ticker;
} mapper_cuda_search_t;

/*
 * Error report
 */
mapper_cuda_search_t* g_mapper_cuda_search; // Global searches on going
pthread_mutex_t mapper_cuda_error_report_mutex = PTHREAD_MUTEX_INITIALIZER;
void mapper_cuda_error_report(FILE* stream) {
//  // Display thread info
//  const uint64_t threads_id = gem_thread_get_thread_id();
//  if (threads_id==0) {
//    fprintf(stream,"GEM::Running-Thread (threadID = MASTER)\n");
//  }
//  // Display Threads-Info
//  MUTEX_BEGIN_SECTION(mapper_cuda_error_report_mutex) {
//    const uint64_t num_threads = g_mapper_cuda_search->mapper_parameters->system.num_threads;
//    uint64_t i;
//    for (i=0;i<num_threads;++i) {
//      mapper_cuda_search_t* const mapper_cuda_search = g_mapper_cuda_search + i; // Thread
//      fprintf(stream,"GEM::Running-Thread (threadID = %lu)\n",mapper_cuda_search->thread_id);
//      // Display Input State
//      const sequence_t* const sequence = archive_search_get_sequence(mapper_cuda_search->);
//      tab_global_inc();
//      mapper_display_input_state(stream,mapper_cuda_search->buffered_fasta_input,sequence);
//      tab_global_dec();
//      // Display Output State (TODO?)
//      // Display Search State (TODO?)
//    }
//  // Display stats until now (if possible) (TODO?)
//  } MUTEX_END_SECTION(mapper_cuda_error_report_mutex);
}
/*
 * SE CUDA Mapper
 */
void* mapper_SE_CUDA_thread(mapper_cuda_search_t* const mapper_search) {
  // GEM-thread error handler
  gem_thread_register_id(mapper_search->thread_id+1);
  PROF_START(GP_MAPPER_CUDA_THREAD);
  // Parameters
  mapper_parameters_t* const parameters = mapper_search->mapper_parameters;
  const mapper_parameters_cuda_t* const cuda_parameters = &parameters->cuda;
  // Create new buffered reader/writer
  mapper_search->buffered_fasta_input = buffered_input_file_new(parameters->input_file,cuda_parameters->input_buffer_lines);
  buffered_output_file_t* const buffered_output_file = buffered_output_file_new(parameters->output_file);
  buffered_input_file_attach_buffered_output(mapper_search->buffered_fasta_input,buffered_output_file);
  // Create search-group & matches
  archive_search_group_t* const search_group = mapper_search->search_group;
  archive_search_group_init_bpm_buffers(search_group); // Init BPM-buffers
  text_collection_t* const text_collection = archive_search_group_get_text_collection(search_group);
  matches_t* const matches = matches_new();
  matches_configure(matches,text_collection);
  // FASTA/FASTQ reading loop
  archive_search_t* archive_search_generate = NULL;
  bpm_gpu_buffer_t* bpm_gpu_buffer;
  archive_search_t* archive_search_select;
  uint64_t reads_processed = 0;
  while (true) {
    PROF_START(GP_MAPPER_CUDA_THREAD_GENERATING);
    // Check the end_of_block (We cannot reload input-buffer until all the searches of the previous block are solved)
    if (buffered_input_file_eob(mapper_search->buffered_fasta_input) && archive_search_group_is_empty(search_group)) {
      if (buffered_input_file_reload__dump_attached(mapper_search->buffered_fasta_input)==INPUT_STATUS_EOF) {
        PROF_STOP(GP_MAPPER_CUDA_THREAD_GENERATING);
        break;
      }
    }
    /*
     * First-stage: Generation. Keep processing the current input-block
     */
    while (!buffered_input_file_eob(mapper_search->buffered_fasta_input)) {
      // Request a clean archive-search
      archive_search_generate = archive_search_group_allocate(search_group);
      // Parse Sequence
      const error_code_t error_code = input_fasta_parse_sequence(
          mapper_search->buffered_fasta_input,archive_search_get_sequence(archive_search_generate),
          parameters->io.fastq_strictly_normalized,parameters->io.fastq_try_recovery,false);
      gem_cond_fatal_error(error_code==INPUT_STATUS_FAIL,MAPPER_CUDA_ERROR_PARSING);
      // Generate Candidates (Search into the archive)
      archive_search_generate_candidates(archive_search_generate);
      // Add archive-search to group (Put candidates in buffer)
      if (!archive_search_group_add_search(search_group,archive_search_generate)) break; // Go to select-candidates
      archive_search_generate = NULL; // Last archive-search is in BPM-buffer
    }
    PROF_STOP(GP_MAPPER_CUDA_THREAD_GENERATING);
    /*
     * Second-stage: Retrieval (Select & Output)
     */
    PROF_START(GP_MAPPER_CUDA_THREAD_SELECTING);
    // Process all search-groups generated
    archive_search_group_retrieve_begin(search_group); // Start retrieving
    while (archive_search_group_get_search(search_group,&archive_search_select,&bpm_gpu_buffer)) {
      // Retrieve candidates
      matches_clear(matches); // Clear Matches
      text_collection_clear(text_collection); // Clear text-collection
      archive_search_retrieve_candidates(archive_search_select,bpm_gpu_buffer,matches);
      // Finish Search
      archive_search_finish_search(archive_search_select,matches);
      // Output matches
      mapper_SE_output_matches(parameters,buffered_output_file,archive_search_select,matches);
      // Update processed
      if (++reads_processed == MAPPER_TICKER_STEP) {
        ticker_update_mutex(mapper_search->ticker,reads_processed);
        reads_processed=0;
      }
    }
    archive_search_group_clear(search_group); // Reset search-group
    // Check if the last archive-search couldn't fit into the BPM-buffer
    if (archive_search_generate!=NULL) {
      archive_search_generate_candidates(archive_search_generate);
      archive_search_group_add_search(search_group,archive_search_generate);
      archive_search_generate = NULL;
    }
    PROF_STOP(GP_MAPPER_CUDA_THREAD_SELECTING);
  }
  // Clean up
  ticker_update_mutex(mapper_search->ticker,reads_processed); // Update processed
  buffered_input_file_close(mapper_search->buffered_fasta_input);
  buffered_output_file_close(buffered_output_file);
  matches_delete(matches);
  PROF_STOP(GP_MAPPER_CUDA_THREAD);
  pthread_exit(0);
}
/*
 * PE CUDA Mapper
 */
void* mapper_PE_CUDA_thread(mapper_cuda_search_t* const mapper_search) {
  // GEM-thread error handler
  gem_thread_register_id(mapper_search->thread_id+1);
  PROF_START(GP_MAPPER_CUDA_THREAD);
  // Parameters
  mapper_parameters_t* const parameters = mapper_search->mapper_parameters;
  const mapper_parameters_cuda_t* const cuda_parameters = &parameters->cuda;
  // Create new buffered reader/writer
  buffered_output_file_t* const buffered_output_file = buffered_output_file_new(parameters->output_file);
  mapper_PE_prepare_io_buffers(parameters,cuda_parameters->input_buffer_lines,
      &mapper_search->buffered_fasta_input_end1,&mapper_search->buffered_fasta_input_end2,buffered_output_file);
  buffered_input_file_t* const buffered_fasta_input_end1 = mapper_search->buffered_fasta_input_end1;
  buffered_input_file_t* const buffered_fasta_input_end2 = mapper_search->buffered_fasta_input_end2;
  // Create search-group, archive_search pointers & paired-end matches
  archive_search_group_t* const search_group = mapper_search->search_group;
  archive_search_group_init_bpm_buffers(search_group); // Init BPM-buffers
  text_collection_t* const text_collection = archive_search_group_get_text_collection(search_group);
  paired_matches_t* const paired_matches = paired_matches_new();
  paired_matches_configure(paired_matches,text_collection);
  archive_search_t *archive_search_generate_end1 = NULL, *archive_search_generate_end2 = NULL;
  archive_search_t *archive_search_select_end1, *archive_search_select_end2;
  bpm_gpu_buffer_t *bpm_gpu_buffer_end1, *bpm_gpu_buffer_end2;
  // FASTA/FASTQ reading loop
  error_code_t error_code = 0;
  uint64_t reads_processed = 0;
  while (true) {
    PROF_START(GP_MAPPER_CUDA_THREAD_GENERATING);
    // Check the end_of_block (We cannot reload input-buffer until all the searches of the previous block are solved)
    if (buffered_input_file_eob(mapper_search->buffered_fasta_input_end1) && archive_search_group_is_empty(search_group)) {
      error_code = mapper_PE_reload_buffers(parameters,buffered_fasta_input_end1,buffered_fasta_input_end2);
      if (error_code==INPUT_STATUS_EOF) break;
    }
    /*
     * First-stage: Generation. Keep processing the current input-block
     */
    while (!buffered_input_file_eob(mapper_search->buffered_fasta_input_end1)) {
      if (buffered_input_file_eob(mapper_search->buffered_fasta_input_end2)) {
        MAPPER_ERROR_PE_PARSE_UNSYNCH_INPUT_FILES(parameters);
      }
      // Request a clean archive-search
      archive_search_group_allocate_pe(search_group,&archive_search_generate_end1,&archive_search_generate_end2);
      // Parse Sequence
      error_code = mapper_PE_parse_paired_sequences(parameters,
          buffered_fasta_input_end1,buffered_fasta_input_end2,
          archive_search_generate_end1,archive_search_generate_end2);
      gem_cond_fatal_error(error_code==INPUT_STATUS_FAIL,MAPPER_CUDA_ERROR_PARSING);
      // Begin Search
      archive_search_pe_generate_candidates(archive_search_generate_end1,archive_search_generate_end2,paired_matches);
      // Add archive-search to group (Put candidates in buffer)
      if (!archive_search_group_add_paired_search(search_group,
          archive_search_generate_end1,archive_search_generate_end2)) break; // Go to select-candidates
      archive_search_generate_end1 = NULL; // Last archive-search is in BPM-buffer
    }
    PROF_STOP(GP_MAPPER_CUDA_THREAD_GENERATING);
    /*
     * Second-stage: Retrieval (Select & Output)
     */
    PROF_START(GP_MAPPER_CUDA_THREAD_SELECTING);
    // Start retrieving
    archive_search_group_retrieve_begin(search_group);
    // Process all search-groups generated
    while (archive_search_group_get_paired_search(search_group,
        &archive_search_select_end1,&bpm_gpu_buffer_end1,
        &archive_search_select_end2,&bpm_gpu_buffer_end2)) {
      // Retrieve candidates
      paired_matches_clear(paired_matches); // Clear Paired Matches
      text_collection_clear(text_collection); // Clear text-collection
      archive_search_retrieve_candidates(archive_search_select_end1,bpm_gpu_buffer_end1,paired_matches->matches_end1);
      archive_search_retrieve_candidates(archive_search_select_end2,bpm_gpu_buffer_end2,paired_matches->matches_end2);
      // Finish Search
      archive_search_pe_finish_search(archive_search_select_end1,archive_search_select_end2,paired_matches);
      // Output matches
      mapper_PE_output_matches(parameters,buffered_output_file,
          archive_search_select_end1,archive_search_select_end2,paired_matches);
      // Update processed
      if (++reads_processed == MAPPER_TICKER_STEP) {
        ticker_update_mutex(mapper_search->ticker,reads_processed);
        reads_processed=0;
      }
    }
    archive_search_group_clear(search_group); // Reset search-group
    // Check if the last archive-search couldn't fit into the BPM-buffer
    PROF_START(GP_MAPPER_CUDA_THREAD_RESTART_UNFIT);
    if (archive_search_generate_end1!=NULL) {
      archive_search_pe_generate_candidates(archive_search_generate_end1,archive_search_generate_end2,paired_matches);
      archive_search_group_add_paired_search(search_group,archive_search_generate_end1,archive_search_generate_end2);
      archive_search_generate_end1 = NULL;
    }
    PROF_STOP(GP_MAPPER_CUDA_THREAD_RESTART_UNFIT);
    PROF_STOP(GP_MAPPER_CUDA_THREAD_SELECTING);
  }
  // Clean up & Quit
  ticker_update_mutex(mapper_search->ticker,reads_processed); // Update processed
  buffered_input_file_close(mapper_search->buffered_fasta_input_end1);
  if (parameters->io.separated_input_files) buffered_input_file_close(mapper_search->buffered_fasta_input_end2);
  buffered_output_file_close(buffered_output_file);
  paired_matches_delete(paired_matches);
  PROF_STOP(GP_MAPPER_CUDA_THREAD);
  pthread_exit(0);
}
/*
 * SE/PE runnable
 */
GEM_INLINE void mapper_CUDA_run(mapper_parameters_t* const mapper_parameters,const bool paired_end) {
  // Check CUDA-Support & parameters compliance
  if (!bpm_gpu_support()) GEM_CUDA_NOT_SUPPORTED();
  mapper_parameters_cuda_t* const cuda_parameters = &mapper_parameters->cuda;
  const uint64_t num_threads = mapper_parameters->system.num_threads;
  const uint64_t num_search_groups_per_thread = cuda_parameters->num_search_groups_per_thread;
  const uint64_t num_search_groups = num_search_groups_per_thread * num_threads;
  // Load GEM-Index
  mapper_load_index(mapper_parameters);
  // Prepare BPM-GPU Buffer
  bpm_gpu_buffer_collection_t* const bpm_gpu_buffer_collection =
      bpm_gpu_init(mapper_parameters->archive->text,num_search_groups,cuda_parameters->bpm_buffer_size,
          mapper_parameters->hints.avg_read_length,mapper_parameters->hints.candidates_per_query,
          mapper_parameters->misc.verbose_dev);
  // Prepare output file (SAM headers)
  if (mapper_parameters->io.output_format==SAM) {
    output_sam_print_header(mapper_parameters->output_file,mapper_parameters->archive,
        &mapper_parameters->io.sam_parameters,mapper_parameters->argc,mapper_parameters->argv);
  }
  // Ticker
  ticker_t ticker;
  ticker_count_reset(&ticker,mapper_parameters->misc.verbose_user,
      paired_end ? "PE::Mapping Sequences" : "SE::Mapping Sequences",0,MAPPER_TICKER_STEP,false);
  ticker_add_process_label(&ticker,"#","sequences processed");
  ticker_add_finish_label(&ticker,"Total","sequences processed");
  ticker_mutex_enable(&ticker);
  // Prepare Mapper searches
  mapper_cuda_search_t* const mapper_search = mm_malloc(num_threads*sizeof(mapper_cuda_search_t));
  // Set error-report function
  g_mapper_cuda_search = mapper_search;
  gem_error_set_report_function(mapper_cuda_error_report);
  /*
   * Launch 'Generating Candidates' threads
   */
  bpm_gpu_buffer_t* const bpm_gpu_buffers = bpm_gpu_buffer_collection->bpm_gpu_buffers;
  uint64_t i, bpm_gpu_buffer_pos=0;
  for (i=0;i<num_threads;++i) {
    // Setup Thread
    mapper_search[i].paired_end = paired_end;
    mapper_search[i].thread_id = i;
    mapper_search[i].thread_data = mm_alloc(pthread_t);
    mapper_search[i].mapper_parameters = mapper_parameters;
    mapper_search[i].search_group = archive_search_group_new(
        mapper_parameters,bpm_gpu_buffers+bpm_gpu_buffer_pos,num_search_groups_per_thread);
    bpm_gpu_buffer_pos += num_search_groups_per_thread;
    mapper_search[i].ticker = &ticker;
    // Launch thread
    gem_cond_fatal_error__perror(
        pthread_create(mapper_search[i].thread_data,0,
            (void* (*)(void*) )(paired_end ? mapper_PE_CUDA_thread : mapper_SE_CUDA_thread),
            (void* )(mapper_search + i)),SYS_THREAD_CREATE);
  }
  // Join all threads
  for (i=0;i<num_threads;++i) {
    gem_cond_fatal_error__perror(pthread_join(*(mapper_search[i].thread_data),0),SYS_THREAD_JOIN);
    archive_search_group_delete(mapper_search[i].search_group);
    mm_free(mapper_search[i].thread_data);
  }
  // Clean up
  ticker_finish(&ticker);
  ticker_mutex_cleanup(&ticker);
  mm_free(mapper_search); // Delete mapper-CUDA searches
  bpm_gpu_destroy(bpm_gpu_buffer_collection); // Delete GPU-buffer collection
}
/*
 * SE-CUDA runnable
 */
GEM_INLINE void mapper_SE_CUDA_run(mapper_parameters_t* const mapper_parameters) {
  mapper_CUDA_run(mapper_parameters,false);
}
/*
 * PE-CUDA Mapper
 */
GEM_INLINE void mapper_PE_CUDA_run(mapper_parameters_t* const mapper_parameters) {
  mapper_CUDA_run(mapper_parameters,true);
}
