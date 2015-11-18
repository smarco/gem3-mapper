/*
 * PROJECT: GEMMapper
 * FILE: mapper_cuda.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "mapper_cuda.h"
#include "mapper_cuda_se.h"
#include "mapper_cuda_pe.h"
#include "archive_search.h"
#include "search_pipeline.h"
#include "gpu_buffer_collection.h"
#include "report_stats.h"

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
//      fprintf(stream,"GEM::Running-Thread (threadID = %"PRIu64")\n",mapper_cuda_search->thread_id);
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
 * SE/PE runnable
 */
GEM_INLINE void mapper_cuda_prepare_ticker(
    ticker_t* const ticker,const bool paired_end,const bool verbose_user) {
  ticker_count_reset(ticker,verbose_user,
      paired_end ? "PE::Mapping Sequences" : "SE::Mapping Sequences",0,MAPPER_TICKER_STEP,false);
  ticker_add_process_label(ticker,"#","sequences processed");
  ticker_add_finish_label(ticker,"Total","sequences processed");
  ticker_mutex_enable(ticker);
}
GEM_INLINE void mapper_cuda_setup_thread(
    mapper_cuda_search_t* const mapper_search,const uint64_t thread_id,
    mapper_parameters_t* const mapper_parameters,gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t gpu_buffers_offset,mapping_stats_t* const mstats,ticker_t* const ticker) {
  // Setup Thread
//  mapper_search[i].paired_end = paired_end;
  mapper_search->thread_id = thread_id;
  mapper_search->thread_data = mm_alloc(pthread_t);
  mapper_search->mapper_parameters = mapper_parameters;
  mapper_search->gpu_buffer_collection = gpu_buffer_collection;
  mapper_search->gpu_buffers_offset = gpu_buffers_offset;
  if (mstats) {
    mapper_search->mapping_stats = mstats + thread_id;
    init_mapping_stats(mapper_search->mapping_stats);
  } else {
    mapper_search->mapping_stats = NULL;
  }
  mapper_search->ticker = ticker;
}
GEM_INLINE void mapper_cuda_run(mapper_parameters_t* const mapper_parameters,const bool paired_end) {
  // Parameters
  mapper_parameters_cuda_t* const cuda_parameters = &mapper_parameters->cuda;
  // Check CUDA-Support
  if (!gpu_supported()) GEM_CUDA_NOT_SUPPORTED();
  // Load GEM-Index
  mapper_load_index(mapper_parameters);
  // Prepare GPU Buffers
  const uint64_t num_threads = mapper_parameters->system.num_threads;
  const uint64_t num_gpu_buffers_per_thread =
      cuda_parameters->num_fmi_bsearch_buffers +
      cuda_parameters->num_fmi_decode_buffers +
      cuda_parameters->num_bpm_buffers;
  const uint64_t num_gpu_buffers = num_gpu_buffers_per_thread * num_threads;
  gpu_buffer_collection_t* const gpu_buffer_collection = gpu_buffer_collection_new(mapper_parameters->archive,
      num_gpu_buffers,cuda_parameters->gpu_buffer_size,mapper_parameters->misc.verbose_dev);
  // I/O (SAM headers)
  archive_t* const archive = mapper_parameters->archive;
  const bool bisulfite_index = (archive->type == archive_dna_bisulfite);
  if (mapper_parameters->io.output_format==SAM) {
    output_sam_print_header(mapper_parameters->output_file,archive,
        &mapper_parameters->io.sam_parameters,mapper_parameters->argc,mapper_parameters->argv);
    mapper_parameters->io.sam_parameters.bisulfite_output = bisulfite_index;
  }
  // Ticker
  ticker_t ticker;
  mapper_cuda_prepare_ticker(&ticker,paired_end,mapper_parameters->misc.verbose_user);
  // Mapping stats
  mapping_stats_t* const global_mapping_stats = mapper_parameters->global_mapping_stats;
  mapping_stats_t* const mstats = global_mapping_stats ? mm_calloc(num_threads,mapping_stats_t,false) : NULL;
  // Mapper Searches (threads)
  mapper_cuda_search_t* const mapper_search = mm_malloc(num_threads*sizeof(mapper_cuda_search_t));
  // Error-report function
  g_mapper_cuda_search = mapper_search;
  gem_error_set_report_function(mapper_cuda_error_report);
  /*
   * Launch threads
   */
  uint64_t i, gpu_buffers_offset = 0;
  for (i=0;i<num_threads;++i) {
    // Setup Thread
  	mapper_cuda_setup_thread(mapper_search+i,i,mapper_parameters,
  	    gpu_buffer_collection,gpu_buffers_offset,mstats,&ticker);
  	gpu_buffers_offset += num_gpu_buffers_per_thread;
    // Launch thread
    gem_cond_fatal_error__perror(
        pthread_create(mapper_search[i].thread_data,0,
            (void* (*)(void*) )(paired_end ? mapper_cuda_pe_thread : mapper_cuda_se_thread),
            (void* )(mapper_search + i)),SYS_THREAD_CREATE);
  }
  // Join all threads
  for (i=0;i<num_threads;++i) {
    gem_cond_fatal_error__perror(pthread_join(*(mapper_search[i].thread_data),0),SYS_THREAD_JOIN);
    search_pipeline_delete(mapper_search[i].search_pipeline);
    mm_free(mapper_search[i].thread_data);
  }
  // Clean up
  ticker_finish(&ticker);
  ticker_mutex_cleanup(&ticker);
	// Merge report stats
	if (mstats) {
    merge_mapping_stats(global_mapping_stats,mstats,num_threads);
    mm_free(mstats);
	}
  mm_free(mapper_search); // Delete mapper-CUDA searches
  gpu_buffer_collection_delete(gpu_buffer_collection); // Delete GPU-buffer collection
}
/*
 * SE-CUDA runnable
 */
GEM_INLINE void mapper_cuda_se_run(mapper_parameters_t* const mapper_parameters) {
  mapper_cuda_run(mapper_parameters,false);
}
/*
 * PE-CUDA Mapper
 */
GEM_INLINE void mapper_cuda_pe_run(mapper_parameters_t* const mapper_parameters) {
  mapper_cuda_run(mapper_parameters,true);
}
