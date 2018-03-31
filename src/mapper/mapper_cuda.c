/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   Mapper-CUDA module implements the high-level mapping workflow using GPU(s)
 */

#include "mapper/mapper_cuda.h"
#include "mapper/mapper_cuda_se.h"
#include "mapper/mapper_cuda_pe.h"
#include "mapper/mapper_io.h"
#include "archive/search/archive_search.h"
#include "gpu/gpu_buffer_collection.h"
#include "stats/report_stats.h"
#include "profiler/profiler.h"

/*
 * Error report
 */
mapper_cuda_search_t* g_mapper_cuda_search; // Global searches on going
pthread_mutex_t mapper_cuda_error_report_mutex = PTHREAD_MUTEX_INITIALIZER;
void mapper_cuda_error_report(FILE* stream) {
  fprintf(stream,"GEM::Unexpected error occurred. Sorry for the inconvenience\n"
                 "     Please re-run in CPU mode to isolate the problem\n"
                 "     Feedback and bug reporting it's highly appreciated,\n"
                 "     => Please report or email (gem.mapper.dev@gmail.com)\n");
}
/*
 * Index loader
 */
void mapper_load_gpu_index(mapper_parameters_t* const parameters) {
  TIMER_START(&parameters->loading_time);
  // Parameters
  mapper_parameters_cuda_t* const cuda_parameters = &parameters->cuda;
  const uint64_t num_threads = parameters->system.num_threads;
  const uint64_t num_gpu_buffers_per_thread =
      cuda_parameters->num_fmi_bsearch_buffers +
      cuda_parameters->num_fmi_decode_buffers +
      cuda_parameters->num_kmer_filter_buffers +
      cuda_parameters->num_bpm_distance_buffers +
      cuda_parameters->num_bpm_align_buffers;
  const uint64_t num_gpu_buffers = num_gpu_buffers_per_thread * num_threads;
  // Prepare/Check index path
  char* const gpu_index_name = gem_strcat(parameters->io.index_file_name,".gpu");
  gem_cond_fatal_error_msg(!gem_access(gpu_index_name,FM_READ),"Couldn't load gpu-index '%s'",gpu_index_name);
  gem_cond_log(parameters->misc.verbose_user,"[Loading GPU Index '%s']",gpu_index_name);
  // Load GPU Index & Prepare Buffers
  parameters->gpu_buffer_collection =
      gpu_buffer_collection_new(
        gpu_index_name,num_gpu_buffers,
        cuda_parameters->gpu_buffer_size,
        cuda_parameters->gpu_devices,
        parameters->misc.verbose_dev);
  free(gpu_index_name);
  TIMER_STOP(&parameters->loading_time);
}
/*
 * SE/PE runnable
 */
void mapper_cuda_prepare_ticker(
    ticker_t* const ticker,
    const uint64_t ticker_step,
    const bool paired_end,
    const bool verbose_user) {
  ticker_count_reset(ticker,verbose_user,
      paired_end ? "PE::Mapping Sequences" : "SE::Mapping Sequences",0,ticker_step,false);
  ticker_add_process_label(ticker,"#","sequences processed");
  ticker_add_finish_label(ticker,"Total","sequences processed");
  ticker_mutex_enable(ticker);
}
void mapper_cuda_setup_thread(
    mapper_cuda_search_t* const mapper_search,
    const uint64_t thread_id,
    mapper_parameters_t* const mapper_parameters,
    const uint64_t gpu_buffers_offset,
    mapping_stats_t* const mstats,
    ticker_t* const ticker) {
  // Setup Thread
  mapper_search->thread_id = thread_id;
  mapper_search->thread_data = mm_alloc(pthread_t);
  mapper_search->mapper_parameters = mapper_parameters;
  mapper_search->gpu_buffer_collection = mapper_parameters->gpu_buffer_collection;
  mapper_search->gpu_buffers_offset = gpu_buffers_offset;
  if (mstats) {
    mapper_search->mapping_stats = mstats + thread_id;
    init_mapping_stats(mapper_search->mapping_stats);
  } else {
    mapper_search->mapping_stats = NULL;
  }
  mapper_search->ticker = ticker;
}
void mapper_cuda_run(mapper_parameters_t* const mapper_parameters,const bool paired_end) {
  // Check CUDA-Support
  if (!gpu_supported()) GEM_CUDA_NOT_SUPPORTED();
  // Parameters
  // Parameters
  mapper_parameters_cuda_t* const cuda_parameters = &mapper_parameters->cuda;
  const uint64_t num_threads = mapper_parameters->system.num_threads;
  const uint64_t num_gpu_buffers_per_thread =
      cuda_parameters->num_fmi_bsearch_buffers +
      cuda_parameters->num_fmi_decode_buffers +
      cuda_parameters->num_kmer_filter_buffers +
      cuda_parameters->num_bpm_distance_buffers +
      cuda_parameters->num_bpm_align_buffers;
  // Load GPU GEM-Index & Prepare Buffers
  mapper_load_gpu_index(mapper_parameters);
  // Load GEM-Index
  mapper_load_index(mapper_parameters);
  gem_cond_fatal_error_msg(
      mapper_parameters->archive->text->run_length,
      "Archive RL-Text not supported for CUDA (use standard index)");
  gem_cond_fatal_error_msg(
      mapper_parameters->archive->type == archive_dna_forward,
      "Archive no-complement not supported for paired-end mode (use standard index)");
  gem_cond_fatal_error_msg(
      !mapper_parameters->archive->gpu_index,
      "Archive was not built with CUDA support (please re-run gem-indexer for GPU)");
  // I/O (SAM headers)
  archive_t* const archive = mapper_parameters->archive;
  const bool bisulfite_index = (archive->type == archive_dna_bisulfite);
  if (mapper_parameters->io.output_format==SAM) {
    output_sam_print_header(
        mapper_parameters->output_file,archive,&mapper_parameters->io.sam_parameters,
        mapper_parameters->argc,mapper_parameters->argv,mapper_parameters->gem_version);
    mapper_parameters->io.sam_parameters.bisulfite_output = bisulfite_index;
  }
  // Ticker
  ticker_t ticker;
  mapper_cuda_prepare_ticker(
      &ticker,mapper_parameters->io.mapper_ticker_step,
      paired_end,mapper_parameters->misc.verbose_user);
  // Mapping stats
  mapping_stats_t* const mstats = mapper_parameters->global_mapping_stats ?
      mm_calloc(num_threads,mapping_stats_t,false) : NULL;
  // Mapper Searches (threads)
  mapper_cuda_search_t* const mapper_search = mm_malloc(num_threads*sizeof(mapper_cuda_search_t));
  // Error-report function
  g_mapper_cuda_search = mapper_search;
  gem_error_set_report_function(mapper_cuda_error_report);
  /*
   * Launch threads
   */
  PROFILE_VTUNE_START(); // Vtune
  uint64_t i, gpu_buffers_offset = 0;
  for (i=0;i<num_threads;++i) {
    // Setup Thread
  	mapper_cuda_setup_thread(mapper_search+i,i,mapper_parameters,gpu_buffers_offset,mstats,&ticker);
  	gpu_buffers_offset += num_gpu_buffers_per_thread;
    // Launch thread
    gem_cond_fatal_error(
        pthread_create(mapper_search[i].thread_data,0,
            (void* (*)(void*) )(paired_end ? mapper_cuda_pe_thread : mapper_cuda_se_thread),
            (void* )(mapper_search + i)),SYS_THREAD_CREATE);
  }
  // Join all threads
  for (i=0;i<num_threads;++i) {
    gem_cond_fatal_error(pthread_join(*(mapper_search[i].thread_data),0),SYS_THREAD_JOIN);
    mm_free(mapper_search[i].thread_data);
  }
  PROFILE_VTUNE_STOP(); // Vtune
  // Clean up
  ticker_finish(&ticker);
  ticker_mutex_cleanup(&ticker);
  // Merge report stats
  if (mstats) {
    merge_mapping_stats(mapper_parameters->global_mapping_stats,mstats,num_threads);
    mm_free(mstats);
  }
  mm_free(mapper_search); // Delete mapper-CUDA searches
  gpu_buffer_collection_delete(mapper_parameters->gpu_buffer_collection); // Delete GPU-buffer collection
}
/*
 * SE-CUDA runnable
 */
void mapper_cuda_se_run(mapper_parameters_t* const mapper_parameters) {
  mapper_cuda_run(mapper_parameters,false);
}
/*
 * PE-CUDA Mapper
 */
void mapper_cuda_pe_run(mapper_parameters_t* const mapper_parameters) {
  mapper_cuda_run(mapper_parameters,true);
}
