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
 *   Mapper module provides high-level functions to run
 *   the standard mapper workflow (SE/PE)
 */

#include "mapper/mapper.h"
#include "mapper/mapper_io.h"
#include "archive/search/archive_search_se.h"
#include "archive/search/archive_search_pe.h"
#include "text/sequence_bisulfite.h"
#include "io/output_sam.h"
#include "profiler/profiler.h"

/*
 * Profile Level
 */
#define PROFILE_LEVEL PHIGH

/*
 * Debug
 */
//#define DEBUG_MAPPER_DISPLAY_EACH_READ_TIME

/*
 * Error Report
 */
mapper_search_t* g_mapper_searches; // Global searches on going
pthread_mutex_t mapper_error_report_mutex = PTHREAD_MUTEX_INITIALIZER;
void mapper_error_report_cmd(
    FILE* stream,
    mapper_parameters_t* const mapper_parameters) {
  // Display header
  fprintf(stream,"GEM::Version %s\n",mapper_parameters->gem_version);
  fprintf(stream,"GEM::CMD");
  // Print CMD line used
  uint64_t i;
  for (i=0;i<mapper_parameters->argc;++i) {
    fprintf(stream," %s",mapper_parameters->argv[i]);
  }
  fprintf(stream,"\n");
}
void mapper_error_report_input_state(
    FILE* stream,
    buffered_input_file_t* const buffered_fasta_input,
    const sequence_t* const sequence) {
  // Display header
  fprintf(stream,"GEM::Input.State\n");
  // Check NULL
  if (sequence==NULL) { fprintf(stream,"Sequence is NULL\n"); return; }
  if (buffered_fasta_input==NULL) { fprintf(stream,"Buffered_fasta_input is NULL\n"); return; }
  // Dump FASTA/FASTQ read
  if (!string_is_null(&sequence->tag) && !string_is_null(&sequence->read)) {
    const bool has_qualities = sequence_has_qualities(sequence);
    char* const end_tag =
        (sequence->end_info == paired_end1) ? "/1" :
      ( (sequence->end_info == paired_end2) ? "/2" : " " );
    fprintf(stream,"Sequence (File '%s' Line '%"PRIu64"')\n",
        buffered_input_file_get_file_name(buffered_fasta_input),
        buffered_fasta_input->current_buffer_line_no - (has_qualities ? 4 : 2));
    if (has_qualities) {
      if (!string_is_null(&sequence->qualities)) {
        fprintf(stream,"@%"PRIs"%s\n%"PRIs"\n+\n%"PRIs"\n",
            PRIs_content(&sequence->tag),end_tag,
            PRIs_content(&sequence->read),
            PRIs_content(&sequence->qualities));
      } else {
        fprintf(stream,"@%"PRIs"%s\n%"PRIs"\n+\n<<Null Qualities>>\n",
            PRIs_content(&sequence->tag),end_tag,
            PRIs_content(&sequence->read));
      }
    } else {
      fprintf(stream,">%"PRIs"%s\n%"PRIs"\n",
          PRIs_content(&sequence->tag),end_tag,
          PRIs_content(&sequence->read));
    }
  } else {
    fprintf(stream,"Current sequence is <<Empty>>\n");
  }
}
void mapper_error_report(FILE* stream) {
  // Select thread
  const uint64_t thread_id = gem_thread_get_thread_id();
  if (thread_id==0) {
    mapper_parameters_t* const mapper_parameters = g_mapper_searches->mapper_parameters;
    MUTEX_BEGIN_SECTION(mapper_parameters->error_report_mutex) {
      fprintf(stream,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
      fprintf(stream,"GEM::Unexpected error occurred. Sorry for the inconvenience\n"
                     "     Feedback and bug reporting it's highly appreciated,\n"
                     "     => Please report or email (gem.mapper.dev@gmail.com)\n");
      fprintf(stream,"GEM::Running-Thread (threadID = MASTER)\n");
      fprintf(stream,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
      mapper_error_report_cmd(stream,mapper_parameters); // Display CMD used
    } MUTEX_END_SECTION(mapper_parameters->error_report_mutex);
  } else {
    mapper_search_t* const mapper_search = g_mapper_searches + (thread_id-1); // Thread
    mapper_parameters_t* const mapper_parameters = mapper_search->mapper_parameters;
    MUTEX_BEGIN_SECTION(mapper_parameters->error_report_mutex) {
      fprintf(stream,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
      fprintf(stream,"GEM::Unexpected error occurred. Sorry for the inconvenience\n"
                     "     Feedback and bug reporting it's highly appreciated,\n"
                     "     => Please report or email (gem.mapper.dev@gmail.com)\n");
      fprintf(stream,"GEM::Running-Thread (threadID = %"PRIu64")\n",thread_id);
      fprintf(stream,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
      mapper_error_report_cmd(stream,mapper_parameters); // Display CMD used
      fprintf(stream,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
      // Display Input State
      mapper_io_handler_t* const mapper_io_handler = mapper_search->mapper_io_handler;
      if (!mapper_search->mapper_io_handler->paired_end) {
        const sequence_t* const sequence = *mapper_search->sequence_end1;
        mapper_error_report_input_state(stream,mapper_io_handler->buffered_fasta_input_end1,sequence);
      } else {
        const sequence_t* const sequence_end1 = *mapper_search->sequence_end1;
        const sequence_t* const sequence_end2 = *mapper_search->sequence_end2;
        mapper_error_report_input_state(stream,mapper_io_handler->buffered_fasta_input_end1,sequence_end1);
        mapper_error_report_input_state(stream,mapper_io_handler->buffered_fasta_input_end2,sequence_end2);
      }
    } MUTEX_END_SECTION(mapper_parameters->error_report_mutex);
  }
}
/*
 * SE Mapper
 */
void* mapper_se_thread(mapper_search_t* const mapper_search) {
  // GEM-thread error handler
  gem_thread_register_id(mapper_search->thread_id+1);

  // Parameters
  mapper_parameters_t* const parameters = mapper_search->mapper_parameters;
  search_parameters_t* const search_parameters = &mapper_search->mapper_parameters->search_parameters;
  archive_t* const archive = parameters->archive;
  mapper_io_handler_t* mapper_io_handler;
  archive_search_handlers_t* archive_search_handlers;
  archive_search_t* archive_search;
  matches_t* matches;
  sequence_t* sequence;

  // Init search structures
  archive_search_handlers = archive_search_handlers_new(archive);
  archive_search_se_new(search_parameters,false,&archive_search);
  matches = matches_new();

  // Init I/O handler
  mapper_io_handler = mapper_io_handler_new_se(
					       parameters,parameters->io.input_buffer_size,
					       archive_search_handlers->mm_allocator);

  // Set structures pointers (DEBUG)
  mapper_search->mapper_io_handler = mapper_io_handler;
  mapper_search->sequence_end1 = &sequence;

  // Initialize bisulfite handling
  bool bisulfite_index = mapper_io_handler->archive_bisulfite;
  bool non_stranded = false;
  bisulfite_conversion_t bisulfite_conversion = C2T_conversion;
  if(bisulfite_index) {
      switch(mapper_io_handler->bisulfite_read) {
        case bisulfite_G2A:
          bisulfite_conversion = G2A_conversion;
          break;
        case bisulfite_disabled:
          bisulfite_conversion = no_conversion;
          bisulfite_index = false;
          break;
        case bisulfite_non_stranded:
          non_stranded = true;
          break;
        default:
          break;
      }
  } else bisulfite_conversion = no_conversion;
  // FASTA/FASTQ reading loop
  uint64_t reads_processed = 0;

  while (mapper_read_sequence(mapper_io_handler,true,&sequence)) {
    if(bisulfite_index) {
		 if(mapper_io_handler->bisulfite_read == bisulfite_inferred_C2T_G2A)
				bisulfite_conversion = sequence->end_info == paired_end1 ? C2T_conversion : G2A_conversion;
		 else if(mapper_io_handler->bisulfite_read == bisulfite_inferred_G2A_C2T)
				bisulfite_conversion = sequence->end_info == paired_end2 ? C2T_conversion : G2A_conversion;
		 else if(non_stranded) {
			 bisulfite_conversion = sequence_bisulfite_check_cg_depletion_se(sequence) ? G2A_conversion : C2T_conversion;
		 }
	 }
    // Prepare search
    archive_search_handlers_prepare_se(archive_search,sequence,bisulfite_conversion,archive_search_handlers);

    // Search into the archive
#ifdef DEBUG_MAPPER_DISPLAY_EACH_READ_TIME
    gem_timer_t timer;
    TIMER_RESTART(&timer); archive_search_se(archive_search,matches); TIMER_STOP(&timer);
    fprintf(stderr,"Done %s in %2.4f ms.\n",sequence->tag.buffer,TIMER_GET_TOTAL_MS(&timer));
#else
    archive_search_se(archive_search,matches);
#endif
    if(non_stranded) {

      // For the non-stranded case we re-do the search with the other conversion
      // if we haven't already found a good mapping

      bool remap = false;
      if(!matches_get_num_match_traces(matches)) remap = true;
      else {
          match_trace_t* primary_match = (matches_get_match_traces(matches))[0];
          if(primary_match->edit_distance>1) remap = true;
      }
      if(remap) {

        // Prepare sequence
        const bool run_length_pattern = archive_search->archive->text->run_length;
        archive_search->approximate_search.bisulfite_conversion = bisulfite_conversion == C2T_conversion ? G2A_conversion : C2T_conversion;
        approximate_search_prepare(&archive_search->approximate_search,run_length_pattern,sequence);
#ifdef DEBUG_MAPPER_DISPLAY_EACH_READ_TIME
        gem_timer_t timer;
        TIMER_RESTART(&timer); archive_search_se(archive_search,matches); TIMER_STOP(&timer);
        fprintf(stderr,"Done %s in %2.4f ms.\n",sequence->tag.buffer,TIMER_GET_TOTAL_MS(&timer));
#else
        archive_search_se(archive_search,matches);
#endif
      }
    }
    // Output matches
    mapper_io_handler_output_matches(mapper_io_handler,
        archive_search,matches,mapper_search->mapping_stats);

    // Update processed
    if (++reads_processed == parameters->io.mapper_ticker_step) {
      ticker_update_mutex(mapper_search->ticker,reads_processed);
      reads_processed=0;
    }

    // Clear (these are not needed => just wipe clean all mm_allocator)
    // archive_search_destroy(archive_search);
    // filtering_candidates_clear(&archive_search_handlers->filtering_candidates_end1,true);
    // Clear
    archive_search_handlers_clear(archive_search_handlers);
    matches_clear(matches);
  }
  // Update processed
  ticker_update_mutex(mapper_search->ticker,reads_processed);
  // Clean up
  matches_delete(matches);
  archive_search_delete(archive_search);
  archive_search_handlers_delete(archive_search_handlers);
  mapper_io_handler_delete(mapper_io_handler);
  pthread_exit(0);
}
/*
 * PE Mapper
 */
void* mapper_pe_thread(mapper_search_t* const mapper_search) {
  // GEM-thread error handler
  gem_thread_register_id(mapper_search->thread_id+1);

  // Parameters
  mapper_parameters_t* const parameters = mapper_search->mapper_parameters;
  search_parameters_t* const search_parameters = &mapper_search->mapper_parameters->search_parameters;
  archive_t* const archive = parameters->archive;
  mapper_io_handler_t* mapper_io_handler;
  archive_search_handlers_t* archive_search_handlers;
  archive_search_t* archive_search_end1;
  archive_search_t* archive_search_end2;
  paired_matches_t* paired_matches;
  sequence_t* sequence_end1;
  sequence_t* sequence_end2;

  // Init search structures
  archive_search_handlers = archive_search_handlers_new(archive);
  archive_search_pe_new(search_parameters,false,&archive_search_end1,&archive_search_end2);
  paired_matches = paired_matches_new();

  // Init I/O handler
  mapper_io_handler = mapper_io_handler_new_pe(
      parameters,parameters->io.input_buffer_size,
      archive_search_handlers->mapper_stats,
      archive_search_handlers->mm_allocator);

  // Initialize bisulfite handling
  bool bisulfite_index = mapper_io_handler->archive_bisulfite;
  bool non_stranded = bisulfite_index && (mapper_io_handler->bisulfite_read == bisulfite_non_stranded);
  bisulfite_conversion_t bisulfite_conversion1, bisulfite_conversion2;
  if(bisulfite_index) {
    if(mapper_io_handler->bisulfite_read == bisulfite_inferred_C2T_G2A) {
			bisulfite_conversion1 = C2T_conversion;
			bisulfite_conversion2 = G2A_conversion;
		} else {
			bisulfite_conversion2 = C2T_conversion;
			bisulfite_conversion1 = G2A_conversion;
		}
  } else {
    non_stranded = false;
    bisulfite_conversion1 = bisulfite_conversion2 = no_conversion;
  }

  // Set structures pointers (DEBUG)
  mapper_search->mapper_io_handler = mapper_io_handler;
  mapper_search->sequence_end1 = &sequence_end1;
  mapper_search->sequence_end2 = &sequence_end2;

  // FASTA/FASTQ reading loop
  uint64_t reads_processed = 0;
  while (mapper_read_paired_sequence(mapper_io_handler,true,&sequence_end1,&sequence_end2)) {
//    // DEBUG
//    if (gem_streq(sequence_end1->tag.buffer,"Sim.Illumina.l100.0000491385/1")) {
//      printf("HERE\n");
//    }
    if(non_stranded) {
      if(sequence_bisulfite_check_cg_depletion_pe(sequence_end1,sequence_end2)) {
        bisulfite_conversion1 = G2A_conversion;
        bisulfite_conversion2 = C2T_conversion;
      } else {
        bisulfite_conversion1 = C2T_conversion;
        bisulfite_conversion2 = G2A_conversion;
      }
    }
    // Prepare Search
    archive_search_handlers_prepare_pe(
       archive_search_end1,archive_search_end2,
			 sequence_end1,sequence_end2,
       bisulfite_conversion1, bisulfite_conversion2,
		archive_search_handlers);

    // Search into the archive
#ifdef DEBUG_MAPPER_DISPLAY_EACH_READ_TIME
    gem_timer_t timer;
    TIMER_RESTART(&timer); archive_search_pe(archive_search_end1,archive_search_end2,paired_matches); TIMER_STOP(&timer);
    fprintf(stderr,"Done %s in %2.4f ms.\n",sequence_end1->tag.buffer,TIMER_GET_TOTAL_MS(&timer));
#else
    archive_search_pe(archive_search_end1,archive_search_end2,paired_matches);
#endif
    if(non_stranded) {
      bool remap = false;
      if(!paired_matches_is_mapped(paired_matches)) remap = true;
      else {
          paired_map_t* primary_map = (paired_matches_get_maps(paired_matches))[0];
          if(primary_map->edit_distance>1) remap = true;
      }
      if(remap) {
        const bool run_length_pattern = archive_search_end1->archive->text->run_length;
        archive_search_end1->approximate_search.bisulfite_conversion = bisulfite_conversion2;
        archive_search_end2->approximate_search.bisulfite_conversion = bisulfite_conversion1;
        approximate_search_prepare(&archive_search_end1->approximate_search,run_length_pattern,sequence_end1);
        approximate_search_prepare(&archive_search_end2->approximate_search,run_length_pattern,sequence_end2);

        // Search into the archive
#ifdef DEBUG_MAPPER_DISPLAY_EACH_READ_TIME
        gem_timer_t timer;
        TIMER_RESTART(&timer); archive_search_pe(archive_search_end1,archive_search_end2,paired_matches); TIMER_STOP(&timer);
        fprintf(stderr,"Done %s in %2.4f ms.\n",sequence_end1->tag.buffer,TIMER_GET_TOTAL_MS(&timer));
#else
        archive_search_pe(archive_search_end1,archive_search_end2,paired_matches);
#endif
      }
    }
    // Output matches
    mapper_io_handler_output_paired_matches(
					    mapper_io_handler,archive_search_end1,archive_search_end2,
					    paired_matches,mapper_search->mapping_stats);
    //output_fastq(mapper_search->buffered_output_file,&archive_search_end1->sequence);
    //output_fastq(mapper_search->buffered_output_file,&archive_search_end2->sequence);

    // Update processed
    if (++reads_processed == parameters->io.mapper_ticker_step) {
      ticker_update_mutex(mapper_search->ticker,reads_processed);
      reads_processed=0;
    }

    // Clear (these are not needed => just wipe clean all mm_allocator)
//    archive_search_destroy(archive_search_end1);
//    archive_search_destroy(archive_search_end2);
//    filtering_candidates_clear(&archive_search_handlers->filtering_candidates_end1,true);
//    filtering_candidates_clear(&archive_search_handlers->filtering_candidates_end2,true);
    // Clear
    archive_search_handlers_clear(archive_search_handlers);
    paired_matches_clear(paired_matches,true);
  }
  // Update processed
  ticker_update_mutex(mapper_search->ticker,reads_processed);

  // Clean up
  archive_search_delete(archive_search_end1);
  archive_search_delete(archive_search_end2);
  paired_matches_delete(paired_matches);
  archive_search_handlers_delete(archive_search_handlers);
  mapper_io_handler_delete(mapper_io_handler);
  pthread_exit(0);
}
/*
 * SE/PE runnable
 */
void mapper_run(mapper_parameters_t* const mapper_parameters,const bool paired_end) {
  // Load GEM-Index
  mapper_load_index(mapper_parameters);
  gem_cond_fatal_error_msg(
      paired_end && mapper_parameters->archive->text->run_length,
      "Archive RL-text not supported for paired-end mode (use standard index)");
  gem_cond_fatal_error_msg(
      paired_end && (mapper_parameters->archive->type == archive_dna_forward),
      "Archive no-complement not supported for paired-end mode (use standard index)");
  if (mapper_parameters->archive->gpu_index) {
    gem_warn_msg(
        "Running CPU-mode using GPU-index; "
        "CPU-index can achieve better performance "
        "(gem-indexer --gpu-index=false)");
  }
  // Setup threads
  const uint64_t num_threads = mapper_parameters->system.num_threads;
  mapper_search_t* const mapper_search = mm_calloc(num_threads,mapper_search_t,false);
  // Set error-report function
  g_mapper_searches = mapper_search;
  MUTEX_INIT(mapper_parameters->error_report_mutex);
  gem_error_set_report_function(mapper_error_report);
  // Prepare output file/parameters (SAM headers)
  archive_t* const archive = mapper_parameters->archive;
  const bool bisulfite_index = (archive->type == archive_dna_bisulfite);
  if (mapper_parameters->io.output_format==SAM) {
    output_sam_print_header(
        mapper_parameters->output_file,archive,&mapper_parameters->io.sam_parameters,
        mapper_parameters->argc,mapper_parameters->argv,mapper_parameters->gem_version);
    mapper_parameters->io.sam_parameters.bisulfite_output = bisulfite_index;
  }
  // Setup Ticker
  ticker_t ticker;
  ticker_count_reset(&ticker,mapper_parameters->misc.verbose_user,
      paired_end ? "PE::Mapping Sequences" : "SE::Mapping Sequences",0,
      mapper_parameters->io.mapper_ticker_step,false);
  ticker_add_process_label(&ticker,"#","sequences processed");
  ticker_add_finish_label(&ticker,"Total","sequences processed");
  ticker_mutex_enable(&ticker);
	// Allocate per thread mapping stats
  mapping_stats_t* const mstats = mapper_parameters->global_mapping_stats ?
      mm_calloc(num_threads,mapping_stats_t,false) : NULL;
  // Launch threads
  PROF_START(GP_MAPPER_MAPPING);
  pthread_handler_t mapper_thread;
  if (paired_end) {
    mapper_thread = (pthread_handler_t) mapper_pe_thread;
  } else {
    mapper_thread = (pthread_handler_t) mapper_se_thread;
  }
  uint64_t i;
  PROFILE_VTUNE_START(); // Vtune
  for (i=0;i<num_threads;++i) {
    // Setup thread
    mapper_search[i].thread_id = i;
    mapper_search[i].thread_data = mm_alloc(pthread_t);
    mapper_search[i].mapper_parameters = mapper_parameters;
    mapper_search[i].ticker = &ticker;
		if (mstats)	{
			 mapper_search[i].mapping_stats = mstats + i;
			 init_mapping_stats(mstats + i);
		} else {
		  mapper_search[i].mapping_stats = NULL;
		}
    // Launch thread
    gem_cond_fatal_error(
        pthread_create(mapper_search[i].thread_data,0,
            mapper_thread,(void*)(mapper_search+i)),SYS_THREAD_CREATE);
  }
  // Join all threads
  for (i=0;i<num_threads;++i) {
    gem_cond_fatal_error(pthread_join(*(mapper_search[i].thread_data),0),SYS_THREAD_JOIN);
    mm_free(mapper_search[i].thread_data);
  }
  PROFILE_VTUNE_STOP(); // Vtune
  ticker_finish(&ticker);
  ticker_mutex_cleanup(&ticker);
  PROF_STOP(GP_MAPPER_MAPPING);
	// Merge report stats
	if (mstats) {
		 merge_mapping_stats(mapper_parameters->global_mapping_stats,mstats,num_threads);
		 mm_free(mstats);
	}
  // Clean up
	MUTEX_DESTROY(mapper_parameters->error_report_mutex);
  mm_free(mapper_search);
}
void mapper_se_run(mapper_parameters_t* const mapper_parameters) {
  mapper_run(mapper_parameters,false);
}
void mapper_pe_run(mapper_parameters_t* const mapper_parameters) {
  mapper_run(mapper_parameters,true);
}
