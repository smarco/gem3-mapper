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
#include "mapper/mapper_bisulfite.h"
#include "archive/search/archive_search_se.h"
#include "archive/search/archive_search_pe.h"
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
      fprintf(stream,"GEM::Running-Thread (threadID = %lu)\n",thread_id);
      fprintf(stream,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
      mapper_error_report_cmd(stream,mapper_parameters); // Display CMD used
      fprintf(stream,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
      // Display Input State
      if (!mapper_search->paired_end) {
        const sequence_t* const sequence = archive_search_get_sequence(mapper_search->archive_search);
        mapper_error_report_input_state(stream,mapper_search->buffered_fasta_input,sequence);
      } else {
        const sequence_t* const sequence_end1 = archive_search_get_sequence(mapper_search->archive_search_end1);
        const sequence_t* const sequence_end2 = archive_search_get_sequence(mapper_search->archive_search_end2);
        mapper_error_report_input_state(stream,mapper_search->buffered_fasta_input_end1,sequence_end1);
        mapper_error_report_input_state(stream,mapper_search->buffered_fasta_input_end2,sequence_end2);
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

  // Create new buffered reader/writer
  mapper_parameters_t* const parameters = mapper_search->mapper_parameters;
  mapper_se_prepare_io_buffers(parameters,parameters->io.input_buffer_size,
      &mapper_search->buffered_fasta_input,&mapper_search->buffered_output_file);

  // Create an Archive-Search
  archive_t* const archive = parameters->archive;
  const bool bisulfite_index = (archive->type == archive_dna_bisulfite); // Bisulfite
  search_parameters_t* const search_parameters = &mapper_search->mapper_parameters->search_parameters;
  archive_search_handlers_t* const archive_search_handlers = archive_search_handlers_new();
  archive_search_se_new(archive,search_parameters,false,NULL,&mapper_search->archive_search);
  archive_search_t* const archive_search = mapper_search->archive_search;
  archive_search_handlers_inject_se(archive_search,archive_search_handlers);
  mapper_search->matches = matches_new();
  matches_t* const matches = mapper_search->matches;

  // FASTA/FASTQ reading loop
  uint64_t reads_processed = 0;
  while (mapper_se_read_single_sequence(mapper_search->archive_search,mapper_search->buffered_fasta_input)) {
//    if (gem_streq(mapper_search->archive_search->sequence.tag.buffer,"H.Sapiens.1M.Illumina.l100.low.000001140/1")) {
//      printf("HERE\n");
//    }

    // Bisulfite: Fully convert reads before searching into archive, making a copy of the original
    if (bisulfite_index) mapper_bisulfite_process_sequence_se(archive_search,search_parameters);

    // Search into the archive
#ifdef DEBUG_MAPPER_DISPLAY_EACH_READ_TIME
    gem_timer_t timer;
    TIMER_RESTART(&timer);
    archive_search_se(archive_search,matches);
    TIMER_STOP(&timer);
    fprintf(stderr,"Done %s in %2.4f ms.\n",archive_search->sequence.tag.buffer,TIMER_GET_TOTAL_MS(&timer));
#else
    archive_search_se(archive_search,matches);
#endif

    // Bisulfite: Copy back original read
    if (bisulfite_index) mapper_bisulfite_restore_sequence_se(archive_search);

    // Output matches
    mapper_se_output_matches(parameters,mapper_search->buffered_output_file,
        archive_search,matches,mapper_search->mapping_stats);
    // output_fastq(mapper_search->buffered_output_file,&mapper_search->archive_search->sequence);

    // Update processed
    if (++reads_processed == parameters->io.mapper_ticker_step) {
      ticker_update_mutex(mapper_search->ticker,reads_processed);
      reads_processed=0;
    }

    // Clear
    archive_search_handlers_clear(archive_search_handlers);
    matches_clear(matches);
  }
  // Update processed
  ticker_update_mutex(mapper_search->ticker,reads_processed);

  // Clean up
  buffered_input_file_close(mapper_search->buffered_fasta_input);
  buffered_output_file_close(mapper_search->buffered_output_file);
  archive_search_delete(archive_search);
  matches_delete(matches);
  archive_search_handlers_delete(archive_search_handlers);

  pthread_exit(0);
}
/*
 * PE Mapper
 */
void* mapper_pe_thread(mapper_search_t* const mapper_search) {
  // GEM-thread error handler
  gem_thread_register_id(mapper_search->thread_id+1);

  // Create new buffered reader/writer
  mapper_parameters_t* const parameters = mapper_search->mapper_parameters;
  mapper_pe_prepare_io_buffers(
      parameters,parameters->io.input_buffer_size,&mapper_search->buffered_fasta_input_end1,
      &mapper_search->buffered_fasta_input_end2,&mapper_search->buffered_output_file);

  // Create an Archive-Search
  archive_t* const archive = parameters->archive;
  const bool bisulfite_index = (archive->type == archive_dna_bisulfite);
  search_parameters_t* const search_parameters = &mapper_search->mapper_parameters->search_parameters;
  mapper_search->archive_search_handlers = archive_search_handlers_new();
  archive_search_pe_new(parameters->archive,search_parameters,false,NULL,
      &mapper_search->archive_search_end1,&mapper_search->archive_search_end2);
  archive_search_handlers_inject_pe(mapper_search->archive_search_end1,
      mapper_search->archive_search_end2,mapper_search->archive_search_handlers);
  mapper_search->paired_matches = paired_matches_new();
  archive_search_t* const archive_search_end1 = mapper_search->archive_search_end1;
  archive_search_t* const archive_search_end2 = mapper_search->archive_search_end2;
  paired_matches_t* const paired_matches = mapper_search->paired_matches;

  // FASTA/FASTQ reading loop
  uint64_t reads_processed = 0;
  while (mapper_pe_read_paired_sequences(mapper_search->mapper_parameters,
      mapper_search->buffered_fasta_input_end1,mapper_search->buffered_fasta_input_end2,
      mapper_search->archive_search_handlers,mapper_search->archive_search_end1,
      mapper_search->archive_search_end2)) {
//    if (gem_streq(mapper_search->archive_search_end1->sequence.tag.buffer,"H.Sapiens.1M.Illumina.l100.low.000001140/1")) {
//      printf("HERE\n");
//    }

    // Bisulfite: Fully convert reads before searching into archive, making a copy of the original
    if (bisulfite_index) mapper_bisulfite_process_sequence_pe(archive_search_end1,archive_search_end2);

    // Search into the archive
#ifdef DEBUG_MAPPER_DISPLAY_EACH_READ_TIME
    gem_timer_t timer;
    TIMER_RESTART(&timer);
    archive_search_pe(archive_search_end1,archive_search_end2,paired_matches);
    TIMER_STOP(&timer);
    fprintf(stderr,"Done %s in %2.4f ms.\n",
      archive_search_end1->sequence.tag.buffer,TIMER_GET_TOTAL_MS(&timer));
#else
    archive_search_pe(archive_search_end1,archive_search_end2,paired_matches);
#endif

    // Bisulfite: Copy back original read
    if (bisulfite_index) mapper_bisulfite_restore_sequence_pe(archive_search_end1,archive_search_end2);

    // Output matches
    mapper_pe_output_matches(parameters,mapper_search->buffered_output_file,
        archive_search_end1,archive_search_end2,paired_matches,mapper_search->mapping_stats);
    //output_fastq(mapper_search->buffered_output_file,&archive_search_end1->sequence);
    //output_fastq(mapper_search->buffered_output_file,&archive_search_end2->sequence);

    // Update processed
    if (++reads_processed == parameters->io.mapper_ticker_step) {
      ticker_update_mutex(mapper_search->ticker,reads_processed);
      reads_processed=0;
    }

    // Clear
    archive_search_handlers_clear(mapper_search->archive_search_handlers);
    paired_matches_clear(paired_matches,true);
  }
  // Update processed
  ticker_update_mutex(mapper_search->ticker,reads_processed);

  // Clean up
  buffered_input_file_close(mapper_search->buffered_fasta_input_end1);
  if (parameters->io.separated_input_files) {
    buffered_input_file_close(mapper_search->buffered_fasta_input_end2);
  }
  buffered_output_file_close(mapper_search->buffered_output_file);
  archive_search_delete(mapper_search->archive_search_end1);
  archive_search_delete(mapper_search->archive_search_end2);
  paired_matches_delete(mapper_search->paired_matches);
  archive_search_handlers_delete(mapper_search->archive_search_handlers);

  pthread_exit(0);
}
/*
 * SE/PE runnable
 */
void mapper_run(mapper_parameters_t* const mapper_parameters,const bool paired_end) {
  // Load GEM-Index
  mapper_load_index(mapper_parameters);
  gem_cond_fatal_error_msg(paired_end && mapper_parameters->archive->text->run_length,
      "Archive.RL-Text not supported for Paired-End Mode (yet...)");
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
    mapper_search[i].paired_end = paired_end;
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
    gem_cond_fatal_error__perror(
        pthread_create(mapper_search[i].thread_data,0,
            mapper_thread,(void*)(mapper_search+i)),SYS_THREAD_CREATE);
  }
  // Join all threads
  for (i=0;i<num_threads;++i) {
    gem_cond_fatal_error__perror(pthread_join(*(mapper_search[i].thread_data),0),SYS_THREAD_JOIN);
    mm_free(mapper_search[i].thread_data);
  }
  PROFILE_VTUNE_STOP(); // Vtune
  ticker_finish(&ticker);
  ticker_mutex_cleanup(&ticker);
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
