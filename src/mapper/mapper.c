/*
 * PROJECT: GEMMapper
 * FILE: mapper.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "mapper/mapper.h"
#include "mapper/mapper_bisulfite.h"
#include "archive/archive_search.h"
#include "archive/archive_search_se.h"
#include "archive/archive_search_pe.h"
#include "stats/report_stats.h"
//#include "/opt/intel/vtune_amplifier_xe/include/ittnotify.h"
//#include "/usr/local/software/intel/vtune_amplifier_xe_2016.3.0.463186/include/ittnotify.h"

/*
 * Debug/Profile
 */
#define MAPPER_OUTPUT

/*
 * Error Messages
 */
#define GEM_ERROR_MAPPER_PE_PARSE_UNSYNCH_INPUT_FILES_PAIR "Parsing Input Files. Files '%s,%s' doesn't contain the same number of reads (cannot pair)"
#define GEM_ERROR_MAPPER_PE_PARSE_UNSYNCH_INPUT_FILES_SINGLE "Parsing Input File. File '%s' doesn't contain a pair number of reads (cannot pair)"
#define GEM_ERROR_MAPPER_PE_PARSE_UNSYNCH_INPUT_FILES_EOF "Parsing Input Files. File '%s' could not read second end (unexpected end-of-file)"
#define GEM_ERROR_MAPPER_PE_PARSE_UNSYNCH_INPUT_FILES_NOT_EOF "Parsing Input Files. File '%s' has too many reads (expected end-of-file)"


/*
 * Error Macros
 */
#define MAPPER_ERROR_PE_PARSE_UNSYNCH_INPUT_FILES(mapper_parameters) \
  if (mapper_parameters->io.separated_input_files) { \
    gem_fatal_error(MAPPER_PE_PARSE_UNSYNCH_INPUT_FILES_PAIR, \
        input_file_sliced_get_file_name(mapper_parameters->input_file_end1), \
        input_file_sliced_get_file_name(mapper_parameters->input_file_end2)); \
  } else { \
    gem_fatal_error(MAPPER_PE_PARSE_UNSYNCH_INPUT_FILES_SINGLE, \
        input_file_sliced_get_file_name(mapper_parameters->input_file)); \
  }

/*
 * Profile Level
 */
#define PROFILE_LEVEL PHIGH

/*
 * Report
 */
void mapper_display_input_state(
    FILE* stream,
    buffered_input_file_t* const buffered_fasta_input,
    const sequence_t* const sequence) {
  // Check NULL
  if (sequence==NULL) { tab_fprintf(stream,"Sequence is NULL\n"); return; }
  if (buffered_fasta_input==NULL) { tab_fprintf(stream,"Buffered_fasta_input is NULL\n"); return; }
  // Dump FASTA/FASTQ read
  if (!string_is_null(&sequence->tag) && !string_is_null(&sequence->read)) {
    const bool has_qualities = sequence_has_qualities(sequence);
    char* const end_tag =
        (sequence->end_info == paired_end1) ? "/1" :
      ( (sequence->end_info == paired_end2) ? "/2" : " " );
    tab_fprintf(stream,"Sequence (File '%s' Line '%"PRIu64"')\n",
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
    tab_fprintf(stream,"Current sequence is <<Empty>>\n");
  }
}
/*
 * Error Report
 */
mapper_search_t* g_mapper_searches; // Global searches on going
pthread_mutex_t mapper_error_report_mutex = PTHREAD_MUTEX_INITIALIZER;
void mapper_error_report(FILE* stream) {
  // Display thread info
  const uint64_t threads_id = gem_thread_get_thread_id();
  if (threads_id==0) {
    fprintf(stream,"GEM::Running-Thread (threadID = MASTER)\n");
  }
  // Display Threads-Info
  MUTEX_BEGIN_SECTION(mapper_error_report_mutex) {
    const uint64_t num_threads = g_mapper_searches->mapper_parameters->system.num_threads;
    uint64_t i;
    for (i=0;i<num_threads;++i) {
      mapper_search_t* const mapper_search = g_mapper_searches + i; // Thread
      if (mapper_search->paired_end) {
        if (mapper_search->buffered_fasta_input == NULL) continue;
        fprintf(stream,"GEM::Running-Thread (threadID = %"PRIu64") currently processing\n",mapper_search->thread_id);
        // Display Input State
        const sequence_t* const sequence = archive_search_get_sequence(mapper_search->archive_search);
        tab_global_inc();
        mapper_display_input_state(stream,mapper_search->buffered_fasta_input,sequence);
        tab_global_dec();
      } else {
        if (mapper_search->buffered_fasta_input_end1 == NULL || mapper_search->buffered_fasta_input_end2) continue;
        fprintf(stream,"GEM::Running-Thread (threadID = %"PRIu64") currently processing\n",mapper_search->thread_id);
        // Display Input State
        const sequence_t* const sequence_end1 = archive_search_get_sequence(mapper_search->archive_search_end1);
        const sequence_t* const sequence_end2 = archive_search_get_sequence(mapper_search->archive_search_end2);
        tab_global_inc();
        mapper_display_input_state(stream,mapper_search->buffered_fasta_input_end1,sequence_end1);
        mapper_display_input_state(stream,mapper_search->buffered_fasta_input_end2,sequence_end2);
        tab_global_dec();
      }
      // Display Output State (TODO?)
      // Display Search State (TODO?)
    }
  // Display stats until now (if possible) (TODO?)
  } MUTEX_END_SECTION(mapper_error_report_mutex);
}
/*
 * Mapper Parameters
 */
void mapper_parameters_set_defaults_io(mapper_parameters_io_t* const io) {
  const uint64_t num_processors = system_get_num_processors();
  /* Input */
  io->index_file_name=NULL;
  io->separated_input_files=false;
  io->input_file_name=NULL;
  io->input_file_name_end1=NULL;
  io->input_file_name_end2=NULL;
  io->input_compression=FM_REGULAR_FILE;
  /* I/O */
  io->input_block_size = BUFFER_SIZE_32M;
  io->input_num_blocks = num_processors;
  io->input_buffer_size = BUFFER_SIZE_4M;
  io->output_file_name=NULL;
  io->output_compression=FM_REGULAR_FILE;
  /* I/O Attributes */
  io->fastq_strictly_normalized = false;
  /* Output */
  io->output_format = SAM;
  output_sam_parameters_set_defaults(&io->sam_parameters);
  output_map_parameters_set_defaults(&io->map_parameters);
  io->output_buffer_size = BUFFER_SIZE_4M;
  io->output_num_buffers = 10*num_processors; // Lazy allocation
  io->report_file_name = NULL;
}
void mapper_parameters_set_defaults_system(mapper_parameters_system_t* const system) {
  /* System */
  const uint64_t num_processors = system_get_num_processors();
  system->num_threads=num_processors;
  system->max_memory=0;
  system->tmp_folder=NULL;
}
void mapper_parameters_set_defaults_cuda(mapper_parameters_cuda_t* const cuda) {
  /* CUDA settings */
  const uint64_t num_processors = system_get_num_processors();
  /* CUDA */
  cuda->cuda_enabled=false;
  /* I/O */
  cuda->input_block_size = BUFFER_SIZE_32M;
  cuda->input_buffer_size = BUFFER_SIZE_4M;
  cuda->output_buffer_size = BUFFER_SIZE_4M;
  cuda->output_num_buffers = 10*num_processors; // Lazy allocation
  /* BPM Buffers */
  cuda->gpu_buffer_size = BUFFER_SIZE_1M;
  cuda->num_fmi_bsearch_buffers = 2;
  cuda->num_fmi_decode_buffers = 3;
  cuda->num_bpm_buffers = 3;
  /* Stages Configuration */
  cuda->cpu_emulation=false;
}
//void mapper_parameters_set_defaults_hints(mapper_parameters_hints_t* const hints) {
//  /* Hints */
//}
void mapper_parameters_set_defaults_misc(mapper_parameters_misc_t* const misc) {
  /* QC */
  misc->quality_control = false;
  misc->profile = false;
  misc->profile_reduce_type = reduce_sample;
  /* Verbose */
  misc->verbose_user=true;
  misc->verbose_dev=false;
}
void mapper_parameters_set_defaults(mapper_parameters_t* const mapper_parameters) {
  /* CMD line */
  mapper_parameters->argc = 0;
  mapper_parameters->argv = NULL;
  /* GEM Structures */
  mapper_parameters->archive = NULL;
  mapper_parameters->input_file = NULL;
  mapper_parameters->input_file_end1 = NULL;
  mapper_parameters->input_file_end2 = NULL;
  MUTEX_INIT(mapper_parameters->input_file_mutex);
  mapper_parameters->output_stream = NULL;
  mapper_parameters->output_file = NULL;
  /* Mapper Type */
  mapper_parameters->mapper_type = mapper_se;
	/* Stats Report */
	mapper_parameters->global_mapping_stats = NULL;
  /* I/O Parameters */
  mapper_parameters_set_defaults_io(&mapper_parameters->io);
  /* Search Parameters (single-end/paired-end) */
  search_parameters_init(&mapper_parameters->search_parameters);
  /* System */
  mapper_parameters_set_defaults_system(&mapper_parameters->system);
  /* CUDA settings */
  mapper_parameters_set_defaults_cuda(&mapper_parameters->cuda);
  /* Hints */
  // mapper_parameters_set_defaults_hints(&mapper_parameters->hints);
  /* Miscellaneous */
  mapper_parameters_set_defaults_misc(&mapper_parameters->misc);
}
/*
 * Mapper parameters display
 */
void mapper_parameters_print(FILE* const stream,mapper_parameters_t* const parameters) {
  tab_fprintf(stream,"[GEM]>Mapper.parameters\n");
  /* CMD line */
  uint64_t i;
  tab_fprintf(stream,"  => Application %s\n",parameters->argv[0]);
  tab_fprintf(stream,"  => Arguments   ");
  for (i=1;i<parameters->argc;++i) {
    fprintf(stream,"%s ",parameters->argv[i]);
  }
  fprintf(stream,"\n");
}
/*
 * Index loader
 */
void mapper_load_index(mapper_parameters_t* const parameters) {
  TIMER_START(&parameters->loading_time);
  PROFILE_START(GP_MAPPER_LOAD_INDEX,PROFILE_LEVEL);
  // Load archive
  gem_cond_log(parameters->misc.verbose_user,"[Loading GEM index '%s']",parameters->io.index_file_name);
  parameters->archive = archive_read(parameters->io.index_file_name,false);
  if (parameters->misc.verbose_dev) archive_print(gem_error_get_stream(),parameters->archive);
  PROFILE_STOP(GP_MAPPER_LOAD_INDEX,PROFILE_LEVEL);
  TIMER_STOP(&parameters->loading_time);
}
/*
 * Input (Low-level)
 */
void mapper_se_prepare_io_buffers(
    const mapper_parameters_t* const parameters,
    const uint64_t input_buffer_size,
    buffered_input_file_t** const buffered_fasta_input,
    buffered_output_file_t** const buffered_output_file) {
  *buffered_fasta_input = buffered_input_file_new(parameters->input_file,input_buffer_size);
  *buffered_output_file = buffered_output_file_new(parameters->output_file);
  buffered_input_file_attach_buffered_output(*buffered_fasta_input,*buffered_output_file);
}
uint64_t mapper_pe_reload_buffers(
    mapper_parameters_t* const parameters,
    buffered_input_file_t* const buffered_fasta_input_end1,
    buffered_input_file_t* const buffered_fasta_input_end2) {
  /// Check end-of-block
  error_code_t error_code;
  if (buffered_input_file_eob(buffered_fasta_input_end1)) {
    if (!parameters->io.separated_input_files) {
      error_code = buffered_input_file_reload(buffered_fasta_input_end1,0);
      if (error_code==INPUT_STATUS_EOF) return INPUT_STATUS_EOF;
      const uint64_t buffered_lines = buffered_fasta_input_end1->num_lines;
      if ((buffered_lines % 8) != 0) {
        MAPPER_ERROR_PE_PARSE_UNSYNCH_INPUT_FILES(parameters);
      }
    } else {
      // Dump buffer (explicitly before synch-reload)
      buffered_output_file_dump_buffer(buffered_fasta_input_end1->attached_buffered_output_file);
      // Reload buffers (in synch)
      MUTEX_BEGIN_SECTION(parameters->input_file_mutex) {
        // Check synch
        if (!buffered_input_file_eob(buffered_fasta_input_end2)) {
          MAPPER_ERROR_PE_PARSE_UNSYNCH_INPUT_FILES(parameters);
        }
        // Reload end1
        error_code = buffered_input_file_reload(buffered_fasta_input_end1,0);
        // Reload end2
        const uint64_t num_lines_read = buffered_fasta_input_end1->num_lines;
        error_code = buffered_input_file_reload(buffered_fasta_input_end2,num_lines_read);
        if (buffered_fasta_input_end2->num_lines != num_lines_read) {
          MAPPER_ERROR_PE_PARSE_UNSYNCH_INPUT_FILES(parameters);
        }
        if (num_lines_read==0) {
          MUTEX_END_SECTION(parameters->input_file_mutex);
          return INPUT_STATUS_EOF;
        }
      } MUTEX_END_SECTION(parameters->input_file_mutex);
    }
  }
  // OK
  return INPUT_STATUS_OK;
}
void mapper_pe_prepare_io_buffers(
    const mapper_parameters_t* const parameters,
    const uint64_t input_buffer_size,
    buffered_input_file_t** const buffered_fasta_input_end1,
    buffered_input_file_t** const buffered_fasta_input_end2,
    buffered_output_file_t** const buffered_output_file) {
  if (parameters->io.separated_input_files) {
    *buffered_fasta_input_end1 = buffered_input_file_new(parameters->input_file_end1,input_buffer_size);
    *buffered_fasta_input_end2 = buffered_input_file_new(parameters->input_file_end2,input_buffer_size);
  } else {
    *buffered_fasta_input_end1 = buffered_input_file_new(parameters->input_file,input_buffer_size);
    *buffered_fasta_input_end2 = *buffered_fasta_input_end1;
  }
  *buffered_output_file = buffered_output_file_new(parameters->output_file);
  buffered_input_file_attach_buffered_output(*buffered_fasta_input_end1,*buffered_output_file);
}
error_code_t mapper_pe_parse_paired_sequences(
    const mapper_parameters_t* const parameters,
    buffered_input_file_t* const buffered_fasta_input_end1,
    buffered_input_file_t* const buffered_fasta_input_end2,
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2) {
  error_code_t error_code;
  // Read end1
  error_code = input_fasta_parse_sequence(buffered_fasta_input_end1,
      archive_search_get_sequence(archive_search_end1),
      parameters->io.fastq_strictly_normalized,false);
  if (gem_expect_false(error_code!=INPUT_STATUS_OK)) return error_code;
  // Read end2
  error_code = input_fasta_parse_sequence(buffered_fasta_input_end2,
      archive_search_get_sequence(archive_search_end2),
      parameters->io.fastq_strictly_normalized,false);
  if (gem_expect_false(error_code!=INPUT_STATUS_OK)) return error_code;
  // OK
  PROF_ADD_COUNTER(GP_MAPPER_NUM_READS,2);
  return INPUT_STATUS_OK;
}
/*
 * Input (High-level)
 */
error_code_t mapper_se_read_single_sequence(mapper_search_t* const mapper_search) {
  const mapper_parameters_io_t* const parameters_io = &mapper_search->mapper_parameters->io;
  sequence_t* const sequence = archive_search_get_sequence(mapper_search->archive_search);
  const error_code_t error_code = input_fasta_parse_sequence(
      mapper_search->buffered_fasta_input,sequence,
      parameters_io->fastq_strictly_normalized,true);
  if (gem_expect_false(error_code==INPUT_STATUS_FAIL)) pthread_exit(0); // Abort
  PROF_INC_COUNTER(GP_MAPPER_NUM_READS);
  return error_code; // OK
}
error_code_t mapper_pe_read_paired_sequences(mapper_search_t* const mapper_search) {
  error_code_t error_code;
  // Check/Reload buffers manually
  error_code = mapper_pe_reload_buffers(mapper_search->mapper_parameters,
      mapper_search->buffered_fasta_input_end1,mapper_search->buffered_fasta_input_end2);
  if (error_code==INPUT_STATUS_EOF) return INPUT_STATUS_EOF;
  if (gem_expect_false(error_code==INPUT_STATUS_FAIL)) pthread_exit(0); // Abort
  // Read end1 & end2
  error_code = mapper_pe_parse_paired_sequences(mapper_search->mapper_parameters,
      mapper_search->buffered_fasta_input_end1,mapper_search->buffered_fasta_input_end2,
      mapper_search->archive_search_end1,mapper_search->archive_search_end2);
  if (gem_expect_false(error_code==INPUT_STATUS_FAIL)) pthread_exit(0); // Abort
  // OK
  return INPUT_STATUS_OK;
}
/*
 * Output
 */
void mapper_se_output_matches(
    mapper_parameters_t* const parameters,
    buffered_output_file_t* const buffered_output_file,
    archive_search_t* const archive_search,
    matches_t* const matches,
    mapping_stats_t *mstats) {
#ifdef MAPPER_OUTPUT
  switch (parameters->io.output_format) {
    case MAP:
      output_map_single_end_matches(buffered_output_file,archive_search,matches,&parameters->io.map_parameters);
      break;
    case SAM:
      output_sam_single_end_matches(buffered_output_file,archive_search,matches,&parameters->io.sam_parameters);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
	if (mstats) collect_se_mapping_stats(archive_search,matches,mstats);
#endif
}
void mapper_pe_output_matches(
    mapper_parameters_t* const parameters,
    buffered_output_file_t* const buffered_output_file,
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches,
    mapping_stats_t* const mstats) {
#ifdef MAPPER_OUTPUT
  switch (parameters->io.output_format) {
    case MAP:
      output_map_paired_end_matches(buffered_output_file,archive_search_end1,
          archive_search_end2,paired_matches,&parameters->io.map_parameters);
      break;
    case SAM:
      output_sam_paired_end_matches(buffered_output_file,
          archive_search_end1,archive_search_end2,
          paired_matches,&parameters->io.sam_parameters);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
	if (mstats) collect_pe_mapping_stats(archive_search_end1,archive_search_end2,paired_matches,mstats);
#endif
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
  search_parameters_t* const search_parameters = &mapper_search->mapper_parameters->search_parameters;
  mm_search_t* const mm_search = mm_search_new();
  archive_search_se_new(parameters->archive,search_parameters,false,NULL,&mapper_search->archive_search);
  archive_search_se_inject_mm(mapper_search->archive_search,mm_search);
  matches_t* const matches = matches_new(mm_search->mm_stack);
  matches_configure(matches,mapper_search->archive_search->text_collection);
  archive_t* const archive = parameters->archive;
  archive_search_t* const archive_search = mapper_search->archive_search;
  const bool bisulfite_index = (archive->type == archive_dna_bisulfite); // Bisulfite

  // FASTA/FASTQ reading loop
  uint64_t reads_processed = 0;
  while (mapper_se_read_single_sequence(mapper_search)) {
//    if (gem_streq(mapper_search->archive_search->sequence.tag.buffer,"H.Sapiens.1M.Illumina.l100.low.000000460")) {
//      printf("HERE\n");
//    }

    // Bisulfite: Fully convert reads before searching into archive, making a copy of the original
    if (bisulfite_index) mapper_bisulfite_process_sequence_se(archive_search,search_parameters);

    // Search into the archive
    archive_search_se(archive_search,matches);

    // Bisulfite: Copy back original read
    if (bisulfite_index) mapper_bisulfite_restore_sequence_se(archive_search);

    // Output matches
    mapper_se_output_matches(parameters,mapper_search->buffered_output_file,
        archive_search,matches,mapper_search->mapping_stats);
    // output_fastq(mapper_search->buffered_output_file,&mapper_search->archive_search->sequence);

    // Update processed
    if (++reads_processed == MAPPER_TICKER_STEP) {
      ticker_update_mutex(mapper_search->ticker,reads_processed);
      reads_processed=0;
    }

    // Clear
    mm_search_clear(mm_search);
    matches_clear(matches);
  }
  // Update processed
  ticker_update_mutex(mapper_search->ticker,reads_processed);

  // Clean up
  buffered_input_file_close(mapper_search->buffered_fasta_input);
  buffered_output_file_close(mapper_search->buffered_output_file);
  archive_search_delete(archive_search);
  matches_delete(matches);
  mm_search_delete(mm_search);

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
  mm_search_t* const mm_search = mm_search_new();
  search_parameters_t* const base_search_parameters = &mapper_search->mapper_parameters->search_parameters;
  mapper_stats_template_init(mm_search->mapper_stats,
      base_search_parameters->search_paired_parameters.min_initial_template_estimation,
      base_search_parameters->search_paired_parameters.max_initial_template_estimation);
  archive_search_pe_new(parameters->archive,base_search_parameters,false,NULL,
      &mapper_search->archive_search_end1,&mapper_search->archive_search_end2);
  archive_search_pe_inject_mm(mapper_search->archive_search_end1,mapper_search->archive_search_end2,mm_search);
  mapper_search->paired_matches = paired_matches_new(mm_search->mm_stack);
  paired_matches_configure(mapper_search->paired_matches,mapper_search->archive_search_end1->text_collection);
  archive_t* const archive = parameters->archive;
  archive_search_t* const archive_search_end1 = mapper_search->archive_search_end1;
  archive_search_t* const archive_search_end2 = mapper_search->archive_search_end2;
  paired_matches_t* const paired_matches = mapper_search->paired_matches;
  const bool bisulfite_index = (archive->type == archive_dna_bisulfite);

  // FASTA/FASTQ reading loop
  uint64_t reads_processed = 0;
  while (mapper_pe_read_paired_sequences(mapper_search)) {
//    if (gem_streq(mapper_search->archive_search_end1->sequence.tag.buffer,"H.Sapiens.1M.Illumina.l100.low.000216715")) {
//      printf("HERE\n");
//    }

    // Bisulfite: Fully convert reads before searching into archive, making a copy of the original
    if (bisulfite_index) mapper_bisulfite_process_sequence_pe(archive_search_end1,archive_search_end2);

    // Search into the archive
    archive_search_pe(archive_search_end1,archive_search_end2,paired_matches);

    // Bisulfite: Copy back original read
    if (bisulfite_index) mapper_bisulfite_restore_sequence_pe(archive_search_end1,archive_search_end2);

    // Output matches
    mapper_pe_output_matches(parameters,mapper_search->buffered_output_file,
        archive_search_end1,archive_search_end2,paired_matches,mapper_search->mapping_stats);
    //output_fastq(mapper_search->buffered_output_file,&archive_search_end1->sequence);
    //output_fastq(mapper_search->buffered_output_file,&archive_search_end2->sequence);

    // Update processed
    if (++reads_processed == MAPPER_TICKER_STEP) {
      ticker_update_mutex(mapper_search->ticker,reads_processed);
      reads_processed=0;
    }

    // Clear
    mm_search_clear(mm_search);
    paired_matches_clear(paired_matches,true);
  }
  // Update processed
  ticker_update_mutex(mapper_search->ticker,reads_processed);

  // Clean up
  buffered_input_file_close(mapper_search->buffered_fasta_input_end1);
  if (parameters->io.separated_input_files) buffered_input_file_close(mapper_search->buffered_fasta_input_end2);
  buffered_output_file_close(mapper_search->buffered_output_file);
  archive_search_delete(mapper_search->archive_search_end1);
  archive_search_delete(mapper_search->archive_search_end2);
  paired_matches_delete(mapper_search->paired_matches);
  mm_search_delete(mm_search);

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
  mapper_search_t* const mapper_search = mm_calloc(num_threads,mapper_search_t,false); // Allocate mapper searches
  // Set error-report function
  g_mapper_searches = mapper_search;
  gem_error_set_report_function(mapper_error_report);
  // Prepare output file/parameters (SAM headers)
  archive_t* const archive = mapper_parameters->archive;
  const bool bisulfite_index = (archive->type == archive_dna_bisulfite);
  if (mapper_parameters->io.output_format==SAM) {
    output_sam_print_header(mapper_parameters->output_file,archive,
        &mapper_parameters->io.sam_parameters,mapper_parameters->argc,mapper_parameters->argv);
    mapper_parameters->io.sam_parameters.bisulfite_output = bisulfite_index;
  }
  // Setup Ticker
  ticker_t ticker;
  ticker_count_reset(&ticker,mapper_parameters->misc.verbose_user,
      paired_end ? "PE::Mapping Sequences" : "SE::Mapping Sequences",0,MAPPER_TICKER_STEP,false);
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
  //__itt_resume();
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
        pthread_create(mapper_search[i].thread_data,0,mapper_thread,(void*)(mapper_search+i)),SYS_THREAD_CREATE);
  }
  // Join all threads
  for (i=0;i<num_threads;++i) {
    gem_cond_fatal_error__perror(pthread_join(*(mapper_search[i].thread_data),0),SYS_THREAD_JOIN);
    mm_free(mapper_search[i].thread_data);
  }
  //__itt_pause();
  ticker_finish(&ticker);
  ticker_mutex_cleanup(&ticker);
	// Merge report stats
	if (mstats) {
		 merge_mapping_stats(mapper_parameters->global_mapping_stats,mstats,num_threads);
		 mm_free(mstats);
	}
  // Clean up
  mm_free(mapper_search);
}
void mapper_se_run(mapper_parameters_t* const mapper_parameters) {
  mapper_run(mapper_parameters,false);
}
void mapper_pe_run(mapper_parameters_t* const mapper_parameters) {
  mapper_run(mapper_parameters,true);
}
