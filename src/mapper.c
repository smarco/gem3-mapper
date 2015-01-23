/*
 * PROJECT: GEMMapper
 * FILE: mapper.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "mapper.h"
#include "archive_search.h"

/*
 * Report
 */
void mapper_display_input_state(
    FILE* stream,buffered_input_file_t* const buffered_fasta_input,const sequence_t* const sequence) {
  // Check NULL
  if (sequence==NULL) { tab_fprintf(stream,"Sequence is NULL\n"); return; }
  if (buffered_fasta_input==NULL) { tab_fprintf(stream,"Buffered_fasta_input is NULL\n"); return; }
  // Dump FASTA/FASTQ read
  if (!string_is_null(&sequence->tag) && !string_is_null(&sequence->read)) {
    const bool has_qualities = sequence_has_qualities(sequence);
    char* const end_tag =
        (sequence->attributes.end_info == PAIRED_END1) ? "/1" :
      ( (sequence->attributes.end_info == PAIRED_END2) ? "/2" : " " );
    tab_fprintf(stream,"Sequence (File '%s' Line '%lu')\n",
        input_file_get_file_name(buffered_fasta_input->input_file),
        buffered_fasta_input->input_buffer->current_line_num - (has_qualities ? 4 : 2));
    if (has_qualities) {
      fprintf(stream,"@%"PRIs"%s\n%"PRIs"\n+\n%"PRIs"\n",
          PRIs_content(&sequence->tag),end_tag,
          PRIs_content(&sequence->read),
          PRIs_content(&sequence->qualities));
    } else {
      fprintf(stream,">%"PRIs"%s\n%"PRIs"\n",
          PRIs_content(&sequence->tag),end_tag,
          PRIs_content(&sequence->read));
    }
  } else {
    tab_fprintf(stream,"Sequence <<Empty>>\n");
  }
}
/*
 * Error Report
 */
mapper_search_t* g_mapper_searches; // Global searches on going
void mapper_error_report(FILE* stream) {
//  // Display thread info
//  const uint64_t threads_id = gem_thread_get_thread_id();
//  if (threads_id==0) {
//    fprintf(stream,"GEM::Running-Thread (threadID = MASTER)\n");
//  } else {
//    uint64_t i;
//    // Display Threads-Info
//    const uint64_t num_threads = g_mapper_searches->mapper_parameters->system.num_threads;
//    for (i=0;i<num_threads;++i) {
//      mapper_search_t* const mapper_search = g_mapper_searches + i;
//      // Thread
//      fprintf(stream,"GEM::Running-Thread (threadID = %lu)\n",mapper_search->thread_id);
////      // Display Input State
////      const sequence_t* const sequence = archive_search_get_sequence(mapper_search->archive_search);
////      tab_global_inc();
////      mapper_display_input_state(stream,mapper_search->buffered_fasta_input,sequence);
////      tab_global_dec();
//      // Display Output State
//      // Display Search State
//      // TODO ¿More useful info?
//    }
////    // Display stats until now (if possible)
////    PROF_BLOCK() {
////
////    }
//  }
}
/*
 * Mapper Parameters
 */
GEM_INLINE void mapper_parameters_set_defaults_io(mapper_parameters_io_t* const io) {
  /* Input */
  io->index_file_name=NULL;
  io->check_index=false;
  io->separated_input_files=false;
  io->input_file_name=NULL;
  io->input_file_name_end1=NULL;
  io->input_file_name_end2=NULL;
  io->input_compression=FM_REGULAR_FILE;
  io->input_block_size = BUFFER_SIZE_64M;
  io->input_buffer_lines = (2*4*NUM_LINES_5K); // 2l-Paired x 4l-FASTQRecord x 5K-BufferSize
  io->output_file_name=NULL;
  io->output_compression=FM_REGULAR_FILE;
  /* I/O Attributes */
  io->fastq_strictly_normalized = false;
  io->fastq_try_recovery = false;
  /* Output */
  io->output_format = MAP;
  output_sam_parameters_set_defaults(&io->sam_parameters);
  const uint64_t num_processors = system_get_num_processors();
  io->output_buffer_size = BUFFER_SIZE_4M;
  io->output_num_buffers = 5*num_processors; // Lazy allocation
}
GEM_INLINE void mapper_parameters_set_defaults_system(mapper_parameters_system_t* const system) {
  /* System */
  const uint64_t num_processors = system_get_num_processors();
  system->num_threads=num_processors;
  system->max_memory=0;
  system->tmp_folder=NULL;
}
GEM_INLINE void mapper_parameters_set_defaults_cuda(mapper_parameters_cuda_t* const cuda) {
  /* CUDA settings */
  const uint64_t num_processors = system_get_num_processors();
  /* CUDA */
  cuda->cuda_enabled=false;
  /* I/O */
  cuda->input_block_size = BUFFER_SIZE_64M;
  cuda->input_num_buffers = 2*num_processors;
  cuda->input_buffer_lines = (2*4*NUM_LINES_5K); // 2l-Paired x 4l-FASTQRecord x 5K-BufferSize
  cuda->output_buffer_size = BUFFER_SIZE_4M;
  cuda->output_num_buffers = 5*num_processors; // Lazy allocation
  /* BPM Buffers */
  cuda->num_search_groups_per_thread = 3;
  cuda->bpm_buffer_size = BUFFER_SIZE_4M;
}
GEM_INLINE void mapper_parameters_set_defaults_hints(mapper_parameters_hints_t* const hints) {
  /* Hints */
  hints->avg_read_length=150;
  hints->std_read_length=50;
  hints->candidates_per_query=20;
}
GEM_INLINE void mapper_parameters_set_defaults_misc(mapper_parameters_misc_t* const misc) {
  /* QC */
  misc->quality_control = false;
  misc->profile = false;
  misc->profile_reduce_type = reduce_sample;
  /* Verbose */
  misc->verbose_user=true;
  misc->verbose_dev=false;
}
GEM_INLINE void mapper_parameters_set_defaults(mapper_parameters_t* const mapper_parameters) {
  /* CMD line */
  mapper_parameters->argc = 0;
  mapper_parameters->argv = NULL;
  /* GEM Structures */
  mapper_parameters->archive = NULL;
  mapper_parameters->input_file = NULL;
  mapper_parameters->input_file_end1 = NULL;
  mapper_parameters->input_file_end2 = NULL;
  mapper_parameters->output_stream = NULL;
  mapper_parameters->output_file = NULL;
  /* Mapper Type */
  mapper_parameters->mapper_type = mapper_se;
  /* I/O Parameters */
  mapper_parameters_set_defaults_io(&mapper_parameters->io);
  /* Search Parameters (single-end/paired-end) */
  approximate_search_parameters_init(&mapper_parameters->search_parameters);
  /* Select Parameters */
  archive_select_parameters_init(&mapper_parameters->select_parameters);
  /* System */
  mapper_parameters_set_defaults_system(&mapper_parameters->system);
  /* CUDA settings */
  mapper_parameters_set_defaults_cuda(&mapper_parameters->cuda);
  /* Hints */
  mapper_parameters_set_defaults_hints(&mapper_parameters->hints);
  /* Miscellaneous */
  mapper_parameters_set_defaults_misc(&mapper_parameters->misc);
}
/*
 * Mapper parameters display
 */
GEM_INLINE void mapper_parameters_print_mode(FILE* const stream,const mapper_type mapper_type,const bool cuda_enabled) {
  /* Mapper Mode */
  switch (mapper_type) {
    case mapper_se:
      tab_fprintf(stream,"  => Mapper.type\tMapperSE");
      break;
    case mapper_pe:
      tab_fprintf(stream,"  => Mapper.type\ttMapperPE");
      break;
    case mapper_graph:
      tab_fprintf(stream,"  => Mapper.type\tMapper-Graph");
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  // Print CUDA enabled
  fprintf(stream,"%s\n",cuda_enabled?"\t[CUDA]":" ");
}
GEM_INLINE void mapper_parameters_print_io_input_file(
    FILE* const stream,const char* const label,const char* const input_file_name,
    const bool fastq_strictly_normalized,const bool fastq_try_recovery) {
  if (input_file_name!=NULL) {
    tab_fprintf(stream,"      => %s\t%s (%luMB)%s%s\n",
        label,input_file_name,CONVERT_B_TO_MB(gem_file_size(input_file_name)),
        (fastq_strictly_normalized) ? "\t[Normalize=ON]" : " ",
        (fastq_try_recovery) ? "\t[Recovery=ON]" : " ");
  } else {
    tab_fprintf(stream,"      => %s\t<<stdin>>%s%s\n",label,
        (fastq_strictly_normalized) ? "\t[Normalize=ON]" : " ",
        (fastq_try_recovery) ? "\t[Recovery=ON]" : " ");
  }
}
GEM_INLINE void mapper_parameters_print_io(FILE* const stream,mapper_parameters_io_t* const io) {
  /* Index */
  tab_fprintf(stream,"    => Index\n");
  tab_fprintf(stream,"      => File\t%s (%luMB) %s\n",
      io->index_file_name,CONVERT_B_TO_MB(gem_file_size(io->index_file_name)),io->check_index?"[checked]":" ");
  /* Input */
  tab_fprintf(stream,"    => Input\n");
  if (!io->separated_input_files) {
    mapper_parameters_print_io_input_file(stream,"File",io->input_file_name,
        io->fastq_strictly_normalized,io->fastq_try_recovery); // Single-File
  } else {
    mapper_parameters_print_io_input_file(stream,"File.End1",io->input_file_name,
        io->fastq_strictly_normalized,io->fastq_try_recovery); // End-1
    mapper_parameters_print_io_input_file(stream,"File.End2",io->input_file_name,
        io->fastq_strictly_normalized,io->fastq_try_recovery); // End-2
  }
  tab_fprintf(stream,"      => Block.size\t%luMB\n",CONVERT_B_TO_MB(io->input_block_size));
//  tab_fprintf(stream,"    => Num.buffers\t%lu\n",io->input_num_buffers);
  tab_fprintf(stream,"      => Buffer.lines\t%luKlines\n",io->input_buffer_lines/1000);
  /* Output */
  tab_fprintf(stream,"  => Output\n");
  tab_fprintf(stream,"    => File\t%s",(io->output_file_name!=NULL) ? io->output_file_name : "<<stdout>>");
  switch (io->output_compression) {
    case FM_GZIPPED_FILE: fprintf(stream,"\t\[GZIP]n"); break;
    case FM_BZIPPED_FILE: fprintf(stream,"\t[BZIP]\n"); break;
    default: fprintf(stream,"\n"); break;
  }
  switch (io->output_format) {
    case MAP:
      tab_fprintf(stream,"    => Format\tMAP\n");
      break;
    case SAM:
      tab_fprintf(stream,"    => Format\tSAM\n");
//      output_sam_parameters_t sam_parameters; // TODO
      break;
    default:
      tab_fprintf(stream,"    => Format\tUNKNOWN?\n");
      break;
  }
  tab_fprintf(stream,"    => Buffer.size\t%luMB\n",CONVERT_B_TO_MB(io->output_buffer_size));
  tab_fprintf(stream,"    => Num.buffers\t%lu\n",io->output_num_buffers);
}
GEM_INLINE void mapper_parameters_print_search_parameters(FILE* const stream,search_parameters_t* const search_parameters) {
//  /*
//   * Search parameters
//   */
//  /* Mapping strategy (Mapping mode + properties) */
//  mapping_mode_t mapping_mode;
//  float mapping_degree;
//  /* Qualities */
//  quality_model_t quality_model;
//  quality_format_t quality_format;
//  uint64_t quality_threshold;
//  /* Error Model (Regulates the number of Mismatch/Indels) */
//  float max_search_error;
//  float max_filtering_error;
//  float complete_strata_after_best;
//  float min_matching_length;
//  /* Matches search (Regulates the number of matches) */
//  uint64_t max_search_matches;
//  /* Replacements (Regulates the bases that can be replaced/mismatched) */
//  char mismatch_alphabet[DNA_EXT_RANGE];
//  uint64_t mismatch_alphabet_length;
//  bool allowed_chars[256];
//  bool allowed_enc[DNA_EXT_RANGE];
//  /* Alignment Model/Score */
//  alignment_model_t alignment_model;
//  uint64_t matching_score;
//  uint64_t mismatch_penalty;
//  uint64_t gap_open_penalty;
//  uint64_t gap_extension_penalty;
//  /*
//   * Internals
//   */
//  /* Soft Region Profile Parameters */
//  uint64_t srp_region_th; // Max. number of candidates allowed per region
//  uint64_t srp_max_steps; // Max. number of characters to explore to improve the region
//  uint64_t srp_dec_factor; // Decreasing factor per step in region exploration
//  uint64_t srp_region_type_th; // Threshold to classify regions {ZERO,NON_ZERO}
//  /* Hard Region Profile Parameters */
//  uint64_t hrp_region_th;
//  uint64_t hrp_max_steps;
//  uint64_t hrp_dec_factor;
//  uint64_t hrp_region_type_th;
//  /* Read recovery parameters */
//  uint64_t rrp_region_th;
//  uint64_t rrp_max_steps;
//  uint64_t rrp_dec_factor;
//  uint64_t rrp_region_type_th;
//  /* Filtering parameters */
//  uint64_t filtering_threshold;
//  float filtering_region_factor;
//  uint64_t pa_filtering_threshold;
//  /* Checkers */
//  check_matches_t check_matches;
}
GEM_INLINE void mapper_parameters_print_paired_end(FILE* const stream,search_parameters_t* const search_parameters) {
//  /* Paired-end mode/alg */
//  bool paired_end;
//  bool map_both_ends;
//  uint64_t max_extendable_candidates;
//  uint64_t max_matches_per_extension;
//  /* Template allowed length */
//  uint64_t min_template_length;
//  uint64_t max_template_length;
//  /* Concordant Orientation */
//  bool pair_orientation_FR;
//  bool pair_orientation_RF;
//  bool pair_orientation_FF;
//  bool pair_orientation_RR;
//  /* Discordant Orientation */
//  bool discordant_pair_orientation_FR;
//  bool discordant_pair_orientation_RF;
//  bool discordant_pair_orientation_FF;
//  bool discordant_pair_orientation_RR;
//  /* Pair allowed lay-outs */
//  bool pair_layout_separate;
//  bool pair_layout_overlap;
//  bool pair_layout_contain;
//  bool pair_layout_dovetail;
}
GEM_INLINE void mapper_parameters_print_select_parameters(FILE* const stream,select_parameters_t* const select_parameters) {
  /* Select Parameters */
}
GEM_INLINE void mapper_parameters_print_system(FILE* const stream,mapper_parameters_system_t* const system) {
  /* System */
  tab_fprintf(stream,"  => Num.Threads\t%lu\n",system->num_threads);
  if (system->max_memory!=0) {
    tab_fprintf(stream,"  => Max.Memory\t%lu\n",system->max_memory);
  }
  if (system->tmp_folder!=NULL) {
    tab_fprintf(stream,"  => Tmp.path\t%s\n",system->tmp_folder);
  }
}
GEM_INLINE void mapper_parameters_print_cuda(FILE* const stream,mapper_parameters_cuda_t* const cuda) {
  /* CUDA settings */
}
GEM_INLINE void mapper_parameters_print_hints(FILE* const stream,mapper_parameters_hints_t* const hints) {
  /* Hints */
}
GEM_INLINE void mapper_parameters_print_misc(FILE* const stream,mapper_parameters_misc_t* const misc) {
  /* Miscellaneous */

}
GEM_INLINE void mapper_parameters_print(FILE* const stream,mapper_parameters_t* const parameters) {
  tab_fprintf(stream,"[GEM]>Mapper.parameters\n");
  /* CMD line */
  uint64_t i;
  tab_fprintf(stream,"  => Application %s\n",parameters->argv[0]);
  tab_fprintf(stream,"  => Arguments   ");
  for (i=1;i<parameters->argc;++i) {
    fprintf(stream,"%s ",parameters->argv[i]);
  }
  fprintf(stream,"\n");
  /* Mapper Mode */
  tab_fprintf(stream,"  => Mapper.Mode\n");
  tab_global_inc();
  mapper_parameters_print_mode(stream,parameters->mapper_type,parameters->cuda.cuda_enabled);
  tab_global_dec();
  /* I/O Parameters */
  tab_fprintf(stream,"  => I/O\n");
  tab_global_inc();
  mapper_parameters_print_io(stream,&parameters->io);
  tab_global_dec();
  /* Single-end Search Parameters */
  tab_fprintf(stream,"  => Single.end.Parameters\n");
  tab_global_inc();
  mapper_parameters_print_search_parameters(stream,&parameters->search_parameters);
  tab_global_dec();
  /* Paired-end Search Parameters */
  tab_fprintf(stream,"  => Paired.end.Parameters\n");
  tab_global_inc();
  mapper_parameters_print_paired_end(stream,&parameters->search_parameters);
  tab_global_dec();
  /* Select Parameters */
  tab_fprintf(stream,"  => Select.Parameters\n");
  tab_global_inc();
  mapper_parameters_print_select_parameters(stream,&parameters->select_parameters);
  tab_global_dec();
  /* System */
  tab_fprintf(stream,"  => System\n");
  tab_global_inc();
  mapper_parameters_print_system(stream,&parameters->system);
  tab_global_dec();
  /* CUDA settings */
  tab_fprintf(stream,"  => CUDA.Setttings\n");
  tab_global_inc();
  mapper_parameters_print_cuda(stream,&parameters->cuda);
  tab_global_dec();
  /* Hints */
  tab_fprintf(stream,"  => Hints\n");
  tab_global_inc();
  mapper_parameters_print_hints(stream,&parameters->hints);
  tab_global_dec();
  /* Miscellaneous */
  tab_fprintf(stream,"  => Misc\n");
  tab_global_inc();
  mapper_parameters_print_misc(stream,&parameters->misc);
  tab_global_dec();
}
/*
 * Index loader
 */
GEM_INLINE void mapper_load_index(mapper_parameters_t* const parameters) {
  PROF_START_TIMER(GP_MAPPER_LOAD_INDEX);
  // Load archive
  gem_cond_log(parameters->misc.verbose_user,"[Loading GEM index '%s']",parameters->io.index_file_name);
  parameters->archive = archive_read(parameters->io.index_file_name,
      parameters->io.check_index,parameters->misc.verbose_dev);
  PROF_STOP_TIMER(GP_MAPPER_LOAD_INDEX);
}
/*
 * I/O
 */
GEM_INLINE void mapper_SE_output_matches(
    const mapper_parameters_t* const parameters,
    buffered_output_file_t* const buffered_output_file,
    sequence_t* const seq_read,matches_t* const matches) {
  switch (parameters->io.output_format) {
    case MAP:
      output_map_single_end_matches(buffered_output_file,seq_read,matches);
      break;
    case SAM:
      output_sam_single_end_matches(buffered_output_file,seq_read,matches,&parameters->io.sam_parameters);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
GEM_INLINE void mapper_PE_output_matches(
    const mapper_parameters_t* const parameters,
    buffered_output_file_t* const buffered_output_file,
    sequence_t* const seq_read_end1,sequence_t* const seq_read_end2,
    paired_matches_t* const paired_matches) {
  switch (parameters->io.output_format) {
    case MAP:
      output_map_paired_end_matches(buffered_output_file,
          seq_read_end1,seq_read_end2,paired_matches);
      break;
    case SAM:
      output_sam_paired_end_matches(buffered_output_file,
          seq_read_end1,seq_read_end2,paired_matches,&parameters->io.sam_parameters);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
GEM_INLINE error_code_t mapper_SE_parse_sequence(mapper_search_t* const mapper_search) {
  const mapper_parameters_t* const parameters = mapper_search->mapper_parameters;
  const error_code_t error_code = input_fasta_parse_sequence(
        mapper_search->buffered_fasta_input,archive_search_get_sequence(mapper_search->archive_search),
        parameters->io.fastq_strictly_normalized,parameters->io.fastq_try_recovery,true);
  if (gem_expect_false(error_code==INPUT_STATUS_FAIL)) pthread_exit(0); // Abort
  return error_code;
}
GEM_INLINE error_code_t mapper_PE_parse_sequence(mapper_search_t* const mapper_search) {
  const mapper_parameters_t* const parameters = mapper_search->mapper_parameters;
  error_code_t error_code;
  // Read end1
  error_code = input_fasta_parse_sequence(
      (parameters->io.separated_input_files) ? mapper_search->buffered_fasta_input_end1 : mapper_search->buffered_fasta_input,
      archive_search_get_sequence(mapper_search->archive_search_end1),
      parameters->io.fastq_strictly_normalized,parameters->io.fastq_try_recovery,true);
  if (gem_expect_false(error_code==INPUT_STATUS_FAIL)) pthread_exit(0); // Abort
  if (gem_expect_false(error_code==INPUT_STATUS_EOF)) return INPUT_STATUS_EOF;
  // Read end2
  error_code = input_fasta_parse_sequence(
      (parameters->io.separated_input_files) ? mapper_search->buffered_fasta_input_end2 : mapper_search->buffered_fasta_input,
      archive_search_get_sequence(mapper_search->archive_search_end2),
      parameters->io.fastq_strictly_normalized,parameters->io.fastq_try_recovery,true);
  if (gem_expect_false(error_code==INPUT_STATUS_FAIL)) pthread_exit(0); // Abort
  if (gem_expect_false(error_code==INPUT_STATUS_EOF)) {
    if (parameters->io.separated_input_files) {
      gem_fatal_error_msg("Parsing Input Files. Files '%s,%s' doesn't contain the same number of reads (cannot pair)",
          input_file_get_file_name(parameters->input_file_end1),input_file_get_file_name(parameters->input_file_end2));
    } else {
      gem_fatal_error_msg("Parsing Input File. File '%s' doesn't contain a pair number of reads (cannot pair)",
          input_file_get_file_name(parameters->input_file));
    }
    pthread_exit(0); // Abort
  }
  return INPUT_STATUS_OK;
}
/*
 * SE Mapper
 */
void* mapper_SE_thread(mapper_search_t* const mapper_search) {
  // GEM-thread error handler
  gem_thread_register_id(mapper_search->thread_id+1);

  // Create new buffered reader/writer
  const mapper_parameters_t* const parameters = mapper_search->mapper_parameters;
  mapper_search->buffered_fasta_input = buffered_input_file_new(parameters->input_file,parameters->io.input_buffer_lines);
  buffered_output_file_t* const buffered_output_file = buffered_output_file_new(parameters->output_file);
  buffered_input_file_attach_buffered_output(mapper_search->buffered_fasta_input,buffered_output_file);

  // Create an Archive-Search
  mm_search_t* const mm_search = mm_search_new(mm_pool_get_slab(mm_pool_32MB));
  mapper_search->archive_search = archive_search_new(
      parameters->archive,&mapper_search->mapper_parameters->search_parameters,
      &mapper_search->mapper_parameters->select_parameters);
  archive_search_configure(mapper_search->archive_search,mm_search);
  matches_t* const matches = matches_new();
  matches_configure(matches,mapper_search->archive_search->text_collection);

  // FASTA/FASTQ reading loop
  uint64_t reads_processed = 0;
  while (mapper_SE_parse_sequence(mapper_search)) {
    // Search into the archive
    archive_search_single_end(mapper_search->archive_search,matches);

    // Output matches
    mapper_SE_output_matches(parameters,buffered_output_file,archive_search_get_sequence(mapper_search->archive_search),matches);

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
  buffered_output_file_close(buffered_output_file);
  archive_search_delete(mapper_search->archive_search);
  matches_delete(matches);
  mm_search_delete(mm_search);

  pthread_exit(0);
}
/*
 * PE Mapper
 */
void* mapper_PE_thread(mapper_search_t* const mapper_search) {
  // GEM-thread error handler
  gem_thread_register_id(mapper_search->thread_id+1);

  // Create new buffered reader/writer
  const mapper_parameters_t* const parameters = mapper_search->mapper_parameters;
  buffered_output_file_t* const buffered_output_file = buffered_output_file_new(parameters->output_file);
  if (parameters->io.separated_input_files) {
    mapper_search->buffered_fasta_input_end1 = buffered_input_file_new(parameters->input_file_end1,parameters->io.input_buffer_lines);
    mapper_search->buffered_fasta_input_end2 = buffered_input_file_new(parameters->input_file_end2,parameters->io.input_buffer_lines);
    buffered_input_file_attach_buffered_output(mapper_search->buffered_fasta_input_end1,buffered_output_file); // Just end1
  } else {
    mapper_search->buffered_fasta_input = buffered_input_file_new(parameters->input_file,parameters->io.input_buffer_lines);
    buffered_input_file_attach_buffered_output(mapper_search->buffered_fasta_input,buffered_output_file);
  }

  // Create an Archive-Search
  mm_search_t* const mm_search = mm_search_new(mm_pool_get_slab(mm_pool_32MB));
  mapper_search->archive_search_end1 = archive_search_new(
      parameters->archive,&mapper_search->mapper_parameters->search_parameters,
      &mapper_search->mapper_parameters->select_parameters);
  mapper_search->archive_search_end2 = archive_search_new(
      parameters->archive,&mapper_search->mapper_parameters->search_parameters,
      &mapper_search->mapper_parameters->select_parameters);
  archive_search_configure(mapper_search->archive_search_end1,mm_search);
  archive_search_configure(mapper_search->archive_search_end2,mm_search);
  mapper_search->paired_matches = paired_matches_new(mm_search->text_collection);
  paired_matches_configure(mapper_search->paired_matches,mapper_search->archive_search_end1->text_collection);

  // FASTA/FASTQ reading loop
  uint64_t reads_processed = 0;
  while (mapper_PE_parse_sequence(mapper_search)) {
    // Search into the archive
    archive_search_paired_end(mapper_search->archive_search_end1,
        mapper_search->archive_search_end2,mapper_search->paired_matches);

    // Output matches
    mapper_PE_output_matches(parameters,buffered_output_file,
        archive_search_get_sequence(mapper_search->archive_search_end1),
        archive_search_get_sequence(mapper_search->archive_search_end2),mapper_search->paired_matches);

    // Update processed
    if (++reads_processed == MAPPER_TICKER_STEP) {
      ticker_update_mutex(mapper_search->ticker,reads_processed);
      reads_processed=0;
    }

    // Clear
    mm_search_clear(mm_search);
    paired_matches_clear(mapper_search->paired_matches);
  }
  // Update processed
  ticker_update_mutex(mapper_search->ticker,reads_processed);

  // Clean up
  if (parameters->io.separated_input_files) {
    buffered_input_file_close(mapper_search->buffered_fasta_input_end1);
    buffered_input_file_close(mapper_search->buffered_fasta_input_end2);
  } else {
    buffered_input_file_close(mapper_search->buffered_fasta_input);
  }
  buffered_output_file_close(buffered_output_file);
  archive_search_delete(mapper_search->archive_search_end1);
  archive_search_delete(mapper_search->archive_search_end2);
  paired_matches_delete(mapper_search->paired_matches);
  mm_search_delete(mm_search);

  pthread_exit(0);
}
/*
 * SE/PE runnable
 */
GEM_INLINE void mapper_run(mapper_parameters_t* const mapper_parameters,const bool paired_end) {
  // Load GEM-Index
  mapper_load_index(mapper_parameters);
  // Setup threads
  const uint64_t num_threads = mapper_parameters->system.num_threads;
  mapper_search_t* const mapper_search = mm_calloc(num_threads,mapper_search_t,false); // Allocate mapper searches
  // Set error-report function
  g_mapper_searches = mapper_search;
  gem_error_set_report_function(mapper_error_report);
  // Prepare output file (SAM headers)
  if (mapper_parameters->io.output_format==SAM) {
    output_sam_print_header(mapper_parameters->output_file,
        mapper_parameters->archive,mapper_parameters->argc,mapper_parameters->argv);
  }
  // Setup Ticker
  ticker_t ticker;
  ticker_count_reset(&ticker,mapper_parameters->misc.verbose_user,
      paired_end ? "PE::Mapping Sequences" : "SE::Mapping Sequences",0,MAPPER_TICKER_STEP,false);
  ticker_add_process_label(&ticker,"#","sequences processed");
  ticker_add_finish_label(&ticker,"Total","sequences processed");
  ticker_mutex_enable(&ticker);
  // Launch threads
  uint64_t i;
  for (i=0;i<num_threads;++i) {
    // Setup thread
    mapper_search[i].thread_id = i;
    mapper_search[i].thread_data = mm_alloc(pthread_t);
    mapper_search[i].mapper_parameters = mapper_parameters;
    mapper_search[i].ticker = &ticker;
    // Launch thread
    gem_cond_fatal_error__perror(
        pthread_create(mapper_search[i].thread_data,0,
            (void* (*)(void*))(paired_end ? mapper_PE_thread : mapper_SE_thread),(void*)(mapper_search+i)),
        SYS_THREAD_CREATE);
  }
  // Join all threads
  for (i=0;i<num_threads;++i) {
    gem_cond_fatal_error__perror(pthread_join(*(mapper_search[i].thread_data),0),SYS_THREAD_JOIN);
    mm_free(mapper_search[i].thread_data);
  }
  ticker_finish(&ticker);
  ticker_mutex_cleanup(&ticker);
  // Clean up
  mm_free(mapper_search);
}
GEM_INLINE void mapper_SE_run(mapper_parameters_t* const mapper_parameters) {
  mapper_run(mapper_parameters,false);
}
GEM_INLINE void mapper_PE_run(mapper_parameters_t* const mapper_parameters) {
  mapper_run(mapper_parameters,true);
}
