/*
 * PROJECT: GEMMapper
 * FILE: mapper.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "mapper.h"
#include "archive_search.h"

/*
 * Mapper Search
 */
typedef struct {
  // Thread Info
  uint64_t thread_id;
  pthread_t* thread_data;
  /* I/O */
  buffered_input_file_t* buffered_fasta_input;
  buffered_output_file_t* buffered_output_file;
  /* Mapper parameters */
  mapper_parameters_t* mapper_parameters;
  /* Archive search */
  archive_search_t* archive_search;
  matches_t* matches;
  mm_search_t* mm_search;
  /* Ticker */
  ticker_t* ticker;
} mapper_search_t;

/*
 * Mapper error-report function
 */
mapper_search_t* g_mapper_searches; // Global searches on going
void mapper_error_report(FILE* stream) {

  // TODO Implement 1 only call of this function, otherwise kill
  // TODO Patch this with crazy checkers that nothing is null

  const uint64_t threads_id = gem_thread_get_thread_id();
  if (threads_id==0) {
    fprintf(stream,"GEM::RunnigThread (threadID = MASTER)\n");
  } else {
    mapper_search_t* const mapper_search = g_mapper_searches + (threads_id-1);
    fprintf(stream,"GEM::RunnigThread (threadID = %lu)\n",mapper_search->thread_id);
    // Dump FASTA/FASTQ read
    const sequence_t* const sequence = archive_search_get_sequence(mapper_search->archive_search);
    if (!string_is_null(&sequence->tag) && !string_is_null(&sequence->read)) {
      const bool has_qualities = sequence_has_qualities(sequence);
      char* const end_tag =
          (sequence->attributes.end_info == PAIRED_END1) ? "/1" :
        ( (sequence->attributes.end_info == PAIRED_END2) ? "/2" : "" );
      fprintf(stream,"GEM::Sequence (File '%s' Line '%lu')\n",
          input_file_get_file_name(mapper_search->buffered_fasta_input->input_file),
          mapper_search->buffered_fasta_input->current_line_num - (has_qualities ? 4 : 2));
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
      fprintf(stream,"GEM::Sequence <<Empty>>\n");
    }
    // TODO ... More useful info
  }
}
/*
 * Mapper Parameters
 */
GEM_INLINE void mapper_parameters_set_defaults(mapper_parameters_t* const mapper_parameters) {
  /* CMD line */
  mapper_parameters->argc = 0;
  mapper_parameters->argv = NULL;
  /* Mapper Type */
  mapper_parameters->mapper_type=mapper_se;
  /* I/O */
  mapper_parameters->index_file_name=NULL;
  mapper_parameters->check_index=false;
  mapper_parameters->input_file_name=NULL;
  mapper_parameters->input_compression=FM_REGULAR_FILE;
  mapper_parameters->output_file_name=NULL;
  mapper_parameters->output_compression=FM_REGULAR_FILE;
  /* I/O */
  mapper_parameters->archive = NULL;
  mapper_parameters->input_file = NULL;
  mapper_parameters->output_stream = NULL;
  mapper_parameters->output_file = NULL;
  const uint64_t num_processors = system_get_num_processors();
  mapper_parameters->max_output_buffers = num_processors + (num_processors/2);
  mapper_parameters->output_format = MAP;
  output_sam_parameters_set_defaults(&mapper_parameters->sam_parameters);
  /* Search Parameters */
  approximate_search_parameters_init(&mapper_parameters->search_parameters);
  /* Select Parameters */
  archive_select_parameters_init(&mapper_parameters->select_parameters);
  /* System */
  mapper_parameters->num_threads=num_processors;
  mapper_parameters->max_memory=0;
  mapper_parameters->tmp_folder=NULL;
  /* Miscellaneous */
  mapper_parameters->stats = false;
  mapper_parameters->verbose_user=true;
  mapper_parameters->verbose_dev=false;
  /* Extras */
}
GEM_INLINE void mapper_parameters_print(
    FILE* const stream,mapper_parameters_t* const parameters,const bool dump_index_info) {
  tab_fprintf(stream,"[GEM]>Mapper.parameters\n");
  /*
   * Mapper Type
   */
  switch (parameters->mapper_type) {
    case mapper_se:
      tab_fprintf(stream,"  => Mapper.type\tMapperSE\n");
      break;
    case mapper_pe:
      tab_fprintf(stream,"  => Mapper.type\ttMapperPE\n");
      break;
    case mapper_se_cuda:
      tab_fprintf(stream,"  => Mapper.type\ttMapperSE-CUDA\n");
      break;
    case mapper_graph:
      tab_fprintf(stream,"  => Mapper.type\tMapper-Graph\n");
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  /*
   * I/O
   */
  tab_fprintf(stream,"  => I/O\n");
  // Index
  if (dump_index_info) {
    tab_global_inc();
    archive_print(stream,parameters->archive);
    tab_global_dec();
  } else {
    tab_fprintf(stream,"    => Index\t%s (size=%2.1fMB,length=%lu)\n",
        parameters->index_file_name,CONVERT_B_TO_MB(archive_get_size(parameters->archive)),
        archive_get_index_length(parameters->archive));
  }

  // Input // TODO PE
  if (parameters->input_file_name!=NULL) {
    tab_fprintf(stream,"    => Input\t%s (size=%2.1fMB)",
        parameters->input_file_name,CONVERT_B_TO_MB(input_file_get_size(parameters->input_file)));
  } else {
    tab_fprintf(stream,"    => Input\t<<stdin>>");
  }
  if (parameters->fastq_strictly_normalized) fprintf(stream,"\t[Normalize=ON]");
  if (parameters->fastq_try_recovery) fprintf(stream,"\t[Recovery=ON]");
  fprintf(stream,"\n");
  // Output
  tab_fprintf(stream,"    => Output\t%s",
      (parameters->output_file_name!=NULL) ? parameters->output_file_name : "<<stdout>>");
  switch (parameters->output_compression) {
    case FM_GZIPPED_FILE:
      fprintf(stream,"\t\[GZIP]n");
      break;
    case FM_BZIPPED_FILE:
      fprintf(stream,"\t[BZIP]\n");
      break;
    default:
      fprintf(stream,"\n");
      break;
  }
  /*
   * Search Parameters
   */
  tab_fprintf(stream,"  => Search Parameters\n");
  // Mapping Strategy

  // Quality Model

  // Error Model

  // Replacements

  // Reporting

  // Alignment Model


  /*
   * Select Parameters
   */
  tab_fprintf(stream,"  => Select Parameters\n");

  /*
   * System
   */
  /* System */
  tab_fprintf(stream,"  => System\n");
  tab_fprintf(stream,"    => Num.Threads\t%lu\n",parameters->num_threads);
  if (parameters->max_memory!=0) {
    tab_fprintf(stream,"    => Max.Memory\t%lu\n",parameters->max_memory);
  }
  if (parameters->tmp_folder!=NULL) {
    tab_fprintf(stream,"    => Tmp.path\t%s\n",parameters->tmp_folder);
  }

  // TODO All


  /* Select Parameters */
  /* Parameters */

  /*
  #define GEM_STATS_SHOW_SE_PARAMETERS(FILE,MAPPER_PARAMS) \
    fprintf(FILE, "PARAMETERS.SE.Mapping\n"); \
    if (MAPPER_PARAMS.percentage_mismatches<0.0) { \
    fprintf(FILE, "  Mismatches          %ld\n", MAPPER_PARAMS.max_mismatches); \
    } else { \
    fprintf(FILE, "  Mismatches          %1.2f\n", MAPPER_PARAMS.percentage_mismatches); \
    } \
    if (MAPPER_PARAMS.percentage_differences<0.0) { \
    fprintf(FILE, "  Differences         %ld\n", MAPPER_PARAMS.max_differences); \
    } else { \
    fprintf(FILE, "  Differences         %1.2f\n", MAPPER_PARAMS.percentage_differences); \
    } \
    fprintf(FILE, "  Delta               %lu\n", MAPPER_PARAMS.delta); \
    fprintf(FILE, "  Max indel length    %lu\n", MAPPER_PARAMS.max_indel_length); \
    fprintf(FILE, "  Replacements        {%s}\n", MAPPER_PARAMS.mismatch_alphabet); \
    if (MAPPER_PARAMS.percentage_incomplete_strata<0.0) { \
    fprintf(FILE, "  Incomplete Strata   %ld\n", MAPPER_PARAMS.incomplete_strata); \
    } else { \
    fprintf(FILE, "  Incomplete Strata   %1.2f\n", MAPPER_PARAMS.percentage_incomplete_strata); \
    }

  #define GEM_STATS_SHOW_PE_PARAMETERS(FILE,PAIR_PARAMS) \
    fprintf(FILE, "PARAMETERS.PE.Mapping\n"); \
    fprintf(FILE, "  Map.Both.Ends          %s\n", GEM_PRINT_BOOL(PAIR_PARAMS.map_both_ends)); \
    if (PAIR_PARAMS.percentage_differences<0.0) \
    fprintf(FILE, "  Mismatches             %ld\n", PAIR_PARAMS.max_differences); \
    else \
    fprintf(FILE, "  Mismatches             %1.2f\n", PAIR_PARAMS.percentage_differences); \
    fprintf(FILE, "  Max.Insertion.Size     %lu\n", PAIR_PARAMS.max_ins_size); \
    fprintf(FILE, "  Min.Insertion.Size     %lu\n", PAIR_PARAMS.min_ins_size); \
    fprintf(FILE, "  Max.Number.Pairs       %lu\n", PAIR_PARAMS.max_extended_matches)
    */


}
/*
 * I/O
 */
GEM_INLINE void mapper_SE_output_matches(
    const mapper_parameters_t* const parameters,
    buffered_output_file_t* const buffered_output_file,
    sequence_t* const seq_read,matches_t* const matches) {
  switch (parameters->output_format) {
    case MAP:
      output_map_single_end_matches(buffered_output_file,seq_read,matches);
      break;
    case SAM:
      output_sam_single_end_matches(buffered_output_file,seq_read,matches,
          &parameters->sam_parameters);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
void* mapper_SE_thread(mapper_search_t* const mapper_search) {
  // GEM-thread error handler
  gem_thread_register_id(mapper_search->thread_id+1);

  // Create new buffered reader/writer
  const mapper_parameters_t* const parameters = mapper_search->mapper_parameters;
  mapper_search->buffered_fasta_input = buffered_input_file_new(parameters->input_file);
  mapper_search->buffered_output_file = buffered_output_file_new(parameters->output_file);
  buffered_input_file_attach_buffered_output(mapper_search->buffered_fasta_input,mapper_search->buffered_output_file);

  // Create an Archive-Search
  mapper_search->mm_search = mm_search_new(mm_pool_get_slab(mm_pool_2MB));
  mapper_search->archive_search = archive_search_new(
      parameters->archive,&mapper_search->mapper_parameters->search_parameters,
      &mapper_search->mapper_parameters->select_parameters);
  archive_search_configure(mapper_search->archive_search,mapper_search->mm_search);
  mapper_search->matches = matches_new(mapper_search->mm_search->text_collection);
  matches_configure(mapper_search->matches,mapper_search->archive_search->text_collection);

  // FASTA/FASTQ reading loop
  uint64_t reads_processed = 0;
  error_code_t error_code;
  while ((error_code=input_fasta_parse_sequence(
      mapper_search->buffered_fasta_input,archive_search_get_sequence(mapper_search->archive_search),
      parameters->fastq_strictly_normalized,parameters->fastq_try_recovery))) { // TODO: Check INPUT_STATUS_FAIL & check skip empty lines

    // Search into the archive
    archive_search_single_end(mapper_search->archive_search,mapper_search->matches);

    // Select matches
    archive_select_matches(mapper_search->archive_search,mapper_search->matches);
    // archive_search_score_matches(archive_search); // TODO

    // Output matches
    mapper_SE_output_matches(
        parameters,mapper_search->buffered_output_file,
        archive_search_get_sequence(mapper_search->archive_search),mapper_search->matches);

    // Update processed
    if (++reads_processed == MAPPER_TICKER_STEP) {
      ticker_update_mutex(mapper_search->ticker,reads_processed);
      reads_processed=0;
    }

    // Clear
    mm_search_clear(mapper_search->mm_search);
    matches_clear(mapper_search->matches);
  }
  // Update processed
  ticker_update_mutex(mapper_search->ticker,reads_processed);

  // Clean up
  buffered_input_file_close(mapper_search->buffered_fasta_input);
  buffered_output_file_close(mapper_search->buffered_output_file);
  archive_search_delete(mapper_search->archive_search);
  matches_delete(mapper_search->matches);
  mm_search_delete(mapper_search->mm_search);

  pthread_exit(0);
}
GEM_INLINE void mapper_SE_run(mapper_parameters_t* const mapper_parameters) {
  // Setup threads
  const uint64_t num_threads = mapper_parameters->num_threads;
  mapper_search_t* const mapper_search = mm_calloc(num_threads,mapper_search_t,false); // Allocate mapper searches
  g_mapper_searches = mapper_search; // Set global searches for error reporting
  // Prepare output file (SAM headers)
  if (mapper_parameters->output_format==SAM) {
    output_sam_print_header(mapper_parameters->output_file,
        mapper_parameters->archive,mapper_parameters->argc,mapper_parameters->argv);
  }
  // Setup Ticker
  ticker_t ticker;
  ticker_count_reset(&ticker,mapper_parameters->verbose_user,"Mapping Sequences",0,MAPPER_TICKER_STEP,true);
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
        pthread_create(mapper_search[i].thread_data,0,(void* (*)(void*))mapper_SE_thread,(void*)(mapper_search+i)),
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
/*
 * PE Mapper
 */
GEM_INLINE void mapper_PE_run(const mapper_parameters_t* const mapper_parameters) {
  GEM_NOT_IMPLEMENTED(); // TODO
}
