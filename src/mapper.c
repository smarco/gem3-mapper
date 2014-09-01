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
} mapper_search_t;

/*
 * Mapper error-report function
 */
mapper_search_t* g_mapper_searches; // Global searches on going
void mapper_error_report(FILE* stream) {
  const uint64_t threads_id = gem_thread_get_thread_id();
  if (threads_id==0) {
    fprintf(stream,"GEM::RunnigThread (threadID = MASTER)\n");
  } else {
    mapper_search_t* const mapper_search = g_mapper_searches + (threads_id-1);
    fprintf(stream,"GEM::RunnigThread (threadID = %lu)\n",mapper_search->thread_id);
    // Dump FASTA/FASTQ read
    const sequence_t* const sequence = archive_search_get_sequence(mapper_search->archive_search);
    if (!string_is_null(sequence->tag) && !string_is_null(sequence->read)) {
      const bool has_qualities = sequence_has_qualities(sequence);
      char* const end_tag =
          (sequence->attributes.end_info == PAIRED_END1) ? "/1" :
        ( (sequence->attributes.end_info == PAIRED_END2) ? "/2" : "" );
      fprintf(stream,"GEM::Sequence (File '%s' Line '%lu')\n",
          input_file_get_file_name(mapper_search->buffered_fasta_input->input_file),
          mapper_search->buffered_fasta_input->current_line_num - (has_qualities ? 4 : 2));
      if (has_qualities) {
        fprintf(stream,"@%"PRIs"%s\n%"PRIs"\n+\n%"PRIs"\n",
            PRIs_content(sequence->tag),end_tag,
            PRIs_content(sequence->read),
            PRIs_content(sequence->qualities));
      } else {
        fprintf(stream,">%"PRIs"%s\n%"PRIs"\n",
            PRIs_content(sequence->tag),end_tag,
            PRIs_content(sequence->read));
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
  /* I/O Attributes */
  mapper_parameters->fastq_strictly_normalized=true;
  mapper_parameters->fastq_try_recovery=false;
  /* Qualities */
  mapper_parameters->quality_format=qualities_ignore;
  mapper_parameters->quality_model=quality_model_type_gem;
  mapper_parameters->quality_threshold=26;
  /* Single-end Alignment */
  mapper_parameters->mapping_mode=mapping_incremental_mapping;
  mapper_parameters->mapping_degree=0;
  mapper_parameters->max_search_error=0.04;
  mapper_parameters->max_filtering_error=0.20;
  mapper_parameters->complete_strata_after_best=0.0;
  mapper_parameters->min_matching_length=0.20;
  mapper_parameters->max_search_matches=ALL;
  mapper_parameters->mismatch_alphabet="ACGT";
  /* Paired-end Alignment */
  /* Reporting */
  mapper_parameters->output_format=MAP;
  mapper_parameters->min_decoded_strata=0;
  mapper_parameters->max_decoded_matches=20;
  mapper_parameters->min_reported_matches=1;
  mapper_parameters->max_reported_matches=100;
  /* System */
  mapper_parameters->num_threads=1; // FIXME
  mapper_parameters->max_memory=0;
  mapper_parameters->tmp_folder=NULL;
  /* Miscellaneous */
  mapper_parameters->user_verbose=false;
  mapper_parameters->dev_verbose=false;
  /* Extras */
}
/*
 * Mapper SE Mapping
 */
GEM_INLINE void mapper_SE_configure_archive_search(
    archive_search_t* const archive_search,const mapper_parameters_t* const parameters) {
  approximate_search_parameters_t* const search_parameters = archive_search_get_search_parameters(archive_search);
  // Configure search
  approximate_search_configure_mapping_strategy(search_parameters,
      parameters->mapping_mode,parameters->mapping_degree);
  approximate_search_configure_quality_model(search_parameters,
      parameters->quality_model,parameters->quality_format,
      parameters->quality_threshold);
  approximate_search_configure_error_model(search_parameters,
      parameters->max_search_error,parameters->max_filtering_error,
      parameters->complete_strata_after_best,parameters->min_matching_length);
  approximate_search_configure_replacements(search_parameters,
      parameters->mismatch_alphabet,gem_strlen(parameters->mismatch_alphabet));
  approximate_search_configure_matches(search_parameters,parameters->max_search_matches);
}
GEM_INLINE void mapper_SE_output_matches(
    const mapper_parameters_t* const parameters,
    buffered_output_file_t* const buffered_output_file,
    const sequence_t* const seq_read,matches_t* const matches) {
  switch (parameters->output_format) {
    case MAP:
      output_map_single_end_matches(buffered_output_file,seq_read,matches);
      break;
    case SAM:
      GEM_NOT_IMPLEMENTED();
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
  mm_stack_t* const mm_stack = mm_stack_new(mm_pool_get_slab(mm_pool_2MB));
  mapper_search->archive_search = archive_search_new(parameters->archive,mm_stack);
  mapper_SE_configure_archive_search(mapper_search->archive_search,parameters);

  // FASTA/FASTQ reading loop
  error_code_t error_code;
  while ((error_code=input_fasta_parse_sequence(
      mapper_search->buffered_fasta_input,archive_search_get_sequence(mapper_search->archive_search),
      parameters->fastq_strictly_normalized,parameters->fastq_try_recovery))) {
    // TODO: Check INPUT_STATUS_FAIL

    // Search into the archive
    mm_stack_free(mm_stack); // Free stack for new use
    archive_search_single_end(mapper_search->archive_search);

    // Select matches
    archive_search_select_matches(mapper_search->archive_search,
        parameters->max_decoded_matches,parameters->min_decoded_strata,
        parameters->min_reported_matches,parameters->max_reported_matches);
    // archive_search_score_matches(archive_search); // TODO

    // Output matches
    mapper_SE_output_matches(
        parameters,mapper_search->buffered_output_file,
        archive_search_get_sequence(mapper_search->archive_search),
        archive_search_get_matches(mapper_search->archive_search));

  }

  // Clean up
  buffered_input_file_close(mapper_search->buffered_fasta_input);
  buffered_output_file_close(mapper_search->buffered_output_file);
  archive_search_delete(mapper_search->archive_search);

  pthread_exit(0);
}
GEM_INLINE void mapper_SE_run(mapper_parameters_t* const mapper_parameters) {
  // Setup threads & launch
  const uint64_t num_threads = mapper_parameters->num_threads;
  PROF_NEW(num_threads+1); // Setup Profiling (Add master thread)
  mapper_search_t* const mapper_search = mm_malloc(num_threads*sizeof(mapper_search_t)); // Allocate mapper searches
  g_mapper_searches = mapper_search; // Set global searches for error reporting
  uint64_t i;
  for (i=0;i<num_threads;++i) {
    // Setup thread
    mapper_search[i].thread_id = i;
    mapper_search[i].thread_data = mm_alloc(pthread_t);
    mapper_search[i].mapper_parameters = mapper_parameters;
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
  // Clean up
  mm_free(mapper_search);
  PROF_DELETE(); // Clean-up Profiler
}
/*
 * PE Mapper
 */
GEM_INLINE void mapper_PE_run(const mapper_parameters_t* const mapper_parameters) {
  GEM_NOT_IMPLEMENTED(); // TODO
}
