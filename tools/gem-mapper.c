/*
 * PROJECT: GEMMapper
 * FILE: gem-mapper.c
 * DATE:5/12/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Genomic Reads Mapper
 */

#include "gem_core.h"

/*
 * GEM-indexer Errors
 */
#define GEM_ERROR_MAPPER_PARSING_MAX_MEMORY "Error parsing parameter --max-memory. '%s' not a valid measure (Eg 2GB)"

/*
 * GEM-Mapper Threads info
 */
typedef struct {
  // Thread Info
  uint64_t id;
  pthread_t* thread;
  // I/O
  buffered_input_file_t* buffered_fasta_input;
  buffered_output_file_t* buffered_output_file;
  // Sequence
  sequence_t* seq_read;
} gem_mapper_thread_t;
// Thread Global Info
gem_mapper_thread_t* gem_mapper_threads = NULL;
/*
 * GEM-Mapper error report function
 */
void gem_mapper_error_report(FILE* stream) {
  const uint64_t threads_id = gem_thread_get_thread_id();
  if (threads_id==0) {
    fprintf(stream,"GEM::RunnigThread (threadID = MASTER)\n");
  } else {
    gem_mapper_thread_t* const thread = gem_mapper_threads + (threads_id-1);
    fprintf(stream,"GEM::RunnigThread (threadID = %lu)\n",thread->id);
    // Dump FASTA/FASTQ read
    if (!string_is_null(thread->seq_read->tag) && !string_is_null(thread->seq_read->read)) {
      const bool has_qualities = sequence_has_qualities(thread->seq_read);
      char* const end_tag = (thread->seq_read->attributes.end_info == PAIRED_END1) ? "/1" :
          ( (thread->seq_read->attributes.end_info == PAIRED_END2) ? "/2" : "" );
      fprintf(stream,"GEM::Sequence (File '%s' Line '%lu')\n",
          input_file_get_file_name(thread->buffered_fasta_input->input_file),
          thread->buffered_fasta_input->current_line_num - (has_qualities ? 4 : 2));
      if (has_qualities) {
        fprintf(stream,"@%"PRIs"%s\n%"PRIs"\n+\n%"PRIs"\n",
            PRIs_content(thread->seq_read->tag),end_tag,
            PRIs_content(thread->seq_read->read),
            PRIs_content(thread->seq_read->qualities));
      } else {
        fprintf(stream,">%"PRIs"%s\n%"PRIs"\n",
            PRIs_content(thread->seq_read->tag),end_tag,
            PRIs_content(thread->seq_read->read));
      }
    } else {
      fprintf(stream,"GEM::Sequence <<Empty>>\n");
    }
    // TODO ... More useful info
  }
}
/*
 * GEM-Mapper Parameters
 */
typedef struct {
  /* I/O */
  char *index_file_name;
  bool check_index;
  char *input_file_name;
  fm_type input_compression;
  char *output_file_name;
  fm_type output_compression;
  /* I/O Attributes */
  bool fastq_strictly_normalized;
  bool fastq_try_recovery;
  /* Qualities */
  quality_format_t quality_format;
  quality_model_t quality_model;
  uint64_t quality_threshold;
    // TODO quality_score
  /* Single-end Alignment */
  mapping_mode_t mapping_mode;
  float mapping_degree;
  float max_search_error;
  float max_filtering_error;
  float complete_strata_after_best;
  float min_matching_length;
  uint64_t max_search_matches;
  char* mismatch_alphabet;
  /* Paired-end Alignment */
  /* Reporting */
  file_format_t output_format;
  uint64_t min_decoded_strata;
  uint64_t max_decoded_matches;
  uint64_t min_reported_matches;
  uint64_t max_reported_matches;
  /* System */
  uint64_t num_threads;
  uint64_t max_memory;
  char* tmp_folder;
  /* Debug/Temporal */
  /* Miscellaneous */
  bool verbose;
  /* Extras */
} gem_mapper_parameters_t;
// Defaults
gem_mapper_parameters_t parameters = {
    /* I/O */
    .index_file_name=NULL,
    .check_index=false,
    .input_file_name=NULL,
    .input_compression=FM_REGULAR_FILE,
    .output_file_name=NULL,
    .output_compression=FM_REGULAR_FILE,
    /* I/O Attributes */
    .fastq_strictly_normalized=true,
    .fastq_try_recovery=false,
    /* Qualities */
    .quality_format=qualities_ignore,
    .quality_model=quality_model_type_gem,
    .quality_threshold=26,
       // TODO quality_score
    /* Single-end Alignment */
    .mapping_mode=mapping_incremental_mapping,
    .mapping_degree=0,
    .max_search_error=0.04,
    .max_filtering_error=0.20,
    .complete_strata_after_best=0.0,
    .min_matching_length=0.20,
    .max_search_matches=ALL,
    .mismatch_alphabet="ACGT",
    /* Paired-end Alignment */
    /* Reporting */
    .output_format=MAP,
    .min_decoded_strata=0,
    .max_decoded_matches=20,
    .min_reported_matches=1,
    .max_reported_matches=ALL,
    /* System */
    .num_threads=1,
    .max_memory=0,
    .tmp_folder=NULL,
    /* Miscellaneous */
    .verbose=false,
    /* Extras */
};
/*
 * GEM-Mapper Data
 */
typedef struct {
  /* I/O */
  archive_t* gem_archive;
  input_file_t* input_file;
  FILE* output_stream;
  output_file_t* output_file;
  /* Single-end Alignment */
  /* Paired-end Alignment */
  /* Reporting */
  /* Miscellaneous */
  /* Extras */
} gem_mapper_data_t;
// Defaults
gem_mapper_data_t mapper_data = {
    /* I/O */
    .gem_archive = NULL,
    .input_file = NULL,
    .output_stream = NULL,
    .output_file = NULL,
    /* Single-end Alignment */
    /* Paired-end Alignment */
    /* Reporting */
    /* Miscellaneous */
    /* Extras */
};
/*
 * Configure search parameters
 */
GEM_INLINE void gem_mapper_configure_search_parameters(approximate_search_parameters_t* const search_parameters) {
  // Configure search
  approximate_search_configure_mapping_strategy(search_parameters,
      parameters.mapping_mode,parameters.mapping_degree);
  approximate_search_configure_quality_model(search_parameters,
      parameters.quality_model,parameters.quality_format,parameters.quality_threshold);
  approximate_search_configure_error_model(search_parameters,
      parameters.max_search_error,parameters.max_filtering_error,
      parameters.complete_strata_after_best,parameters.min_matching_length);
  approximate_search_configure_replacements(search_parameters,
      parameters.mismatch_alphabet,gem_strlen(parameters.mismatch_alphabet));
  approximate_search_configure_matches(search_parameters,parameters.max_search_matches);
}
/*
 * GEM-Mapper SE Mapping
 */
GEM_INLINE void gem_mapper_output_matches(
    buffered_output_file_t* const buffered_output_file,
    sequence_t* const seq_read,matches_t* const matches) {
  switch (parameters.output_format) {
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
void* gem_mapper_thread(uint64_t thread_id) {
  // GEM-thread error handler
  gem_thread_register_id(thread_id+1);
  gem_mapper_thread_t* const mapper_thread = gem_mapper_threads+thread_id;

  // Create new buffered reader/writer
  mapper_thread->seq_read = sequence_new(); // Allocate Sequence
  mapper_thread->buffered_fasta_input = buffered_input_file_new(mapper_data.input_file);
  mapper_thread->buffered_output_file = buffered_output_file_new(mapper_data.output_file);
  buffered_input_file_attach_buffered_output(mapper_thread->buffered_fasta_input,mapper_thread->buffered_output_file);

  // Create an Archive-Search
  archive_search_t* const archive_search = archive_search_new(mapper_data.gem_archive);
  gem_mapper_configure_search_parameters(archive_search_get_search_parameters(archive_search));

  // FASTA/FASTQ reading loop
  error_code_t error_code;
  while ((error_code=input_fasta_parse_sequence(
      mapper_thread->buffered_fasta_input,mapper_thread->seq_read,
      parameters.fastq_strictly_normalized,parameters.fastq_try_recovery))) {

    // Search into the archive
    archive_search_single_end(archive_search,mapper_thread->seq_read);

    // Decode matches
    archive_search_decode_matches(archive_search,
        parameters.max_decoded_matches,parameters.min_decoded_strata,
        parameters.min_reported_matches,parameters.max_reported_matches);
    // archive_search_score_matches(archive_search); // TODO

    // Output matches
    gem_mapper_output_matches(mapper_thread->buffered_output_file,
        mapper_thread->seq_read,archive_search_get_matches(archive_search));
  }

  // Clean up
  buffered_input_file_close(mapper_thread->buffered_fasta_input);
  buffered_output_file_close(mapper_thread->buffered_output_file);
  sequence_delete(mapper_thread->seq_read);

  pthread_exit(0);
}
GEM_INLINE void gem_mapper_launch_mapping_threads() {
  uint64_t i;
  // Allocate threads info
  gem_mapper_threads = mm_malloc(parameters.num_threads*sizeof(gem_mapper_thread_t));
  // Setup threads & launch them
  for (i=0;i<parameters.num_threads;++i) {
    // Allocate thread info
    gem_mapper_threads[i].id = i;
    gem_mapper_threads[i].thread = mm_alloc(pthread_t);
    // Launch thread
    gem_cond_fatal_error__perror(
        pthread_create(gem_mapper_threads[i].thread,0,(void* (*)(void*))gem_mapper_thread,(void*)i),
        SYS_THREAD_CREATE);
  }
  // Join all threads
  for (i=0;i<parameters.num_threads;++i) {
    gem_cond_fatal_error__perror(pthread_join(*(gem_mapper_threads[i].thread),0),SYS_THREAD_JOIN);
    mm_free(gem_mapper_threads[i].thread);
  }
  // Free & return
  mm_free(gem_mapper_threads);
}
/*
 * GEM-mapper I/O related functions
 */
GEM_INLINE void gem_mapper_load_index() {
  // Load archive
  mapper_data.gem_archive = archive_read(parameters.index_file_name,parameters.check_index,parameters.verbose);
}
GEM_INLINE void gem_mapper_open_input() {
  // Open input file
  if (parameters.input_file_name==NULL) {
    switch (parameters.input_compression) {
      case FM_GZIPPED_FILE:
        mapper_data.input_file = input_gzip_stream_open(stdin);
        break;
      case FM_BZIPPED_FILE:
        mapper_data.input_file = input_bzip_stream_open(stdin);
        break;
      default:
        mapper_data.input_file = input_stream_open(stdin);
        break;
    }
  } else {
     mapper_data.input_file = input_file_open(parameters.input_file_name,false);
  }
  // TODO: Checks
}
GEM_INLINE void gem_mapper_close_input() {
  input_file_close(mapper_data.input_file);
}
GEM_INLINE void gem_mapper_open_output() {
  // Open output stream
  mapper_data.output_stream = (parameters.output_file_name==NULL) ?
      stdout : fopen(parameters.output_file_name,"w");
  // Open output file // TODO Unsorted
  switch (parameters.output_compression) {
    case FM_GZIPPED_FILE:
      mapper_data.output_file = output_gzip_stream_new(mapper_data.output_stream,SORTED_FILE);
      break;
    case FM_BZIPPED_FILE:
      mapper_data.output_file = output_bzip_stream_new(mapper_data.output_stream,SORTED_FILE);
      break;
    default:
      mapper_data.output_file = output_stream_new(mapper_data.output_stream,SORTED_FILE);
      break;
  }
}
GEM_INLINE void gem_mapper_close_output() {
  output_file_close(mapper_data.output_file);
  if (parameters.output_file_name!=NULL) fclose(mapper_data.output_stream);
}
/*
 * GEM-Mapper options Menu
 */
option_t gem_mapper_options[] = {
  /* I/O */
  { 'I', "index", REQUIRED, TYPE_STRING, 2 , true, "<index_file.gem>" , "" },
  { 204, "check-index", NO_ARGUMENT, TYPE_NONE, 2 , false, "" , "" },
  { 'i', "input", REQUIRED, TYPE_STRING, 2 , true, "<file>" , "(FASTA/FASTQ, default=stdin)" },
  { 200, "gzip-input", NO_ARGUMENT, TYPE_NONE, 2, true, "", "(gzip input)" },
  { 201, "bzip-input", NO_ARGUMENT, TYPE_NONE, 2, true, "", "(bzip input)" },
  { 'o', "output", REQUIRED, TYPE_STRING, 2 , true, "<output_prefix>" , "(default=stdout)" },
  { 202, "gzip-output", NO_ARGUMENT, TYPE_NONE, 2, true, "", "(gzip output)" },
  { 203, "bzip-output", NO_ARGUMENT, TYPE_NONE, 2, true, "", "(bzip output)" },
  /* Qualities */
  { 'q', "quality-format", REQUIRED, TYPE_STRING, 3, true, "'ignore'|'offset-33'|'offset-64'", "(default=offset-33)" },
  { 'Q', "quality-model", REQUIRED, TYPE_STRING, 3, false, "'gem'|'flat'", "(default=gem)" },
  { 300, "gem-quality-threshold", REQUIRED, TYPE_INT, 3, true, "<number>", "(default=26, that is e<=2e-3)" },
  /* Single-end Alignment */
  { 400, "mapping-mode", REQUIRED, TYPE_INT, 4 , false, "'incremental'|'adaptive'|'fixed'|'fast'|'brute-force'" , "(default=fast)" },
  { 401, "mapping-degree", REQUIRED, TYPE_INT, 4 , false, "<number|percentage>" , "(default=0)" },
  { 'e', "max-search-error", REQUIRED, TYPE_INT, 4 , true, "<number|percentage>" , "(default=0.04, 4%)" },
  { 'E', "max-filtering-error", REQUIRED, TYPE_INT, 4 , false, "<number|percentage>" , "(default=0.2, 20%)" },
  { 's', "complete-strata-after-best", REQUIRED, TYPE_INT, 4 , true, "<number|percentage>" , "(default=0)" },
  { 402, "min-matching-length", REQUIRED, TYPE_INT, 4 , false, "<number|percentage>" , "(default=0.20, 20%)" },
  { 403, "max-search-matches", REQUIRED, TYPE_INT, 4 , true, "<number>" , "(not limited by default)" },
  { 404, "mismatch-alphabet", REQUIRED, TYPE_STRING, 4 , false, "<symbols>" , "(default=\"ACGT\")" },
  /* Paired-end Alignment */
  /* Reporting */
  { 'F', "output-format", REQUIRED, TYPE_STRING, 6 , false, "'MAP'|'SAM'" , "(default=MAP)" },
  { 'D', "min-decoded-strata", REQUIRED, TYPE_INT, 6 , false, "<number>" , "(stratum-wise, default=0)" },
  { 'd', "max-decoded-matches", REQUIRED, TYPE_INT, 6 , false, "<number>|'all'" , "(stratum-wise, default=20)" },
  { 'm', "min-reported-matches", REQUIRED, TYPE_INT, 6 , false, "<number>|'all'" , "(default=1)" },
  { 'M', "max-reported-matches", REQUIRED, TYPE_INT, 6 , true, "<number>|'all'" , "(default=all)" },
  /* System */
  { 't', "threads", REQUIRED, TYPE_STRING, 7 , true, "<number>" , "(default=1)" },
  { 700, "max-memory", REQUIRED, TYPE_STRING, 7 , true, "<maximum-memory>" , "(Eg 2GB)" },
  { 701, "tmp-folder", REQUIRED, TYPE_STRING, 7 , true, "<temporal_dir_path>" , "(/tmp/)" },
  /* Debug */
  /* Miscellaneous */
  { 'v', "verbose", NO_ARGUMENT, TYPE_NONE, 9 , true, "" , "(disabled)" },
  { 'h', "help", NO_ARGUMENT, TYPE_NONE, 9 , true, "" , "(print usage)" },
  { 'H', "help", NO_ARGUMENT, TYPE_NONE, 9 , false, "" , "(print usage + extras)" },
  /* Extras */
  {  0, "", 0, 0, 0, false, "", ""}
};
char* gem_mapper_groups[] = {
  /*  0 */ "Null",
  /*  1 */ "Unclassified",
  /*  2 */ "I/O",
  /*  3 */ "Qualities",
  /*  4 */ "Single-end Alignment",
  /*  5 */ "Paired-end Alignment",
  /*  6 */ "Reporting",
  /*  7 */ "System",
  /*  8 */ "Debug",
  /*  9 */ "Miscellaneous",
  /* 10 */ "Extras"
};
void usage(const bool print_inactive) {
  fprintf(stderr, "USAGE: ./gem-mapper [ARGS]...\n");
  options_fprint_menu(stderr,gem_mapper_options,gem_mapper_groups,false,print_inactive);
}
void parse_arguments(int argc,char** argv) {
  struct option* getopt_options = options_adaptor_getopt(gem_mapper_options);
  string_t* const getopt_short_string = options_adaptor_getopt_short(gem_mapper_options);
  char* const getopt_short = string_get_buffer(getopt_short_string);
  int option, option_index;
  while (true) {
    // Get option &  Select case
    if ((option=getopt_long(argc,argv,getopt_short,getopt_options,&option_index))==-1) break;
    switch (option) {
    /* I/O */
    case 'I': // --index
      parameters.index_file_name = optarg;
      break;
    case 204: // --check-index
      parameters.check_index = true;
      break;
    case 'i': // --input
      parameters.input_file_name = optarg;
      break;
    case 200: // --gzip-input
      parameters.input_compression = FM_GZIPPED_FILE;
      break;
    case 201: // --bzip-input
      parameters.input_compression = FM_BZIPPED_FILE;
      break;
    case 'o': // --output
      parameters.output_file_name = optarg;
      break;
    case 202: // --gzip-output
      parameters.output_compression = FM_GZIPPED_FILE;
      break;
    case 203: // --bzip-output
      parameters.output_compression = FM_BZIPPED_FILE;
      break;
    /* Qualities */
    case 'q': // --quality-format
      if (gem_strcaseeq(optarg,"ignore")) {
        parameters.quality_format=qualities_ignore;
      } else if (gem_strcaseeq(optarg,"offset-33")) {
        parameters.quality_format=qualities_offset_33;
      } else if (gem_strcaseeq(optarg,"offset-64")) {
        parameters.quality_format=qualities_offset_64;
      } else {
        gem_fatal_error_msg("Option '-q|--quality-format' must be 'ignore', 'offset-33' or 'offset-64'");
      }
      break;
    case 'Q': // --quality-model
      if (gem_strcaseeq(optarg,"gem")) {
        parameters.quality_model=quality_model_type_gem;
        // TODO parameters.quality_score=gem_quality_score;
      } else if (gem_strcaseeq(optarg,"flat")) {
        parameters.quality_model=quality_model_type_flat;
        // TODO parameters.quality_score=flat_quality_score;
      } else
      gem_fatal_error_msg("Option '-Q|--quality-model' must be 'gem' or 'flat'");
      break;
    case 300: // --gem-quality-threshold
      parameters.quality_threshold = atoll(optarg);
      break;
    /* Single-end Alignment */
    case 400: // --mapping-mode
      if (gem_strcaseeq(optarg,"incremental")) {
        parameters.mapping_mode = mapping_incremental_mapping;
        break;
      }
      if (gem_strcaseeq(optarg,"adaptive")) {
        parameters.mapping_mode = mapping_adaptive_filtering;
        break;
      }
      if (gem_strcaseeq(optarg,"fixed")) {
        parameters.mapping_mode = mapping_fixed_filtering;
        break;
      }
      if (gem_strcaseeq(optarg,"fast")) {
        parameters.mapping_mode = mapping_fast;
        break;
      }
      if (gem_strcaseeq(optarg,"brute-force")) {
        parameters.mapping_mode = mapping_neighborhood_search;
        break;
      }
      // TODO error
      break;
    case 401: // --mapping-degree
      parameters.mapping_degree = atof(optarg);
      break;
    case 'e': // --max-search-error
      parameters.max_search_error = atof(optarg);
      break;
    case 'E': // --max-filtering-error
      parameters.max_filtering_error = atof(optarg);
      break;
    case 's': // --complete-strata-after-best
      parameters.complete_strata_after_best = atof(optarg);
      break;
    case 402: // --min-matching-length
      parameters.min_matching_length = atof(optarg);
      break;
    case 403: // --max-search-matches
      parameters.max_search_matches = atol(optarg);
      break;
    case 404: // --mismatch-alphabet
      parameters.mismatch_alphabet = optarg;
      break;
    /* Paired-end Alignment */
    /* Reporting */
    case 'F': // --output-format
      if (gem_strcaseeq(optarg,"MAP")) {
        parameters.output_format = MAP;
      } else if (gem_strcaseeq(optarg,"SAM")) {
        parameters.output_format = SAM;
      } else {
        gem_fatal_error_msg("Option '-F|--output-format' must be 'MAP' or 'SAM'");
      }
      break;
      break;
    case 'D': // --min-decoded-strata
      parameters.min_decoded_strata = (gem_strcaseeq(optarg,"all")) ? ALL : atol(optarg);
      break;
    case 'd': // --max-decoded-matches
      parameters.max_decoded_matches = (gem_strcaseeq(optarg,"all")) ? ALL : atol(optarg);
      break;
    case 'm': // --min-reported-matches
      parameters.min_reported_matches = (gem_strcaseeq(optarg,"all")) ? ALL : atol(optarg);
      break;
    case 'M': // --max-reported-matches
      parameters.max_reported_matches = (gem_strcaseeq(optarg,"all")) ? ALL : atol(optarg);
      break;
    /* System */
    case 't': // --threads
      parameters.num_threads = atol(optarg);
      break;
    case 700: // --max-memory
      gem_cond_fatal_error(input_text_parse_size(optarg,&(parameters.max_memory)),PARSING_SIZE,"--max-memory",optarg);
      break;
    case 701: // --tmp-folder
      parameters.tmp_folder = optarg;
      break;
    /* Debug */
    /* Misc */
    case 'v':
      parameters.verbose = true;
      break;
    case 'h':
      usage(false);
      exit(1);
    case 'H':
      usage(true);
      exit(1);
    case '?':
    default:
      gem_fatal_error_msg("Option not recognized");
    }
  }
  /*
   * Parameters Check
   */
  // System
  if (parameters.max_memory==0) {
    parameters.max_memory = mm_get_available_mem();
  }
  // Free
  string_delete(getopt_short_string);
  mm_free(getopt_options);
}

int main(int argc,char** argv) {
  // Parsing command-line options
  parse_arguments(argc,argv); // Parse cmd-line

  // GEM Runtime setup
  gem_runtime_init(parameters.num_threads,parameters.max_memory,parameters.tmp_folder,gem_mapper_error_report);

  // Open Input/Output File(s)
  gem_mapper_open_input();
  gem_mapper_open_output();

  // Load GEM-Index
  gem_mapper_load_index();

  // Launch mapping threads
  gem_mapper_launch_mapping_threads();

  // CleanUP
  // TODO Clean index
  // TODO clean parameters
  gem_mapper_close_input();  // Close I/O files
  gem_mapper_close_output(); // Close I/O files

  return 0;
}

