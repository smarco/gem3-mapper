/*
 * PROJECT: GEMMapper
 * FILE: gem-mapper.c
 * DATE:5/12/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Genomic Read Mapper
 */

#include "gem_core.h"

/*
 * GEM-mapper I/O related functions
 */
GEM_INLINE void gem_mapper_load_index(mapper_parameters_t* const parameters) {
  // Load archive
  gem_cond_log(parameters->user_verbose,"[Loading GEM index '%s']",parameters->index_file_name);
  parameters->archive = archive_read(parameters->index_file_name,parameters->check_index,parameters->dev_verbose);
}
GEM_INLINE void gem_mapper_open_input(mapper_parameters_t* const parameters) {
  // Open input file
  if (parameters->input_file_name==NULL) {
    gem_cond_log(parameters->user_verbose,"[Reading input file from stdin]");
    switch (parameters->input_compression) {
      case FM_GZIPPED_FILE:
        parameters->input_file = input_gzip_stream_open(stdin);
        break;
      case FM_BZIPPED_FILE:
        parameters->input_file = input_bzip_stream_open(stdin);
        break;
      default:
        parameters->input_file = input_stream_open(stdin);
        break;
    }
  } else {
    gem_cond_log(parameters->user_verbose,"[Opening input file '%s']",parameters->input_file_name);
    parameters->input_file = input_file_open(parameters->input_file_name,false);
  }
  // TODO: Checks
}
GEM_INLINE void gem_mapper_close_input(mapper_parameters_t* const parameters) {
  input_file_close(parameters->input_file);
}
GEM_INLINE void gem_mapper_open_output(mapper_parameters_t* const parameters) {
  // Open output stream
  if (parameters->output_file_name==NULL) {
    gem_cond_log(parameters->user_verbose,"[Outputting to stdout]");
    parameters->output_stream = stdout;
  } else {
    gem_cond_log(parameters->user_verbose,"[Outputting to '%s']",parameters->output_file_name);
    parameters->output_stream = fopen(parameters->output_file_name,"w");
  }
  // Open output file // TODO Unsorted
  switch (parameters->output_compression) {
    case FM_GZIPPED_FILE:
      parameters->output_file = output_gzip_stream_new(parameters->output_stream,SORTED_FILE);
      break;
    case FM_BZIPPED_FILE:
      parameters->output_file = output_bzip_stream_new(parameters->output_stream,SORTED_FILE);
      break;
    default:
      parameters->output_file = output_stream_new(parameters->output_stream,SORTED_FILE);
      break;
  }
}
GEM_INLINE void gem_mapper_close_output(mapper_parameters_t* const parameters) {
  output_file_close(parameters->output_file);
  if (parameters->output_file_name!=NULL) fclose(parameters->output_stream);
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
  { 400, "mapping-mode", REQUIRED, TYPE_INT, 4 , false, "'incremental'|'adaptive'|'fixed'|'fast'|'brute-force'|'massive'" , "(default=fast)" },
  { 401, "mapping-degree", REQUIRED, TYPE_INT, 4 , false, "<number|percentage>" , "(default=0)" },
  { 402, "cuda", NO_ARGUMENT, TYPE_NONE, 4, true, "", ""},
  { 'e', "max-search-error", REQUIRED, TYPE_INT, 4 , true, "<number|percentage>" , "(default=0.04, 4%)" },
  { 'E', "max-filtering-error", REQUIRED, TYPE_INT, 4 , false, "<number|percentage>" , "(default=0.2, 20%)" },
  { 's', "complete-strata-after-best", REQUIRED, TYPE_INT, 4 , true, "<number|percentage>" , "(default=0)" },
  { 403, "min-matching-length", REQUIRED, TYPE_INT, 4 , false, "<number|percentage>" , "(default=0.20, 20%)" },
  { 404, "max-search-matches", REQUIRED, TYPE_INT, 4 , true, "<number>" , "(not limited by default)" },
  { 405, "mismatch-alphabet", REQUIRED, TYPE_STRING, 4 , false, "<symbols>" , "(default=\"ACGT\")" },
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
  { 'V', "developer-verbose", NO_ARGUMENT, TYPE_NONE, 9 , false, "" , "(disabled)" },
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
void parse_arguments(int argc,char** argv,mapper_parameters_t* const parameters) {
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
      parameters->index_file_name = optarg;
      break;
    case 204: // --check-index
      parameters->check_index = true;
      break;
    case 'i': // --input
      parameters->input_file_name = optarg;
      break;
    case 200: // --gzip-input
      parameters->input_compression = FM_GZIPPED_FILE;
      break;
    case 201: // --bzip-input
      parameters->input_compression = FM_BZIPPED_FILE;
      break;
    case 'o': // --output
      parameters->output_file_name = optarg;
      break;
    case 202: // --gzip-output
      parameters->output_compression = FM_GZIPPED_FILE;
      break;
    case 203: // --bzip-output
      parameters->output_compression = FM_BZIPPED_FILE;
      break;
    /* Qualities */
    case 'q': // --quality-format
      if (gem_strcaseeq(optarg,"ignore")) {
        parameters->quality_format=qualities_ignore;
      } else if (gem_strcaseeq(optarg,"offset-33")) {
        parameters->quality_format=qualities_offset_33;
      } else if (gem_strcaseeq(optarg,"offset-64")) {
        parameters->quality_format=qualities_offset_64;
      } else {
        gem_fatal_error_msg("Option '-q|--quality-format' must be 'ignore', 'offset-33' or 'offset-64'");
      }
      break;
    case 'Q': // --quality-model
      if (gem_strcaseeq(optarg,"gem")) {
        parameters->quality_model=quality_model_type_gem;
        // TODO parameters->quality_score=gem_quality_score;
      } else if (gem_strcaseeq(optarg,"flat")) {
        parameters->quality_model=quality_model_type_flat;
        // TODO parameters->quality_score=flat_quality_score;
      } else
      gem_fatal_error_msg("Option '-Q|--quality-model' must be 'gem' or 'flat'");
      break;
    case 300: // --gem-quality-threshold
      parameters->quality_threshold = atoll(optarg);
      break;
    /* Single-end Alignment */
    case 400: // --mapping-mode
      if (gem_strcaseeq(optarg,"incremental")) {
        parameters->mapping_mode = mapping_incremental_mapping;
        break;
      }
      if (gem_strcaseeq(optarg,"adaptive")) {
        parameters->mapping_mode = mapping_adaptive_filtering;
        break;
      }
      if (gem_strcaseeq(optarg,"fixed")) {
        parameters->mapping_mode = mapping_fixed_filtering;
        break;
      }
      if (gem_strcaseeq(optarg,"fast")) {
        parameters->mapping_mode = mapping_fast;
        break;
      }
      if (gem_strcaseeq(optarg,"brute-force")) {
        parameters->mapping_mode = mapping_neighborhood_search;
        break;
      }
      if (gem_strcaseeq(optarg,"massive")) {
        parameters->mapping_mode = mapping_massive_filtering;
        break;
      }
      gem_fatal_error_msg("Option '--mapping-mode' must be 'incremental'|'adaptive'|'fixed'|'fast'|'brute-force'|'massive'");
      break;
    case 401: // --mapping-degree
      parameters->mapping_degree = atof(optarg);
      break;
    case 402: // --cuda
      parameters->mapper_type = mapper_se_cuda;
      break;
    case 'e': // --max-search-error
      parameters->max_search_error = atof(optarg);
      break;
    case 'E': // --max-filtering-error
      parameters->max_filtering_error = atof(optarg);
      break;
    case 's': // --complete-strata-after-best
      parameters->complete_strata_after_best = atof(optarg);
      break;
    case 403: // --min-matching-length
      parameters->min_matching_length = atof(optarg);
      break;
    case 404: // --max-search-matches
      parameters->max_search_matches = atol(optarg);
      break;
    case 405: // --mismatch-alphabet
      parameters->mismatch_alphabet = optarg;
      break;
    /* Paired-end Alignment */
    /* Reporting */
    case 'F': // --output-format
      if (gem_strcaseeq(optarg,"MAP")) {
        parameters->output_format = MAP;
      } else if (gem_strcaseeq(optarg,"SAM")) {
        parameters->output_format = SAM;
      } else {
        gem_fatal_error_msg("Option '-F|--output-format' must be 'MAP' or 'SAM'");
      }
      break;
      break;
    case 'D': // --min-decoded-strata
      parameters->min_decoded_strata = (gem_strcaseeq(optarg,"all")) ? ALL : atol(optarg);
      break;
    case 'd': // --max-decoded-matches
      parameters->max_decoded_matches = (gem_strcaseeq(optarg,"all")) ? ALL : atol(optarg);
      break;
    case 'm': // --min-reported-matches
      parameters->min_reported_matches = (gem_strcaseeq(optarg,"all")) ? ALL : atol(optarg);
      break;
    case 'M': // --max-reported-matches
      parameters->max_reported_matches = (gem_strcaseeq(optarg,"all")) ? ALL : atol(optarg);
      break;
    /* System */
    case 't': // --threads
      parameters->num_threads = atol(optarg);
      break;
    case 700: // --max-memory
      gem_cond_fatal_error(input_text_parse_size(optarg,&(parameters->max_memory)),PARSING_SIZE,"--max-memory",optarg);
      break;
    case 701: // --tmp-folder
      parameters->tmp_folder = optarg;
      break;
    /* Debug */
    /* Misc */
    case 'v':
      parameters->user_verbose = true;
      break;
    case 'V':
      parameters->dev_verbose = true;
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
  if (parameters->max_memory==0) {
    parameters->max_memory = mm_get_available_mem();
  }
  // Check min_reported_matches <= max_reported_matches
  // TODO
  // Free
  string_delete(getopt_short_string);
  mm_free(getopt_options);
}
/*
 * Main
 */
int main(int argc,char** argv) {
  // Parsing command-line options
  mapper_parameters_t parameters;
  mapper_cuda_parameters_t cuda_parameters;
  mapper_parameters_set_defaults(&parameters); // Set defaults
  mapper_cuda_parameters_set_defaults(&cuda_parameters); // Set defaults
  parse_arguments(argc,argv,&parameters); // Parse cmd-line

  // Runtime setup
  gem_runtime_init(parameters.max_memory,parameters.tmp_folder,mapper_error_report);

  // Open Input/Output File(s)
  gem_mapper_open_input(&parameters);
  gem_mapper_open_output(&parameters);

  // Load GEM-Index
  gem_mapper_load_index(&parameters);

  // Launch threads
  switch (parameters.mapper_type) {
    case mapper_se:
      mapper_SE_run(&parameters);
      break;
    case mapper_pe:
      GEM_NOT_IMPLEMENTED(); // TODO
      break;
    case mapper_se_cuda:
      // Force Massive-filtering mapping mode
      parameters.mapping_mode = mapping_massive_filtering;
      mapper_SE_CUDA_run(&parameters,&cuda_parameters); // SE-CUDA mapping threads (Producer-Consumer)
      GEM_NOT_IMPLEMENTED(); // TODO
      break;
    case mapper_graph:
      GEM_NOT_IMPLEMENTED(); // TODO
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }

  // CleanUP
  archive_delete(parameters.archive); // Delete archive
  gem_mapper_close_input(&parameters); // Close I/O files
  gem_mapper_close_output(&parameters); // Close I/O files

  return 0;
}
