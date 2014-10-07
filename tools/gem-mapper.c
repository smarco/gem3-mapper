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
  PROF_START_TIMER(GP_MAPPER_LOAD_INDEX);
  // Load archive
  gem_cond_log(parameters->verbose_user,"[Loading GEM index '%s']",parameters->index_file_name);
  parameters->archive = archive_read(parameters->index_file_name,parameters->check_index,parameters->verbose_dev);
  PROF_STOP_TIMER(GP_MAPPER_LOAD_INDEX);
}
GEM_INLINE void gem_mapper_open_input(mapper_parameters_t* const parameters) {
  // Open input file
  if (parameters->input_file_name==NULL) {
    gem_cond_log(parameters->verbose_user,"[Reading input file from stdin]");
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
    gem_cond_log(parameters->verbose_user,"[Opening input file '%s']",parameters->input_file_name);
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
    gem_cond_log(parameters->verbose_user,"[Outputting to stdout]");
    parameters->output_stream = stdout;
  } else {
    gem_cond_log(parameters->verbose_user,"[Outputting to '%s']",parameters->output_file_name);
    parameters->output_stream = fopen(parameters->output_file_name,"w");
  }
  // Open output file
  switch (parameters->output_compression) {
    case FM_GZIPPED_FILE:
      parameters->output_file = output_gzip_stream_new(parameters->output_stream,parameters->max_output_buffers);
      break;
    case FM_BZIPPED_FILE:
      parameters->output_file = output_bzip_stream_new(parameters->output_stream,parameters->max_output_buffers);
      break;
    default:
      parameters->output_file = output_stream_new(parameters->output_stream,parameters->max_output_buffers);
      break;
  }
}
GEM_INLINE void gem_mapper_close_output(mapper_parameters_t* const parameters) {
  output_file_close(parameters->output_file);
}
GEM_INLINE void gem_mapper_print_stats(
    mapper_parameters_t* const parameters,mapper_cuda_parameters_t* const cuda_parameters) {
  PROF_SUM_REDUCE();
  // Sys
  // TODO
  mapper_profile_print_io(gem_info_get_stream());
  // Parameters
  // TODO
  // Mapper
  if (parameters->mapper_type==mapper_se) {
    // Mapper-SE
    switch (parameters->search_parameters.mapping_mode) {
      case mapping_adaptive_filtering:
        mapper_profile_print_mapper_adaptive(gem_info_get_stream());
        break;
      case mapping_incremental_mapping:
      case mapping_fixed_filtering:
      case mapping_fast:
      case mapping_neighborhood_search:
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
  } else if (parameters->mapper_type==mapper_pe) {
  } else if (parameters->mapper_type==mapper_se_cuda) {
    switch (parameters->search_parameters.mapping_mode) {
      case mapping_adaptive_filtering:
        mapper_profile_print_mapper_adaptive(gem_info_get_stream());
        mapper_profile_print_archive_search_group(gem_info_get_stream());
        break;
      case mapping_incremental_mapping:
      case mapping_fixed_filtering:
      case mapping_fast:
      case mapping_neighborhood_search:
      default:
        GEM_INVALID_CASE();
        break;
    }
  } else if (parameters->mapper_type==mapper_graph) {
  }
  // Others
  // TODO
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
#ifdef HAVE_CUDA
  { 402, "cuda", NO_ARGUMENT, TYPE_NONE, 4, true, "", ""},
#endif
  { 'e', "max-search-error", REQUIRED, TYPE_INT, 4 , true, "<number|percentage>" , "(default=0.04, 4%)" },
  { 'E', "max-filtering-error", REQUIRED, TYPE_INT, 4 , false, "<number|percentage>" , "(default=0.2, 20%)" },
  { 's', "complete-strata-after-best", REQUIRED, TYPE_INT, 4 , true, "<number|percentage>" , "(default=0)" },
  { 403, "min-matching-length", REQUIRED, TYPE_INT, 4 , false, "<number|percentage>" , "(default=0.20, 20%)" },
  { 404, "max-search-matches", REQUIRED, TYPE_INT, 4 , true, "<number>" , "(not limited by default)" },
  { 405, "mismatch-alphabet", REQUIRED, TYPE_STRING, 4 , false, "<symbols>" , "(default=\"ACGT\")" },
  /* Paired-end Alignment */

  /* Alignment Score */
  { 600, "alignment-model", REQUIRED, TYPE_STRING, 8 , false, "'none'|'hamming'|'edit'|'gap-affine'" , "(default=edit)" },
  { 'A', "matching-score", REQUIRED, TYPE_INT, 8 , false, "" , "(default=)" },
  { 'B', "mismatch-penalty", REQUIRED, TYPE_INT, 8 , false, "" , "(default=)" },
  { 'O', "gap-open-penalty", REQUIRED, TYPE_INT, 8 , false, "" , "(default=)" },
  { 'X', "gap-extension-penalty", REQUIRED, TYPE_INT, 8 , false, "" , "(default=)" },
  /* Mapping Quality */

  /* Reporting */
  { 'F', "output-format", REQUIRED, TYPE_STRING, 8 , true, "'MAP'|'SAM'" , "(default=MAP)" },
  { 'D', "min-decoded-strata", REQUIRED, TYPE_INT, 8 , false, "<number>" , "(stratum-wise, default=0)" },
  { 'd', "max-decoded-matches", REQUIRED, TYPE_INT, 8 , false, "<number>|'all'" , "(stratum-wise, default=20)" },
  { 'm', "min-reported-matches", REQUIRED, TYPE_INT, 8 , false, "<number>|'all'" , "(default=1)" },
  { 'M', "max-reported-matches", REQUIRED, TYPE_INT, 8 , true, "<number>|'all'" , "(default=all)" },
  /* System */
  { 't', "threads", REQUIRED, TYPE_STRING, 9 , true, "<number>" , "(default=1)" },
  { 900, "max-memory", REQUIRED, TYPE_STRING, 9 , true, "<maximum-memory>" , "(Eg 2GB)" },
  { 901, "tmp-folder", REQUIRED, TYPE_STRING, 9 , true, "<temporal_dir_path>" , "(/tmp/)" },
  /* Debug */
  /* Miscellaneous */
  { 1100, "stats", NO_ARGUMENT, TYPE_NONE, 11 , true, "" , "(disabled)" },
  { 'v', "verbose", NO_ARGUMENT, TYPE_NONE, 11 , true, "" , "(disabled)" },
  { 'V', "developer-verbose", NO_ARGUMENT, TYPE_NONE, 11 , false, "" , "(disabled)" },
  { 'h', "help", NO_ARGUMENT, TYPE_NONE, 11 , true, "" , "(print usage)" },
  { 'H', "help", NO_ARGUMENT, TYPE_NONE, 11 , false, "" , "(print usage + extras)" },
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
  /*  6 */ "Alignment Score",
  /*  7 */ "Mapping Quality",
  /*  8 */ "Reporting",
  /*  9 */ "System",
  /* 10 */ "Debug",
  /* 11 */ "Miscellaneous",
  /* 12 */ "Extras",
};
void usage(const bool print_inactive) {
  fprintf(stderr, "USAGE: ./gem-mapper [ARGS]...\n");
  options_fprint_menu(stderr,gem_mapper_options,gem_mapper_groups,false,print_inactive);
}
void parse_arguments(
    int argc,char** argv,mapper_parameters_t* const parameters,
    mapper_cuda_parameters_t* const cuda_parameters) {
  // Check number of parameters (quick usage & exit)
  if (argc<=1) {
    usage(false);
    exit(0);
  }
  // Set CMD line
  parameters->argc = argc;
  parameters->argv = argv;
  // Parse parameters
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
        parameters->search_parameters.quality_format=qualities_ignore;
      } else if (gem_strcaseeq(optarg,"offset-33")) {
        parameters->search_parameters.quality_format=qualities_offset_33;
      } else if (gem_strcaseeq(optarg,"offset-64")) {
        parameters->search_parameters.quality_format=qualities_offset_64;
      } else {
        gem_fatal_error_msg("Option '-q|--quality-format' must be 'ignore', 'offset-33' or 'offset-64'");
      }
      break;
    case 'Q': // --quality-model
      if (gem_strcaseeq(optarg,"gem")) {
        parameters->search_parameters.quality_model=quality_model_type_gem;
        // TODO parameters->quality_score=gem_quality_score;
      } else if (gem_strcaseeq(optarg,"flat")) {
        parameters->search_parameters.quality_model=quality_model_type_flat;
        // TODO parameters->quality_score=flat_quality_score;
      } else
      gem_fatal_error_msg("Option '-Q|--quality-model' must be 'gem' or 'flat'");
      break;
    case 300: // --gem-quality-threshold
      parameters->search_parameters.quality_threshold = atoll(optarg);
      break;
    /* Single-end Alignment */
    case 400: // --mapping-mode
      if (gem_strcaseeq(optarg,"incremental")) {
        parameters->search_parameters.mapping_mode = mapping_incremental_mapping;
        break;
      }
      if (gem_strcaseeq(optarg,"adaptive")) {
        parameters->search_parameters.mapping_mode = mapping_adaptive_filtering;
        break;
      }
      if (gem_strcaseeq(optarg,"fixed")) {
        parameters->search_parameters.mapping_mode = mapping_fixed_filtering;
        break;
      }
      if (gem_strcaseeq(optarg,"fast")) {
        parameters->search_parameters.mapping_mode = mapping_fast;
        break;
      }
      if (gem_strcaseeq(optarg,"brute-force")) {
        parameters->search_parameters.mapping_mode = mapping_neighborhood_search;
        break;
      }
      gem_fatal_error_msg("Option '--mapping-mode' must be 'incremental'|'adaptive'|'fixed'|'fast'|'brute-force'");
      break;
    case 401: // --mapping-degree
      parameters->search_parameters.mapping_degree = atof(optarg);
      break;
    case 402: // --cuda
      if (!bpm_gpu_support()) GEM_CUDA_NOT_SUPPORTED();
      parameters->mapper_type = mapper_se_cuda;
      break;
    case 'e': // --max-search-error
      parameters->search_parameters.max_search_error = atof(optarg);
      break;
    case 'E': // --max-filtering-error
      parameters->search_parameters.max_filtering_error = atof(optarg);
      break;
    case 's': // --complete-strata-after-best
      parameters->search_parameters.complete_strata_after_best = atof(optarg);
      break;
    case 403: // --min-matching-length
      parameters->search_parameters.min_matching_length = atof(optarg);
      break;
    case 404: // --max-search-matches
      parameters->search_parameters.max_search_matches = atol(optarg);
      break;
    case 405: { // --mismatch-alphabet
      const char* const mismatch_alphabet = optarg;
      approximate_search_configure_replacements(
          &parameters->search_parameters,mismatch_alphabet,gem_strlen(mismatch_alphabet));
      break;
    }
    /* Paired-end Alignment */
    /* Alignment Score */
    case 600: // --alignment-model
      if (gem_strcaseeq(optarg,"none")) {
        parameters->select_parameters.alignment_model = alignment_model_none;
      } else if (gem_strcaseeq(optarg,"hamming")) {
        parameters->select_parameters.alignment_model = alignment_model_hamming;
      } else if (gem_strcaseeq(optarg,"edit") || gem_strcaseeq(optarg,"levenshtein") ) {
        parameters->select_parameters.alignment_model = alignment_model_levenshtein;
      } else if (gem_strcaseeq(optarg,"gap-affine'")) {
        parameters->select_parameters.alignment_model = alignment_model_gap_affine;
      } else {
        gem_fatal_error_msg("Option '--alignment-model' must be 'none'|'hamming'|'edit'|'gap-affine'");
      }
      break;
    /* Mapping Quality */
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
    case 'D': // --min-decoded-strata
      parameters->select_parameters.min_decoded_strata = (gem_strcaseeq(optarg,"all")) ? ALL : atol(optarg);
      break;
    case 'd': // --max-decoded-matches
      parameters->select_parameters.max_decoded_matches = (gem_strcaseeq(optarg,"all")) ? ALL : atol(optarg);
      break;
    case 'm': // --min-reported-matches
      parameters->select_parameters.min_reported_matches = (gem_strcaseeq(optarg,"all")) ? ALL : atol(optarg);
      break;
    case 'M': // --max-reported-matches
      parameters->select_parameters.max_reported_matches = (gem_strcaseeq(optarg,"all")) ? ALL : atol(optarg);
      break;
    /* System */
    case 't': // --threads
      parameters->num_threads = atol(optarg);
      break;
    case 900: // --max-memory
      gem_cond_fatal_error(input_text_parse_size(optarg,&(parameters->max_memory)),PARSING_SIZE,"--max-memory",optarg);
      break;
    case 901: // --tmp-folder
      parameters->tmp_folder = optarg;
      break;
    /* Debug */
    /* Misc */
    case 1100: // --stats
      parameters->stats = true;
      break;
    case 'v':
      parameters->verbose_user = true;
      break;
    case 'V':
      parameters->verbose_user = true;
      parameters->verbose_dev = true;
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
  // I/O Parameters
  gem_cond_fatal_error_msg(parameters->index_file_name==NULL,"Index file required");
  gem_cond_fatal_error_msg(parameters->input_file_name==NULL,"Input file required");
  gem_cond_fatal_error_msg(parameters->output_file_name!=NULL &&
      gem_streq(parameters->input_file_name,parameters->output_file_name),
      "Input-file and output-file must be different");
  gem_cond_fatal_error_msg(gem_streq(parameters->index_file_name,parameters->input_file_name),
      "Index-file and input-file must be different");
  gem_cond_fatal_error_msg(parameters->output_file_name!=NULL &&
      gem_streq(parameters->index_file_name,parameters->output_file_name),
      "Index-file and output-file must be different");
  /*
   * Search Parameters
   */
  /* Mapping strategy (Mapping mode + properties) */
  /* Qualities */
//  uint64_t quality_threshold;
  gem_cond_fatal_error_msg(parameters->search_parameters.quality_threshold > 94,
      "Quality threshold is too high (please consider lowering it or run ignoring qualities)");
  gem_cond_fatal_error_msg(parameters->search_parameters.quality_threshold == 0,
      "Quality threshold is zero (all base-calls will be considered wrong))");
  /* Error Model */
//  gem_cond_fatal_error_msg(parameters->search_parameters.max_search_error >
//      parameters->search_parameters.max_filtering_error,
//      "Option '--max-search-error' is greater than 'max-filtering-error' [Inefficiency]"); // FIXME
  /* Matches search */
  /* Replacements (Regulates the bases that can be replaced/mismatched) */
  /* Alignment Score */
  /*
   * Select Parameters
   */
  /* Single-end Alignment */
  /* Paired-end Alignment */
  /* Alignment Score */
  /* Reporting */
  gem_cond_fatal_error_msg(parameters->select_parameters.min_reported_matches >
      parameters->select_parameters.max_reported_matches,
      "Option '--max-reported-matches' must be greater or equal than 'min-reported-matches'");
  // System
  /*
   * TODO
   *   check path (gem_access(char* const path,const fm_mode mode))
   *   check number of cpus (force to increase it)
   */
  // Free
  string_destroy(getopt_short_string);
  mm_free(getopt_short_string);
  mm_free(getopt_options);
}
/*
 * Main
 */
int main(int argc,char** argv) {
  // Mapper timer
  gem_timer_t mapper_timer;
  TIMER_RESTART(&mapper_timer);
  // Parsing command-line options
  mapper_parameters_t parameters;
  mapper_cuda_parameters_t cuda_parameters;
  mapper_parameters_set_defaults(&parameters); // Set defaults
  mapper_cuda_parameters_set_defaults(&cuda_parameters);
  parse_arguments(argc,argv,&parameters,&cuda_parameters); // Parse cmd-line

  // Runtime setup
  const uint64_t total_threads = (parameters.mapper_type==mapper_se_cuda) ?
      cuda_parameters.num_generating_threads + cuda_parameters.num_selecting_threads + 1:
      parameters.num_threads + 1;
  gem_runtime_init(total_threads,parameters.max_memory,parameters.tmp_folder,mapper_error_report);
  PROF_START_TIMER(GP_MAPPER_ALL);

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
      mapper_SE_CUDA_run(&parameters,&cuda_parameters); // SE-CUDA mapping threads (Producer-Consumer)
      break;
    case mapper_graph:
      GEM_NOT_IMPLEMENTED(); // TODO
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  PROF_STOP_TIMER(GP_MAPPER_ALL);

  // Stats
  if (parameters.stats) gem_mapper_print_stats(&parameters,&cuda_parameters);

  // CleanUP
  archive_delete(parameters.archive); // Delete archive
  gem_mapper_close_input(&parameters); // Close I/O files
  gem_mapper_close_output(&parameters); // Close I/O files
  gem_runtime_destroy();

  // Display end banner
  TIMER_STOP(&mapper_timer);
  gem_cond_log(parameters.verbose_user,
      "[GEMMapper terminated successfully in %2.3f min.]\n",TIMER_GET_TOTAL_M(&mapper_timer));

  return 0;
}
